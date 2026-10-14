"""``prepare.*`` step implementations: genotypes, references, datasets, mappings.

Each function reads only its declared inputs and params, writes every declared output
atomically, and turns :class:`anqchor.io.InputError` into :class:`StepError`.
File formats are documented on :class:`anqchor.pipelines.layout.Layout`.
"""

from __future__ import annotations

import functools
import time
from collections import Counter
from typing import Any

import pandas as pd

from anqchor.io.cov import read_and_encode, write_cov_raw
from anqchor.io.genotypes import (
    check_list_style,
    chrom_style,
    entries_from_pfile,
    first_pvar_chrom,
    link_standard_prefix,
    missing_files,
    parse_pgen_list,
    read_prepared_pgen_list,
    read_psam_iids,
)
from anqchor.io.maps import (
    best_rule,
    filter_sample_map,
    map_traits_by_file,
    map_traits_by_rule,
    read_sample_map,
    read_trait_map,
)
from anqchor.io.pairs import (
    Variant,
    check_presence,
    chrom_order,
    read_pairs_table,
    read_reference_pairs,
    read_variant_list,
    sort_variants,
    write_pairs,
)
from anqchor.io.pheno import (
    annotate_traits,
    convert_to_regenie,
    read_trait_annotation,
    write_traits,
)
from anqchor.io.tables import InputError, duplicates, examples, read_lines, write_lines
from anqchor.pipelines.layout import Layout
from anqchor.steps import StepContext, StepError, register
from anqchor.util.files import read_header, read_json, write_json

SCHEMA_VERSION = 1


def input_errors(fn: Any) -> Any:
    """Re-raise :class:`InputError` as :class:`StepError` (no traceback for data problems)."""

    @functools.wraps(fn)
    def wrapper(ctx: StepContext) -> None:
        try:
            fn(ctx)
        except InputError as e:
            raise StepError(str(e)) from e

    return wrapper


def _log(ctx: StepContext, msg: str, *args: object) -> None:
    if ctx.log is not None:
        ctx.log.info(msg, *args)


# ------------------------------------------------------------------- genotypes
@register("prepare_genotypes")
@input_errors
def prepare_genotypes(ctx: StepContext) -> None:
    """Validate the genotype files; write the absolute pgen list, samples and summary."""
    p = ctx.params
    if p.get("pgen_list"):
        entries = parse_pgen_list(p["pgen_list"])
        source = {"pgen_list": p["pgen_list"]}
    elif p.get("pfile"):
        entries = entries_from_pfile(p["pfile"])
        source = {"pfile": p["pfile"]}
    else:
        raise StepError("prepare.genotypes needs 'pgen_list' or 'pfile'")
    missing = [str(f) for e in entries for f in missing_files(e)]
    if missing:
        hint = ""
        if any(m.endswith(".pvar.idx.gz") or m.endswith(".tbi") for m in missing):
            hint = " Create indexed pvar files with 'anqchor prepare index-pvar --pfile PREFIX'."
        raise StepError(f"genotype files not found: {examples(missing, 6)}.{hint}")
    first = first_pvar_chrom(entries[0].pivar)
    style = chrom_style(first)
    check_list_style(entries, style)
    samples: list[str] | None = None
    psam_of_samples = ""
    for psam in dict.fromkeys(e.psam for e in entries):
        iids = read_psam_iids(psam)
        if samples is None:
            dup = duplicates(iids)
            if dup:
                raise StepError(f"{psam}: duplicate sample IIDs {examples(dup)}")
            samples, psam_of_samples = iids, str(psam)
        elif iids != samples:
            only_a = [s for s in samples if s not in set(iids)]
            only_b = [s for s in iids if s not in set(samples)]
            detail = (
                f"only in {psam_of_samples}: {examples(only_a, 3)}; only in {psam}: {examples(only_b, 3)}"
                if only_a or only_b
                else "same IDs in a different order"
            )
            raise StepError(
                f"genotype files must share one sample list, but {psam} differs from "
                f"{psam_of_samples} ({detail})"
            )
    assert samples is not None
    lines = []
    for i, e in enumerate(entries):
        prefix = e.standard_prefix()
        if prefix is None:  # non-standard names: link them under one prefix (4-column list)
            prefix = link_standard_prefix(
                e, ctx.path(Layout.geno_dir) / "links" / f"{i + 1}.{e.chrom}"
            )
        lines.append(e.line(prefix))
    write_lines(ctx.path(Layout.pgen_list), lines)
    write_lines(ctx.path(Layout.geno_samples), samples)
    write_json(
        ctx.path(Layout.geno_json),
        {
            "schema_version": SCHEMA_VERSION,
            **source,
            "n_samples": len(samples),
            "chrom_style": style,
            "chromosomes": list(dict.fromkeys(e.chrom for e in entries)),
            "entries": [e.as_json() for e in entries],
            "n_entries": len(entries),
        },
    )
    _log(
        ctx,
        "%d genotype file set(s), %d samples, chromosome style %s",
        len(entries),
        len(samples),
        style,
    )


# ------------------------------------------------------------------ references
@register("prepare_reference")
@input_errors
def prepare_reference(ctx: StepContext) -> None:
    """Normalise a reference and check which of its variants exist in the genotypes."""
    p = ctx.params
    ref = p["reference"]
    geno = read_json(ctx.path(Layout.geno_json))
    style = geno["chrom_style"]
    entries = read_prepared_pgen_list(ctx.path(Layout.pgen_list))
    columns: dict[str, Any] = dict(p.get("columns") or {})
    counts: dict[str, Any] = {"n_pairs_input": 0, "n_pairs": 0, "n_invalid": 0, "n_duplicates": 0}
    variants: dict[str, Variant] = {}
    renamed = 0
    table = pd.DataFrame(columns=["trait", "variant", "beta", "se", "pval"])
    invalid_examples: list[str] = []
    if p.get("pairs"):
        rp = read_reference_pairs(p["pairs"], columns, style)
        table = rp.table
        variants.update(rp.variants)
        renamed += rp.n_renamed
        counts.update(
            n_pairs_input=rp.n_input,
            n_pairs=len(rp.table),
            n_invalid=rp.n_invalid,
            n_duplicates=rp.n_duplicates,
        )
        invalid_examples = rp.invalid_examples
        if rp.n_input and not len(rp.table):
            raise StepError(
                f"reference {ref}: none of the {rp.n_input} pairs is valid (variant IDs must be "
                f"CHROM:POS:REF:ALT and beta/se finite numbers), e.g. {examples(rp.invalid_examples, 3)}"
            )
    n_list = n_list_invalid = 0
    if p.get("variants"):
        lv, n_list_invalid, n_ren = read_variant_list(p["variants"], style)
        n_list = len(lv)
        for k, v in lv.items():
            variants.setdefault(k, v)
        renamed += n_ren
    if not variants:
        raise StepError(f"reference {ref}: no valid variants")
    present, method = check_presence(list(variants.values()), entries)
    order = chrom_order(entries)
    ordered = sort_variants(variants.values(), order)
    rows = [
        f"{v.cpra}\t{v.chrom}\t{v.pos}\t{v.ref}\t{v.alt}\t{int(v.is_indel)}\t{int(v in present)}"
        for v in ordered
    ]
    n_present = len(present)
    absent = [v for v in ordered if v not in present]
    if n_present == 0:
        ref_chroms = sorted({v.chrom for v in variants.values()})
        raise StepError(
            f"reference {ref}: none of its {len(variants)} variants exists in the genotypes. "
            f"Reference chromosomes (after normalisation): {examples(ref_chroms, 5)}; genotype "
            f"chromosomes: {examples(list(order), 5)}. Check the chromosome naming, the genome "
            f"build (e.g. GRCh37 vs GRCh38) and the allele order; first variants: "
            f"{examples([v.cpra for v in ordered], 3)}"
        )
    per_chrom: dict[str, dict[str, int]] = {}
    pair_chrom = Counter(variants[v].chrom for v in table["variant"]) if len(table) else Counter()
    for v in ordered:
        d = per_chrom.setdefault(
            v.chrom, {"n_pairs": pair_chrom.get(v.chrom, 0), "n_variants": 0, "n_present": 0}
        )
        d["n_variants"] += 1
        d["n_present"] += int(v in present)
    write_pairs(ctx.path(Layout.ref_pairs(ref)), table, ["trait", "variant", "beta", "se", "pval"])
    write_lines(
        ctx.path(Layout.ref_variants(ref)),
        ["variant\tchrom\tpos\tref\talt\tis_indel\tpresent", *rows],
    )
    write_lines(
        ctx.path(Layout.ref_variants_present(ref)), [v.cpra for v in ordered if v in present]
    )
    write_json(
        ctx.path(Layout.ref_json(ref)),
        {
            "schema_version": SCHEMA_VERSION,
            "reference": ref,
            "kind": p.get("kind", "cis"),
            **counts,
            "n_traits": int(table["trait"].nunique()) if len(table) else 0,
            "n_variants": len(variants),
            "n_variants_list": n_list,
            "n_variants_list_invalid": n_list_invalid,
            "n_present": n_present,
            "n_absent": len(absent),
            "n_absent_snv": sum(1 for v in absent if not v.is_indel),
            "n_absent_indel": sum(1 for v in absent if v.is_indel),
            "absent_examples": [v.cpra for v in absent[:20]],
            "invalid_examples": invalid_examples,
            "per_chrom": per_chrom,
            "chrom_renamed": renamed,
            "chrom_style": style,
            "presence_method": method,
        },
    )
    _log(
        ctx,
        "reference %s: %d pairs, %d variants, %d present (%s)",
        ref,
        len(table),
        len(variants),
        n_present,
        method,
    )


# -------------------------------------------------------------------- datasets
@register("prepare_dataset")
@input_errors
def prepare_dataset(ctx: StepContext) -> None:
    """Convert a dataset to regenie layout, annotate traits, encode covariates, map samples."""
    p = ctx.params
    ds = p["dataset"]
    t0 = time.time()
    geno = read_json(ctx.path(Layout.geno_json))
    style = geno["chrom_style"]
    geno_samples = read_lines(ctx.path(Layout.geno_samples))
    gset = set(geno_samples)
    missing = list(p.get("missing_str") or ["NA"])
    fmt = p.get("pheno_format", "regenie")
    conv = convert_to_regenie(p["pheno"], fmt, ctx.path(Layout.ds_pheno(ds)), missing)
    m = conv.matrix
    samples, traits = m.samples, m.features
    t_pheno = time.time() - t0
    # trait loci
    annot = None
    annot_source = None
    if p.get("trait_annot"):
        annot = read_trait_annotation(p["trait_annot"])
        annot_source = "trait_annot"
    elif m.loci is not None:
        annot = m.loci.dropna(subset=["start", "end"])
        annot_source = "pheno"
    loci = annotate_traits(traits, annot, style, p.get("trait_id_rule") or "as_is")
    write_traits(ctx.path(Layout.ds_traits(ds)), loci.table)
    # covariates
    cov_info: dict[str, Any] = {"has_cov": False, "n_cov_columns": 0}
    if p.get("cov"):
        enc = read_and_encode(
            p["cov"],
            p.get("cov_format", "regenie"),
            p.get("cov_categorical") or [],
            p.get("cov_missing", "error"),
            missing,
        )
        if enc.values.shape[1] == 0:
            raise StepError(
                f"dataset {ds}: no usable covariate column left in {p['cov']} (all constant?)"
            )
        assert enc.raw is not None
        write_cov_raw(ctx.path(Layout.ds_cov(ds)), enc.raw)
        write_json(ctx.path(Layout.ds_cov_json(ds)), enc.as_json())
        cov_samples = set(enc.values.index)
        n_without = sum(1 for s in samples if s not in cov_samples)
        if n_without == len(samples):
            raise StepError(
                f"dataset {ds}: no phenotype sample has covariates; covariate sample IDs look like "
                f"{examples(list(enc.values.index), 3)}, phenotype IDs like {examples(samples, 3)}"
            )
        cov_info = {
            "has_cov": True,
            "n_cov_columns": enc.values.shape[1],
            "cov_columns": enc.columns,
            "n_cov_samples": enc.values.shape[0],
            "n_pheno_without_cov": n_without,
            "cov_categorical": enc.categorical,
            "cov_dropped_constant": enc.dropped_constant,
            "cov_missing_policy": enc.missing_policy,
            "cov_n_missing_handled": enc.n_missing_handled,
        }
    ds_warnings: list[str] = []
    if cov_info.get("n_pheno_without_cov"):
        ds_warnings.append(
            f"{cov_info['n_pheno_without_cov']} of {len(samples)} phenotyped sample(s) have no covariates "
            "and are left out of every association with covariates"
        )
    if cov_info.get("cov_n_missing_handled"):
        how = (
            "dropped"
            if cov_info.get("cov_missing_policy") == "drop_samples"
            else "imputed with the column mean"
        )
        what = (
            "sample(s) with a missing covariate"
            if how == "dropped"
            else "missing covariate value(s)"
        )
        ds_warnings.append(f"{cov_info['cov_n_missing_handled']} {what} {how} (cov_missing)")
    # sample map
    smap_info: dict[str, Any]
    if p.get("sample_map"):
        rows = read_sample_map(p["sample_map"])
        res = filter_sample_map(rows, gset, samples)
        pairs = res.pairs
        smap_info = {
            "sample_map": p["sample_map"],
            "n_sample_map_rows": res.n_rows,
            "n_sample_map_rows_geno_missing": res.n_rows_geno_missing,
            "n_sample_map_rows_pheno_missing": res.n_rows_pheno_missing,
            "n_geno_with_several_pheno": res.n_geno_with_several_pheno,
        }
    else:
        pairs = [(s, s) for s in samples if s in gset]
        smap_info = {"sample_map": None}
    mapped = {ph for _, ph in pairs}
    unmapped = [s for s in samples if s not in mapped]
    write_lines(ctx.path(Layout.ds_sample_map(ds)), (f"{g}\t{ph}" for g, ph in pairs))
    write_lines(ctx.path(Layout.ds_samples(ds)), samples)
    n_overlap = sum(1 for s in samples if s in gset)
    cov_set = set(enc.values.index) if p.get("cov") else None
    n_analysable = sum(1 for s in samples if s in gset and (cov_set is None or s in cov_set))
    if p.get("sample_map") and n_overlap < len(pairs):
        ds_warnings.append(
            f"sample_map is used by match only: cis and trans match phenotype and genotype IDs "
            f"directly, and only {n_overlap} of the {len(pairs)} mapped sample(s) have identical IDs"
        )
    write_json(
        ctx.path(Layout.ds_json(ds)),
        {
            "schema_version": SCHEMA_VERSION,
            "dataset": ds,
            "omics": p.get("omics", "omics"),
            "pheno": p["pheno"],
            "pheno_format": fmt,
            "pheno_fast_path": conv.fast_path,
            "n_samples_pheno": len(samples),
            "n_traits": len(traits),
            "n_missing_values": conv.n_missing,
            "n_samples_geno": len(geno_samples),
            "n_samples_geno_overlap": n_overlap,
            "n_samples_mapped": len(pairs),
            "n_samples_analysable": n_analysable,
            "warnings": ds_warnings,
            "n_pheno_unmapped": len(unmapped),
            "pheno_unmapped_examples": unmapped[:20],
            **smap_info,
            **cov_info,
            "has_loci": loci.n_with_loci > 0,
            "n_traits_with_loci": loci.n_with_loci,
            "trait_annot_source": annot_source,
            "trait_annot_rule": loci.rule,
            "trait_id_rule": p.get("trait_id_rule", "as_is"),
            "seconds_pheno": round(t_pheno, 2),
            "seconds_total": round(time.time() - t0, 2),
        },
    )
    if not pairs:
        _log(ctx, "WARNING dataset %s: no phenotype sample matches a genotyped sample", ds)
    for w in ds_warnings:
        _log(ctx, "WARNING dataset %s: %s", ds, w)
    _log(
        ctx,
        "dataset %s: %d samples x %d traits (%s), %d mapped to genotypes, %d analysable by cis/trans "
        "(genotypes, phenotypes and covariates), %d traits with loci, %.1fs",
        ds,
        len(samples),
        len(traits),
        "copied" if conv.fast_path else "converted",
        len(pairs),
        n_analysable,
        loci.n_with_loci,
        time.time() - t0,
    )


# -------------------------------------------------------------------- mappings
@register("prepare_mapping")
@input_errors
def prepare_mapping(ctx: StepContext) -> None:
    """Map a reference's traits onto a dataset's traits and rewrite its pairs."""
    p = ctx.params
    ds, ref = p["dataset"], p["reference"]
    ds_traits = read_header(ctx.path(Layout.ds_pheno(ds)))[2:]
    pairs = read_pairs_table(ctx.path(Layout.ref_pairs(ref)))
    variants = read_pairs_table(ctx.path(Layout.ref_variants(ref)))
    present = dict(zip(variants["variant"], variants["present"], strict=True))
    ref_traits = list(dict.fromkeys(pairs["trait"].tolist()))
    rule = p.get("trait_id_rule") or "as_is"
    if p.get("trait_map"):
        tm = map_traits_by_file(
            ref_traits, ds_traits, read_trait_map(p["trait_map"]), p["trait_map"]
        )
    else:
        tm = map_traits_by_rule(ref_traits, ds_traits, rule)
    if not tm.mapping:
        suggestion, n = best_rule(ref_traits, ds_traits)
        hint = (
            f" trait_id_rule {suggestion!r} would map {n} trait(s)."
            if n > 0 and suggestion != rule
            else " No trait_id_rule matches; provide a trait_map file (ref_trait pheno_trait)."
        )
        raise StepError(
            f"no trait of reference {ref} maps onto dataset {ds} "
            f"({'trait_map ' + str(p['trait_map']) if p.get('trait_map') else 'trait_id_rule ' + rule}). "
            f"Reference traits look like {examples(ref_traits, 3)}; dataset traits look like "
            f"{examples(ds_traits, 3)}.{hint}"
        )
    sub = pairs[pairs["trait"].isin(tm.mapping.keys())].copy()
    sub["ref_trait"] = sub["trait"]
    sub["trait"] = sub["ref_trait"].map(tm.mapping)
    sub["present"] = sub["variant"].map(lambda v: present.get(v, "0"))
    write_lines(
        ctx.path(Layout.dsref_trait_map(ds, ref)), (f"{r}\t{d}" for r, d in tm.mapping.items())
    )
    write_pairs(
        ctx.path(Layout.dsref_pairs(ds, ref)),
        sub,
        ["trait", "variant", "beta", "se", "pval", "ref_trait", "present"],
    )
    used = set(sub.loc[sub["present"] == "1", "variant"])
    write_lines(
        ctx.path(Layout.dsref_variants(ds, ref)), [v for v in variants["variant"] if v in used]
    )
    mapped_ds = set(tm.mapping.values())
    write_lines(ctx.path(Layout.dsref_traits(ds, ref)), [t for t in ds_traits if t in mapped_ds])
    write_json(
        ctx.path(Layout.dsref_json(ds, ref)),
        {
            "schema_version": SCHEMA_VERSION,
            "dataset": ds,
            "reference": ref,
            "rule": tm.rule,
            "trait_map": p.get("trait_map"),
            "n_ref_traits": len(ref_traits),
            "n_dataset_traits": len(ds_traits),
            "n_mapped_traits": len(tm.mapping),
            "n_unmapped_traits": len(tm.unmapped),
            "unmapped_examples": tm.unmapped[:20],
            "n_pairs": len(pairs),
            "n_pairs_mapped": len(sub),
            "n_pairs_present": int((sub["present"] == "1").sum()),
            "n_variants_present": len(used),
            "n_traits_with_present_variant": int(sub.loc[sub["present"] == "1", "trait"].nunique()),
        },
    )
    _log(
        ctx,
        "%s -> %s: %d of %d reference traits mapped (%s)",
        ref,
        ds,
        len(tm.mapping),
        len(ref_traits),
        tm.rule,
    )
