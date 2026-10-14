"""Steps shared by several pipelines: input preparation, strategy inputs and PRS.

Steps are identified by name; :class:`~anqchor.engine.StepGraph` de-duplicates
identical definitions, so every pipeline can simply ask for the steps it needs.
"""

from __future__ import annotations

import functools
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

from anqchor.config.model import InputSetConfig
from anqchor.engine.step import Step
from anqchor.pipelines.context import BuildContext
from anqchor.qpgen import commands as qc
from anqchor.util.files import open_text, resolve_near


def read_pgen_list_files(pgen_list: Path) -> list[Path]:
    """Genotype files (``.pgen``, ``.psam``, ``.pvar.idx.gz``) named in a pgen list,
    resolved like ``prepare.genotypes`` does (:func:`anqchor.util.files.resolve_near`).
    Used only to declare them as make prerequisites, so replacing any of them re-runs
    every step that reads genotypes; missing files are reported by the step."""
    files: list[Path] = []
    try:
        with open_text(pgen_list) as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                f = s.split()
                if len(f) == 4:
                    prefix = str(resolve_near(f[3], pgen_list))
                    paths = [Path(prefix + suf) for suf in (".pgen", ".psam", ".pvar.idx.gz")]
                elif len(f) == 6:
                    paths = [resolve_near(n, pgen_list) for n in f[3:6]]
                else:
                    continue
                for p in paths:
                    if p.exists() and p not in files:
                        files.append(p)
    except OSError:
        pass
    return files


def prepare_genotypes(ctx: BuildContext) -> Step:
    g = ctx.cfg.genotypes
    L = ctx.layout
    inputs: list[Path] = []
    if g.pgen_list is not None:
        inputs.append(g.pgen_list)
        inputs += read_pgen_list_files(g.pgen_list)
    else:
        assert g.pfile is not None
        inputs += [Path(str(g.pfile) + s) for s in (".pgen", ".psam", ".pvar.idx.gz")]
    return Step(
        name="prepare.genotypes",
        kind="prepare_genotypes",
        inputs=tuple(inputs),
        outputs=(L.pgen_list, L.geno_samples, L.geno_json),
        params={
            "pgen_list": str(g.pgen_list) if g.pgen_list else None,
            "pfile": str(g.pfile) if g.pfile else None,
        },
        description="validate genotype files; write an absolute pgen list and the sample list",
    )


def prepare_reference(ctx: BuildContext, ref: str) -> Step:
    r = ctx.cfg.references[ref]
    L = ctx.layout
    inputs = [p for p in (r.pairs, r.variants) if p is not None]
    return Step(
        name=f"prepare.reference.{ref}",
        kind="prepare_reference",
        inputs=(*inputs, L.pgen_list, L.geno_json),
        outputs=(
            L.ref_pairs(ref),
            L.ref_variants(ref),
            L.ref_variants_present(ref),
            L.ref_json(ref),
        ),
        params={"reference": ref, **r.model_dump(mode="json")},
        description=f"normalise reference {ref} and check which variants exist in the genotypes",
    )


def prepare_dataset(ctx: BuildContext, ds: str) -> Step:
    d = ctx.cfg.datasets[ds]
    L = ctx.layout
    inputs = [p for p in (d.pheno, d.cov, d.trait_annot, d.sample_map) if p is not None]
    outputs = [
        L.ds_pheno(ds),
        L.ds_traits(ds),
        L.ds_samples(ds),
        L.ds_sample_map(ds),
        L.ds_json(ds),
    ]
    if d.cov is not None:
        outputs += [L.ds_cov(ds), L.ds_cov_json(ds)]
    return Step(
        name=f"prepare.dataset.{ds}",
        kind="prepare_dataset",
        inputs=(*inputs, L.geno_samples, L.geno_json),
        outputs=tuple(outputs),
        params={"dataset": ds, **d.model_dump(mode="json")},
        description=f"convert dataset {ds} to regenie layout, check covariates, annotate traits",
    )


def prepare_mapping(ctx: BuildContext, ds: str, ref: str) -> Step:
    d = ctx.cfg.datasets[ds]
    L = ctx.layout
    inputs = [L.ds_pheno(ds), L.ds_traits(ds), L.ds_json(ds), L.ref_pairs(ref), L.ref_variants(ref)]
    if d.trait_map is not None:
        inputs.append(d.trait_map)
    return Step(
        name=f"prepare.mapping.{ds}.{ref}",
        kind="prepare_mapping",
        inputs=tuple(inputs),
        outputs=(
            L.dsref_trait_map(ds, ref),
            L.dsref_pairs(ds, ref),
            L.dsref_variants(ds, ref),
            L.dsref_traits(ds, ref),
            L.dsref_json(ds, ref),
        ),
        params={
            "dataset": ds,
            "reference": ref,
            "trait_id_rule": d.trait_id_rule,
            "trait_map": str(d.trait_map) if d.trait_map else None,
        },
        description=f"map reference {ref} traits onto dataset {ds}",
    )


def prs(ctx: BuildContext, ref: str) -> Step:
    L = ctx.layout
    m = ctx.cfg.pipelines.match
    inputs: list[Path] = [L.pgen_list, L.ref_pairs(ref)]
    if m.prs_samples is not None:
        inputs.append(m.prs_samples)
    argv = qc.pair_prs(
        ctx.qpgen,
        genotypes=ctx.genotypes(),
        pairs=L.ref_pairs(ref),
        out_prefix=L.prs_prefix(ref),
        sample=m.prs_samples,
    )
    return Step(
        name=f"prs.{ref}",
        commands=(qc.to_shell(argv),),
        inputs=tuple(inputs),
        outputs=(L.prs_matrix(ref), L.prs_se(ref)),
        pipeline="match",
        description=f"polygenic scores from reference {ref} (qpgentools pair-prs)",
    )


def strategy_inputs(
    ctx: BuildContext, pipeline: Literal["cis", "trans"], ds: str, strategy: InputSetConfig
) -> Step:
    """Covariates (and optional phenotype override) for one strategy."""
    d = ctx.cfg.datasets[ds]
    L = ctx.layout
    inputs: list[Path] = [L.ds_json(ds), L.ds_samples(ds)]
    outputs: list[Path] = [L.strategy_json(pipeline, ds, strategy.name)]
    cov = strategy.cov
    if cov == "default":
        if d.cov is not None:
            inputs += [L.ds_cov(ds), L.ds_cov_json(ds)]
            outputs.append(L.strategy_cov(pipeline, ds, strategy.name))
    elif cov is not None:
        if cov.include is not None:
            if d.cov is None:
                raise ValueError(
                    f"strategy {strategy.name}: cov.include needs covariates in dataset {ds}"
                )
            inputs += [L.ds_cov(ds), L.ds_cov_json(ds)]
        else:
            assert cov.file is not None
            inputs.append(cov.file)
        outputs.append(L.strategy_cov(pipeline, ds, strategy.name))
    if strategy.pheno is not None:
        inputs.append(strategy.pheno)
        outputs.append(L.strategy_pheno(pipeline, ds, strategy.name))
    return Step(
        name=f"{pipeline}.{ds}.{strategy.name}.inputs",
        kind="strategy_inputs",
        inputs=tuple(inputs),
        outputs=tuple(outputs),
        params={
            "pipeline": pipeline,
            "dataset": ds,
            "strategy": strategy.model_dump(mode="json"),
            "dataset_cov_categorical": d.cov_categorical,
            "dataset_cov_missing": d.cov_missing,
            "missing_str": d.missing_str,
        },
        pipeline=pipeline,
        description=f"covariates/phenotypes for {pipeline} strategy {strategy.name} on {ds}",
    )


def strategy_files(
    ctx: BuildContext, pipeline: Literal["cis", "trans"], ds: str, strategy: InputSetConfig
) -> tuple[Path, Path | None]:
    """(phenotype matrix, covariate file or None) a strategy's association runs use."""
    step = strategy_inputs(ctx, pipeline, ds, strategy)
    L = ctx.layout
    pheno = L.strategy_pheno(pipeline, ds, strategy.name) if strategy.pheno else L.ds_pheno(ds)
    cov_path = L.strategy_cov(pipeline, ds, strategy.name)
    cov = cov_path if cov_path in step.outputs else None
    return pheno, cov


def shared_prepare_steps(ctx: BuildContext, datasets: list[str], refs: list[str]) -> list[Step]:
    """Genotype, reference, dataset and mapping preparation for the given names."""
    steps = [prepare_genotypes(ctx)]
    steps += [prepare_reference(ctx, r) for r in refs]
    steps += [prepare_dataset(ctx, d) for d in datasets]
    steps += [prepare_mapping(ctx, d, r) for d in datasets for r in refs]
    return steps


@functools.lru_cache(maxsize=32)
def _categorical_levels(
    path: str, fmt: str, missing: tuple[str, ...], declared: tuple[str, ...]
) -> tuple[tuple[str, int], ...]:
    from anqchor.io.cov import categorical_levels
    from anqchor.io.tables import InputError

    try:
        return tuple(categorical_levels(path, fmt, list(missing), list(declared)).items())
    except (InputError, OSError):
        return tuple((n, 0) for n in declared)  # the prepare step reports the problem with the file


def _dataset_categorical(
    path: str, fmt: str, missing: tuple[str, ...], declared: tuple[str, ...]
) -> tuple[str, ...]:
    return tuple(n for n, _ in _categorical_levels(path, fmt, missing, declared))


QPGEN_MAX_CAT_LEVELS = 10
"""qpgen's default ``--max-cat-levels``."""


def max_cat_levels(
    ctx: BuildContext, ds: str, cats: Sequence[str], strategy: InputSetConfig | None = None
) -> int | None:
    """``--max-cat-levels`` for a qpgen call using categorical covariates ``cats``.

    anqchor validates categorical covariates itself (prepare step), so qpgen must accept
    every level count anqchor accepts. Returns the largest level count among ``cats``
    when it exceeds qpgen's default of 10, else ``None`` (flag omitted)."""
    if not cats:
        return None
    d = ctx.cfg.datasets[ds]
    cov = strategy.cov if strategy is not None else "default"
    if cov is not None and cov != "default" and cov.file is not None:
        src = (str(cov.file), cov.format)
    elif d.cov is not None:
        src = (str(d.cov), d.cov_format)
    else:
        return None
    levels = dict(
        _categorical_levels(src[0], src[1], tuple(d.missing_str), tuple(d.cov_categorical))
    )
    n = max((levels.get(c, 0) for c in cats), default=0)
    return n if n > QPGEN_MAX_CAT_LEVELS else None


def strategy_cat_cov_names(ctx: BuildContext, ds: str, strategy: InputSetConfig) -> list[str]:
    """Categorical covariates a strategy passes to qpgen as ``--cat-cov-names``.

    They are the dataset's categorical covariates (declared with ``cov_categorical`` /
    ``--cat-cov-names``, or any column with a non-numeric value) among the covariates the
    strategy uses. Read from the covariate file when the steps are built, so the qpgen
    command line is fixed in the Makefile. A strategy with its own covariate file uses
    its own detected categorical covariates."""
    d = ctx.cfg.datasets[ds]
    cov = strategy.cov
    if cov is None:
        return []
    if cov != "default" and cov.file is not None:
        return list(
            _dataset_categorical(
                str(cov.file), cov.format, tuple(d.missing_str), tuple(d.cov_categorical)
            )
        )
    if d.cov is None:
        return []
    cats = dataset_cat_cov_names(ctx, ds)
    if cov == "default" or not cats:
        return cats
    from anqchor.io.cov import expand_name

    wanted: set[str] = set()
    for token in cov.include or []:
        wanted.update(expand_name(token.strip(), cats))
    return [c for c in cats if c in wanted]


def dataset_cat_cov_names(ctx: BuildContext, ds: str) -> list[str]:
    """Categorical covariates of a dataset's ``cov`` (declared or non-numeric), passed to
    qpgen as ``--cat-cov-names``; empty without covariates."""
    d = ctx.cfg.datasets[ds]
    if d.cov is None:
        return []
    return list(
        _dataset_categorical(
            str(d.cov), d.cov_format, tuple(d.missing_str), tuple(d.cov_categorical)
        )
    )
