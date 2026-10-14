"""``anqchor trans compare``: compare finished trans runs.

A trans run (``anqchor trans run --out DIR``) tests one set of inputs (genotypes,
phenotypes, covariates, phenotype transform, reference) over a variant × gene rectangle
and summarises it (DESIGN 11.2-Q). This module reads the summaries of several runs,
labels each run (its directory name by default) and writes the comparison tables that
``anqchor trans report`` renders (DESIGN 11.2-R):

* ``runs.tsv``: one row per run: label, run directory, the inputs (genotypes, phenotypes,
  covariates, ``rint``, reference), the rectangle and its cells per class, replication of
  the reference pairs per class (over each run's own tested pairs and over the pairs tested
  in every run), λ_GC and observed vs expected counts of the cross-chromosome cells, and the
  hotspot count at the run's default threshold and minimum.
* ``hotspot_grid.tsv``: hotspot counts per run, threshold and minimum genes (long format).
* ``hotspots.tsv``, ``hotspot_genes.tsv``, ``qq.tsv``, ``reference_pairs.tsv``: the runs'
  tables stacked with a ``run`` column (the run label); ``reference_pairs.tsv`` adds
  ``tested_in_all``.
* ``summary.json`` (``kind: comparison``): the runs, the union of the hotspot thresholds,
  each run's summary sections (``per_run``), the comparability checks and the warnings.
* ``compare.json`` and ``settings/<label>.anqchor.yaml`` (shared with ``cis compare``).

Comparability: λ_GC, hotspots and the observed counts are statistics of each run's whole
rectangle; they compare directly only when the runs tested the same cells. The comparison
checks the genotypes, the reference, the rectangle (tested variants and genes), the cell
counts per class and the window and bin settings, and warns about every difference. The
reference pairs are matched across runs by ``(ref_trait, variant)``; the ``*_common``
counts are over the pairs tested in every run.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from anqchor.config.load import ConfigError, load_resolved_config
from anqchor.config.model import Config, CovSpec
from anqchor.io.tables import read_lines
from anqchor.pipelines.layout import Layout
from anqchor.steps import StepError
from anqchor.steps.compare import (
    COMPARISON_DATASET,
    check_unique,
    clean_label,
    file_sha256,
    prepare_out_dir,
    run_labels,
    write_settings_and_manifest,
)
from anqchor.steps.trans import REPLICATION_LEVELS, SCHEMA_VERSION, write_tsv
from anqchor.util.files import read_json, sha256_text, write_json

TABLES = (
    "runs.tsv",
    "hotspot_grid.tsv",
    "hotspots.tsv",
    "hotspot_genes.tsv",
    "qq.tsv",
    "reference_pairs.tsv",
    "summary.json",
)
"""Tables of a trans comparison directory (all read by the report)."""

LEVELS = ("sign", *REPLICATION_LEVELS)
"""Replication levels of the reference pairs: same sign, p < 0.05, p < 1e-5."""

CLASSES = ("cis", "trans")

_STR = {
    "variant": str,
    "trait": str,
    "chrom": str,
    "trait_chrom": str,
    "top_trait": str,
    "top_variant": str,
    "variant_tested": str,
    "ref_trait": str,
    "run": str,
}


@dataclass
class TransRun:
    """One compared unit: a trans run directory and one of its datasets."""

    label: str
    run_dir: Path
    dataset: str
    cfg: Config
    summary: dict[str, Any]
    pairs: pd.DataFrame
    bins: pd.DataFrame
    genes: pd.DataFrame
    qq: pd.DataFrame
    variants: list[str]
    traits: list[str]


RUN_TABLES = ("reference_pairs.tsv", "hotspots.tsv", "hotspot_genes.tsv", "qq.tsv", "variants.tsv")
"""Summary tables of a trans run that the comparison reads (besides ``summary.json``)."""


def _read(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", dtype=_STR)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    except (OSError, ValueError) as e:
        raise StepError(f"cannot read {path}: {e}") from e


def load_run(run_dir: Path, label: str) -> list[TransRun]:
    """The trans summaries of every dataset in ``run_dir`` (one entry per dataset)."""
    run_dir = Path(run_dir)
    if not (run_dir / Layout.resolved_config).exists():
        raise StepError(f"{run_dir} is not an anqchor run directory (no {Layout.resolved_config})")
    try:
        cfg = load_resolved_config(run_dir)
    except ConfigError as e:
        raise StepError(f"{run_dir}: {e}") from e
    datasets = [
        ds
        for ds in cfg.datasets_for("trans")
        if (run_dir / Layout.trans_summary(ds, "summary.json")).exists()
    ]
    if not datasets:
        raise StepError(
            f"{run_dir} is not a trans run directory; give directories written by "
            "'anqchor trans run'"
        )
    out: list[TransRun] = []
    for ds in datasets:

        def p(name: str, ds: str = ds) -> Path:
            return run_dir / Layout.trans_summary(ds, name)

        summary = read_json(p("summary.json"))
        if summary.get("schema_version") != SCHEMA_VERSION:
            raise StepError(
                f"{run_dir}: trans summaries of schema {summary.get('schema_version')}, expected "
                f"{SCHEMA_VERSION} (written by another anqchor version); re-run "
                f"'anqchor trans run --out {run_dir}'"
            )
        for name in RUN_TABLES:
            if not p(name).exists():
                raise StepError(
                    f"{run_dir}: missing {Layout.trans_summary(ds, name)}; re-run "
                    f"'anqchor trans run --out {run_dir}'"
                )
        vt = _read(p("variants.tsv"))
        traits_file = run_dir / Layout.trans_traits(ds)
        out.append(
            TransRun(
                label=clean_label(label if len(datasets) == 1 else f"{label}:{ds}"),
                run_dir=run_dir.resolve(),
                dataset=ds,
                cfg=cfg,
                summary=summary,
                pairs=_read(p("reference_pairs.tsv")),
                bins=_read(p("hotspots.tsv")),
                genes=_read(p("hotspot_genes.tsv")),
                qq=_read(p("qq.tsv")),
                variants=sorted(vt["variant"].astype(str)) if "variant" in vt else [],
                traits=sorted(read_lines(traits_file)) if traits_file.exists() else [],
            )
        )
    return out


def _with_run(df: pd.DataFrame, label: str) -> pd.DataFrame:
    out = df.copy()
    if "run" in out.columns:
        out = out.drop(columns=["run"])
    out.insert(0, "run", label)
    return out


def describe(e: TransRun) -> dict[str, Any]:
    """The inputs of one run (the first columns of ``runs.tsv``)."""
    cfg, ds = e.cfg, e.dataset
    d = cfg.datasets[ds]
    ref = cfg.reference_for("trans")
    rc = cfg.references[ref]
    s0 = cfg.pipelines.trans.strategies[0]
    include = None
    if isinstance(s0.cov, CovSpec) and s0.cov.include is not None:
        include = ",".join(s0.cov.include)
    gj = read_json(e.run_dir / Layout.geno_json) if (e.run_dir / Layout.geno_json).exists() else {}
    s = e.summary
    return {
        "label": e.label,
        "run_dir": str(e.run_dir),
        "dataset": ds,
        "genotypes": str(cfg.genotypes.pgen_list or cfg.genotypes.pfile),
        "n_geno_samples": gj.get("n_samples"),
        "pheno": str(d.pheno),
        "n_samples": s.get("n_samples"),
        "cov": str(d.cov) if d.cov else None,
        "cov_include": include if include is not None else ("all" if d.cov else None),
        "n_cov": s.get("n_cov") or 0,
        "cov_columns": ",".join(s.get("cov_columns") or []),
        "rint": bool(s.get("rint")),
        "reference": ref,
        "reference_pairs": str(rc.pairs if rc.pairs is not None else rc.variants),
    }


def _replication(e: TransRun, common: set[tuple[str, str]]) -> dict[str, Any]:
    """Per class: pairs, tested, replicated per level; over the run's own pairs and over the
    pairs tested in every run (``_common``)."""
    out: dict[str, Any] = {}
    df = e.pairs
    for cls in CLASSES:
        sub = df[df["pair_class"] == cls] if len(df) else df
        tested = sub[sub["status"] == "tested"] if len(sub) else sub
        out[f"{cls}_n_pairs"] = len(sub)
        out[f"{cls}_n_tested"] = len(tested)
        for lv in LEVELS:
            out[f"{cls}_n_rep_{lv}"] = int(sub[f"rep_{lv}"].astype(bool).sum()) if len(sub) else 0
        if len(tested):
            keys = list(zip(tested["ref_trait"], tested["variant"], strict=True))
            on = np.array([k in common for k in keys], dtype=bool)
            t2 = tested[on]
        else:
            t2 = tested
        out[f"{cls}_n_common"] = len(t2)
        for lv in LEVELS:
            out[f"{cls}_n_rep_{lv}_common"] = (
                int(t2[f"rep_{lv}"].astype(bool).sum()) if len(t2) else 0
            )
    return out


def _alpha_key(a: float) -> str:
    """``0.001`` → ``1e-3``, ``1e-05`` → ``1e-5`` (column suffixes of ``runs.tsv``)."""
    m, x = f"{a:.0e}".split("e")
    return f"{m}e{int(x)}"


def _metrics(e: TransRun) -> dict[str, Any]:
    """Cells, inflation and default hotspot count of one run (``runs.tsv`` columns)."""
    s = e.summary
    cells = s.get("cells") or {}
    inf = s.get("inflation") or {}
    hot = s.get("hotspots") or {}
    row: dict[str, Any] = {
        "n_variants_tested": s.get("n_variants_tested"),
        "n_traits_tested": s.get("n_traits_tested"),
        "n_cells": s.get("n_cells"),
        **{f"n_cells_{k}": cells.get(k) for k in ("cis", "cross_chrom", "same_chrom_distal")},
        "n_cells_trans": s.get("n_cells_trans"),
        "cis_window_bp": (s.get("params") or {}).get("cis_window_bp"),
        "lambda_gc": inf.get("lambda_gc"),
        "n_inflation_cells": inf.get("n_cells"),
        "n_inflation_cells_excluded": inf.get("n_cells_excluded"),
    }
    for c in inf.get("counts") or []:
        k = _alpha_key(float(c["alpha"]))
        row[f"observed_p{k}"] = c.get("observed")
        row[f"expected_p{k}"] = c.get("expected")
    d = hot.get("default") or {}
    row.update(
        {
            "hotspot_bin_bp": hot.get("bin_bp"),
            "n_bins": hot.get("n_bins"),
            "bonferroni_p": hot.get("bonferroni_p"),
            "hotspot_threshold": d.get("threshold"),
            "hotspot_min_genes": d.get("k"),
            "n_hotspots": d.get("n_hotspots"),
        }
    )
    return row


def grid_rows(e: TransRun) -> list[dict[str, Any]]:
    """``hotspot_grid.tsv`` rows of one run: threshold × minimum genes."""
    hot = e.summary.get("hotspots") or {}
    ths, ks, counts = hot.get("thresholds") or [], hot.get("k") or [], hot.get("counts") or []
    rows = []
    for i, t in enumerate(ths):
        for j, k in enumerate(ks):
            rows.append(
                {
                    "run": e.label,
                    "threshold": t["key"],
                    "label": t["label"],
                    "alpha": t["alpha"],
                    "log10p": t["log10p"],
                    "k": k,
                    "n_hotspots": counts[i][j],
                }
            )
    return rows


def union_thresholds(entries: list[TransRun]) -> list[dict[str, Any]]:
    """Every threshold key of the runs' grids, loosest first (Bonferroni by its first alpha)."""
    seen: dict[str, dict[str, Any]] = {}
    for e in entries:
        for t in (e.summary.get("hotspots") or {}).get("thresholds") or []:
            seen.setdefault(t["key"], {k: t[k] for k in ("key", "label", "alpha", "log10p")})
    return sorted(seen.values(), key=lambda t: -float(t["alpha"]))


def _hash_lines(lines: list[str]) -> str:
    return sha256_text("\n".join(lines))


def comparability(entries: list[TransRun]) -> tuple[dict[str, Any], list[str]]:
    """Which inputs and settings the runs share, and a warning for each difference."""
    warnings: list[str] = []

    def per_run(fn: Any) -> dict[str, Any]:
        return {e.label: fn(e) for e in entries}

    def ref_file(e: TransRun) -> Path:
        ref = e.cfg.reference_for("trans")
        f = e.run_dir / Layout.ref_pairs(ref)
        return f if f.exists() else e.run_dir / Layout.ref_variants(ref)

    geno = per_run(lambda e: file_sha256(e.run_dir / Layout.geno_samples))
    pheno_samples = per_run(lambda e: file_sha256(e.run_dir / Layout.ds_samples(e.dataset)))
    ref = per_run(lambda e: file_sha256(ref_file(e)))
    rect = per_run(lambda e: _hash_lines([*e.variants, "--", *e.traits]))
    cells = per_run(lambda e: e.summary.get("cells") or {})
    samples = per_run(lambda e: e.summary.get("n_samples"))
    window = per_run(lambda e: (e.summary.get("params") or {}).get("cis_window_bp"))
    bin_bp = per_run(lambda e: (e.summary.get("hotspots") or {}).get("bin_bp"))
    excluded = per_run(lambda e: (e.summary.get("inflation") or {}).get("n_cells_excluded") or 0)

    def same(d: dict[str, Any]) -> bool:
        vals = [repr(v) for v in d.values()]
        return len(set(vals)) <= 1

    def listing(d: dict[str, Any]) -> str:
        return ", ".join(f"{k} {v:,}" if isinstance(v, int) else f"{k} {v}" for k, v in d.items())

    out = {
        "same_genotype_samples": same(geno),
        "same_reference": same(ref),
        "same_rectangle": same(rect),
        "same_cells": same(cells),
        "same_samples": same(samples) and same(pheno_samples),
        "same_cis_window": same(window),
        "same_bin_width": same(bin_bp),
        "reference_sha256": ref,
        "rectangle_sha256": rect,
        "n_samples": samples,
        "rint": per_run(lambda e: bool(e.summary.get("rint"))),
    }
    if len(entries) > 1:
        if not out["same_genotype_samples"]:
            warnings.append("the runs used genotype files with different samples")
        if not out["same_reference"]:
            warnings.append("the runs used different references")
        n = per_run(lambda e: int(e.summary.get("n_cells") or 0))
        if not out["same_rectangle"] or not same(n):
            warnings.append(
                "the runs tested different variant–gene tests (" + listing(n) + "): λGC, "
                "observed counts and hotspots are each over the run's own tests"
            )
        elif not out["same_cells"] and out["same_cis_window"]:

            def split(c: dict[str, Any]) -> str:
                return " / ".join(
                    f"{int(c.get(k) or 0):,}" for k in ("cis", "cross_chrom", "same_chrom_distal")
                )

            warnings.append(
                "the runs split their tests differently into cis / different chromosome / same "
                "chromosome (" + ", ".join(f"{k} {split(v)}" for k, v in cells.items()) + ")"
            )
        if not same(samples):
            warnings.append(
                "the runs tested different numbers of samples (" + listing(samples) + ")"
            )
        elif not same(pheno_samples):
            warnings.append("the runs tested different samples (same number)")
        if not out["same_cis_window"]:
            warnings.append("the runs used different cis windows (" + listing(window) + ")")
        if not out["same_bin_width"]:
            warnings.append(
                "the runs used different hotspot bin widths ("
                + listing(bin_bp)
                + "); hotspot counts are not comparable"
            )
        if len({repr(v) for v in excluded.values()}) > 1:
            warnings.append(
                "the runs left different numbers of tests out of the inflation check ("
                + listing(excluded)
                + ")"
            )
    return out, warnings


def compare_runs(
    run_dirs: list[Path],
    out_dir: Path,
    labels: list[str] | None = None,
    argv: list[str] | None = None,
) -> dict[str, Any]:
    """Compare trans runs and write the comparison tables into ``out_dir``. Returns the summary."""
    labs = run_labels(run_dirs, labels)
    prepare_out_dir(out_dir, "trans", run_dirs)
    entries: list[TransRun] = []
    for run, lab in zip(run_dirs, labs, strict=True):
        entries += load_run(run, lab)
    check_unique([e.label for e in entries])

    comp, warnings = comparability(entries)

    # reference pairs tested in every run (matched by reference trait and variant)
    common: set[tuple[str, str]] | None = None
    for e in entries:
        t = e.pairs[e.pairs["status"] == "tested"] if len(e.pairs) else e.pairs
        keys = set(zip(t["ref_trait"], t["variant"], strict=True)) if len(t) else set()
        common = keys if common is None else common & keys
    common = common or set()

    run_rows: list[dict[str, Any]] = []
    for e in entries:
        run_rows.append({**describe(e), **_replication(e, common), **_metrics(e)})
    runs = pd.DataFrame(run_rows)
    if len(entries) > 1:
        for cls in CLASSES:
            n_c = int(runs[f"{cls}_n_common"].min())
            n_m = int(runs[f"{cls}_n_tested"].max())
            if 0 < n_c < n_m:
                warnings.append(
                    f"the runs tested different {cls} reference pairs: {n_c:,} tested in every "
                    f"run, up to {n_m:,} in one run; replicated {cls} pairs are counted over the "
                    f"{n_c:,}"
                )
            elif n_c < n_m:
                warnings.append(
                    f"no {cls} reference pair was tested in every run (up to {n_m:,} in one "
                    f"run); replicated {cls} pairs are counted over each run's own pairs"
                )

    grid = pd.DataFrame([r for e in entries for r in grid_rows(e)])
    # hotspots of every run at one setting (the first run's opening threshold and minimum)
    d0 = (entries[0].summary.get("hotspots") or {}).get("default") or {}
    t0, k0 = d0.get("threshold"), d0.get("k")
    shared: list[Any] = []
    for e in entries:
        hit = (
            grid[(grid["run"] == e.label) & (grid["threshold"] == t0) & (grid["k"] == k0)]
            if len(grid)
            else grid
        )
        shared.append(int(hit["n_hotspots"].iloc[0]) if len(hit) else None)
    runs["n_hotspots_shared"] = pd.array(shared, dtype="Int64")
    opening = {(r["hotspot_threshold"], r["hotspot_min_genes"]) for r in run_rows}
    if len(opening) > 1:
        warnings.append(
            "the runs open their hotspots at different thresholds or minimum genes; "
            "n_hotspots_shared counts every run at the first run's setting"
        )
    bins = pd.concat([_with_run(e.bins, e.label) for e in entries], ignore_index=True)
    genes = pd.concat([_with_run(e.genes, e.label) for e in entries], ignore_index=True)
    qq = pd.concat([_with_run(e.qq, e.label) for e in entries], ignore_index=True)
    pairs_parts = []
    for e in entries:
        pdf = _with_run(e.pairs, e.label)
        if len(pdf):
            pdf["tested_in_all"] = [
                (str(t), str(v)) in common and s == "tested"
                for t, v, s in zip(pdf["ref_trait"], pdf["variant"], pdf["status"], strict=True)
            ]
        pairs_parts.append(pdf)
    pairs = pd.concat(pairs_parts, ignore_index=True)

    for e in entries:
        warnings += [f"{e.label}: {w}" for w in e.summary.get("warnings") or []]
    first = entries[0]
    per_run = {
        e.label: {
            k: e.summary.get(k)
            for k in (
                "run",
                "dataset",
                "rint",
                "n_cov",
                "n_samples",
                "n_variants",
                "n_traits",
                "n_variants_tested",
                "n_traits_tested",
                "n_cells",
                "n_cells_trans",
                "cells",
                "n_reference_pairs",
                "replication",
                "hotspots",
                "inflation",
                "params",
                "warnings",
            )
        }
        for e in entries
    }
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "kind": "comparison",
        "pipeline": "trans",
        "dataset": COMPARISON_DATASET,
        "labels": [e.label for e in entries],
        "reference": ", ".join(dict.fromkeys(e.cfg.reference_for("trans") for e in entries)),
        "runs": [describe(e) for e in entries],
        "replication_levels": {"sign": None, **REPLICATION_LEVELS},
        "n_common_pairs": {
            cls: int(runs[f"{cls}_n_common"].min()) if len(runs) else 0 for cls in CLASSES
        },
        "thresholds": union_thresholds(entries),
        "comparability": comp,
        "per_run": per_run,
        "definitions": first.summary.get("definitions"),
        "warnings": warnings,
    }

    out_dir = Path(out_dir)
    prepare_out_dir(out_dir, "trans", run_dirs)
    out_dir.mkdir(parents=True, exist_ok=True)
    T = Layout.compare_table
    write_tsv(out_dir / T("runs.tsv"), runs)
    write_tsv(out_dir / T("hotspot_grid.tsv"), grid)
    write_tsv(out_dir / T("hotspots.tsv"), bins)
    write_tsv(out_dir / T("hotspot_genes.tsv"), genes)
    write_tsv(out_dir / T("qq.tsv"), qq)
    write_tsv(out_dir / T("reference_pairs.tsv"), pairs)
    write_json(out_dir / T("summary.json"), summary)
    write_settings_and_manifest(out_dir, "trans", run_rows, argv)
    return summary
