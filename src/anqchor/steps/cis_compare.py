"""``anqchor cis compare``: compare finished cis runs.

A cis run (``anqchor cis run --out DIR``) evaluates replication of a reference's known
cis pairs for one set of inputs: genotypes, phenotypes, covariates (optionally a subset
of their columns) and reference. Runs to be compared may differ in any of these. This
module reads the summaries of several runs, labels each run (directory name by default)
and writes comparison tables in the same layout as a run's summaries, so the cis report
renders them with no special casing:

* ``runs.tsv``: one row per compared run: label, run directory, dataset, the inputs used
  (genotypes, phenotypes, covariates and the covariate columns, ``rint``, reference) and
  counts.
* ``strategies.tsv``: one row per run (``run`` = ``label`` = the run label), with
  replication counts over the run's own tested pairs and over the pairs tested in every
  compared run (``n_common_tested``, ``n_rep_<level>_common``), and the effective sample
  size relative to the first run (``mean_z2_common``, ``rel_neff``). ``z`` and ``mean_z2``
  are re-derived from each run's stored ``beta`` and ``se`` (``z = beta / se``).
* ``pairs.long.tsv.gz``, ``traits.tsv``, ``power_curve.tsv``: as in a run.
* ``summary.json``: as in a run, plus ``runs`` and ``comparability`` notes.
* ``compare.json``: the manifest (runs and labels) that ``anqchor cis report`` uses to
  recognise a comparison directory.

Pairs are matched across runs by ``(ref_trait, variant)``, so runs whose phenotype files
name traits differently still line up through the reference. Runs may differ in the
phenotype transform (``rint``); each run's rows carry its own statistics.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from anqchor.config.load import ConfigError, load_resolved_config
from anqchor.config.model import Config, CovSpec
from anqchor.pipelines.layout import Layout
from anqchor.steps import StepError
from anqchor.steps.cis import (
    SCHEMA_VERSION,
    active_levels,
    add_levels,
    common_keys,
    level_cutoffs,
    level_description,
    neff_columns,
    power_curve,
    recommend,
    stringency_order,
    trait_table,
    wald_z,
    write_tsv,
)
from anqchor.steps.compare import (
    COMPARISON_DATASET,
    MANIFEST,
    check_unique,
    clean_label,
    default_labels,
    file_sha256,
    is_comparison_dir,
    prepare_out_dir,
    run_labels,
    write_settings_and_manifest,
)
from anqchor.util.files import read_json, write_json

SUMMARY_FILES = (
    "pairs.long.tsv.gz",
    "traits.tsv",
    "strategies.tsv",
    "power_curve.tsv",
    "summary.json",
)
"""Tables of a cis comparison directory read by the report (with ``runs.tsv``)."""

__all__ = [
    "COMPARISON_DATASET",
    "MANIFEST",
    "SUMMARY_FILES",
    "compare_runs",
    "default_labels",
    "describe",
    "is_comparison_dir",
    "load_run",
]


@dataclass
class RunEntry:
    """One compared unit: a cis run directory and one of its datasets."""

    label: str
    run_dir: Path
    dataset: str
    cfg: Config
    summary: dict[str, Any]
    strategies: pd.DataFrame
    long: pd.DataFrame


_clean_label = clean_label


def load_run(run_dir: Path, label: str) -> list[RunEntry]:
    """The cis summaries of every dataset in ``run_dir`` (one entry per dataset)."""
    run_dir = Path(run_dir)
    if not (run_dir / Layout.resolved_config).exists():
        raise StepError(f"{run_dir} is not an anqchor run directory (no {Layout.resolved_config})")
    try:
        cfg = load_resolved_config(run_dir)
    except ConfigError as e:
        raise StepError(f"{run_dir}: {e}") from e
    entries: list[RunEntry] = []
    datasets = [
        ds
        for ds in cfg.datasets_for("cis")
        if (run_dir / Layout.cis_summary(ds, "summary.json")).exists()
    ]
    if not datasets:
        raise StepError(
            f"{run_dir} has no cis results; run 'anqchor cis run --out {run_dir}' first"
        )
    for ds in datasets:

        def p(name: str, ds: str = ds) -> Path:
            return run_dir / Layout.cis_summary(ds, name)

        version = read_json(p("summary.json")).get("schema_version")
        if version != SCHEMA_VERSION:
            raise StepError(
                f"{run_dir}: cis summaries of schema {version}, expected {SCHEMA_VERSION} "
                f"(written by another anqchor version); re-run 'anqchor cis run --out {run_dir}'"
            )

        ent_label = label if len(datasets) == 1 else f"{label}:{ds}"
        entries.append(
            RunEntry(
                label=_clean_label(ent_label),
                run_dir=run_dir.resolve(),
                dataset=ds,
                cfg=cfg,
                summary=read_json(p("summary.json")),
                strategies=pd.read_csv(p("strategies.tsv"), sep="\t"),
                long=pd.read_csv(
                    p("pairs.long.tsv.gz"),
                    sep="\t",
                    dtype={"trait": str, "ref_trait": str, "variant": str, "run": str},
                ),
            )
        )
    return entries


def _relabel(entry: RunEntry) -> dict[str, str]:
    """Run-internal label → comparison label (the run label when the run has a single
    strategy, ``<label>:<strategy label>`` otherwise)."""
    names = list(dict.fromkeys(entry.strategies["run"].astype(str)))
    if len(names) == 1:
        return {names[0]: entry.label}
    return {s: _clean_label(f"{entry.label}:{s}") for s in names}


def describe(entry: RunEntry) -> dict[str, Any]:
    """One ``runs.tsv`` row: what went into the run."""
    cfg, ds = entry.cfg, entry.dataset
    d = cfg.datasets[ds]
    ref = cfg.reference_for("cis")
    strategies = cfg.pipelines.cis.strategies
    cov_cols: list[str] = []
    for s in strategies:
        sj = entry.run_dir / Layout.strategy_json("cis", ds, s.name)
        if sj.exists():
            cov_cols = list(read_json(sj).get("cov_columns") or [])
    geno_json = entry.run_dir / Layout.geno_json
    ds_json = entry.run_dir / Layout.ds_json(ds)
    gj = read_json(geno_json) if geno_json.exists() else {}
    dj = read_json(ds_json) if ds_json.exists() else {}
    s0 = strategies[0] if len(strategies) == 1 else None
    include = None
    if s0 is not None and isinstance(s0.cov, CovSpec) and s0.cov.include is not None:
        include = ",".join(s0.cov.include)
    return {
        "label": entry.label,
        "run_dir": str(entry.run_dir),
        "dataset": ds,
        "genotypes": str(cfg.genotypes.pgen_list or cfg.genotypes.pfile),
        "n_geno_samples": gj.get("n_samples"),
        "pheno": str(d.pheno),
        "n_pheno_samples": dj.get("n_samples_pheno"),
        "n_samples_analysable": dj.get("n_samples_analysable", dj.get("n_samples_mapped")),
        "cov": str(d.cov) if d.cov else None,
        "cov_include": include if include is not None else ("all" if d.cov else None),
        "n_cov": len(cov_cols),
        "cov_columns": ",".join(cov_cols),
        "reference": ref,
        "reference_pairs": str(cfg.references[ref].pairs),
        "rint": bool(entry.summary.get("rint", False)),
        "engine": cfg.pipelines.cis.assoc_engine,
        "n_pairs_mapped": entry.summary.get("n_pairs_mapped"),
    }


def compare_runs(
    run_dirs: list[Path],
    out_dir: Path,
    labels: list[str] | None = None,
    argv: list[str] | None = None,
) -> dict[str, Any]:
    """Compare cis runs and write the comparison tables into ``out_dir``. Returns the summary."""
    labels = run_labels(run_dirs, labels)
    prepare_out_dir(out_dir, "cis", run_dirs)
    entries: list[RunEntry] = []
    for run, lab in zip(run_dirs, labels, strict=True):
        entries += load_run(run, lab)
    check_unique([e.label for e in entries])

    level_cfg = entries[0].cfg.pipelines.cis.levels.model_dump(mode="json")
    levels = active_levels(level_cfg)
    warnings: list[str] = []

    # --------------------------------------------------------------- comparability
    ref_hashes: dict[str, str] = {}
    for e in entries:
        ref = e.cfg.reference_for("cis")
        ref_hashes[e.label] = file_sha256(e.run_dir / Layout.ref_pairs(ref))
    same_reference = len(set(ref_hashes.values())) == 1
    if not same_reference:
        warnings.append(
            "the runs used different reference pairs; compare the common-pair columns "
            "(n_rep_<level>_common), which count only pairs tested in every run"
        )
    lv_sets = {
        e.label: tuple(active_levels(e.cfg.pipelines.cis.levels.model_dump(mode="json")))
        for e in entries
    }
    if len(set(lv_sets.values())) > 1:
        warnings.append(
            "the runs used different replication levels; the levels of the first run are used"
        )

    # --------------------------------------------------------------- long + rows
    longs, rows, run_rows = [], [], []
    for e in entries:
        mapping = _relabel(e)
        lg = e.long.copy()
        lg["run"] = lg["run"].astype(str).map(mapping)
        # z from the stored statistics, never trusted from the run (one definition for all)
        tested = (lg["status"] == "tested").to_numpy()
        z = wald_z(
            lg["beta"].to_numpy(dtype=float),
            lg["se"].to_numpy(dtype=float),
            lg["log10p"].to_numpy(dtype=float),
        )
        lg["z"] = np.where(tested, z, np.nan)
        # Re-derive the levels with the comparison's level configuration (per run).
        parts = [
            add_levels(
                g.drop(
                    columns=[c for c in g.columns if c.startswith("rep_") or c == "sign_concordant"]
                ),
                level_cfg,
            )
            for _, g in lg.groupby("run", sort=False)
        ]
        lg = pd.concat(parts, ignore_index=True)
        longs.append(lg)
        st = e.strategies.copy()
        st["run"] = st["run"].astype(str).map(mapping)
        st["label"] = st["run"]
        rows.append(st)
        run_rows.append(describe(e))
    long = pd.concat(longs, ignore_index=True)
    strategies = pd.concat(rows, ignore_index=True)
    comp_labels = list(dict.fromkeys(strategies["run"]))
    # The comparison's levels may differ from a run's own: recount from the long table.
    for i, r in strategies.iterrows():
        sub = long[long["run"] == r["run"]]
        n_t = int((sub["status"] == "tested").sum())
        z = sub.loc[sub["status"] == "tested", "z"].to_numpy(dtype=float)
        z = z[np.isfinite(z)]
        strategies.loc[i, "mean_z2"] = float(np.mean(z * z)) if z.size else np.nan
        for lv in levels:
            strategies.loc[i, f"n_rep_{lv}"] = int(sub[f"rep_{lv}"].sum())
            strategies.loc[i, f"frac_rep_{lv}"] = (
                (int(sub[f"rep_{lv}"].sum()) / n_t) if n_t else np.nan
            )

    # common pairs: tested in every compared run, levels recomputed on them
    keys = common_keys(long, comp_labels)
    for i, r in strategies.iterrows():
        sub = long[long["run"] == r["run"]]
        mask = [k in keys for k in zip(sub["ref_trait"], sub["variant"], strict=True)]
        common = add_levels(
            sub[mask].drop(
                columns=[c for c in sub.columns if c.startswith("rep_") or c == "sign_concordant"]
            ),
            level_cfg,
        )
        strategies.loc[i, "n_common_tested"] = len(keys)
        for lv in levels:
            strategies.loc[i, f"n_rep_{lv}_common"] = (
                int(common[f"rep_{lv}"].sum()) if len(common) else 0
            )
    common_cols = ["n_common_tested", *(f"n_rep_{lv}_common" for lv in levels)]
    strategies[common_cols] = strategies[common_cols].astype(
        int
    )  # counts, not floats, in summary.json
    neff = neff_columns(long, comp_labels)
    strategies["mean_z2_common"] = [neff[lab]["mean_z2_common"] for lab in strategies["run"]]
    strategies["rel_neff"] = [neff[lab]["rel_neff"] for lab in strategies["run"]]
    if len(entries) > 1:
        n_common = strategies["n_common_tested"].min()
        n_max = strategies["n_tested"].max()
        if n_common < n_max:
            warnings.append(
                f"the runs tested different pairs: {int(n_common)} pairs were tested in every run "
                f"(up to {int(n_max)} in one run); the *_common columns compare those only"
            )

    # --------------------------------------------------------------- tables
    rows_dicts = strategies.to_dict(orient="records")
    longest = max(
        comp_labels,
        key=lambda lab: int(((long["run"] == lab) & (long["status"] == "tested")).sum()),
    )
    order = stringency_order(level_cutoffs(long[long["run"] == longest], levels, level_cfg))
    rec = recommend(rows_dicts, levels, level_cfg, order)
    traits = trait_table(long, comp_labels, levels, level_cfg)
    curve = power_curve(long, comp_labels)
    for e in entries:
        warnings += [f"{e.label}: {w}" for w in e.summary.get("warnings") or []]

    first = entries[0].summary
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "kind": "comparison",
        "dataset": COMPARISON_DATASET,
        "reference": ", ".join(dict.fromkeys(e.cfg.reference_for("cis") for e in entries)),
        "engine": ", ".join(dict.fromkeys(e.cfg.pipelines.cis.assoc_engine for e in entries)),
        "levels": level_description(level_cfg),
        "level_order": order,
        "transform_semantics": first.get("transform_semantics", {}),
        "baseline": comp_labels[0],
        "runs": run_rows,
        "comparability": {
            "same_reference": same_reference,
            "reference_sha256": ref_hashes,
            "rint": {r["label"]: r["rint"] for r in run_rows},
        },
        "n_ref_pairs": first.get("n_ref_pairs"),
        "n_pairs_mapped": first.get("n_pairs_mapped"),
        "n_pairs_unmapped": first.get("n_pairs_unmapped"),
        "n_traits_mapped": first.get("n_traits_mapped"),
        "n_variant_absent": first.get("n_variant_absent"),
        "n_absent_variants": first.get("n_absent_variants"),
        "n_absent_variants_indel": first.get("n_absent_variants_indel"),
        "n_absent_variants_snv": first.get("n_absent_variants_snv"),
        "absent_variant_examples": first.get("absent_variant_examples", []),
        "unresolved_fraction": max(
            float(e.summary.get("unresolved_fraction") or 0.0) for e in entries
        ),
        "max_unresolved_fraction": first.get("max_unresolved_fraction"),
        "per_run": {r["label"]: {k: v for k, v in r.items() if k != "label"} for r in rows_dicts},
        "ranking": {
            lv: [r["label"] for r in sorted(rows_dicts, key=lambda r: -r[f"n_rep_{lv}"])]
            for lv in levels
        },
        "recommended": rec,
        "warnings": warnings,
        "sentinels": None,
    }

    out_dir = Path(out_dir)
    prepare_out_dir(out_dir, "cis", run_dirs)
    out_dir.mkdir(parents=True, exist_ok=True)
    T = Layout.compare_table
    write_tsv(out_dir / T("runs.tsv"), pd.DataFrame(run_rows))
    write_tsv(out_dir / T("pairs.long.tsv.gz"), long)
    write_tsv(out_dir / T("traits.tsv"), traits)
    write_tsv(out_dir / T("strategies.tsv"), strategies)
    write_tsv(out_dir / T("power_curve.tsv"), curve)
    write_json(out_dir / T("summary.json"), summary)
    # Settings needed to restore every compared run exactly, and the compare command.
    write_settings_and_manifest(out_dir, "cis", run_rows, argv)
    return summary
