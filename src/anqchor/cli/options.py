"""Reusable Typer option declarations and the flag → config-override mapping."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Any

import typer

Config = Annotated[
    Path | None,
    typer.Option("--config", "-c", help="YAML configuration file.", exists=True, dir_okay=False),
]
Out = Annotated[Path | None, typer.Option("--out", "-o", help="Run directory (overrides out_dir).")]
Jobs = Annotated[
    int, typer.Option("--jobs", "-j", min=1, help="Steps to run in parallel (make -j).")
]
Threads = Annotated[
    int | None,
    typer.Option(
        "--threads",
        min=1,
        help="Threads for qpgentools match-prs-pheno, which only identity matching runs (cis and trans use single-threaded qpgen commands). Changing it re-runs nothing.",
    ),
]
Dataset = Annotated[
    list[str] | None, typer.Option("--dataset", help="Restrict to this dataset (repeatable).")
]
Strategy = Annotated[
    list[str] | None, typer.Option("--strategy", help="Restrict to this strategy (repeatable).")
]
DryRun = Annotated[
    bool,
    typer.Option("--dry-run", help="Write config and Makefile, show what would run, run nothing."),
]
Until = Annotated[str | None, typer.Option("--until", help="Run up to and including this step.")]
Only = Annotated[
    str | None,
    typer.Option("--only", help="Re-run only this step (its upstream steps must be up to date)."),
]
Force = Annotated[
    bool, typer.Option("--force", help="Re-run the selected steps and everything upstream.")
]
KeepGoing = Annotated[
    bool, typer.Option("--keep-going", help="Continue with independent steps after a failure.")
]
PlotlyJs = Annotated[
    str | None,
    typer.Option("--plotlyjs", help="'embed' (offline, default) or 'cdn' (small files)."),
]

# Inputs for flag-only invocations (no --config).
GenoList = Annotated[
    Path | None,
    typer.Option(
        "--geno-list", help="qpgen pgen list (CHROM BEG END PREFIX).", exists=True, dir_okay=False
    ),
]
PFile = Annotated[
    str | None, typer.Option("--pfile", help="Single PLINK2 prefix with .pgen/.psam/.pvar.idx.gz.")
]
Pheno = Annotated[
    Path | None, typer.Option("--pheno", help="Phenotype matrix.", exists=True, dir_okay=False)
]
PhenoFormat = Annotated[
    str | None,
    typer.Option("--pheno-format", help="regenie | tensorqtl | tsv-sample-row | tsv-sample-col"),
]
Cov = Annotated[
    Path | None, typer.Option("--cov", help="Covariate matrix.", exists=True, dir_okay=False)
]
CovFormat = Annotated[
    str | None,
    typer.Option("--cov-format", help="regenie | tensorqtl | tsv-sample-row | tsv-sample-col"),
]
CatCovNames = Annotated[
    str | None,
    typer.Option(
        "--cat-cov-names",
        help="Comma-separated categorical covariates of --cov (e.g. sex,batch); their text levels "
        "are encoded by qpgen. Non-numeric covariate columns are treated as categorical automatically.",
    ),
]
TraitAnnot = Annotated[
    Path | None,
    typer.Option(
        "--trait-annot",
        help="Trait loci (BED-like or tensorQTL matrix).",
        exists=True,
        dir_okay=False,
    ),
]
TraitIdRule = Annotated[
    str | None,
    typer.Option(
        "--trait-id-rule", help="as_is | before_colon | strip_version | before_colon_strip_version"
    ),
]
Pairs = Annotated[
    Path | None,
    typer.Option(
        "--pairs",
        help="Reference pairs: trait variant beta se \\[pval].",
        exists=True,
        dir_okay=False,
    ),
]
SampleMap = Annotated[
    Path | None,
    typer.Option(
        "--sample-map",
        help="geno_id/pheno_id map for identity matching.",
        exists=True,
        dir_okay=False,
    ),
]
DatasetName = Annotated[
    str, typer.Option("--dataset-name", help="Dataset name used with flag-only inputs.")
]
RefName = Annotated[
    str, typer.Option("--reference-name", help="Reference name used with flag-only inputs.")
]


def flag_overrides(
    *,
    geno_list: Path | None = None,
    pfile: str | None = None,
    pheno: Path | None = None,
    pheno_format: str | None = None,
    cov: Path | None = None,
    cov_format: str | None = None,
    cat_cov_names: str | None = None,
    trait_annot: Path | None = None,
    trait_id_rule: str | None = None,
    sample_map: Path | None = None,
    pairs: Path | None = None,
    variants: Path | None = None,
    ref_kind: str = "cis",
    dataset_name: str = "dataset",
    reference_name: str = "reference",
) -> dict[str, Any]:
    """Nested override dict from input flags; empty when no input flag was given."""
    o: dict[str, Any] = {}
    if geno_list is not None or pfile is not None:
        o["genotypes"] = {"pgen_list": str(geno_list) if geno_list else None, "pfile": pfile}
        o["genotypes"] = {k: v for k, v in o["genotypes"].items() if v is not None}
    ds = {
        "pheno": str(pheno) if pheno else None,
        "pheno_format": pheno_format,
        "cov": str(cov) if cov else None,
        "cov_format": cov_format,
        "cov_categorical": [c.strip() for c in cat_cov_names.split(",") if c.strip()]
        if cat_cov_names
        else None,
        "trait_annot": str(trait_annot) if trait_annot else None,
        "trait_id_rule": trait_id_rule,
        "sample_map": str(sample_map) if sample_map else None,
    }
    ds = {k: v for k, v in ds.items() if v is not None}
    if ds:
        o["datasets"] = {dataset_name: ds}
    if pairs is not None or variants is not None:
        ref: dict[str, Any] = {"kind": ref_kind}
        if pairs is not None:
            ref["pairs"] = str(pairs)
        if variants is not None:
            ref["variants"] = str(variants)
        o["references"] = {reference_name: ref}
    return o


def set_nested(d: dict[str, Any], dotted: str, value: Any) -> None:
    """``set_nested(d, "pipelines.match.rint", False)``; ``None`` values are skipped."""
    if value is None:
        return
    keys = dotted.split(".")
    cur = d
    for k in keys[:-1]:
        cur = cur.setdefault(k, {})
    cur[keys[-1]] = value
