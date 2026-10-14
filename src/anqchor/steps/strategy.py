"""``<pipeline>.<dataset>.<strategy>.inputs``: covariates and phenotypes of one strategy.

Covariates come from the dataset's covariates (``default`` or a named subset via
``include``), written with categorical covariates as text so that qpgen encodes them
(``--cat-cov-names``); from a separate file (read and validated like the dataset's, categorical levels kept as
text); or are absent (``cov: null``). A strategy may also override the phenotype matrix, which is
converted to regenie layout.
"""

from __future__ import annotations

from typing import Any

from anqchor.io.cov import (
    original_names,
    raw_as_text,
    read_and_encode,
    read_cov_raw,
    select_columns,
)
from anqchor.io.pheno import convert_to_regenie, write_regenie
from anqchor.io.tables import atomic_copy_gz, examples, read_lines
from anqchor.pipelines.layout import Layout
from anqchor.steps import StepContext, StepError, register
from anqchor.steps.prepare import input_errors
from anqchor.util.files import read_json, write_json


def _check_overlap(ds: str, what: str, ids: list[str], ds_samples: list[str]) -> int:
    n = len(set(ids) & set(ds_samples))
    if n == 0:
        raise StepError(
            f"{what}: no sample in common with dataset {ds} (IDs look like {examples(ids, 3)}; "
            f"dataset samples look like {examples(ds_samples, 3)})"
        )
    return n


@register("strategy_inputs")
@input_errors
def strategy_inputs(ctx: StepContext) -> None:
    """Write the strategy's covariate file (if any), phenotype override and summary."""
    p = ctx.params
    pipeline, ds = p["pipeline"], p["dataset"]
    strat: dict[str, Any] = p["strategy"]
    name = strat["name"]
    cov = strat.get("cov", "default")
    out_cov = Layout.strategy_cov(pipeline, ds, name)
    out_pheno = Layout.strategy_pheno(pipeline, ds, name)
    ds_samples = read_lines(ctx.path(Layout.ds_samples(ds)))
    missing = list(p.get("missing_str") or ["NA"])
    info: dict[str, Any] = {"cov_source": None, "cov_columns": [], "n_cov": 0}
    label = f"{pipeline} strategy {name} on {ds}"
    if cov == "default":
        if out_cov in ctx.outputs:
            atomic_copy_gz(ctx.path(Layout.ds_cov(ds)), ctx.path(out_cov))
            meta = read_json(ctx.path(Layout.ds_cov_json(ds)))
            info = {
                "cov_source": "default",
                "cov_columns": meta["raw_columns"],
                "n_cov": len(meta["columns"]),
                "cat_cov_names": meta["categorical_names"],
            }
        else:
            info["cov_source"] = "default (dataset has no covariates)"
    elif cov is not None and cov.get("include") is not None:
        meta = read_json(ctx.path(Layout.ds_cov_json(ds)))
        cols, empty = select_columns(cov["include"], meta["groups"], meta["columns"])
        if empty and ctx.log is not None:
            ctx.log.warning(
                "%s: covariate(s) %s were dropped as constant; nothing to include",
                label,
                ", ".join(empty),
            )
        if not cols:
            raise StepError(f"{label}: cov.include {cov['include']} selects no covariate column")
        raw = read_cov_raw(ctx.path(Layout.ds_cov(ds)))
        names = [
            n
            for n in original_names(cov["include"], meta["groups"], meta["columns"])
            if n in raw.columns
        ]
        write_regenie(ctx.path(out_cov), raw[names])
        cats = [n for n in meta.get("categorical_names", []) if n in names]
        info = {
            "cov_source": "include",
            "cov_include": cov["include"],
            "cov_columns": names,
            "n_cov": len(cols),
            "cov_empty": empty,
            "cat_cov_names": cats,
        }
    elif cov is not None and cov.get("file"):
        enc = read_and_encode(
            cov["file"],
            cov.get("format", "regenie"),
            p.get("dataset_cov_categorical") or [],
            p.get("dataset_cov_missing", "error"),
            missing,
            ignore_unknown_categorical=True,
        )
        if enc.values.shape[1] == 0:
            raise StepError(f"{label}: no usable covariate column in {cov['file']}")
        _check_overlap(
            ds, f"{label}: covariate file {cov['file']}", list(enc.values.index), ds_samples
        )
        assert enc.raw is not None
        write_regenie(ctx.path(out_cov), raw_as_text(enc.raw))
        info = {
            "cov_source": "file",
            "cov_file": cov["file"],
            "cov_columns": [str(c) for c in enc.raw.columns],
            "n_cov": len(enc.columns),
            "cat_cov_names": enc.categorical_names,
            "cov_encoding": enc.as_json(),
        }
    pheno_rel = Layout.ds_pheno(ds)
    pheno_info: dict[str, Any] = {}
    if strat.get("pheno"):
        conv = convert_to_regenie(
            strat["pheno"], strat.get("pheno_format", "regenie"), ctx.path(out_pheno), missing
        )
        n = _check_overlap(
            ds, f"{label}: phenotype file {strat['pheno']}", conv.matrix.samples, ds_samples
        )
        pheno_rel = out_pheno
        pheno_info = {
            "pheno_source": strat["pheno"],
            "n_samples_pheno": len(conv.matrix.samples),
            "n_traits": len(conv.matrix.features),
            "n_samples_in_dataset": n,
            "pheno_fast_path": conv.fast_path,
        }
    write_json(
        ctx.path(Layout.strategy_json(pipeline, ds, name)),
        {
            "schema_version": 1,
            "pipeline": pipeline,
            "dataset": ds,
            "strategy": name,
            **info,
            "pheno": str(pheno_rel),
            **pheno_info,
            # trans strategies list their transforms; a cis run has one (pipelines.cis.rint)
            **({"transforms": strat["transforms"]} if "transforms" in strat else {}),
        },
    )
    if ctx.log is not None:
        ctx.log.info("%s: %d covariate(s), phenotypes %s", label, info["n_cov"], pheno_rel)
