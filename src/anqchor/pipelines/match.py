"""The ``match`` pipeline: sample identity checks between genotypes and omics data.

Steps (DESIGN.md section 6.1), besides the shared prepare and PRS steps:

* ``match.<ds>.qpgen`` (kind ``match_qpgen``): ``qpgentools match-prs-pheno``.
* ``match.<ds>.summarize`` (kind ``summarize_match``): per-sample calls and summaries.
* ``report.match`` (kind ``report_match``): ``report/match.html``.

The step implementations are in :mod:`anqchor.steps.match` and
:mod:`anqchor.report.match`.
"""

from __future__ import annotations

from pathlib import Path

from anqchor.engine.step import Step
from anqchor.pipelines import common
from anqchor.pipelines.context import BuildContext
from anqchor.steps.match import SUMMARY_FILES

REPORT_STEP = "report.match"


def qpgen_step(ctx: BuildContext, ds: str, ref: str) -> Step:
    """``match.<ds>.qpgen``: run ``match-prs-pheno`` for one dataset."""
    L = ctx.layout
    m = ctx.cfg.pipelines.match
    d = ctx.cfg.datasets[ds]
    use_cov = bool(m.use_cov and d.cov is not None)
    inputs: list[Path] = [L.prs_matrix(ref), L.ds_pheno(ds)]
    if use_cov:
        inputs.append(L.ds_cov(ds))
    inputs += [
        L.ds_samples(ds),
        L.ds_sample_map(ds),
        L.dsref_trait_map(ds, ref),
        L.ds_json(ds),
        L.dsref_json(ds, ref),
    ]
    cats = common.dataset_cat_cov_names(ctx, ds) if use_cov else []
    # Only when needed, so that the step signature of ordinary covariates is unchanged.
    extra = {"max_cat_levels": n} if (n := common.max_cat_levels(ctx, ds, cats)) else {}
    return Step(
        name=f"match.{ds}.qpgen",
        kind="match_qpgen",
        inputs=tuple(inputs),
        outputs=(*L.match_qpgen_outputs(ds), L.match_trait_map(ds), L.match_params(ds)),
        params={
            "qpgen": str(ctx.qpgen),
            "ref": ref,
            "dataset": ds,
            "has_sample_map": d.sample_map is not None,
            "cat_cov_names": cats,
            "use_cov": use_cov,
            "rint": m.rint,
            "mahalanobis": m.mahalanobis,
            "lambda": m.lambda_,
            "min_self_for_auto": m.min_self_for_auto,
            "z_threshold": m.z_threshold,
            "z_diff": m.z_diff,
            "min_weight": m.min_weight,
            "missing_str": list(d.missing_str),
            **extra,
        },
        pipeline="match",
        description=f"match omics samples of {ds} to genotypes (qpgentools match-prs-pheno)",
    )


def summarize_step(ctx: BuildContext, ds: str) -> Step:
    """``match.<ds>.summarize``: samples, swaps, rescue candidates, duplicates, relatives,
    weights."""
    L = ctx.layout
    return Step(
        name=f"match.{ds}.summarize",
        kind="summarize_match",
        inputs=(
            *L.match_qpgen_outputs(ds),
            L.match_params(ds),
            L.match_trait_map(ds),
            L.ds_json(ds),
            L.ds_samples(ds),
            L.ds_sample_map(ds),
        ),
        outputs=tuple(L.match_summary(ds, n) for n in SUMMARY_FILES),
        params={"dataset": ds},
        pipeline="match",
        description=f"summarise identity calls for {ds}",
    )


def report_step(ctx: BuildContext, datasets: list[str]) -> Step:
    """``report.match``: the HTML report over every dataset."""
    L = ctx.layout
    m = ctx.cfg.pipelines.match
    return Step(
        name=REPORT_STEP,
        kind="report_match",
        inputs=tuple(L.match_summary(ds, n) for ds in datasets for n in SUMMARY_FILES),
        outputs=(L.report("match"),),
        params={
            "datasets": list(datasets),
            "omics": {ds: ctx.cfg.datasets[ds].omics for ds in datasets},
            "project": ctx.cfg.project,
            "z_threshold": m.z_threshold,
            "z_diff": m.z_diff,
            "plotlyjs": ctx.cfg.report.plotlyjs,
        },
        pipeline="match",
        description="identity-matching report (report/match.html)",
    )


def build(ctx: BuildContext) -> list[Step]:
    """All steps of the ``match`` pipeline, shared prepare and PRS steps included."""
    datasets = ctx.cfg.datasets_for("match")
    ref = ctx.cfg.reference_for("match")
    steps = common.shared_prepare_steps(ctx, datasets, [ref])
    steps.append(common.prs(ctx, ref))
    for ds in datasets:
        steps.append(qpgen_step(ctx, ds, ref))
        steps.append(summarize_step(ctx, ds))
    steps.append(report_step(ctx, datasets))
    return steps
