"""The ``cis`` pipeline: cis-xQTL power of one set of inputs, measured by replication.

For every dataset the known (trait, variant) pairs of the cis reference are tested with
qpgen, then ``summarize_cis`` writes the summary tables and ``report_cis`` renders
``report/cis.html``. Runs are compared with ``anqchor cis compare``.

Association engines (``pipelines.cis.assoc_engine``):

``pair`` (default)
    One ``qpgentools pair-assoc`` run per strategy over the reference pairs
    (``--skip-rint`` unless ``pipelines.cis.rint``).
``rect``
    One ``qpgentools rect-assoc`` run per strategy over the rectangle of the reference's
    present variants × its mapped traits (``--rint-after-adj`` with ``pipelines.cis.rint``).
    The known cells are extracted. The other cells are summarised as local controls: they
    are in LD with neighbouring eQTLs, so they are not a null.

Phenotype transform: a cis run tests ONE transform (DESIGN 11.2-N). By default the
phenotypes as given; with ``pipelines.cis.rint`` (flag ``--rint``) a rank-based inverse
normal transform after covariate adjustment. Raw and RINT are compared as two runs with
``anqchor cis compare``. Semantics (verified in ``extern/qpgen`` ``cmd_pair_assoc.cpp``,
``cmd_rect_assoc.cpp`` and ``assoc_utils.cpp``):

* as given: the phenotype is residualised on the covariates (with an intercept) by
  linear regression, then regressed on the genotype.
* ``rint``: the *covariate-adjusted* phenotype is rank-based inverse-normal transformed
  (Blom offsets, average ranks for ties), then regressed on the genotype. This is the
  same transformation for both engines: pair-assoc's ``*_RINT`` columns equal
  ``rect-assoc --rint-after-adj`` in what is transformed and when.
* Both engines residualise the *genotypes* on the same covariates (Frisch-Waugh-Lovell),
  so ``BETA`` is the partial effect and the two engines give identical statistics for
  the same (variant, trait) pair. With covariates both engines use ``N - 2 - k`` residual
  degrees of freedom (k = rank of the covariates), ``N - 2`` without. anqchor requires
  qpgen 0.2.0 or later, which has these fixes (older development builds residualised
  only the phenotype in ``pair-assoc`` or used ``N - 2`` with covariates).
"""

from __future__ import annotations

from pathlib import Path

from anqchor.config.model import CisStrategyConfig
from anqchor.engine.step import Step
from anqchor.pipelines import common
from anqchor.pipelines.context import BuildContext
from anqchor.qpgen import commands as qc
from anqchor.steps.cis import SCHEMA_VERSION as SUMMARY_SCHEMA_VERSION

REPORT_STEP = "report.cis"

SUMMARY_FILES = (
    "pairs.long.tsv.gz",
    "traits.tsv",
    "strategies.tsv",
    "power_curve.tsv",
    "summary.json",
)
"""Files written by ``cis.<ds>.summarize`` (under :meth:`Layout.cis_summary`)."""

TRANSFORM_SEMANTICS = {
    "as_given": "phenotype residualised on the covariates, then regressed on the genotype",
    "rint": (
        "rank-based inverse normal transform (Blom) applied after covariate adjustment, "
        "then regressed on the genotype (pair-assoc *_RINT columns = rect-assoc --rint-after-adj)"
    ),
    "engines": (
        "pair-assoc and rect-assoc both residualise phenotype and genotypes on the covariates "
        "(Frisch-Waugh-Lovell) and give identical statistics (qpgen 0.2.0 or later)"
    ),
    "degrees_of_freedom": "N - 2 - rank(covariates) with covariates, N - 2 without",
}
"""Recorded in ``summary.json`` (``transform_semantics``); ``summary.json`` ``rint`` says
which transform the run applied."""


def transform_name(rint: bool) -> str:
    """``rint`` or ``raw`` (the rect-assoc output file suffix)."""
    return "rint" if rint else "raw"


def assoc_output(ctx: BuildContext, ds: str, strategy: CisStrategyConfig) -> Path:
    """Association output file a strategy produces for ``ds``."""
    L = ctx.layout
    c = ctx.cfg.pipelines.cis
    if c.assoc_engine == "pair":
        return L.cis_pair_assoc(ds, strategy.name)
    return L.rect_assoc("cis", ds, strategy.name, transform_name(c.rint))


def run_labels(ctx: BuildContext) -> dict[str, str]:
    """strategy name → label in the summaries: the run directory name for a run with one
    strategy (as ``cis compare`` labels runs by default), the strategy name otherwise."""
    strategies = ctx.cfg.pipelines.cis.strategies
    if len(strategies) == 1:
        name = Path(ctx.cfg.out_dir).name or strategies[0].name
        return {strategies[0].name: name.replace("/", "_")}
    return {s.name: s.name for s in strategies}


def covariate_text(strategy: CisStrategyConfig) -> str:
    """The covariates of a run as shown in the report: ``all``, ``none``, the included
    names (``sex, PC1-5``) or the covariate file name."""
    cov = strategy.cov
    if cov is None:
        return "none"
    if cov == "default":
        return "all"
    if cov.include is not None:
        return ", ".join(cov.include)
    return Path(cov.file).name if cov.file is not None else ""


def assoc_steps(ctx: BuildContext, ds: str, ref: str, strategy: CisStrategyConfig) -> list[Step]:
    """qpgen association step of one strategy on one dataset."""
    L = ctx.layout
    c = ctx.cfg.pipelines.cis
    pheno, cov = common.strategy_files(ctx, "cis", ds, strategy)
    cats = common.strategy_cat_cov_names(ctx, ds, strategy) if cov is not None else []
    base_inputs: list[Path] = [L.pgen_list, pheno]
    if cov is not None:
        base_inputs.append(cov)
    s = strategy.name
    out = assoc_output(ctx, ds, strategy)
    what = "with RINT after covariate adjustment" if c.rint else "phenotypes as given"
    if c.assoc_engine == "pair":
        argv = qc.pair_assoc(
            ctx.qpgen,
            genotypes=ctx.genotypes(),
            pheno=pheno,
            cov=cov,
            cat_cov_names=cats,
            max_cat_levels=common.max_cat_levels(ctx, ds, cats, strategy),
            pairs=L.dsref_pairs(ds, ref),
            out=out,
            skip_rint=not c.rint,
        )
        return [
            Step(
                name=f"cis.{ds}.{s}.assoc",
                commands=(qc.to_shell(argv),),
                inputs=(*base_inputs, L.dsref_pairs(ds, ref)),
                outputs=(out,),
                pipeline="cis",
                description=f"test the {ref} pairs on {ds}, {what} (qpgentools pair-assoc)",
            )
        ]
    t = transform_name(c.rint)
    argv = qc.rect_assoc(
        ctx.qpgen,
        genotypes=ctx.genotypes(),
        pheno=pheno,
        cov=cov,
        cat_cov_names=cats,
        max_cat_levels=common.max_cat_levels(ctx, ds, cats, strategy),
        var_list=L.dsref_variants(ds, ref),
        pheno_list=L.dsref_traits(ds, ref),
        rint_after_adj=c.rint,
        out=out,
    )
    return [
        Step(
            name=f"cis.{ds}.{s}.rect.{t}",
            commands=(qc.to_shell(argv),),
            inputs=(*base_inputs, L.dsref_variants(ds, ref), L.dsref_traits(ds, ref)),
            outputs=(out,),
            pipeline="cis",
            description=f"test the {ref} variant x trait rectangle on {ds}, {what} (qpgentools rect-assoc)",
        )
    ]


def summarize_step(ctx: BuildContext, ds: str, ref: str) -> Step:
    """``cis.<ds>.summarize``: replication of the reference pairs."""
    L = ctx.layout
    c = ctx.cfg.pipelines.cis
    inputs: list[Path] = [
        L.dsref_pairs(ds, ref),
        L.dsref_json(ds, ref),
        L.ref_variants(ref),
        L.ds_json(ds),
    ]
    if c.assoc_engine == "rect":
        inputs.append(L.ds_traits(ds))
    labels = run_labels(ctx)
    strategies = []
    for s in c.strategies:
        inputs.append(L.strategy_json("cis", ds, s.name))
        f = assoc_output(ctx, ds, s)
        if f not in inputs:
            inputs.append(f)
        strategies.append(
            {
                "name": s.name,
                "label": labels[s.name],
                "file": str(f),
                "strategy_json": str(L.strategy_json("cis", ds, s.name)),
            }
        )
    return Step(
        name=f"cis.{ds}.summarize",
        kind="summarize_cis",
        inputs=tuple(inputs),
        outputs=tuple(L.cis_summary(ds, n) for n in SUMMARY_FILES),
        params={
            "dataset": ds,
            "reference": ref,
            "engine": c.assoc_engine,
            "rint": c.rint,
            "strategies": strategies,
            "levels": c.levels.model_dump(mode="json"),
            "max_unresolved_fraction": c.max_unresolved_fraction,
            "power_step": 0.5,
            "transform_semantics": TRANSFORM_SEMANTICS,
            # a new summary schema re-runs the step (and refreshes stale summaries)
            "schema_version": SUMMARY_SCHEMA_VERSION,
        },
        pipeline="cis",
        description=f"replication of {ref} pairs on {ds}",
    )


def report_step(ctx: BuildContext, datasets: list[str], ref: str) -> Step:
    """``report.cis``: the HTML report over every dataset's summaries."""
    L = ctx.layout
    inputs = [L.cis_summary(ds, n) for ds in datasets for n in SUMMARY_FILES]
    pairs = ctx.cfg.references[ref].pairs
    labels = run_labels(ctx)
    return Step(
        name=REPORT_STEP,
        kind="report_cis",
        inputs=tuple(inputs),
        outputs=(L.report("cis"),),
        params={
            "project": ctx.cfg.project,
            "datasets": datasets,
            "reference": ref,
            "reference_file": Path(pairs).name if pairs is not None else None,
            "covariates": {
                labels[s.name]: covariate_text(s) for s in ctx.cfg.pipelines.cis.strategies
            },
            "engine": ctx.cfg.pipelines.cis.assoc_engine,
            "plotlyjs": ctx.cfg.report.plotlyjs,
        },
        pipeline="cis",
        description="cis report (report/cis.html)",
    )


def build(ctx: BuildContext) -> list[Step]:
    """All steps of the cis pipeline, including the shared prepare steps it needs."""
    datasets = ctx.cfg.datasets_for("cis")
    ref = ctx.cfg.reference_for("cis")
    steps = common.shared_prepare_steps(ctx, datasets, [ref])
    for ds in datasets:
        for s in ctx.cfg.pipelines.cis.strategies:
            steps.append(common.strategy_inputs(ctx, "cis", ds, s))
            steps += assoc_steps(ctx, ds, ref, s)
        steps.append(summarize_step(ctx, ds, ref))
    steps.append(report_step(ctx, datasets, ref))
    return steps
