"""The ``trans`` pipeline: trans-xQTL replication, hotspots and calibration (DESIGN 8, 11.2-Q).

A trans run tests one set of inputs and one phenotype transform over a variant × trait
rectangle with ``qpgentools rect-assoc`` and summarises it:

* ``trans.<ds>.lists`` writes the rectangle's variant and trait lists; the traits are
  restricted to the chromosomes chosen by ``trait_chroms`` (default ``genotyped``: the
  chromosomes of the genotype set, so a chr21/22 genotype set tests chr21/22 genes);
* ``trans.<ds>.<s>.inputs`` writes the run's covariates (shared step kind);
* ``trans.<ds>.<s>.rect.<t>`` runs ``rect-assoc``; ``<t>`` is ``raw`` (phenotypes as given)
  or ``rint`` with ``pipelines.trans.rint`` (``--rint-after-adj``);
* ``trans.<ds>.summarize`` classifies every cell (cis / trans) and writes the replication of
  the reference pairs per class, the hotspot bins and grid, and the calibration of the
  cross-chromosome cells;
* ``report.trans`` renders ``report/trans.html`` (``trans all`` / ``trans report`` only).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from anqchor.config.model import StrategyConfig
from anqchor.engine.step import Step
from anqchor.pipelines import common
from anqchor.pipelines.cis import transform_name
from anqchor.pipelines.context import BuildContext
from anqchor.qpgen import commands as qc
from anqchor.steps.trans import SCHEMA_VERSION as SUMMARY_SCHEMA_VERSION

REPORT_STEP = "report.trans"

SUMMARY_FILES = (
    "summary.json",
    "reference_pairs.tsv",
    "hotspots.tsv",
    "hotspot_genes.tsv",
    "qq.tsv",
    "variants.tsv",
    "pairs.sig.tsv.gz",
)
"""Files written by ``trans.<ds>.summarize`` (see :meth:`Layout.trans_summary`)."""


def covariate_text(strategy: StrategyConfig) -> str:
    """The run's covariates as shown in the report: ``all``, ``none``, the included names
    (``sex, PC1-5``) or the covariate file name (as the cis report)."""
    cov = strategy.cov
    if cov is None:
        return "none"
    if cov == "default":
        return "all"
    if cov.include is not None:
        return ", ".join(cov.include)
    return Path(cov.file).name if cov.file is not None else ""


def run_label(ctx: BuildContext) -> str:
    """The run's label in the summaries: the run directory name (as ``cis``)."""
    name = Path(ctx.cfg.out_dir).name or ctx.cfg.pipelines.trans.strategies[0].name
    return name.replace("/", "_")


def trait_chroms_param(value: object) -> str | list[str]:
    """``trait_chroms`` as a JSON-able step parameter."""
    return list(value) if isinstance(value, (list, tuple)) else str(value)


def uses_pheno_list(ctx: BuildContext) -> bool:
    """qpgen needs ``--pheno-list`` whenever the trait set is restricted at all."""
    t = ctx.cfg.pipelines.trans
    return t.traits != "all" or t.trait_chroms != "all"


def _has_pairs(ctx: BuildContext, ref: str) -> bool:
    return ctx.cfg.references[ref].pairs is not None


def lists_step(ctx: BuildContext, ds: str, ref: str) -> Step:
    """``trans.<ds>.lists``: the rectangle's variant list and trait list."""
    t = ctx.cfg.pipelines.trans
    L = ctx.layout
    inputs: list[Path] = [L.geno_json, L.ds_traits(ds)]
    if t.variants == "reference":
        variants: str = "reference"
        inputs.append(L.ref_variants_present(ref))
    else:
        variants = str(t.variants)
        inputs.append(Path(t.variants))
    if t.traits == "reference":
        traits: str = "reference"
        inputs.append(L.dsref_traits(ds, ref))
    elif t.traits == "all":
        traits = "all"
        inputs.append(L.ds_pheno(ds))
    else:
        traits = str(t.traits)
        inputs += [Path(t.traits), L.ds_pheno(ds)]
    return Step(
        name=f"trans.{ds}.lists",
        kind="trans_lists",
        inputs=tuple(inputs),
        outputs=(L.trans_variants(ds), L.trans_traits(ds), L.trans_lists_json(ds)),
        params={
            "dataset": ds,
            "reference": ref,
            "variants": variants,
            "traits": traits,
            "trait_chroms": trait_chroms_param(t.trait_chroms),
        },
        pipeline="trans",
        description=f"variant and trait lists of the trans rectangle for {ds}",
    )


def rect_step(ctx: BuildContext, ds: str, strategy: Any) -> Step:
    """``trans.<ds>.<strategy>.rect.<raw|rint>``: the run's ``rect-assoc``."""
    t = ctx.cfg.pipelines.trans
    transform = transform_name(t.rint)
    L = ctx.layout
    pheno, cov = common.strategy_files(ctx, "trans", ds, strategy)
    cats = common.strategy_cat_cov_names(ctx, ds, strategy) if cov is not None else []
    use_pheno_list = uses_pheno_list(ctx)
    argv = qc.rect_assoc(
        ctx.qpgen,
        genotypes=ctx.genotypes(),
        pheno=pheno,
        var_list=L.trans_variants(ds),
        out=L.rect_assoc("trans", ds, strategy.name, transform),
        cov=cov,
        cat_cov_names=cats,
        max_cat_levels=common.max_cat_levels(ctx, ds, cats, strategy),
        pheno_list=L.trans_traits(ds) if use_pheno_list else None,
        rint_after_adj=t.rint,
        min_maf=t.min_maf,
        min_mac=t.min_mac,
        max_chunk_vars=t.max_chunk_vars,
    )
    inputs: list[Path] = [
        L.pgen_list,
        L.geno_json,
        pheno,
        L.trans_variants(ds),
        L.strategy_json("trans", ds, strategy.name),
    ]
    if cov is not None:
        inputs.append(cov)
    if use_pheno_list:
        inputs.append(L.trans_traits(ds))
    return Step(
        name=f"trans.{ds}.{strategy.name}.rect.{transform}",
        commands=(qc.to_shell(argv),),
        inputs=tuple(inputs),
        outputs=(L.rect_assoc("trans", ds, strategy.name, transform),),
        pipeline="trans",
        description=f"test the variant x trait rectangle on {ds}, "
        + ("with RINT after covariate adjustment" if t.rint else "phenotypes as given")
        + " (qpgentools rect-assoc)",
    )


def summarize_step(ctx: BuildContext, ds: str, ref: str) -> Step:
    """``trans.<ds>.summarize``: stream the rectangle and write the summaries."""
    t = ctx.cfg.pipelines.trans
    L = ctx.layout
    s = t.strategies[0]
    rect = L.rect_assoc("trans", ds, s.name, transform_name(t.rint))
    inputs: list[Path] = [
        L.geno_json,
        L.ds_traits(ds),
        L.trans_variants(ds),
        L.trans_traits(ds),
        L.trans_lists_json(ds),
        L.ds_json(ds),
        L.strategy_json("trans", ds, s.name),
        rect,
    ]
    known = bool(t.known_pairs) and _has_pairs(ctx, ref)
    if known:
        inputs.append(L.dsref_pairs(ds, ref))
    lv = t.levels
    return Step(
        name=f"trans.{ds}.summarize",
        kind="summarize_trans",
        inputs=tuple(inputs),
        outputs=tuple(L.trans_summary(ds, n) for n in SUMMARY_FILES),
        params={
            "dataset": ds,
            "reference": ref,
            "run": run_label(ctx),
            "strategy": s.name,
            "rect": str(rect),
            "rint": t.rint,
            "variants": "reference" if t.variants == "reference" else str(t.variants),
            "traits": t.traits if t.traits in ("reference", "all") else str(t.traits),
            "trait_chroms": trait_chroms_param(t.trait_chroms),
            "known_pairs": known,
            "cis_window_bp": t.cis_window_bp,
            "min_maf": t.min_maf,
            "min_mac": t.min_mac,
            "levels": {"nominal": lv.nominal, "strict": lv.strict, "bonferroni": lv.bonferroni},
            "hotspot_bin_bp": t.hotspot_bin_bp,
            "hotspot_min_traits": t.hotspot_min_traits,
            "max_sig_rows": t.max_sig_rows,
            # a new summary schema re-runs the step (and refreshes stale summaries)
            "schema_version": SUMMARY_SCHEMA_VERSION,
        },
        pipeline="trans",
        description=f"replication, hotspots and calibration of the trans rectangle on {ds}",
    )


def report_step(ctx: BuildContext, datasets: list[str], ref: str) -> Step:
    """``report.trans``: the HTML report from the summaries."""
    L = ctx.layout
    inputs: list[Path] = [L.trans_summary(ds, n) for ds in datasets for n in SUMMARY_FILES]
    rc = ctx.cfg.references[ref]
    src = rc.pairs if rc.pairs is not None else rc.variants
    return Step(
        name=REPORT_STEP,
        kind="report_trans",
        inputs=tuple(inputs),
        outputs=(L.report("trans"),),
        params={
            "project": ctx.cfg.project,
            "datasets": datasets,
            "reference": ref,
            "reference_file": Path(src).name if src is not None else None,
            "covariates": covariate_text(ctx.cfg.pipelines.trans.strategies[0]),
            "plotlyjs": ctx.cfg.report.plotlyjs,
        },
        pipeline="trans",
        description="trans report (report/trans.html)",
    )


def build(ctx: BuildContext) -> list[Step]:
    """All steps of the trans pipeline, including the shared prepare steps it needs."""
    t = ctx.cfg.pipelines.trans
    ref = ctx.cfg.reference_for("trans")
    datasets = ctx.cfg.datasets_for("trans")
    steps = common.shared_prepare_steps(ctx, datasets, [ref])
    if not _has_pairs(ctx, ref):
        # A variants-only trans reference has no traits to map; nothing reads the mapping.
        steps = [s for s in steps if not s.name.startswith("prepare.mapping.")]
    for ds in datasets:
        steps.append(lists_step(ctx, ds, ref))
        for s in t.strategies:
            steps.append(common.strategy_inputs(ctx, "trans", ds, s))
            steps.append(rect_step(ctx, ds, s))
        steps.append(summarize_step(ctx, ds, ref))
    steps.append(report_step(ctx, datasets, ref))
    return steps
