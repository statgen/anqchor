"""The report index step (links every pipeline report of the run)."""

from __future__ import annotations

from pathlib import Path

from anqchor.engine.step import Step
from anqchor.pipelines.context import BuildContext


def index_step(ctx: BuildContext, pipelines: list[str]) -> Step:
    """``report.index``: the landing page linking the reports of ``pipelines``.

    Its inputs are every pipeline report plus each dataset's ``summary.json``, from
    which the page draws one tile per pipeline and dataset.
    """
    L = ctx.layout
    summaries = {"match": L.match_summary, "cis": L.cis_summary, "trans": L.trans_summary}
    inputs: list[Path] = []
    for p in pipelines:
        for ds in ctx.cfg.datasets_for(p):  # type: ignore[arg-type]
            inputs.append(summaries[p](ds, "summary.json"))
        inputs.append(L.report(p))
    return Step(
        name="report.index",
        kind="report_index",
        inputs=tuple(inputs),
        outputs=(L.report("index"),),
        params={
            "project": ctx.cfg.project,
            "pipelines": pipelines,
            "datasets": {p: ctx.cfg.datasets_for(p) for p in pipelines},  # type: ignore[arg-type]
            "plotlyjs": ctx.cfg.report.plotlyjs,
        },
        pipeline="report",
        description="report index linking the pipeline reports",
    )
