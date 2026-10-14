"""Pipelines turn a configuration into a graph of steps.

Each pipeline module exposes ``build(ctx: BuildContext) -> list[Step]`` returning its
own steps plus the shared steps it needs (see :mod:`anqchor.pipelines.common`), and
``REPORT_STEP`` naming the step that renders its HTML report.
"""

from __future__ import annotations

import dataclasses
import importlib
from collections.abc import Sequence

from anqchor.engine.step import Step, StepGraph
from anqchor.pipelines.context import BuildContext
from anqchor.pipelines.layout import Layout

PIPELINES = ("match", "cis", "trans")


def pipeline_module(name: str):  # type: ignore[no-untyped-def]
    if name not in PIPELINES:
        raise ValueError(f"unknown pipeline {name!r}")
    return importlib.import_module(f"anqchor.pipelines.{name}")


QPGEN_ID_PARAM = "qpgen_sha256"
"""Step parameter holding the qpgen binary's content hash (see :func:`with_tool_identity`)."""


def runs_qpgen(step: Step, ctx: BuildContext) -> bool:
    """Whether ``step`` runs ``ctx.qpgen``: a shell command naming the binary, or a
    Python step that receives it as its ``qpgen`` parameter."""
    binary = str(ctx.qpgen)
    return step.params.get("qpgen") == binary or any(binary in c for c in step.commands)


def with_tool_identity(step: Step, ctx: BuildContext) -> Step:
    """``step`` with the qpgen binary's sha256 in its parameters when it runs qpgen.

    The signature (``steps/<name>.json``) otherwise names the binary only by its path, so a
    rebuilt ``qpgentools`` at the same path would leave old results marked up to date."""
    if ctx.qpgen_sha256 is None or not runs_qpgen(step, ctx):
        return step
    return dataclasses.replace(step, params={**step.params, QPGEN_ID_PARAM: ctx.qpgen_sha256})


def is_report_step(step: Step) -> bool:
    """Steps that only render HTML (``report.<pipeline>``, ``report.index``)."""
    return step.name.startswith("report.")


def build_graph(
    ctx: BuildContext, pipelines: Sequence[str], *, reports: bool = True, index: bool = True
) -> tuple[StepGraph, dict[str, list[str]]]:
    """Steps of the requested pipelines, optionally with their HTML reports and the index.

    Args:
        reports: include the ``report.<pipeline>`` steps. Without them the graph ends at
            the summary statistics (``anqchor <pipeline> run``).
        index: add ``report.index`` (only meaningful with ``reports``).

    Returns the graph and the make groups: one per pipeline (its steps and everything
    upstream), ``prepare`` and ``all``.
    """
    graph = StepGraph()
    groups: dict[str, list[str]] = {}
    with_index = bool(pipelines) and reports and index
    for name in pipelines:
        steps: list[Step] = [with_tool_identity(s, ctx) for s in pipeline_module(name).build(ctx)]
        if not reports:
            steps = [s for s in steps if not is_report_step(s)]
        elif with_index:
            # the reports link to the index page (``has_index``; see report.render.index_nav)
            steps = [
                dataclasses.replace(s, params={**s.params, "has_index": True})
                if is_report_step(s)
                else s
                for s in steps
            ]
        graph.extend(steps)
        groups[name] = [s.name for s in steps]
    if with_index:
        from anqchor.pipelines.report import index_step

        idx = index_step(ctx, list(pipelines))
        graph.add(idx)
        for name in pipelines:
            groups[name].append(idx.name)
    groups["prepare"] = [s.name for s in graph if s.name.startswith("prepare.")]
    return graph, groups


__all__ = [
    "PIPELINES",
    "BuildContext",
    "Layout",
    "build_graph",
    "is_report_step",
    "pipeline_module",
]
