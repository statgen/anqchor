"""Python step implementations and their registry.

A Python step is declared by a pipeline as ``Step(kind="<kind>", params={...})`` and
executed by make as ``python -m anqchor _step <name> --run-dir .``. The runner loads
``steps/<name>.json`` and calls the function registered for ``kind`` with a
:class:`StepContext`.

Rules for implementations:

* Read only the step's declared ``inputs`` and ``params``; never read the config file.
  Anything the step depends on must be in ``params`` so that changing it re-runs the step.
* Write every declared output (atomically, see :func:`anqchor.util.files.atomic_write`).
* Raise :class:`StepError` with an actionable message for data problems.
"""

from __future__ import annotations

import importlib
import json
import os
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from anqchor.util.log import get_logger

STEP_MODULES = (
    "anqchor.steps.prepare",
    "anqchor.steps.strategy",
    "anqchor.steps.match",
    "anqchor.steps.cis",
    "anqchor.steps.trans",
    "anqchor.report.match",
    "anqchor.report.cis",
    "anqchor.report.trans",
    "anqchor.report.index",
)

StepFunction = Callable[["StepContext"], None]
REGISTRY: dict[str, StepFunction] = {}


class StepError(Exception):
    """A data or input problem detected by a step (reported without a traceback)."""


def register(kind: str) -> Callable[[StepFunction], StepFunction]:
    """Decorator registering a step implementation under ``kind``."""

    def deco(fn: StepFunction) -> StepFunction:
        if kind in REGISTRY and REGISTRY[kind] is not fn:
            raise ValueError(f"step kind {kind!r} registered twice")
        REGISTRY[kind] = fn
        return fn

    return deco


@dataclass
class StepContext:
    """Everything a step implementation receives."""

    name: str
    run_dir: Path
    params: dict[str, Any]
    inputs: list[Path]
    outputs: list[Path]
    log: Any = field(default=None)

    def path(self, rel: str | Path) -> Path:
        """Absolute path of a run-directory-relative path."""
        p = Path(rel)
        return p if p.is_absolute() else self.run_dir / p


def load_all() -> None:
    """Import every step module so that their ``@register`` decorators run."""
    for mod in STEP_MODULES:
        try:
            importlib.import_module(mod)
        except ModuleNotFoundError as e:  # a module not implemented yet
            if e.name != mod:
                raise


def run_step_now(step: Any, run_dir: Path) -> None:
    """Execute a Python :class:`~anqchor.engine.step.Step` directly, without make.

    Used by ``anqchor <pipeline> report``, which re-renders a report from existing
    summaries every time it is called. The step's declared inputs must exist."""
    run_dir = Path(run_dir).resolve()
    missing = [str(p) for p in step.inputs if not (run_dir / p).exists()]
    if missing:
        raise StepError("missing inputs (run the analysis first):\n  " + "\n  ".join(missing))
    load_all()
    fn = REGISTRY.get(step.kind)
    if fn is None:
        raise StepError(f"no implementation registered for step kind {step.kind!r}")
    ctx = StepContext(
        name=step.name,
        run_dir=run_dir,
        params=json.loads(json.dumps(step.params, default=str)),
        inputs=list(step.inputs),
        outputs=list(step.outputs),
        log=get_logger(f"step.{step.name}"),
    )
    cwd = os.getcwd()
    os.chdir(run_dir)
    try:
        fn(ctx)
    finally:
        os.chdir(cwd)


def run_step(name: str, run_dir: Path) -> None:
    """Execute the Python step ``name`` of the run in ``run_dir``."""
    run_dir = Path(run_dir).resolve()
    sig_path = run_dir / "steps" / f"{name}.json"
    if not sig_path.exists():
        raise StepError(f"unknown step {name}: {sig_path} not found (regenerate the Makefile)")
    sig = json.loads(sig_path.read_text(encoding="utf-8"))
    kind = sig.get("kind")
    if not kind:
        raise StepError(f"step {name} is not a Python step")
    load_all()
    fn = REGISTRY.get(kind)
    if fn is None:
        raise StepError(f"no implementation registered for step kind {kind!r}")
    ctx = StepContext(
        name=name,
        run_dir=run_dir,
        params=sig.get("params", {}),
        inputs=[Path(p) for p in sig.get("inputs", [])],
        outputs=[Path(p) for p in sig.get("outputs", [])],
        log=get_logger(f"step.{name}"),
    )
    cwd = os.getcwd()
    os.chdir(run_dir)
    try:
        fn(ctx)
    finally:
        os.chdir(cwd)
