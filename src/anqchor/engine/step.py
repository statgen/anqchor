"""The :class:`Step` unit of work and the graph of steps that forms a run.

A step declares its input files, output files and shell commands. Steps are linked
by shared paths: when one step's input is another step's output, the Makefile makes
the consumer depend on the producer. Python-implemented steps run as
``python -m anqchor _step <name> --run-dir <dir>`` and receive their parameters from
the step's signature file ``steps/<name>.json`` (see :mod:`anqchor.engine.minimake`).
"""

from __future__ import annotations

import shlex
import sys
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from anqchor import __version__

STEP_NAME_CHARS = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")


@dataclass(frozen=True)
class Step:
    """One unit of work.

    Attributes:
        name: unique dotted name, e.g. ``match.rnaseq.qpgen``. Also a make target.
        inputs: files read by the step. Relative paths are relative to the run directory.
        outputs: files written by the step (relative to the run directory).
        commands: shell commands run in order from the run directory; empty for a pure
            Python step whose ``kind`` is set (the command is then generated).
        kind: name of a registered Python step implementation (``anqchor.steps``), or
            ``None`` for a shell-only step.
        params: JSON-serialisable parameters. They are part of the step signature, so
            changing any of them re-runs the step and everything downstream.
        pipeline: the pipeline the step belongs to (``prepare`` for shared steps).
        description: one line shown by ``anqchor steps``.
    """

    name: str
    inputs: tuple[Path, ...] = ()
    outputs: tuple[Path, ...] = ()
    commands: tuple[str, ...] = ()
    kind: str | None = None
    params: dict[str, Any] = field(default_factory=dict, hash=False, compare=False)
    pipeline: str = "prepare"
    description: str = ""

    def __post_init__(self) -> None:
        if not self.name or set(self.name) - STEP_NAME_CHARS:
            raise ValueError(f"invalid step name {self.name!r}")
        if not self.outputs:
            raise ValueError(f"step {self.name} declares no outputs")
        if isinstance(self.commands, str):  # a bare string would iterate per character
            object.__setattr__(self, "commands", (self.commands,))
        object.__setattr__(self, "commands", tuple(self.commands))
        if not self.commands and self.kind is None:
            raise ValueError(f"step {self.name} has neither commands nor a Python kind")
        for c in self.commands:
            if not isinstance(c, str) or "\n" in c or "\r" in c:
                raise ValueError(
                    f"step {self.name}: each command must be a single-line string "
                    f"(one make recipe line), got {c!r}"
                )
        for attr in ("inputs", "outputs"):
            value = getattr(self, attr)
            if isinstance(value, (str, Path)):
                value = (value,)
            object.__setattr__(self, attr, tuple(Path(p) for p in value))

    def all_commands(self) -> tuple[str, ...]:
        """Shell commands including the generated Python-step invocation."""
        if self.kind is None:
            return self.commands
        return (*self.commands, python_step_command(self.name))

    def signature(self) -> dict[str, Any]:
        """Content of ``steps/<name>.json``: everything that defines the step's result.

        Python steps also record the anqchor version, so upgrading anqchor re-runs the
        steps whose code it implements (qpgen steps carry the binary's hash in
        ``params``, see :func:`anqchor.pipelines.with_tool_identity`)."""
        sig: dict[str, Any] = {
            "name": self.name,
            "kind": self.kind,
            "pipeline": self.pipeline,
            "inputs": [str(p) for p in self.inputs],
            "outputs": [str(p) for p in self.outputs],
            "commands": list(self.commands),
            "params": self.params,
        }
        if self.kind is not None:
            sig["anqchor_version"] = __version__
        return sig


def python_step_command(name: str) -> str:
    """Shell command that runs a Python step with the current interpreter."""
    return f"{shlex.quote(sys.executable)} -m anqchor _step {shlex.quote(name)} --run-dir ."


class StepGraph:
    """An ordered, de-duplicated collection of steps with producer lookup."""

    def __init__(self, steps: Iterable[Step] = ()) -> None:
        self._steps: dict[str, Step] = {}
        self._producer: dict[Path, str] = {}
        for s in steps:
            self.add(s)

    def add(self, step: Step) -> None:
        """Add a step. Re-adding an identical step (same signature) is a no-op, which is
        how pipelines share prepare and PRS steps."""
        prev = self._steps.get(step.name)
        if prev is not None:
            if prev.signature() != step.signature():
                raise ValueError(f"conflicting definitions for step {step.name}")
            return
        for out in step.outputs:
            other = self._producer.get(out)
            if other is not None:
                raise ValueError(f"output {out} is produced by both {other} and {step.name}")
            self._producer[out] = step.name
        self._steps[step.name] = step

    def extend(self, steps: Iterable[Step]) -> None:
        """Add several steps (see :meth:`add`)."""
        for s in steps:
            self.add(s)

    def __iter__(self) -> Iterator[Step]:
        return iter(self._steps.values())

    def __len__(self) -> int:
        return len(self._steps)

    def __contains__(self, name: object) -> bool:
        return name in self._steps

    def get(self, name: str) -> Step:
        """The step called ``name`` (``KeyError`` naming it when unknown)."""
        try:
            return self._steps[name]
        except KeyError:
            raise KeyError(f"unknown step {name!r}") from None

    def names(self) -> list[str]:
        """Step names in insertion order."""
        return list(self._steps)

    def producer_of(self, path: Path) -> str | None:
        """Name of the step that writes ``path`` (``None`` for an external input)."""
        return self._producer.get(Path(path))

    def upstream(self, name: str) -> list[str]:
        """Names of the steps that ``name`` directly depends on."""
        deps: list[str] = []
        for inp in self.get(name).inputs:
            prod = self._producer.get(inp)
            if prod is not None and prod != name and prod not in deps:
                deps.append(prod)
        return deps

    def closure(self, names: Iterable[str]) -> list[str]:
        """``names`` plus everything upstream, in topological order."""
        order: list[str] = []
        seen: set[str] = set()
        visiting: set[str] = set()

        def visit(n: str) -> None:
            if n in seen:
                return
            if n in visiting:
                raise ValueError(f"dependency cycle through step {n}")
            visiting.add(n)
            for d in self.upstream(n):
                visit(d)
            visiting.discard(n)
            seen.add(n)
            order.append(n)

        for n in names:
            visit(n)
        return order

    def downstream(self, name: str) -> list[str]:
        """Steps that (transitively) depend on ``name``."""
        result: list[str] = []
        frontier = [name]
        while frontier:
            cur = frontier.pop()
            for s in self._steps.values():
                if s.name not in result and cur in self.upstream(s.name):
                    result.append(s.name)
                    frontier.append(s.name)
        return result

    def by_pipeline(self, pipeline: str) -> list[Step]:
        """Steps whose ``pipeline`` attribute equals ``pipeline``."""
        return [s for s in self._steps.values() if s.pipeline == pipeline]
