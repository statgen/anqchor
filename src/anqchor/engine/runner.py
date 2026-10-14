"""Invoke make on a run directory and report step status."""

from __future__ import annotations

import enum
import os
import re
import subprocess
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from anqchor.engine.minimake import log_path, signature_path, stamp_path
from anqchor.engine.step import StepGraph

FAILED_RE = re.compile(r"^\[anqchor\] FAILED (\S+) ")
RUNNING_RE = re.compile(r"^\[anqchor\] running (\S+)")
DRY_RUNNING_RE = re.compile(r"^echo '\[anqchor\] running (\S+)'")
"""``make -n`` prints the recipe instead of running it; this finds the steps it would run."""
NO_RULE_RE = re.compile(r"No rule to make target [`'](.+?)'(?:, needed by [`'](.+?)')?")
MAKE_ERROR_RE = re.compile(r"^(?:g?make)(?:\[\d+\])?: \*\*\* (.*)$")
MAKE_NOISE_RE = re.compile(
    r"^(?:g?make)(?:\[\d+\])?: (?:Nothing to be done for|[`'].+?' is up to date"
    r"|Entering directory|Leaving directory)"
)
"""make chatter that is not shown: anqchor prints its own summary of what ran."""


class StepStatus(str, enum.Enum):
    """Status of a step as shown by ``anqchor steps``."""

    DONE = "done"
    STALE = "stale"
    MISSING = "missing"


@dataclass
class RunResult:
    """Outcome of one make invocation.

    Attributes:
        returncode: make's exit status.
        ran: steps that started (in start order); empty for a dry run.
        would_run: for a dry run, the steps make would run.
        failed_step: the first step that failed, if any.
        log: that step's log file.
        message: make's own error when it failed outside a step (for example a
            missing input file: ``No rule to make target``), else ``None``.
    """

    returncode: int
    ran: list[str] = field(default_factory=list)
    failed_step: str | None = None
    log: Path | None = None
    would_run: list[str] = field(default_factory=list)
    message: str | None = None

    @property
    def ok(self) -> bool:
        """True when make succeeded."""
        return self.returncode == 0


def make_executable() -> str:
    """``$ANQCHOR_MAKE`` or ``make``."""
    return os.environ.get("ANQCHOR_MAKE", "make")


def invalidate_missing_outputs(
    graph: StepGraph, run_dir: Path, *, dry_run: bool = False
) -> list[str]:
    """Delete the stamp of every step with a missing output so make re-runs it.

    A stamp only proves the step finished once; if a user deleted an output since,
    the step must run again, and so must everything downstream of it. Returns the
    names of the steps whose outputs are missing. With ``dry_run`` nothing is deleted
    (the caller passes the steps to make as ``what_if`` instead).
    """
    run_dir = Path(run_dir)
    invalidated: list[str] = []
    for step in graph:
        stamp = run_dir / stamp_path(step.name)
        if stamp.exists() and any(not (run_dir / o).exists() for o in step.outputs):
            invalidated.append(step.name)
    if dry_run:
        return invalidated
    # Downstream steps go too: make compares mtimes at one-second resolution (3.81),
    # so a producer re-run in the same second as a consumer's last run would not
    # trigger the consumer.
    stale = set(invalidated)
    for name in invalidated:
        stale.update(graph.downstream(name))
    for name in stale:
        stamp = run_dir / stamp_path(name)
        if stamp.exists():
            stamp.unlink()
    return invalidated


def run_make(
    run_dir: Path,
    targets: Sequence[str],
    *,
    jobs: int = 1,
    dry_run: bool = False,
    force: bool = False,
    keep_going: bool = False,
    assume_old: Iterable[str] = (),
    what_if: Iterable[str] = (),
    stream: bool = True,
    env: dict[str, str] | None = None,
) -> RunResult:
    """Run ``make -C run_dir targets`` and parse anqchor's progress lines.

    Args:
        targets: make targets (pipeline groups, step names, or ``all``).
        jobs: ``-j`` value.
        dry_run: ``-n``; prints the commands that would run.
        force: ``-B``; rebuild the targets and everything they depend on.
        assume_old: step names whose stamps make should treat as up to date (``-o``).
        what_if: step names whose signatures make should treat as just modified (``-W``);
            with ``dry_run`` this shows what a parameter change would re-run without
            touching any file.
        env: extra environment variables for make and the steps (e.g. ``ANQCHOR_THREADS``,
            settings that must not be part of any step signature).
    """
    cmd = [make_executable(), "-C", str(run_dir), "-j", str(max(1, jobs))]
    if dry_run:
        cmd.append("-n")
    if force:
        cmd.append("-B")
    if keep_going:
        cmd.append("-k")
    for name in assume_old:
        cmd += ["-o", stamp_path(name)]
    for name in what_if:
        cmd += ["-W", signature_path(name)]
    cmd += list(targets)
    result = RunResult(returncode=0)
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env={**os.environ, **env} if env else None,
        )
    except OSError as e:
        result.returncode = 127
        result.message = f"cannot run {cmd[0]!r}: {e} (install GNU make or set $ANQCHOR_MAKE)"
        return result
    assert proc.stdout is not None
    make_error: str | None = None
    for line in proc.stdout:
        if stream and not MAKE_NOISE_RE.match(line):
            sys.stderr.write(line)
        m = RUNNING_RE.match(line)
        if m:
            result.ran.append(m.group(1))
        m = DRY_RUNNING_RE.match(line)
        if m:
            result.would_run.append(m.group(1))
        m = FAILED_RE.match(line)
        if m and result.failed_step is None:
            result.failed_step = m.group(1)
            result.log = Path(run_dir) / log_path(m.group(1))
        m = MAKE_ERROR_RE.match(line.rstrip("\n"))
        if m and make_error is None:
            make_error = m.group(1)
        if result.message is None:
            m = NO_RULE_RE.search(line)
            if m:
                needed = f" (needed by {m.group(2)})" if m.group(2) else ""
                result.message = f"missing file or unknown target: {m.group(1)}{needed}"
    result.returncode = proc.wait()
    if result.returncode != 0 and result.failed_step is None and result.message is None:
        result.message = make_error or f"make exited with status {result.returncode}"
    return result


def step_status(
    graph: StepGraph, run_dir: Path, names: Iterable[str] | None = None
) -> dict[str, StepStatus]:
    """Status of each step: ``done`` (up to date), ``stale`` (an input or parameter
    changed), or ``missing`` (never run, or an output was deleted)."""
    run_dir = Path(run_dir)
    out: dict[str, StepStatus] = {}
    for name in names if names is not None else graph.names():
        step = graph.get(name)
        stamp = run_dir / stamp_path(name)
        outputs_exist = all((run_dir / o).exists() for o in step.outputs)
        if not outputs_exist:
            out[name] = StepStatus.MISSING
            continue
        if not stamp.exists():  # outputs from an earlier run, invalidated since
            out[name] = StepStatus.STALE
            continue
        rc = subprocess.run(
            [make_executable(), "-C", str(run_dir), "-q", stamp_path(name)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode
        out[name] = StepStatus.DONE if rc == 0 else StepStatus.STALE
    return out
