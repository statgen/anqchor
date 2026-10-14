"""Run provenance: ``run.json`` with versions, inputs and timings."""

from __future__ import annotations

import contextlib
import datetime as _dt
import platform
import socket
import subprocess
import sys
from pathlib import Path
from typing import Any

from anqchor import __version__
from anqchor.engine.minimake import seconds_path
from anqchor.engine.step import StepGraph
from anqchor.util.files import read_json, sha256_file, sha256_text, write_json


def now_iso() -> str:
    """Current local time, ISO 8601 with offset, to the second."""
    return _dt.datetime.now(_dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def qpgen_provenance(binary: Path | None) -> dict[str, Any]:
    """Binary path and sha256, plus the git commit when the binary lives in a checkout."""
    if binary is None or not Path(binary).exists():
        return {"binary": str(binary) if binary else None}
    info: dict[str, Any] = {"binary": str(binary), "sha256": sha256_file(binary)}
    repo = Path(binary).resolve().parent.parent
    try:
        commit = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()
        info["git_commit"] = commit
    except (OSError, subprocess.SubprocessError):
        pass
    return info


def step_seconds(graph: StepGraph, run_dir: Path) -> dict[str, float]:
    """Wall time of the last successful run of each step that has run."""
    out: dict[str, float] = {}
    for name in graph.names():
        p = Path(run_dir) / seconds_path(name)
        if p.exists():
            with contextlib.suppress(ValueError, OSError):
                out[name] = float(p.read_text().strip())
    return out


def write_run_json(
    run_dir: Path,
    *,
    config_text: str,
    qpgen_binary: Path | None,
    argv: list[str] | None = None,
    started: str | None = None,
    finished: str | None = None,
    returncode: int | None = None,
    graph: StepGraph | None = None,
) -> Path:
    """Write or update ``run.json``. The history of invocations is kept in ``invocations``."""
    path = Path(run_dir) / "run.json"
    prev: dict[str, Any] = read_json(path) if path.exists() else {}
    invocation = {
        "argv": argv if argv is not None else sys.argv,
        "started": started,
        "finished": finished,
        "returncode": returncode,
    }
    data: dict[str, Any] = {
        "schema_version": 1,
        "anqchor_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "host": socket.gethostname(),
        "config_sha256": sha256_text(config_text),
        "qpgen": qpgen_provenance(qpgen_binary),
        "invocations": [*prev.get("invocations", [])[-19:], invocation],
    }
    if graph is not None:
        data["step_seconds"] = step_seconds(graph, run_dir)
    write_json(path, data)
    return path
