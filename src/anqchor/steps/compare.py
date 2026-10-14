"""Shared parts of ``anqchor cis compare`` and ``anqchor trans compare``.

A comparison reads the summaries of finished runs (never qpgen), labels each run (its
directory name by default) and writes its tables into a comparison directory. Both
comparisons share the run labels, the reference hash used for the comparability check,
the copy of each run's settings (``settings/<label>.anqchor.yaml``) and the manifest
(``compare.json``) that ``cis report`` / ``trans report`` use to recognise the directory.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from anqchor.pipelines.layout import Layout
from anqchor.steps import StepError
from anqchor.util.files import open_text, read_json, sha256_text, write_json

MANIFEST = Layout.compare_manifest.as_posix()
COMPARISON_DATASET = "comparison"
"""Name under which a comparison appears in a report (it has no dataset of its own)."""


def comparison_kind(path: Path) -> str | None:
    """``cis`` or ``trans`` for a comparison directory, ``None`` otherwise."""
    f = Path(path) / Layout.compare_manifest
    if not f.exists():
        return None
    try:
        return str(read_json(f).get("kind") or "cis")
    except (OSError, ValueError):
        return None


def is_comparison_dir(path: Path, kind: str | None = None) -> bool:
    """True for a comparison directory (of ``kind`` when given)."""
    k = comparison_kind(path)
    return k is not None and (kind is None or k == kind)


def _as_given(p: Path) -> Path:
    """The absolute path as given (``..`` collapsed, symbolic links kept), so that a label
    names the directory the user typed, not the target of a link."""
    return Path(os.path.abspath(p))


def default_labels(run_dirs: Sequence[Path]) -> list[str]:
    """Directory names; ``<parent>_<name>`` where names collide."""
    names = [_as_given(p).name for p in run_dirs]
    if len(set(names)) == len(names):
        return names
    return [f"{_as_given(p).parent.name}_{_as_given(p).name}" for p in run_dirs]


def clean_label(label: str) -> str:
    """Labels name columns and files (``settings/<label>.anqchor.yaml``); keep them free of ``/``."""
    return label.replace("/", "_").strip() or "run"


def run_labels(run_dirs: Sequence[Path], labels: Sequence[str] | None) -> list[str]:
    """The label of each run: ``labels`` (one per run) or the directory names."""
    if not run_dirs:
        raise StepError("give at least one run directory")
    seen: set[Path] = set()
    for p in run_dirs:
        r = Path(p).resolve()
        if r in seen:
            raise StepError(f"{p} given more than once")
        seen.add(r)
    if labels is not None and len(labels) != len(run_dirs):
        raise StepError(f"{len(labels)} --label value(s) for {len(run_dirs)} run(s)")
    return list(labels) if labels else default_labels(run_dirs)


def check_unique(labels: Sequence[str]) -> None:
    dup = sorted({x for x in labels if list(labels).count(x) > 1})
    if dup:
        raise StepError(f"duplicate run label(s) {', '.join(dup)}; give distinct --label values")


def prepare_out_dir(out_dir: Path, kind: str, run_dirs: Sequence[Path]) -> None:
    """Check that ``out_dir`` may hold a ``kind`` comparison and remove what a previous
    comparison there would leave stale (its report; ``settings/`` is pruned when written).

    Refuses a run directory (one of ``run_dirs`` or any directory holding
    ``anqchor.resolved.yaml``) and a comparison of the other kind."""
    out_dir = Path(out_dir)
    resolved = out_dir.resolve()
    if (
        any(Path(p).resolve() == resolved for p in run_dirs)
        or (out_dir / Layout.resolved_config).exists()
    ):
        raise StepError(
            f"--out {out_dir} is a run directory; choose an empty directory or a previous "
            f"{kind} comparison"
        )
    other = comparison_kind(out_dir)
    if other is not None and other != kind:
        raise StepError(
            f"--out {out_dir} holds a {other} comparison; choose an empty directory or a "
            f"previous {kind} comparison"
        )
    stale = out_dir / Layout.report(kind)
    if stale.exists():
        stale.unlink()


def file_sha256(path: Path) -> str:
    """sha256 of a text file's decompressed content (gzip headers carry a timestamp);
    ``?`` when the file is missing."""
    path = Path(path)
    if not path.exists():
        return "?"
    with open_text(path) as fh:
        return sha256_text(fh.read())


def write_settings_and_manifest(
    out_dir: Path,
    kind: str,
    runs: Sequence[dict[str, Any]],
    argv: Sequence[str] | None,
) -> None:
    """Copy each run's ``anqchor.resolved.yaml`` to ``settings/<label>.anqchor.yaml`` and
    write ``compare.json``. ``runs``: dicts with ``label``, ``run_dir`` and ``dataset``."""
    out_dir = Path(out_dir)
    keep = {Layout.compare_settings(r["label"]).name for r in runs}
    settings_dir = (out_dir / Layout.compare_settings("x")).parent
    if settings_dir.is_dir():
        for f in settings_dir.glob("*.anqchor.yaml"):
            if f.name not in keep:
                f.unlink()
    for r in runs:
        src = Path(r["run_dir"]) / Layout.resolved_config
        dst = out_dir / Layout.compare_settings(r["label"])
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    write_json(
        out_dir / Layout.compare_manifest,
        {
            "schema_version": 1,
            "kind": kind,
            "runs": [
                {
                    "label": r["label"],
                    "run_dir": r["run_dir"],
                    "dataset": r["dataset"],
                    "settings": Layout.compare_settings(r["label"]).as_posix(),
                }
                for r in runs
            ],
            "labels": [r["label"] for r in runs],
            "argv": list(argv) if argv is not None else None,
            "cwd": os.getcwd(),
        },
    )
