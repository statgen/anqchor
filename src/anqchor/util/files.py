"""File helpers: gzip-aware open, atomic writes, JSON, hashing."""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import io
import json
import os
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import IO, Any, cast


def is_gzip(path: str | os.PathLike[str]) -> bool:
    """True if the file starts with the gzip magic bytes (bgzip files included)."""
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_text(path: str | os.PathLike[str]) -> IO[str]:
    """Open a plain or gzip/bgzip-compressed text file for reading."""
    if is_gzip(path):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, encoding="utf-8")


def read_header(path: str | os.PathLike[str], comment: str | None = None) -> list[str]:
    """Return the tab-separated fields of the first line (skipping ``comment``-prefixed lines)."""
    with open_text(path) as fh:
        for line in fh:
            if comment and line.startswith(comment):
                continue
            return line.rstrip("\n").rstrip("\r").split("\t")
    return []


@contextlib.contextmanager
def atomic_write(path: str | os.PathLike[str], mode: str = "w") -> Iterator[IO[Any]]:
    """Write to a temporary file in the same directory and rename on success.

    A path ending in ``.gz`` is gzip-compressed (text mode). Readers never see a
    partially written file, which keeps make's view of outputs consistent.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    os.close(fd)
    try:
        if path.suffix == ".gz":
            raw = gzip.open(tmp, "wb", compresslevel=6)
            fh: IO[Any] = (
                io.TextIOWrapper(raw, encoding="utf-8") if "b" not in mode else cast(IO[Any], raw)
            )
        else:
            fh = open(tmp, mode, encoding=None if "b" in mode else "utf-8")
        with fh:
            yield fh
        # mkstemp creates 0600 files; give outputs the usual umask-based permissions so
        # reports and tables are readable by the people they are written for.
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(tmp)
        raise


def write_json(path: str | os.PathLike[str], obj: Any) -> None:
    """Atomically write ``obj`` as indented, key-sorted JSON (NaN becomes null)."""
    text = json.dumps(_jsonable(obj), indent=2, sort_keys=True, allow_nan=False)
    with atomic_write(path) as fh:
        fh.write(text + "\n")


def read_json(path: str | os.PathLike[str]) -> Any:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def write_text_if_changed(path: str | os.PathLike[str], text: str) -> bool:
    """Write ``text`` only if the file content differs; returns True when written.

    Leaving an unchanged file untouched preserves its mtime, which is how command
    signature files avoid triggering needless re-runs in make.
    """
    path = Path(path)
    try:
        if path.read_text(encoding="utf-8") == text:
            return False
    except (FileNotFoundError, UnicodeDecodeError):
        pass
    with atomic_write(path) as fh:
        fh.write(text)
    return True


def sha256_file(path: str | os.PathLike[str], chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _jsonable(obj: Any) -> Any:
    """Convert numpy scalars, paths, sets and non-finite floats into JSON-safe values."""
    import math

    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    try:
        import numpy as np

        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            v = float(obj)
            return v if math.isfinite(v) else None
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return [_jsonable(v) for v in obj.tolist()]
    except ImportError:  # pragma: no cover
        pass
    return obj


def resolve_near(
    name: str | os.PathLike[str], anchor: str | os.PathLike[str], max_up: int = 3
) -> Path:
    """Resolve a possibly relative path written inside the file ``anchor``.

    Absolute paths are returned as is. A relative path is tried against the anchor's
    directory and then up to ``max_up`` of its ancestors (lists are often written
    relative to a project directory above the list itself); the first existing match
    wins. When nothing exists, the path relative to the anchor's directory is returned
    so that error messages name a sensible location.
    """
    p = Path(name).expanduser()
    if p.is_absolute():
        return p
    base = Path(anchor).resolve().parent
    candidates = [base, *list(base.parents)[:max_up]]
    for d in candidates:
        c = d / p
        if c.exists() or Path(str(c) + ".pgen").exists():
            return c.resolve()
    return (base / p).resolve()
