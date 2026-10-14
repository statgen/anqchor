"""Low-level table helpers shared by the input readers.

Everything here is format-agnostic: separator detection, first-line access, atomic
(b)gzip-compatible writers and the :class:`InputError` raised for bad input data.
"""

from __future__ import annotations

import contextlib
import gzip
import io
import os
import re
import shutil
import tempfile
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import IO, Any

import pandas as pd

from anqchor.util.files import is_gzip, open_text

GZ_LEVEL = 3
"""gzip level for the large matrices anqchor writes (level 3 is ~2x faster than 6)."""

INT_RE = re.compile(r"^[+-]?\d+$")


class InputError(ValueError):
    """An input file is malformed or inconsistent. The message says how to fix it."""


def read_user_csv(path: str | os.PathLike[str], what: str, **kw: Any) -> pd.DataFrame:
    """``pd.read_csv`` for a user input file, with compression detected from the file's
    content (gzip/bgzip magic bytes) like :func:`first_line`, never from its extension:
    a bgzip ``.bgz``, a gzip file named ``.tsv`` or a plain file named ``.gz`` all work.
    Decoding and parse errors become an :class:`InputError` naming the file."""
    path = Path(path)
    kw.setdefault("compression", "gzip" if is_gzip(path) else None)
    try:
        return pd.read_csv(path, **kw)
    except (pd.errors.ParserError, ValueError, UnicodeDecodeError, OSError, EOFError) as e:
        raise InputError(f"{path}: cannot parse as a {what}: {e}") from e


def first_line(path: str | os.PathLike[str], skip_prefix: str | None = None) -> str:
    """The first non-empty line of a (possibly gzipped) text file without its newline.

    Lines starting with ``skip_prefix`` (e.g. ``##``) are skipped. Returns ``""`` for an
    empty file.
    """
    with open_text(path) as fh:
        for line in fh:
            s = line.rstrip("\n").rstrip("\r")
            if not s.strip():
                continue
            if skip_prefix and s.startswith(skip_prefix):
                continue
            return s
    return ""


def detect_sep(line: str) -> str:
    """``"\\t"`` when the line has a tab, otherwise the whitespace regex ``r"\\s+"``."""
    return "\t" if "\t" in line else r"\s+"


def split_fields(line: str, sep: str) -> list[str]:
    """Split a line with a separator returned by :func:`detect_sep`."""
    return line.split("\t") if sep == "\t" else line.split()


def is_int(s: str) -> bool:
    return bool(INT_RE.match(s.strip()))


def duplicates(values: Iterable[str]) -> list[str]:
    """Values that occur more than once, in order of their second occurrence."""
    seen: set[str] = set()
    dup: list[str] = []
    dup_set: set[str] = set()
    for v in values:
        if v in seen and v not in dup_set:
            dup.append(v)
            dup_set.add(v)
        seen.add(v)
    return dup


def examples(values: Sequence[str] | Iterable[str], n: int = 5) -> str:
    """``'a', 'b', 'c' (and 7 more)`` for error messages."""
    vals = list(values)
    head = ", ".join(repr(v) for v in vals[:n])
    return head + (f" (and {len(vals) - n} more)" if len(vals) > n else "")


@contextlib.contextmanager
def atomic_text_writer(path: str | os.PathLike[str], level: int = GZ_LEVEL) -> Iterator[IO[str]]:
    """Atomically write a text file; ``.gz`` paths are gzip-compressed at ``level``.

    Same contract as :func:`anqchor.util.files.atomic_write` (temporary file in the
    target directory, renamed on success) but with a configurable compression level,
    which matters for multi-hundred-megabyte phenotype matrices.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    os.close(fd)
    try:
        if path.suffix == ".gz":
            raw = gzip.open(tmp, "wb", compresslevel=level)
            fh: IO[str] = io.TextIOWrapper(raw, encoding="utf-8", newline="\n")
        else:
            fh = open(tmp, "w", encoding="utf-8", newline="\n")
        with fh:
            yield fh
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(tmp)
        raise


def atomic_copy_gz(src: str | os.PathLike[str], dest: str | os.PathLike[str]) -> None:
    """Copy ``src`` to the ``.gz`` path ``dest`` atomically, byte for byte when ``src`` is
    already gzip/bgzip-compressed and compressing it otherwise."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=dest.parent, prefix=f".{dest.name}.", suffix=".tmp")
    os.close(fd)
    try:
        if is_gzip(src):
            shutil.copyfile(src, tmp)
        else:
            with open(src, "rb") as fi, gzip.open(tmp, "wb", compresslevel=GZ_LEVEL) as fo:
                shutil.copyfileobj(fi, fo, 1 << 20)
        os.replace(tmp, dest)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(tmp)
        raise


def write_lines(path: str | os.PathLike[str], lines: Iterable[str]) -> int:
    """Atomically write one string per line; returns the number of lines."""
    n = 0
    with atomic_text_writer(path) as fh:
        for line in lines:
            fh.write(line)
            fh.write("\n")
            n += 1
    return n


def read_lines(path: str | os.PathLike[str]) -> list[str]:
    """Non-empty lines of a text file, stripped of trailing whitespace."""
    with open_text(path) as fh:
        return [s for s in (line.rstrip() for line in fh) if s]


def fmt_float(v: float) -> str:
    """Shortest round-trip text of a float, ``NA`` for NaN, integers without ``.0``."""
    if v != v:
        return "NA"
    r = repr(v)
    return r[:-2] if r.endswith(".0") else r
