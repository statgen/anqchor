"""Find the ``qpgentools`` binary.

Resolution order: explicit path (``--qpgen`` / config ``qpgen.binary``) →
``$ANQCHOR_QPGEN`` → ``qpgentools`` on ``PATH`` → the submodule build
``extern/qpgen/bin/qpgentools`` of a source checkout.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

REQUIRED_COMMANDS = ("pair-prs", "match-prs-pheno", "pair-assoc", "rect-assoc")
REQUIRED_FLAGS = {
    c: ("--cat-cov-names", "--max-cat-levels")
    for c in ("match-prs-pheno", "pair-assoc", "rect-assoc")
}
"""Flags newer than ``qpgentools --version`` can tell apart (both 76f2fd1 and 5b811e2
print 0.2.0): categorical covariates need qpgen 5b811e2 or later."""
ENV_VAR = "ANQCHOR_QPGEN"


class QpgenNotFound(Exception):
    """No usable ``qpgentools`` binary (exit code 4)."""


def source_root() -> Path | None:
    """Repository root when running from a source checkout (``src/anqchor`` layout)."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "extern" / "qpgen").is_dir() and (parent / "pyproject.toml").exists():
            return parent
    return None


def submodule_dir() -> Path | None:
    """``extern/qpgen`` of the source checkout, if any."""
    root = source_root()
    return root / "extern" / "qpgen" if root else None


@dataclass(frozen=True)
class Located:
    """A discovered binary and how it was found."""

    path: Path
    source: str
    """``explicit``, ``env``, ``path`` or ``submodule``."""


def locate(explicit: str | os.PathLike[str] | None = None) -> Located:
    """Return the first candidate that exists and is executable.

    An explicit path or ``$ANQCHOR_QPGEN`` that does not name an executable file raises
    :class:`QpgenNotFound` immediately; the later candidates are only searched when
    neither is set."""
    tried: list[str] = []
    if explicit:
        p = Path(explicit).expanduser()
        if _usable(p):
            return Located(p.resolve(), "explicit")
        raise QpgenNotFound(f"qpgentools not found or not executable: {p}")
    env = os.environ.get(ENV_VAR)
    if env:
        # An explicitly configured binary that is missing is an error, not a reason
        # to silently run a different qpgentools from PATH.
        p = Path(env).expanduser()
        if _usable(p):
            return Located(p.resolve(), "env")
        raise QpgenNotFound(f"${ENV_VAR} names a missing or non-executable file: {env}")
    on_path = shutil.which("qpgentools")
    if on_path:
        return Located(Path(on_path).resolve(), "path")
    tried.append("qpgentools on PATH")
    sub = submodule_dir()
    if sub is not None:
        p = sub / "bin" / "qpgentools"
        if _usable(p):
            return Located(p, "submodule")
        tried.append(str(p))
    raise QpgenNotFound(
        "qpgentools not found (tried: " + "; ".join(tried) + "). Build it with "
        "'anqchor qpgen build', or pass --qpgen PATH, or set $ANQCHOR_QPGEN."
    )


def find_qpgentools(explicit: str | os.PathLike[str] | None = None) -> Path:
    """Path of the first usable ``qpgentools`` (see :func:`locate`)."""
    return locate(explicit).path


def _usable(p: Path) -> bool:
    return p.is_file() and os.access(p, os.X_OK)


def available_commands(binary: Path) -> list[str]:
    """Parse ``qpgentools --help`` for the list of subcommands."""
    try:
        proc = subprocess.run(
            [str(binary), "--help"], capture_output=True, text=True, timeout=60, check=False
        )
    except (OSError, subprocess.SubprocessError) as e:
        raise QpgenNotFound(f"cannot run {binary} --help: {e}") from e
    text = proc.stdout + proc.stderr
    cmds = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("- "):
            cmds.append(s[2:].split(":", 1)[0].strip())
    return cmds


def check(binary: Path) -> list[str]:
    """Return the required subcommands missing from ``binary``."""
    have = set(available_commands(binary))
    return [c for c in REQUIRED_COMMANDS if c not in have]


def missing_flags(binary: Path) -> list[str]:
    """Required flags (:data:`REQUIRED_FLAGS`) missing from ``binary``, as
    ``"<command> <flag>"``; read from each subcommand's ``--help``."""
    out = []
    for cmd, flags in REQUIRED_FLAGS.items():
        try:
            proc = subprocess.run(
                [str(binary), cmd, "--help"],
                capture_output=True,
                text=True,
                timeout=60,
                check=False,
            )
        except (OSError, subprocess.SubprocessError) as e:
            raise QpgenNotFound(f"cannot run {binary} {cmd} --help: {e}") from e
        text = proc.stdout + proc.stderr
        out += [f"{cmd} {f}" for f in flags if f not in text]
    return out
