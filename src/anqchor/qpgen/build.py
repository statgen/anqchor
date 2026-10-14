"""Build qpgen from the ``extern/qpgen`` submodule.

Mirrors qpgen's documented procedure: htslib (autotools), qgenlib (cmake), then qpgen
(cmake + make). Each stage is skipped when its artifact already exists unless
``force`` is set.
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from anqchor.util.log import get_logger

log = get_logger("qpgen.build")


class BuildError(Exception):
    pass


def _run(cmd: Sequence[str], cwd: Path, logfile: Path) -> None:
    log.info("$ (cd %s && %s)", cwd, " ".join(cmd))
    with open(logfile, "a", encoding="utf-8") as fh:
        fh.write(f"\n$ (cd {cwd} && {' '.join(cmd)})\n")
        fh.flush()
        rc = subprocess.run(
            list(cmd), cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, check=False
        ).returncode
    if rc != 0:
        tail = logfile.read_text(encoding="utf-8", errors="replace").splitlines()[-30:]
        raise BuildError(
            f"command failed ({rc}): {' '.join(cmd)}\n" + "\n".join(tail) + f"\nfull log: {logfile}"
        )


def preflight() -> list[str]:
    """Return human-readable problems with the toolchain (empty when fine)."""
    problems = []
    for tool in ("cmake", "make", "autoreconf"):
        if shutil.which(tool) is None:
            problems.append(f"'{tool}' not found on PATH")
    if platform.system() == "Darwin":
        if not any(Path(p).exists() for p in ("/opt/homebrew/opt/libomp", "/usr/local/opt/libomp")):
            problems.append("OpenMP runtime not found: run 'brew install libomp'")
        brew = Path("/opt/homebrew/opt")
        if not (brew.exists() and any(brew.glob("llvm*"))):
            problems.append(
                "Homebrew LLVM not found (recommended on macOS): run 'brew install llvm'"
            )
    return problems


def build(
    src: Path,
    *,
    jobs: int = 4,
    force: bool = False,
    cmake_args: Sequence[str] = (),
    logfile: Path | None = None,
) -> Path:
    """Build qpgen in ``src`` and return the path of ``bin/qpgentools``."""
    src = Path(src).resolve()
    if not (src / "CMakeLists.txt").exists():
        raise BuildError(
            f"{src} is not a qpgen checkout. Run 'git submodule update --init --recursive' first."
        )
    logfile = logfile or (src / "anqchor-build.log")
    logfile.write_text("", encoding="utf-8")
    j = f"-j{max(1, jobs)}"
    sub = src / "submodules"
    for name in ("htslib", "qgenlib", "eigen"):
        d = sub / name
        if not d.is_dir() or not any(d.iterdir()):
            raise BuildError(
                f"submodule {d} is empty. Run 'git submodule update --init --recursive'."
            )
    hts = sub / "htslib"
    if force or not (hts / "libhts.a").exists():
        _run(["autoreconf", "-i"], hts, logfile)
        _run(["./configure"], hts, logfile)
        _run(["make", j], hts, logfile)
    qgen = sub / "qgenlib"
    if force or not (qgen / "lib" / "libqgen.a").exists():
        (qgen / "build").mkdir(exist_ok=True)
        _run(["cmake", ".."], qgen / "build", logfile)
        _run(["make", j], qgen / "build", logfile)
    bdir = src / "build"
    bdir.mkdir(exist_ok=True)
    if force or not (bdir / "CMakeCache.txt").exists() or cmake_args:
        _run(["cmake", *cmake_args, ".."], bdir, logfile)
    _run(["make", j], bdir, logfile)
    binary = src / "bin" / "qpgentools"
    if not binary.exists() or not os.access(binary, os.X_OK):
        raise BuildError(f"build finished but {binary} is missing; see {logfile}")
    return binary
