"""Shared pytest fixtures.

``mage_data`` is the compact MAGE chr21/22 example shipped in the repository
(``examples/mage/data``, 729 samples), so the integration tests need no download.
``qpgen_bin`` skips the requesting test when no built qpgentools binary is found
(``$ANQCHOR_QPGEN``, ``PATH``, then ``extern/qpgen/bin``), so the unit suite runs anywhere.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
FIXTURES = Path(__file__).resolve().parent / "fixtures"
EXAMPLE_DATA = REPO / "examples" / "mage" / "data"
EXAMPLE_PGEN_LIST = EXAMPLE_DATA / "geno" / "mage.pgen_list.tsv"
SUBMODULE_QPGEN = REPO / "extern" / "qpgen" / "bin" / "qpgentools"


def _qpgen_binary() -> Path | None:
    for cand in (os.environ.get("ANQCHOR_QPGEN"), shutil.which("qpgentools"), SUBMODULE_QPGEN):
        if cand and Path(cand).is_file() and os.access(cand, os.X_OK):
            return Path(cand).resolve()
    return None


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    """``tests/fixtures``."""
    return FIXTURES


@pytest.fixture(scope="session")
def qpgen_fixtures() -> Path:
    """Small excerpts of real qpgentools outputs (see ``tests/fixtures/qpgen/README.md``)."""
    return FIXTURES / "qpgen"


@pytest.fixture(scope="session")
def qpgen_bin() -> Path:
    """A built ``qpgentools``; skips the test when none is available."""
    b = _qpgen_binary()
    if b is None:
        pytest.skip("qpgentools not built (set ANQCHOR_QPGEN or run 'anqchor qpgen build')")
    return b


@pytest.fixture(scope="session")
def mage_data() -> Path:
    """The in-repo MAGE chr21/22 example data (``examples/mage/data``)."""
    if not EXAMPLE_PGEN_LIST.exists():
        pytest.skip(f"example data missing: {EXAMPLE_PGEN_LIST}")
    return EXAMPLE_DATA


@pytest.fixture(scope="session")
def mage_pgen_list(mage_data: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The example pgen list rewritten with absolute prefixes (qpgen resolves relative
    prefixes against its working directory)."""
    src = mage_data / "geno" / "mage.pgen_list.tsv"
    out = tmp_path_factory.mktemp("mage") / "pgen_list.abs.tsv"
    lines = []
    for line in src.read_text().splitlines():
        f = line.split("\t")
        if len(f) == 4 and not Path(f[3]).is_absolute():
            f[3] = str(src.parent / f[3])
        lines.append("\t".join(f))
    out.write_text("\n".join(lines) + "\n")
    return out
