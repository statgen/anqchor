"""The command builders against the real qpgentools, and the readers on its fresh output.

Needs the built binary (skipped otherwise) and the in-repo MAGE example data. Uses the
10-pair subset in ``tests/fixtures/qpgen/pairs.sub.tsv`` so each command takes seconds.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from anqchor.qpgen import commands as qc
from anqchor.qpgen.outputs import (
    iter_rect_assoc,
    read_match_all,
    read_match_weights,
    read_pair_assoc,
    read_prs,
    read_rect_header,
)

pytestmark = pytest.mark.integration

PHENO = "pheno/mage_expr.chr21_22.regenie.tsv.gz"


def _run(argv: list[str], cwd: Path) -> str:
    proc = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, check=False, timeout=600)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    return proc.stdout + proc.stderr


def test_builders_run_and_readers_parse(
    qpgen_bin: Path, mage_data: Path, mage_pgen_list: Path, qpgen_fixtures: Path, tmp_path: Path
) -> None:
    geno = qc.Genotypes(pgen_list=mage_pgen_list)
    pairs = qpgen_fixtures / "pairs.sub.tsv"
    pheno = mage_data / PHENO

    _run(qc.pair_prs(qpgen_bin, genotypes=geno, pairs=pairs, out_prefix=tmp_path / "prs"), tmp_path)
    prs = read_prs(tmp_path / "prs.prs.tsv.gz")
    assert prs.shape[0] == 729
    assert prs.shape[1] == 9  # 10 traits, one lost its only (absent) variant

    log = _run(
        qc.match_prs_pheno(
            qpgen_bin,
            prs=tmp_path / "prs.prs.tsv.gz",
            pheno=pheno,
            out_prefix=tmp_path / "match",
            threads=2,
        ),
        tmp_path,
    )
    assert "Analysis finished" in log
    m = read_match_all(tmp_path / "match.match.all.tsv.gz")
    assert len(m) == 729
    w = read_match_weights(tmp_path / "match.weights.tsv.gz")
    assert len(w) == 5  # the 5 of the 9 PRS traits that are in the example expression

    _run(
        qc.pair_assoc(
            qpgen_bin, genotypes=geno, pheno=pheno, pairs=pairs, out=tmp_path / "pa.tsv.gz"
        ),
        tmp_path,
    )
    pa = read_pair_assoc(tmp_path / "pa.tsv.gz")
    assert len(pa) == 5 and not pa.attrs["recomputed_log10p"]  # trait in the expression

    _run(
        qc.rect_assoc(
            qpgen_bin,
            genotypes=geno,
            pheno=pheno,
            var_list=qpgen_fixtures / "rect_variants.txt",
            pheno_list=qpgen_fixtures / "rect_traits.txt",
            out=tmp_path / "rect.tsv.gz",
            min_maf=0.01,
        ),
        tmp_path,
    )
    _, traits = read_rect_header(tmp_path / "rect.tsv.gz")
    assert len(traits) == 4  # NCAPH2 is not in the example expression
    chunks = list(iter_rect_assoc(tmp_path / "rect.tsv.gz"))
    assert sum(len(c.meta) for c in chunks) == 6
    assert all((c.meta["N"] == 729).all() for c in chunks)
