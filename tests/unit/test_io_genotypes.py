"""Pgen lists, psam IDs, chromosome style and pvar indexing (anqchor.io.genotypes)."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fixtures" / "io"))
from anq_io_fixtures import make_geno, needs_htslib

from anqchor.io.genotypes import (
    chrom_style,
    index_pvar,
    normalize_chrom,
    parse_pgen_list,
    pvar_chroms,
    read_psam_iids,
)
from anqchor.io.tables import InputError


@pytest.mark.parametrize(
    ("name", "style", "expected"),
    [
        ("21", "chr", "chr21"),
        ("chr21", "chr", "chr21"),
        ("chr21", "plain", "21"),
        ("X", "plain", "X"),
        ("chrX", "plain", "X"),
        ("23", "chr", "chr23"),
    ],
)
def test_normalize_chrom(name: str, style: str, expected: str) -> None:
    assert normalize_chrom(name, style) == expected  # type: ignore[arg-type]


def test_chrom_style() -> None:
    assert (
        chrom_style("chr1") == "chr"
        and chrom_style("1") == "plain"
        and chrom_style("CHR1") == "chr"
    )


def test_parse_pgen_list_four_and_six_columns(tmp_path: Path) -> None:
    sub = tmp_path / "proj" / "lists"
    sub.mkdir(parents=True)
    (tmp_path / "proj" / "geno").mkdir()
    (tmp_path / "proj" / "geno" / "a.pgen").write_bytes(b"")
    lst = sub / "l.tsv"
    lst.write_text(
        "# comment\n\n1 1 100 geno/a\n2\t1\t100\t/abs/b.pgen\t/abs/b.psam\t/abs/b.pvar.idx.gz\n"
    )
    e = parse_pgen_list(lst)
    assert (
        e[0].pgen == (tmp_path / "proj" / "geno" / "a.pgen").resolve()
    )  # resolved via an ancestor
    assert e[0].pivar.name == "a.pvar.idx.gz"
    assert e[1].psam == Path("/abs/b.psam") and e[1].beg == 1 and e[1].end == 100
    lst.write_text("1 1 100\n")
    with pytest.raises(InputError, match="expected 4 fields"):
        parse_pgen_list(lst)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("#IID\tSEX\nA\t1\nB\t2\n", ["A", "B"]),
        ("#FID\tIID\tSEX\nf\tA\t1\nf\tB\t2\n", ["A", "B"]),
        ("f A 0 0 1 -9\nf B 0 0 1 -9\n", ["A", "B"]),
        ("A\nB\n", ["A", "B"]),
    ],
)
def test_read_psam(tmp_path: Path, text: str, expected: list[str]) -> None:
    p = tmp_path / "x.psam"
    p.write_text(text)
    assert read_psam_iids(p) == expected


@needs_htslib
def test_pvar_chroms_tabix_and_stream(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    make_geno(tmp_path, style="chr")
    piv = tmp_path / "geno" / "set.chr1.pvar.idx.gz"
    assert pvar_chroms(piv) == ["chr1"]
    monkeypatch.setenv("ANQCHOR_NO_TABIX", "1")
    assert pvar_chroms(piv) == ["chr1"]


@needs_htslib
@pytest.mark.parametrize("header", [True, False])
def test_index_pvar(tmp_path: Path, header: bool) -> None:
    prefix = tmp_path / "x"
    if header:
        text = "##fileformat=PVARv1.0\n#CHROM\tPOS\tID\tREF\tALT\n1\t10\trs1\tA\tG\n1\t20\trs2\tC\tT\n2\t5\trs3\tG\tA\n"
    else:  # .bim order: CHROM ID CM POS ALT REF
        text = "1\trs1\t0\t10\tG\tA\n1\trs2\t0\t20\tT\tC\n2\trs3\t0\t5\tA\tG\n"
    (tmp_path / "x.pvar").write_text(text)
    out, n = index_pvar(prefix)
    assert n == 3 and Path(str(out) + ".tbi").exists()
    lines = gzip.decompress(out.read_bytes()).decode().splitlines()
    assert lines[0] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tINDEX"
    assert lines[1] == "1\t10\trs1\tA\tG\t.\t.\t.\t1"
    assert lines[3] == "2\t5\trs3\tG\tA\t.\t.\t.\t3"
    with pytest.raises(InputError, match="--force"):
        index_pvar(prefix)
    assert index_pvar(prefix, force=True)[1] == 3
