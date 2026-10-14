"""Reference pairs, CPRA parsing and presence checks (anqchor.io.pairs)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fixtures" / "io"))
from anq_io_fixtures import make_geno, needs_htslib, tsv, write_text

from anqchor.io.genotypes import parse_pgen_list
from anqchor.io.pairs import (
    Variant,
    check_presence,
    chrom_order,
    parse_cpra,
    read_reference_pairs,
    read_variant_list,
    sort_variants,
)
from anqchor.io.tables import InputError


@pytest.mark.parametrize("bad", ["1:100:A", "1:x:A:G", "1:0:A:G", ":1:A:G", "1:1::G", "1:1:A:G:T"])
def test_parse_cpra_rejects(bad: str) -> None:
    assert parse_cpra(bad) is None


def test_parse_cpra() -> None:
    v = parse_cpra("chr22:50595411:C:GA")
    assert (
        v == Variant("chr22", 50595411, "C", "GA")
        and v.is_indel
        and v.cpra == "chr22:50595411:C:GA"
    )


def test_read_reference_pairs(tmp_path: Path) -> None:
    rows = [
        ["gene", "snp", "b", "s", "p", "other"],
        ["T1", "1:100:A:G", 0.5, 0.1, 1e-3, "x"],
        ["T1", "1:100:A:G", 0.6, 0.1, 1e-3, "dup"],
        ["T2", "1:200:C:T", "NA", 0.1, 0.5, "bad beta"],
        ["T3", "1-200-C-T", 0.1, 0.1, 0.5, "bad variant"],
        ["T4", "2:50:T:C", -0.2, "inf", 0.5, "bad se"],
        ["T5", "chr2:50:T:C", -0.2, 0.05, "NA", "ok"],
    ]
    p = write_text(tmp_path / "pairs.tsv.gz", tsv(rows))
    cols = {"trait": "gene", "variant": "snp", "beta": "b", "se": "s", "pval": "p"}
    rp = read_reference_pairs(p, cols, "chr")
    assert rp.table.values.tolist() == [
        ["T1", "chr1:100:A:G", "0.5", "0.1", "0.001"],
        ["T5", "chr2:50:T:C", "-0.2", "0.05", "NA"],
    ]
    assert (rp.n_input, rp.n_invalid, rp.n_duplicates, rp.n_renamed) == (6, 3, 1, 1)
    rp = read_reference_pairs(p, {**cols, "pval": None}, "plain")
    assert (
        rp.table["pval"].tolist() == ["NA", "NA"] and rp.table["variant"].tolist()[1] == "2:50:T:C"
    )
    with pytest.raises(InputError, match="missing column"):
        read_reference_pairs(p, {"trait": "trait"}, "plain")


def test_read_variant_list(tmp_path: Path) -> None:
    p = tmp_path / "v.txt"
    p.write_text("variant\n1:100:A:G\nchr1:200:C:T\n# c\nbad:1\n1:100:A:G\n")
    v, n_bad, n_ren = read_variant_list(p, "plain")
    assert list(v) == ["1:100:A:G", "1:200:C:T"] and n_bad == 1 and n_ren == 1


@needs_htslib
@pytest.mark.parametrize("style", ["plain", "chr"])
def test_presence_tabix_and_python_agree(tmp_path: Path, style: str) -> None:
    entries = parse_pgen_list(make_geno(tmp_path, style=style))
    c1, c2 = ("1", "2") if style == "plain" else ("chr1", "chr2")
    want = [
        Variant(c1, 100, "A", "G"),
        Variant(c1, 200, "C", "CA"),
        Variant(c1, 200, "C", "G"),
        Variant(c1, 100, "G", "A"),
        Variant(c2, 150, "AT", "A"),
        Variant(c2, 151, "T", "A"),
        Variant("3", 1, "A", "C"),
        Variant(c1, 300, "G", "A"),
    ]
    a, ma = check_presence(want, entries, "tabix")
    b, mb = check_presence(want, entries, "python")
    assert (ma, mb) == ("tabix", "python")
    assert a == b == {want[0], want[1], want[4], want[7]}


def test_sort_variants_follows_list_order() -> None:
    vs = [
        Variant("2", 5, "A", "C"),
        Variant("X", 1, "A", "C"),
        Variant("10", 9, "A", "C"),
        Variant("2", 1, "A", "C"),
    ]
    order = {"10": 0, "2": 1}
    assert [v.cpra for v in sort_variants(vs, order)] == [
        "10:9:A:C",
        "2:1:A:C",
        "2:5:A:C",
        "X:1:A:C",
    ]
    assert chrom_order([]) == {}


def test_lowercase_alleles_normalised() -> None:
    assert parse_cpra("21:100:a:g") == Variant("21", 100, "A", "G")
    assert parse_cpra("21:100:<del>:a") == Variant("21", 100, "<del>", "A")  # symbolic allele kept


def test_gzip_pairs_without_extension(tmp_path: Path) -> None:
    import gzip

    p = tmp_path / "pairs_noext.tsv"
    with gzip.open(p, "wt") as fh:
        fh.write(tsv([["trait", "variant", "beta", "se"], ["T1", "1:100:a:g", 0.5, 0.1]]))
    rp = read_reference_pairs(p, {}, "plain")
    assert rp.table["variant"].tolist() == ["1:100:A:G"]


def test_pval_column_read_without_mapping(tmp_path: Path) -> None:
    """'--pairs' documents 'trait variant beta se [pval]': a 'pval' column is read even
    when no column map names it (flag-only runs), and stays NA when absent."""
    p = write_text(
        tmp_path / "with.tsv",
        tsv([["trait", "variant", "beta", "se", "pval"], ["T1", "1:100:A:G", 0.5, 0.1, 1e-3]]),
    )
    assert read_reference_pairs(p, {}, "plain").table["pval"].tolist() == ["0.001"]
    assert read_reference_pairs(p, {"pval": None}, "plain").table["pval"].tolist() == ["0.001"]
    q = write_text(
        tmp_path / "without.tsv",
        tsv([["trait", "variant", "beta", "se"], ["T1", "1:100:A:G", 0.5, 0.1]]),
    )
    assert read_reference_pairs(q, {}, "plain").table["pval"].tolist() == ["NA"]
