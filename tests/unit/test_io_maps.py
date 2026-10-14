"""Trait-ID rules, trait maps and sample maps (anqchor.io.maps)."""

from __future__ import annotations

from pathlib import Path

import pytest

from anqchor.io.maps import (
    best_rule,
    filter_sample_map,
    map_traits_by_file,
    map_traits_by_rule,
    read_sample_map,
    read_trait_map,
    trait_key,
)
from anqchor.io.tables import InputError


@pytest.mark.parametrize(
    ("rule", "expected"),
    [
        ("as_is", "ENSG01.5:WASH7P"),
        ("before_colon", "ENSG01.5"),
        ("strip_version", "ENSG01.5:WASH7P"),
        ("before_colon_strip_version", "ENSG01"),
    ],
)
def test_trait_key_rules(rule: str, expected: str) -> None:
    assert trait_key("ENSG01.5:WASH7P", rule) == expected


def test_strip_version_only_trailing_digits() -> None:
    assert trait_key("ENSG01.12", "strip_version") == "ENSG01"
    assert trait_key("ENSG01.x", "strip_version") == "ENSG01.x"
    with pytest.raises(InputError, match="unknown trait_id_rule"):
        trait_key("x", "nope")


def test_map_by_rule() -> None:
    res = map_traits_by_rule(
        ["G1.5:A", "G2.1:B", "G9.1:Z"], ["G1.5", "G2.1", "G3.1"], "before_colon"
    )
    assert res.mapping == {"G1.5:A": "G1.5", "G2.1:B": "G2.1"}
    assert res.unmapped == ["G9.1:Z"]


def test_ambiguous_dataset_key() -> None:
    with pytest.raises(
        InputError, match=r"match several dataset traits.*'G1.5:A' -> 'G1.5', 'G1.6'"
    ):
        map_traits_by_rule(["G1.5:A"], ["G1.5", "G1.6"], "before_colon_strip_version")


def test_several_reference_traits_onto_one() -> None:
    with pytest.raises(InputError, match="receive several reference traits"):
        map_traits_by_rule(["G1.5", "G1.6"], ["G1"], "strip_version")


def test_best_rule() -> None:
    assert best_rule(["G1.5:A", "G2.1:B"], ["G1.5", "G2.1"]) == ("before_colon", 2)
    assert best_rule(["x"], ["y"])[1] == 0


def test_trait_map_file(tmp_path: Path) -> None:
    p = tmp_path / "tm.tsv"
    p.write_text("ref_trait\tpheno_trait\nR1\tP1\nR2\tP2\nR3\tPX\n")
    rows = read_trait_map(p)
    assert rows == [("R1", "P1"), ("R2", "P2"), ("R3", "PX")]
    res = map_traits_by_file(["R1", "R2", "R3"], ["P1", "P2"], rows)
    assert res.mapping == {"R1": "P1", "R2": "P2"} and res.unmapped == ["R3"] and res.rule == "file"
    with pytest.raises(InputError, match="map to several dataset traits"):
        map_traits_by_file(["R1"], ["P1", "P2"], [("R1", "P1"), ("R1", "P2")])
    with pytest.raises(InputError, match="receive several reference traits"):
        map_traits_by_file(["R1", "R2"], ["P1"], [("R1", "P1"), ("R2", "P1")])


@pytest.mark.parametrize("header", ["geno_id\tpheno_id", "NWD_ID\tSAMPLE_ID", "pheno_id\tgeno_id"])
def test_read_sample_map_headers(tmp_path: Path, header: str) -> None:
    p = tmp_path / "m.tsv"
    flip = header.startswith("pheno")
    p.write_text(
        header + "\n" + ("P1\tG1\n" if flip else "G1\tP1\n") + ("P1\tG1\n" if flip else "G1\tP1\n")
    )
    assert read_sample_map(p) == [("G1", "P1")]


def test_read_sample_map_needs_header(tmp_path: Path) -> None:
    p = tmp_path / "m.tsv"
    p.write_text("G1\tP1\n")
    with pytest.raises(InputError, match="geno_id pheno_id"):
        read_sample_map(p)


def test_filter_sample_map() -> None:
    rows = [("G1", "P1"), ("G2", "P2"), ("GX", "P3"), ("G3", "PX"), ("G1", "P4")]
    res = filter_sample_map(rows, ["G1", "G2", "G3"], ["P4", "P1", "P2", "P3"])
    assert res.pairs == [("G1", "P4"), ("G1", "P1"), ("G2", "P2")]
    assert (res.n_rows_geno_missing, res.n_rows_pheno_missing, res.n_geno_with_several_pheno) == (
        1,
        1,
        1,
    )


def test_pheno_sample_mapped_to_two_genotypes() -> None:
    with pytest.raises(InputError, match=r"several genotyped samples.*'P1' -> 'G1', 'G2'"):
        filter_sample_map([("G1", "P1"), ("G2", "P1")], ["G1", "G2"], ["P1"])
    # a second genotype that is not genotyped is not a conflict
    assert filter_sample_map([("G1", "P1"), ("G9", "P1")], ["G1"], ["P1"]).pairs == [("G1", "P1")]
