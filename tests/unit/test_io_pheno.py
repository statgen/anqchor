"""Phenotype matrix readers/writers and trait annotation (anqchor.io.pheno)."""

from __future__ import annotations

import filecmp
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fixtures" / "io"))
from anq_io_fixtures import tsv, write_text

from anqchor.io.pheno import (
    annotate_traits,
    convert_to_regenie,
    read_matrix,
    read_trait_annotation,
    write_regenie,
)
from anqchor.io.tables import InputError

VALUES = {"S1": [0.1, 1.5], "S2": [-2.25, 0.3333333333333333], "S3": [3.0, 1e-07]}
TRAITS = ["T1", "T2"]


def _expected() -> pd.DataFrame:
    return pd.DataFrame.from_dict(VALUES, orient="index", columns=TRAITS)


def _write(tmp: Path, fmt: str, header0: str = "#chr", name: str = "m.tsv.gz") -> Path:
    s = list(VALUES)
    if fmt == "regenie":
        rows = [["FID", "IID", *TRAITS]] + [[x, x, *VALUES[x]] for x in s]
    elif fmt == "tsv-sample-row":
        rows = [["ID", *TRAITS]] + [[x, *VALUES[x]] for x in s]
    elif fmt == "tsv-sample-col":
        rows = [["ID", *s]] + [[t, *(VALUES[x][j] for x in s)] for j, t in enumerate(TRAITS)]
    else:
        rows = [[header0, "start", "end", "phenotype_id", *s]] + [
            ["chr1", 100 * (j + 1), 100 * (j + 1) + 1, t, *(VALUES[x][j] for x in s)]
            for j, t in enumerate(TRAITS)
        ]
    return write_text(tmp / name, tsv(rows))


@pytest.mark.parametrize("fmt", ["regenie", "tensorqtl", "tsv-sample-row", "tsv-sample-col"])
def test_every_format_reads_to_samples_by_traits(tmp_path: Path, fmt: str) -> None:
    m = read_matrix(_write(tmp_path, fmt), fmt)
    assert m.samples == ["S1", "S2", "S3"]
    assert m.features == TRAITS
    pd.testing.assert_frame_equal(m.values, _expected(), check_names=False)
    if fmt == "tensorqtl":
        assert m.loci is not None and m.loci["start"].tolist() == [100, 200]


@pytest.mark.parametrize("h0", ["#chr", "#Chr", "chr", "CHROM", "#chrom"])
def test_tensorqtl_header_variants(tmp_path: Path, h0: str) -> None:
    m = read_matrix(_write(tmp_path, "tensorqtl", header0=h0), "tensorqtl")
    assert m.features == TRAITS


def test_tensorqtl_bad_header(tmp_path: Path) -> None:
    p = write_text(
        tmp_path / "x.bed", tsv([["#chr", "beg", "end", "id", "S1"], ["chr1", 1, 2, "T", 0.5]])
    )
    with pytest.raises(InputError, match="tensorqtl"):
        read_matrix(p, "tensorqtl")


def test_regenie_header_checked(tmp_path: Path) -> None:
    p = _write(tmp_path, "tsv-sample-row")
    with pytest.raises(InputError, match="FID IID"):
        read_matrix(p, "regenie")


def test_missing_strings(tmp_path: Path) -> None:
    p = write_text(tmp_path / "m.tsv", tsv([["ID", "T1", "T2"], ["S1", ".", 1], ["S2", "NA", 2]]))
    m = read_matrix(p, "tsv-sample-row", missing=[".", "NA"])
    assert m.n_missing() == 2
    with pytest.raises(InputError, match=r"'\.' for trait 'T1' in sample 'S1'"):
        read_matrix(p, "tsv-sample-row", missing=["NA"])


@pytest.mark.parametrize("fmt", ["tsv-sample-row", "tsv-sample-col"])
def test_non_numeric_names_trait_and_sample(tmp_path: Path, fmt: str) -> None:
    if fmt == "tsv-sample-row":
        rows = [["ID", "T1", "T2"], ["S1", 1, 2], ["S2", 3, "oops"], ["S3", "bad", 1]]
    else:
        rows = [["ID", "S1", "S2", "S3"], ["T1", 1, 3, "bad"], ["T2", 2, "oops", 1]]
    p = write_text(tmp_path / "m.tsv", tsv(rows))
    with pytest.raises(InputError) as e:
        read_matrix(p, fmt)
    msg = str(e.value)
    # first offending cell in file order
    if fmt == "tsv-sample-row":
        assert "'oops' for trait 'T2' in sample 'S2'" in msg
    else:
        assert "'bad' for trait 'T1' in sample 'S3'" in msg


def test_duplicate_ids(tmp_path: Path) -> None:
    p = write_text(
        tmp_path / "a.tsv", tsv([["FID", "IID", "T1"], ["S1", "S1", 1], ["S1", "S1", 2]])
    )
    with pytest.raises(InputError, match="duplicate sample"):
        read_matrix(p, "regenie")
    p = write_text(tmp_path / "b.tsv", tsv([["FID", "IID", "T1", "T1"], ["S1", "S1", 1, 2]]))
    with pytest.raises(InputError, match="duplicate"):
        read_matrix(p, "regenie")
    p = write_text(tmp_path / "c.tsv", tsv([["ID", "S1", "S1"], ["T1", 1, 2]]))
    with pytest.raises(InputError, match="duplicate sample"):
        read_matrix(p, "tsv-sample-col")


def test_fast_path_copies_bytes(tmp_path: Path) -> None:
    src = _write(tmp_path, "regenie")
    res = convert_to_regenie(src, "regenie", tmp_path / "out.tsv.gz")
    assert res.fast_path
    assert filecmp.cmp(src, tmp_path / "out.tsv.gz", shallow=False)
    # plain input is compressed, content unchanged
    plain = tmp_path / "plain.tsv"
    plain.write_text(gzip.decompress(src.read_bytes()).decode())
    res = convert_to_regenie(plain, "regenie", tmp_path / "out2.tsv.gz")
    assert res.fast_path
    assert gzip.decompress((tmp_path / "out2.tsv.gz").read_bytes()).decode() == plain.read_text()


def test_no_fast_path_when_fid_differs_or_custom_missing(tmp_path: Path) -> None:
    p = write_text(
        tmp_path / "a.tsv", tsv([["FID", "IID", "T1"], ["F1", "S1", 1], ["F2", "S2", 2]])
    )
    res = convert_to_regenie(p, "regenie", tmp_path / "o.tsv.gz")
    assert not res.fast_path
    out = pd.read_csv(tmp_path / "o.tsv.gz", sep="\t", dtype=str)
    assert out["FID"].tolist() == ["S1", "S2"]
    p = write_text(
        tmp_path / "b.tsv", tsv([["FID", "IID", "T1"], ["S1", "S1", "-9"], ["S2", "S2", 2]])
    )
    res = convert_to_regenie(p, "regenie", tmp_path / "o2.tsv.gz", missing=["-9"])
    assert not res.fast_path and res.n_missing == 1
    assert "S1\tS1\tNA" in gzip.decompress((tmp_path / "o2.tsv.gz").read_bytes()).decode()


def test_conversion_preserves_doubles(tmp_path: Path) -> None:
    rng = np.random.default_rng(1)
    vals = rng.normal(size=(3, 50))
    rows = [["ID", "S1", "S2", "S3"]] + [
        [f"T{j}", *(repr(float(v)) for v in vals[:, j])] for j in range(50)
    ]
    p = write_text(tmp_path / "m.tsv", tsv(rows))
    convert_to_regenie(p, "tsv-sample-col", tmp_path / "o.tsv.gz")
    back = pd.read_csv(
        tmp_path / "o.tsv.gz", sep="\t", float_precision="round_trip", index_col=1
    ).drop(columns="FID")
    assert np.array_equal(back.to_numpy(), vals)


def test_write_regenie_strings_and_nan(tmp_path: Path) -> None:
    df = pd.DataFrame({"a": [1.0, np.nan], "b": [2.5, 3.0]}, index=["x", "y"])
    write_regenie(tmp_path / "o.tsv", df, trim_int=True)
    assert (tmp_path / "o.tsv").read_text() == "FID\tIID\ta\tb\nx\tx\t1\t2.5\ny\ty\tNA\t3\n"


def test_trait_annotation_layouts(tmp_path: Path) -> None:
    bed = write_text(
        tmp_path / "a.bed.gz",
        tsv(
            [
                ["#Chr", "start", "end", "ID", "S1"],
                ["chr1", 10, 11, "G1.5", 0.1],
                ["chr2", 20, 21, "G2.1", 0.2],
            ]
        ),
    )
    a = read_trait_annotation(bed)
    assert a.to_dict("list") == {
        "trait": ["G1.5", "G2.1"],
        "chrom": ["chr1", "chr2"],
        "start": [10, 20],
        "end": [11, 21],
    }
    nohdr = write_text(tmp_path / "b.bed", tsv([["1", 10, 11, "G1.5"], ["1", 30, 31, "G1.5"]]))
    b = read_trait_annotation(nohdr)
    assert b["trait"].tolist() == ["G1.5"] and b["start"].tolist() == [10]  # first row kept
    tcse = write_text(
        tmp_path / "c.tsv", tsv([["trait", "chrom", "start", "end"], ["G1.5", "1", 10, 11]])
    )
    assert read_trait_annotation(tcse)["chrom"].tolist() == ["1"]
    tcse2 = write_text(tmp_path / "d.tsv", tsv([["G1.5", "1", 10, 11]]))
    assert read_trait_annotation(tcse2)["trait"].tolist() == ["G1.5"]


def test_annotate_traits_rules_and_style(tmp_path: Path) -> None:
    annot = pd.DataFrame(
        {"trait": ["G1.5", "G2.1"], "chrom": ["chr1", "chr2"], "start": [10, 20], "end": [11, 21]}
    )
    res = annotate_traits(["G1.5:ABC", "G2.1:DEF", "G3.1:XYZ"], annot, "plain", "as_is")
    assert res.rule == "before_colon" and res.n_with_loci == 2
    assert res.table["chrom"].tolist()[:2] == ["1", "2"] and pd.isna(res.table["chrom"].iloc[2])
    res = annotate_traits(["G1.5"], annot, "chr")
    assert res.rule == "as_is" and res.table["chrom"].tolist() == ["chr1"]
    res = annotate_traits(["Q"], annot, "chr")
    assert res.n_with_loci == 0


# ------------------------------------------------- input robustness (review findings)
def _gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(text)
    return path


@pytest.mark.parametrize("value", ["inf", "-inf", "1e400"])
def test_infinite_phenotype_rejected(tmp_path: Path, value: str) -> None:
    p = write_text(tmp_path / "p.tsv", f"FID\tIID\tg1\n1\t1\t{value}\n2\t2\t3\n")
    with pytest.raises(InputError, match=r"infinite value .* 'g1' in sample '1'"):
        convert_to_regenie(p, "regenie", tmp_path / "out.tsv.gz")


@pytest.mark.parametrize("name", ["gz_named.tsv", "gz_named.bgz", "plain_named.tsv.gz"])
def test_compression_detected_from_content(tmp_path: Path, name: str) -> None:
    text = "FID\tIID\tg1\n1\t1\t2\n2\t2\t3\n"
    p = tmp_path / name
    if name.startswith("gz_"):
        _gz(p, text)
    else:
        p.write_text(text)
    res = convert_to_regenie(p, "regenie", tmp_path / "out.tsv.gz")
    assert res.matrix.samples == ["1", "2"]
    ann = tmp_path / ("ann_" + name)
    if name.startswith("gz_"):
        _gz(ann, "#chr\tstart\tend\tid\nchr1\t1\t2\tg1\n")
    else:
        ann.write_text("#chr\tstart\tend\tid\nchr1\t1\t2\tg1\n")
    assert read_trait_annotation(ann)["trait"].tolist() == ["g1"]


def test_whitespace_or_quoted_ids_rejected(tmp_path: Path) -> None:
    p = write_text(tmp_path / "p.tsv", "FID\tIID\tg1\nA B\tA B\t2\nC\tC\t3\n")
    with pytest.raises(InputError, match="whitespace"):
        convert_to_regenie(p, "regenie", tmp_path / "o.tsv.gz")
    q = write_text(tmp_path / "q.tsv", 'FID\tIID\tg1\n"1"\t"1"\t2\n2\t2\t3\n')
    with pytest.raises(InputError, match="quote"):
        convert_to_regenie(q, "regenie", tmp_path / "o.tsv.gz")
    c = write_text(tmp_path / "c.tsv", "ID\tS1\tS2\nage at draw\t1\t2\n")
    with pytest.raises(InputError, match="whitespace"):
        read_matrix(c, "tsv-sample-col", numeric=False, what="covariate")


def test_tensorqtl_annotation_column_message(tmp_path: Path) -> None:
    p = write_text(
        tmp_path / "e.bed",
        "#Chr\tstart\tend\tphenotype_id\tgene_name\tstrand\t001\t002\n"
        "chr1\t1\t2\tg1\tA\t+\t0.1\t0.2\n"
        "chr1\t3\t4\tg2\tB\t-\t0.3\t0.4\n",
    )
    with pytest.raises(InputError, match="annotation column"):
        read_matrix(p, "tensorqtl")


def test_headerless_annotation_with_numeric_trait_ids(tmp_path: Path) -> None:
    bed = write_text(tmp_path / "a.bed", "chr1\t100\t200\t1001\nchr2\t300\t400\t1002\n")
    a = read_trait_annotation(bed)
    assert a.to_dict("list") == {
        "trait": ["1001", "1002"],
        "chrom": ["chr1", "chr2"],
        "start": [100, 300],
        "end": [200, 400],
    }
    tcse = write_text(tmp_path / "t.tsv", "1001\t1\t100\t200\n")
    assert read_trait_annotation(tcse).to_dict("list") == {
        "trait": ["1001"],
        "chrom": ["1"],
        "start": [100],
        "end": [200],
    }
    amb = write_text(tmp_path / "x.tsv", "1001\t2001\t100\t200\n")
    with pytest.raises(InputError, match="add a header"):
        read_trait_annotation(amb)
