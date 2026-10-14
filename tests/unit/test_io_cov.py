"""Covariate encoding and selection (anqchor.io.cov)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fixtures" / "io"))
from anq_io_fixtures import tsv, write_text

from anqchor.io.cov import (
    compact_names,
    encode_covariates,
    expand_name,
    read_and_encode,
    select_columns,
    write_cov,
)
from anqchor.io.tables import InputError


def raw(**cols: list[object]) -> pd.DataFrame:
    n = len(next(iter(cols.values())))
    return pd.DataFrame(
        {k: [np.nan if v is None else str(v) for v in vals] for k, vals in cols.items()},
        index=[f"S{i}" for i in range(1, n + 1)],
        dtype=object,
    )


def test_two_level_categorical_is_one_two() -> None:
    enc = encode_covariates(raw(sex=["XY", "XX", "XY"], age=[30, 40, 50]))
    assert enc.columns == ["sex", "age"]
    assert enc.values["sex"].tolist() == [2.0, 1.0, 2.0]
    assert enc.categorical == {"sex": {"XX": 1, "XY": 2}}
    assert enc.groups == {"sex": ["sex"], "age": ["age"]}


def test_k_level_indicators_with_reference_level() -> None:
    enc = encode_covariates(raw(site=["b", "a", "c", "a"], x=[1, 2, 3, 5]))
    assert enc.columns == ["site_b", "site_c", "x"]
    assert enc.values["site_b"].tolist() == [2.0, 1.0, 1.0, 1.0]
    assert enc.values["site_c"].tolist() == [1.0, 1.0, 2.0, 1.0]
    assert enc.groups["site"] == ["site_b", "site_c"]
    assert enc.categorical["site"] == {"a": 1, "b": 2, "c": 3}


def test_forced_numeric_categorical_sorted_numerically() -> None:
    enc = encode_covariates(raw(batch=[10, 2, 1, 2], x=[1, 2, 3, 4]), categorical=["batch"])
    assert enc.categorical["batch"] == {"1": 1, "2": 2, "10": 3}
    assert enc.columns == ["batch_2", "batch_10", "x"]


def test_unknown_categorical_name() -> None:
    with pytest.raises(InputError, match="unknown covariate"):
        encode_covariates(raw(x=[1, 2]), categorical=["sex"])


def test_constant_columns_dropped() -> None:
    enc = encode_covariates(raw(c=[1, 1, 1], s=["A", "A", "A"], x=[1, 2, 3]))
    assert enc.columns == ["x"]
    assert sorted(enc.dropped_constant) == ["c", "s"]
    assert enc.groups["c"] == [] and enc.groups["s"] == []


def test_missing_error() -> None:
    with pytest.raises(InputError, match=r"missing covariate values.*'x' in sample 'S2'"):
        encode_covariates(raw(x=[1, None, 3], y=[1, 2, 4]))


def test_missing_drop_samples() -> None:
    enc = encode_covariates(
        raw(x=[1, None, 3, 4], sex=["M", "F", None, "F"]), missing_policy="drop_samples"
    )
    assert list(enc.values.index) == ["S1", "S4"]
    assert enc.n_missing_handled == 2 and enc.dropped_samples == ["S2", "S3"]


def test_missing_impute_mean() -> None:
    enc = encode_covariates(raw(x=[1, None, 3], sex=["M", "F", None]), missing_policy="impute_mean")
    assert enc.values["x"].tolist() == [1.0, 2.0, 3.0]
    assert enc.values["sex"].tolist() == [2.0, 1.0, 1.5]
    assert enc.n_missing_handled == 2


def test_read_and_encode_sample_col_layout(tmp_path: Path) -> None:
    p = write_text(
        tmp_path / "cov.tab.gz",
        tsv(
            [
                ["id", "S1", "S2", "S3"],
                ["sex", "XX", "XY", "XX"],
                ["PC1", 0.1, 0.2, -0.3],
                ["one", 1, 1, 1],
            ]
        ),
    )
    enc = read_and_encode(p, "tsv-sample-col")
    assert enc.columns == ["sex", "PC1"] and enc.dropped_constant == ["one"]
    write_cov(tmp_path / "o.tsv", enc)
    assert (tmp_path / "o.tsv").read_text().splitlines() == [
        "FID\tIID\tsex\tPC1",
        "S1\tS1\t1\t0.1",
        "S2\tS2\t2\t0.2",
        "S3\tS3\t1\t-0.3",
    ]


@pytest.mark.parametrize(
    ("token", "expected"),
    [
        ("PEER1-3", ["PEER1", "PEER2", "PEER3"]),
        ("PC1-PC2", ["PC1", "PC2"]),
        ("sex", ["sex"]),
        ("a-b", ["a-b"]),
    ],
)
def test_expand_name(token: str, expected: list[str]) -> None:
    assert expand_name(token, ["sex"]) == expected


def test_expand_name_existing_dash_name_kept() -> None:
    assert expand_name("PC1-5", ["PC1-5"]) == ["PC1-5"]


def test_select_columns_original_and_encoded_names() -> None:
    groups = {
        "sex": ["sex"],
        "site": ["site_b", "site_c"],
        "PC1": ["PC1"],
        "PC2": ["PC2"],
        "PC3": ["PC3"],
        "k": [],
    }
    cols = ["sex", "site_b", "site_c", "PC1", "PC2", "PC3"]
    assert select_columns(["PC1-2", "site", "sex", "PC1"], groups, cols) == (
        ["PC1", "PC2", "site_b", "site_c", "sex"],
        [],
    )
    assert select_columns(["site_c", "k"], groups, cols) == (["site_c"], ["k"])
    with pytest.raises(
        InputError,
        match=r"unknown covariate name\(s\) 'PC4', 'PC5'.*available: sex, site, PC1-3, k",
    ):
        select_columns(["PC1-5"], groups, cols)


def test_compact_names() -> None:
    names = ["sex", *[f"PC{i}" for i in range(1, 6)], *[f"PEER{i}" for i in range(1, 61)], "age"]
    assert compact_names(names) == "sex, PC1-5, PEER1-60, age"


def _col(values: list[object], name: str = "age") -> pd.DataFrame:
    return pd.DataFrame({name: values}, index=[f"s{i}" for i in range(len(values))])


@pytest.mark.parametrize("token", ["", "nan", "NaN", ".", "  "])
def test_numeric_column_with_odd_missing_token_is_an_error(token: str) -> None:
    """One empty/'nan'/'.' cell must not turn a continuous covariate into ~N indicators."""
    vals: list[object] = [f"{i / 10:.1f}" for i in range(50)] + [token]
    if token.strip() == "":
        # an empty cell is missing (cov_missing policy applies), never a level
        with pytest.raises(InputError, match="missing covariate values"):
            encode_covariates(_col(vals))
        enc = encode_covariates(_col(vals), missing_policy="drop_samples")
        assert enc.columns == ["age"] and enc.dropped_samples == ["s50"] and not enc.categorical
    else:
        with pytest.raises(InputError, match=r"non-numeric value .* missing_str"):
            encode_covariates(_col(vals))


def test_many_level_auto_categorical_needs_confirmation() -> None:
    vals = [f"id{i}" for i in range(60)]
    with pytest.raises(InputError, match="cov_categorical"):
        encode_covariates(_col(vals, "batch"))
    enc = encode_covariates(_col(vals, "batch"), categorical=["batch"])
    assert len(enc.columns) == 59


def test_infinite_covariate_rejected() -> None:
    with pytest.raises(InputError, match="infinite"):
        encode_covariates(_col(["1", "2", "inf", "3"]))
    with pytest.raises(InputError, match="infinite"):
        encode_covariates(_col(["1", "2", "1e400", "3"]))


def test_colliding_level_names_rejected() -> None:
    with pytest.raises(InputError, match="both encode"):
        encode_covariates(_col(["0", "a b", "a_b", "0", "a b"], "batch"), categorical=["batch"])
