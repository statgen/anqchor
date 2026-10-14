"""strategy_inputs step: covariate selection, file covariates, phenotype override."""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fixtures" / "io"))
from anq_io_fixtures import step_ctx, tsv, write_text

from anqchor.io.cov import categorical_levels, detect_categorical, encode_covariates, write_cov_raw
from anqchor.pipelines.layout import Layout as L
from anqchor.steps import StepError
from anqchor.steps.strategy import strategy_inputs
from anqchor.util.files import write_json

SAMPLES = ["S1", "S2", "S3", "S4"]


def make_dataset(run: Path) -> None:
    """Dataset 'd' with covariates sex, site (3 levels), PC1-3, PEER1-15."""
    import pandas as pd

    raw = {"sex": ["XX", "XY", "XY", "XX"], "site": ["a", "b", "c", "a"]}
    for i in range(1, 4):
        raw[f"PC{i}"] = [str(0.1 * i * k) for k in range(1, 5)]
    for i in range(1, 16):
        raw[f"PEER{i}"] = [str(i + 0.5 * k * k) for k in range(1, 5)]
    enc = encode_covariates(pd.DataFrame(raw, index=SAMPLES, dtype=object))
    assert enc.raw is not None
    write_cov_raw(run / L.ds_cov("d"), enc.raw)
    write_json(run / L.ds_cov_json("d"), enc.as_json())
    (run / L.ds_samples("d")).write_text("\n".join(SAMPLES) + "\n")
    write_json(run / L.ds_json("d"), {"dataset": "d"})


def params(cov: Any, pheno: str | None = None, fmt: str = "regenie") -> dict[str, Any]:
    return {
        "pipeline": "cis",
        "dataset": "d",
        "strategy": {
            "name": "s",
            "cov": cov,
            "pheno": pheno,
            "pheno_format": fmt,
            "transforms": ["raw", "rint"],
        },
        "dataset_cov_categorical": ["sex"],
        "dataset_cov_missing": "error",
        "missing_str": ["NA"],
    }


def outputs(cov: bool, pheno: bool = False) -> list[Path]:
    out = [L.strategy_json("cis", "d", "s")]
    if cov:
        out.append(L.strategy_cov("cis", "d", "s"))
    if pheno:
        out.append(L.strategy_pheno("cis", "d", "s"))
    return out


def header(p: Path) -> list[str]:
    return gzip.decompress(p.read_bytes()).decode().splitlines()[0].split("\t")


def sj(run: Path) -> dict[str, Any]:
    return json.loads((run / L.strategy_json("cis", "d", "s")).read_text())


@pytest.fixture
def run(tmp_path: Path) -> Path:
    r = tmp_path / "run"
    make_dataset(r)
    return r


def test_default_copies_dataset_cov(run: Path) -> None:
    strategy_inputs(step_ctx(run, "x", params("default"), outputs(True)))
    # the strategy gets the covariates with categorical levels as text (qpgen --cat-cov-names)
    assert (run / L.strategy_cov("cis", "d", "s")).read_bytes() == (
        run / L.ds_cov("d")
    ).read_bytes()
    assert header(run / L.strategy_cov("cis", "d", "s"))[2:4] == ["sex", "site"]
    s = sj(run)
    # n_cov counts encoded columns: sex 1, site 2 indicators, PC 3, PEER 15
    assert (
        s["n_cov"] == 1 + 2 + 3 + 15
        and s["pheno"] == str(L.ds_pheno("d"))
        and s["transforms"] == ["raw", "rint"]
    )
    assert s["cat_cov_names"] == ["sex", "site"]


def test_include_with_ranges_and_groups(run: Path) -> None:
    strategy_inputs(
        step_ctx(
            run,
            "x",
            params({"include": ["PEER1-3", "site", "sex", "PC1-2", "PEER2"]}),
            outputs(True),
        )
    )
    cols = ["PEER1", "PEER2", "PEER3", "site", "sex", "PC1", "PC2"]
    assert header(run / L.strategy_cov("cis", "d", "s")) == ["FID", "IID", *cols]
    assert sj(run)["cov_columns"] == cols
    assert sj(run)["n_cov"] == 3 + 2 + 1 + 2  # encoded column count
    assert sj(run)["cat_cov_names"] == ["sex", "site"]  # file order
    lines = (
        gzip.decompress((run / L.strategy_cov("cis", "d", "s")).read_bytes()).decode().splitlines()
    )
    row = lines[1].split("\t")
    assert row[:2] == ["S1", "S1"] and row[5] == "a" and row[6] == "XX"  # text levels kept


def test_include_unknown_name(run: Path) -> None:
    with pytest.raises(
        StepError,
        match=r"unknown covariate name\(s\) 'PEER16', 'PEER17'; available: sex, site, PC1-3, PEER1-15$",
    ):
        strategy_inputs(step_ctx(run, "x", params({"include": ["PEER1-17"]}), outputs(True)))
    with pytest.raises(StepError, match="'age'"):
        strategy_inputs(step_ctx(run, "x", params({"include": ["age"]}), outputs(True)))


def test_null_cov(run: Path) -> None:
    strategy_inputs(step_ctx(run, "x", params(None), outputs(False)))
    assert not (run / L.strategy_cov("cis", "d", "s")).exists()
    assert sj(run)["n_cov"] == 0 and sj(run)["cov_columns"] == []


def test_cov_file_encoded(run: Path, tmp_path: Path) -> None:
    f = write_text(
        tmp_path / "c.tsv",
        tsv([["ID", "sex", "age"], ["S1", "M", 30], ["S2", "F", 40], ["S9", "F", 50]]),
    )
    strategy_inputs(
        step_ctx(run, "x", params({"file": str(f), "format": "tsv-sample-row"}), outputs(True))
    )
    lines = (
        gzip.decompress((run / L.strategy_cov("cis", "d", "s")).read_bytes()).decode().splitlines()
    )
    assert lines == ["FID\tIID\tsex\tage", "S1\tS1\tM\t30", "S2\tS2\tF\t40", "S9\tS9\tF\t50"]
    assert sj(run)["cov_source"] == "file" and sj(run)["cat_cov_names"] == ["sex"]
    g = write_text(tmp_path / "g.tsv", tsv([["ID", "age"], ["X1", 30], ["X2", 40]]))
    with pytest.raises(StepError, match="no sample in common"):
        strategy_inputs(
            step_ctx(run, "x", params({"file": str(g), "format": "tsv-sample-row"}), outputs(True))
        )


def test_pheno_override(run: Path, tmp_path: Path) -> None:
    p = write_text(tmp_path / "p.tsv", tsv([["ID", "S1", "S2"], ["T1", 1.5, 2.5]]))
    strategy_inputs(
        step_ctx(run, "x", params(None, str(p), "tsv-sample-col"), outputs(False, True))
    )
    out = run / L.strategy_pheno("cis", "d", "s")
    assert gzip.decompress(out.read_bytes()).decode() == "FID\tIID\tT1\nS1\tS1\t1.5\nS2\tS2\t2.5\n"
    s = sj(run)
    assert s["pheno"] == str(L.strategy_pheno("cis", "d", "s")) and s["n_traits"] == 1


def test_detect_categorical(tmp_path: Path) -> None:
    f = write_text(
        tmp_path / "c.tsv",
        tsv(
            [
                ["FID", "IID", "sex", "age", "batch"],
                ["S1", "S1", "XX", 30, 1],
                ["S2", "S2", "XY", "NA", 2],
            ]
        ),
    )
    assert detect_categorical(f, "regenie", ["NA"]) == ["sex"]
    assert detect_categorical(f, "regenie", ["NA"], ["batch"]) == ["sex", "batch"]
    assert categorical_levels(f, "regenie", ["NA"], ["batch"]) == {"sex": 2, "batch": 2}
