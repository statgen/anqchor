"""qpgen output readers against excerpts of real qpgentools output
(``tests/fixtures/qpgen``), including the layouts of unpatched dev builds."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from anqchor.qpgen.outputs import (
    MATCH_STATUSES_ALL,
    MATCH_STATUSES_ASSIGNED,
    QpgenOutputError,
    RectChunk,
    _t_log10p_tail,
    iter_rect_assoc,
    parse_auto_lambda,
    read_match_all,
    read_match_assigned,
    read_match_weights,
    read_pair_assoc,
    read_prs,
    read_prs_traits,
    read_rect_header,
    t_log10p,
)

TRAITS5 = [
    "ENSG00000008735.14:MAPK8IP2",
    "ENSG00000015475.19:BID",
    "ENSG00000025708.14:TYMP",
    "ENSG00000025770.19:NCAPH2",
    "ENSG00000159200.18:RCAN1",
]


def _gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


# ------------------------------------------------------------------- pair-prs
def test_read_prs(qpgen_fixtures: Path) -> None:
    df = read_prs(qpgen_fixtures / "prs.prs.tsv.gz")
    assert df.index.name == "IID"
    assert df.shape == (8, 6)
    assert list(df.index[:2]) == ["HG00096", "HG00097"]
    assert df.columns[0] == "ENSG00000008735.14:MAPK8IP2"
    assert df.loc["HG00096", "ENSG00000015475.19:BID"] == pytest.approx(0.034317)
    assert all(pd.api.types.is_float_dtype(t) for t in df.dtypes)
    se = read_prs(qpgen_fixtures / "prs.se.tsv.gz")
    assert se.shape == (8, 6) and (se.to_numpy() >= 0).all()
    assert read_prs_traits(qpgen_fixtures / "prs.prs.tsv.gz") == list(df.columns)


def test_read_prs_bad_header(tmp_path: Path) -> None:
    p = _gz(tmp_path / "x.tsv.gz", "ID\tt1\nA\t1\n")
    with pytest.raises(QpgenOutputError, match="missing columns"):
        read_prs(p)
    with pytest.raises(QpgenOutputError, match="FID IID"):
        read_prs_traits(p)


# ------------------------------------------------------------ match-prs-pheno
def test_read_match_all(qpgen_fixtures: Path) -> None:
    df = read_match_all(qpgen_fixtures / "match.match.all.tsv.gz")
    assert len(df) == 8
    assert set(df["MatchStatus"]) <= set(MATCH_STATUSES_ALL)
    assert {"SELF_BEST", "SELF_LENIENT", "UNCLEAR"} <= set(df["MatchStatus"])
    unclear = df[df["MatchStatus"] == "UNCLEAR"].iloc[0]
    # sample IDs beginning with 'NA' stay strings, never missing values
    assert unclear["ID.Pheno"] == "NA11919" and unclear["ID.self"] == "NA11919"
    assert unclear["ID.1st"] == "HG00237"
    assert unclear["Rank.self"] == 400
    assert unclear["Z.1st"] == pytest.approx(4.884522)
    assert df["ID.2nd"].isna().sum() == 0
    assert pd.api.types.is_float_dtype(df["Z.self"])


def test_read_match_all_auto(qpgen_fixtures: Path) -> None:
    df = read_match_all(qpgen_fixtures / "matchauto.match.all.tsv.gz")
    row = df[df["MatchStatus"] == "SINGLE_NEW_BEST"].iloc[0]
    assert (row["ID.Pheno"], row["ID.1st"]) == ("NA11919", "HG00237")


def test_read_match_assigned_and_weights(qpgen_fixtures: Path) -> None:
    a = read_match_assigned(qpgen_fixtures / "match.match.assigned.tsv.gz")
    assert set(a["MatchStatus"]) <= set(MATCH_STATUSES_ASSIGNED)
    assert {"BEST_MATCH", "LENIENT_MATCH", "NO_MATCH"} <= set(a["MatchStatus"])
    w = read_match_weights(qpgen_fixtures / "match.weights.tsv.gz")
    assert list(w.columns) == ["Trait", "Weight", "N.Obs"]
    assert len(w) == 7 and (w["N.Obs"] == 729).all()
    assert w["Trait"].iloc[0] == "ENSG00000008735.14:MAPK8IP2"


def test_match_missing_columns(tmp_path: Path) -> None:
    p = _gz(tmp_path / "m.tsv.gz", "ID.Pheno\tN.Traits\nA\t1\n")
    with pytest.raises(QpgenOutputError, match="MatchStatus"):
        read_match_all(p)
    with pytest.raises(QpgenOutputError):
        read_match_weights(p)


def test_parse_auto_lambda(qpgen_fixtures: Path) -> None:
    text = (qpgen_fixtures / "matchauto.log").read_text()
    assert parse_auto_lambda(text) == pytest.approx(0.4549)
    assert parse_auto_lambda("no tuning here") is None
    two = "Using automatically tuned lambda = 0.1\nUsing automatically tuned lambda = 2.5e-01 after"
    assert parse_auto_lambda(two) == pytest.approx(0.25)


# ----------------------------------------------------------------- pair-assoc
def test_read_pair_assoc(qpgen_fixtures: Path) -> None:
    df = read_pair_assoc(qpgen_fixtures / "pair_assoc.tsv.gz")
    assert len(df) == 9
    assert df.columns[0] == "TRAIT"  # '#TRAIT' normalised
    assert df.attrs["recomputed_log10p"] is False
    assert str(df["N"].dtype) == "Int64" and (df["N"] == 729).all()
    assert (df["N_RR"] + df["N_RA"] + df["N_AA"] == df["N"]).all()
    assert "21:9357536:G:GTT" not in set(df["ID"])  # absent from the pvar
    row = df[(df["TRAIT"] == "ENSG00000280071.4:GATD3B")].iloc[0]
    assert row["CHROM"] == "21" and row["GENPOS"] == 5101672
    assert row["BETA"] == pytest.approx(-1.09876)
    assert row["LOG10P"] == pytest.approx(77.0227)
    assert row["LOG10P_RINT"] == pytest.approx(76.9814)
    # written LOG10P agrees with the value recomputed from TSTAT and N - 2
    recomputed = t_log10p(df["TSTAT"], df["N"].astype(float) - 2)
    assert np.allclose(recomputed, df["LOG10P"].to_numpy(), atol=2e-3, rtol=1e-4)


def test_read_pair_assoc_short_rows(qpgen_fixtures: Path) -> None:
    """Unpatched dev builds write 19 data fields under the 20-field header."""
    full = read_pair_assoc(qpgen_fixtures / "pair_assoc.tsv.gz")
    short = read_pair_assoc(qpgen_fixtures / "pair_assoc.short.tsv.gz")
    assert short.attrs["recomputed_log10p"] is True
    assert list(short.columns) == list(full.columns)
    for c in ("BETA", "SE", "TSTAT", "BETA_RINT", "LOG10P_RINT"):
        assert np.allclose(short[c], full[c])
    assert np.allclose(short["LOG10P"], full["LOG10P"], atol=2e-3, rtol=1e-4)


def test_read_pair_assoc_errors(tmp_path: Path, qpgen_fixtures: Path) -> None:
    with pytest.raises(QpgenOutputError, match="empty"):
        read_pair_assoc(_gz(tmp_path / "e.tsv.gz", ""))
    with pytest.raises(QpgenOutputError, match="missing columns"):
        read_pair_assoc(_gz(tmp_path / "h.tsv.gz", "#TRAIT\tCHROM\n"))
    lines = gzip.open(qpgen_fixtures / "pair_assoc.tsv.gz", "rt").read().splitlines()
    bad = lines[0] + "\n" + "\t".join(lines[1].split("\t")[:15]) + "\n"
    with pytest.raises(QpgenOutputError, match="row has 15 fields"):
        read_pair_assoc(_gz(tmp_path / "b.tsv.gz", bad))
    only_header = read_pair_assoc(_gz(tmp_path / "o.tsv.gz", lines[0] + "\n"))
    assert len(only_header) == 0 and "LOG10P" in only_header.columns


def test_t_log10p_extremes() -> None:
    out = t_log10p(np.array([0.0, 2.0, -2.0, 80.0]), np.array([100, 100, 100, 700]))
    assert out[0] == pytest.approx(0.0, abs=1e-12)
    assert out[1] == pytest.approx(out[2])
    assert out[1] == pytest.approx(-np.log10(0.04821), abs=1e-3)
    assert np.isfinite(out[3]) and out[3] > 300
    # the exact tail formula used past scipy's underflow agrees where both are finite
    t = np.array([1.0, 2.0, 10.0, 30.0, 37.0, 5.0])
    d = np.array([10.0, 100.0, 100.0, 727.0, 727.0, 3.0])
    assert np.allclose(_t_log10p_tail(t, d), t_log10p(t, d), rtol=1e-10)
    # monotone in |t| through the switch-over, scalar df broadcasts, inf/NaN propagate
    big = t_log10p(np.array([30.0, 40.0, 60.0, 80.0, 200.0]), 700)
    assert np.all(np.diff(big) > 0) and np.all(np.isfinite(big))
    edge = t_log10p(np.array([np.inf, np.nan]), 10)
    assert np.isinf(edge[0]) and np.isnan(edge[1])


# ----------------------------------------------------------------- rect-assoc
def test_read_rect_header(qpgen_fixtures: Path) -> None:
    meta, traits = read_rect_header(qpgen_fixtures / "rect_assoc.tsv.gz")
    assert meta[:3] == ["#CHROM", "GENPOS", "ID"] and "TEST" in meta
    assert traits == TRAITS5  # trait IDs containing '.' and ':' survive


def _collect(path: Path, **kw: object) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    chunks: list[RectChunk] = list(iter_rect_assoc(path, **kw))  # type: ignore[arg-type]
    meta = pd.concat([c.meta for c in chunks], ignore_index=True)
    stats = {k: np.vstack([c.stats[k] for c in chunks]) for k in chunks[0].stats}
    return meta, stats


def test_iter_rect_assoc(qpgen_fixtures: Path) -> None:
    meta, stats = _collect(
        qpgen_fixtures / "rect_assoc.tsv.gz", stats_wanted=("BETA", "SE", "LOG10P"), chunksize=4
    )
    assert len(meta) == 6
    assert list(meta.columns[:3]) == ["CHROM", "GENPOS", "ID"]
    assert "21:9357536:G:GTT" not in set(meta["ID"])
    assert (meta["N"] == 729).all()
    assert stats["BETA"].shape == (6, 5) and stats["LOG10P"].shape == (6, 5)
    assert set(stats) == {"BETA", "SE", "LOG10P"}
    first = meta.index[meta["ID"] == "21:34640365:G:A"][0]
    assert stats["BETA"][first, 0] == pytest.approx(0.0785813)
    assert stats["LOG10P"][first, 4] == pytest.approx(5.42269)  # RCAN1, its cis pair


def test_iter_rect_assoc_doubled_n(qpgen_fixtures: Path) -> None:
    """Unpatched builds report N as the allele number; N is recomputed from counts."""
    meta, stats = _collect(qpgen_fixtures / "rect_assoc.doubledN.tsv.gz")
    ref_meta, ref_stats = _collect(qpgen_fixtures / "rect_assoc.tsv.gz")
    assert (meta["N"] == 729).all()
    assert np.allclose(stats["BETA"], ref_stats["BETA"])
    pd.testing.assert_frame_equal(meta, ref_meta)


def test_rect_assoc_errors(tmp_path: Path) -> None:
    with pytest.raises(QpgenOutputError, match="header start"):
        read_rect_header(_gz(tmp_path / "a.tsv.gz", "CHROM\tGENPOS\tID\n"))
    hdr = "#CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tN_RR\tN_RA\tN_AA\tTEST"
    with pytest.raises(QpgenOutputError, match="do not match"):
        read_rect_header(_gz(tmp_path / "b.tsv.gz", hdr + "\tBETA.t1\tSE.t1\tTSTAT.t1\n"))
    with pytest.raises(QpgenOutputError, match="do not match"):
        read_rect_header(
            _gz(tmp_path / "c.tsv.gz", hdr + "\tBETA.t1\tSE.t2\tTSTAT.t1\tLOG10P.t1\n")
        )
    with pytest.raises(QpgenOutputError, match="missing columns"):
        read_rect_header(
            _gz(tmp_path / "d.tsv.gz", "#CHROM\tGENPOS\tID\tBETA.t\tSE.t\tTSTAT.t\tLOG10P.t\n")
        )
    ok = _gz(tmp_path / "e.tsv.gz", hdr + "\tBETA.t1\tSE.t1\tTSTAT.t1\tLOG10P.t1\n")
    with pytest.raises(ValueError, match="unknown rect-assoc statistic"):
        list(iter_rect_assoc(ok, stats_wanted=("P",)))
    assert list(iter_rect_assoc(ok)) == []  # header only: no chunks
