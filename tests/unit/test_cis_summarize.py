"""summarize_cis: replication levels, BH, pi1, controls, statuses, recommendation."""

from __future__ import annotations

import gzip
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from anqchor.pipelines.cis import SUMMARY_FILES
from anqchor.pipelines.layout import Layout as L
from anqchor.qpgen.outputs import read_pair_assoc
from anqchor.steps import StepContext, StepError
from anqchor.steps.cis import (
    SCHEMA_VERSION,
    ControlStats,
    active_levels,
    add_levels,
    bh_neglog10q,
    lambda_gc_from_median_log10p,
    pair_assoc_stats,
    power_curve,
    read_mapped_pairs,
    recommend,
    rect_stats,
    storey_pi1,
    summarize_cis,
    trait_table,
)
from anqchor.util.files import read_json, write_json

LEVELS = {"nominal": 0.05, "strict": 1e-5, "bonferroni": True, "fdr": 0.05}
PAIR_HEADER = [
    "#TRAIT",
    "CHROM",
    "GENPOS",
    "ID",
    "ALLELE0",
    "ALLELE1",
    "A1FREQ",
    "N",
    "N_RR",
    "N_RA",
    "N_AA",
    "INFO",
    "BETA",
    "SE",
    "TSTAT",
    "LOG10P",
    "BETA_RINT",
    "SE_RINT",
    "TSTAT_RINT",
    "LOG10P_RINT",
]


# ------------------------------------------------------------------ builders
def ref_pairs() -> list[dict[str, Any]]:
    """Six mapped pairs: T1..T5 with present variants, T6 with an absent indel."""
    return [
        {"trait": "T1", "variant": "1:100:A:G", "beta": 0.5, "se": 0.05, "pval": 1e-20},
        {"trait": "T2", "variant": "1:200:C:T", "beta": -0.4, "se": 0.05, "pval": 1e-12},
        {"trait": "T3", "variant": "1:300:G:A", "beta": 0.3, "se": 0.05, "pval": 1e-8},
        {"trait": "T4", "variant": "1:400:T:C", "beta": 0.2, "se": 0.05, "pval": 1e-5},
        {"trait": "T5", "variant": "2:500:A:C", "beta": 0.1, "se": 0.05, "pval": 0.01},
        {"trait": "T6", "variant": "2:600:A:AT", "beta": 0.3, "se": 0.05, "pval": 1e-6},
    ]


def assoc_row(
    trait: str,
    variant: str,
    beta: float,
    log10p: float,
    rint_beta: float | None = None,
    rint_log10p: float | None = None,
) -> list[str]:
    c, pos, a0, a1 = variant.split(":")
    rb = beta if rint_beta is None else rint_beta
    rl = log10p if rint_log10p is None else rint_log10p
    t = math.copysign(3.0, beta)
    return [
        trait,
        c,
        pos,
        variant,
        a0,
        a1,
        "0.3",
        "100",
        "50",
        "40",
        "10",
        "1",
        f"{beta}",
        "0.05",
        f"{t}",
        f"{log10p}",
        f"{rb}",
        "0.05",
        f"{t}",
        f"{rl}",
    ]


def write_pair_assoc(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "\t".join(PAIR_HEADER) + "\n" + "".join("\t".join(r) + "\n" for r in rows)
    path.write_bytes(gzip.compress(text.encode()))


STRONG = [  # strategy "full": every present pair tested, all concordant
    ("T1", "1:100:A:G", 0.45, 30.0),
    ("T2", "1:200:C:T", -0.35, 15.0),
    ("T3", "1:300:G:A", 0.25, 8.0),
    ("T4", "1:400:T:C", 0.15, 4.5),
    ("T5", "2:500:A:C", 0.05, 1.5),
]
WEAK = [  # strategy "none": T5 discordant, T4 not tested (no row)
    ("T1", "1:100:A:G", 0.40, 12.0),
    ("T2", "1:200:C:T", -0.30, 6.0),
    ("T3", "1:300:G:A", 0.20, 3.0),
    ("T5", "2:500:A:C", -0.05, 2.0),
]


def make_run(
    run: Path,
    engine: str = "pair",
    weak: list[tuple[str, str, float, float]] | None = None,
    n_ref_pairs: int = 7,
    rint: bool = False,
) -> dict[str, Any]:
    """A run directory with the prepare outputs summarize_cis reads and two strategies'
    pair-assoc outputs. Returns the step params."""
    pairs = pd.DataFrame(ref_pairs())
    pairs["ref_trait"] = "ref_" + pairs["trait"]
    pairs["present"] = [1, 1, 1, 1, 1, 0]
    p = run / L.dsref_pairs("d", "r")
    p.parent.mkdir(parents=True, exist_ok=True)
    pairs[["trait", "variant", "beta", "se", "pval", "ref_trait", "present"]].to_csv(
        p, sep="\t", index=False
    )
    write_json(
        run / L.dsref_json("d", "r"),
        {
            "n_pairs": n_ref_pairs,
            "n_pairs_mapped": 6,
            "n_unmapped_traits": n_ref_pairs - 6,
            "rule": "as_is",
        },
    )
    v = run / L.ref_variants("r")
    v.parent.mkdir(parents=True, exist_ok=True)
    v.write_text(
        "variant\tchrom\tpos\tref\talt\tis_indel\tpresent\n"
        + "".join(
            f"{r['variant']}\t{r['variant'].split(':')[0]}\t0\tA\tG\t"
            f"{int(r['variant'].endswith('AT'))}\t1\n"
            for r in ref_pairs()
        )
    )
    strategies = []
    for name, rows, ncov in (("full", STRONG, 10), ("none", WEAK if weak is None else weak, 0)):
        write_json(run / L.strategy_json("cis", "d", name), {"n_cov": ncov})
        f = L.cis_pair_assoc("d", name)
        write_pair_assoc(
            run / f, [assoc_row(t, var, b, lp, rint_log10p=lp * 0.9) for t, var, b, lp in rows]
        )
        strategies.append(
            {
                "name": name,
                "label": name,
                "file": str(f),
                "strategy_json": str(L.strategy_json("cis", "d", name)),
            }
        )
    return {
        "dataset": "d",
        "reference": "r",
        "engine": engine,
        "rint": rint,
        "strategies": strategies,
        "levels": dict(LEVELS),
        "max_unresolved_fraction": 0.5,
        "power_step": 0.5,
        "transform_semantics": {"as_given": "x"},
    }


def run_summarize(run: Path, params: dict[str, Any]) -> None:
    summarize_cis(
        StepContext(
            name="cis.d.summarize",
            run_dir=run,
            params=params,
            inputs=[],
            outputs=[L.cis_summary("d", n) for n in SUMMARY_FILES],
        )
    )


# ---------------------------------------------------------------- pure helpers
def test_active_levels_order_and_disable() -> None:
    assert active_levels(LEVELS) == ["sign", "nominal", "fdr", "bonferroni", "strict"]
    assert active_levels({"nominal": None, "strict": 1e-5, "bonferroni": False, "fdr": None}) == [
        "sign",
        "strict",
    ]


def test_bh_matches_reference_implementation() -> None:
    rng = np.random.default_rng(1)
    p = np.concatenate([rng.uniform(size=200), 10 ** -rng.uniform(2, 30, size=50)])
    q_ref = stats.false_discovery_control(p, method="bh")
    q = 10 ** -bh_neglog10q(-np.log10(p))
    np.testing.assert_allclose(q, q_ref, rtol=1e-9)


def test_bh_handles_underflow_and_nan() -> None:
    L = np.array([400.0, 350.0, np.nan, 0.1])
    nlq = bh_neglog10q(L)
    assert np.isfinite(nlq).all()
    assert nlq[0] == pytest.approx(400 - math.log10(4))
    assert nlq[2] == 0.0  # NaN counts as p = 1 -> q = 1


def test_storey_pi1() -> None:
    # half the p-values are 0.75 (> 0.5): pi0 = 0.5 / 0.5 = 1 -> pi1 = 0
    assert storey_pi1(-np.log10(np.array([0.75, 0.75, 1e-3, 1e-4]))) == pytest.approx(0.0)
    assert storey_pi1(np.array([10.0, 20.0])) == pytest.approx(1.0)
    assert math.isnan(storey_pi1(np.array([])))


def test_lambda_gc_and_control_histogram() -> None:
    assert lambda_gc_from_median_log10p(-math.log10(0.5)) == pytest.approx(1.0, rel=1e-9)
    rng = np.random.default_rng(2)
    L = -np.log10(rng.uniform(size=100_000))
    c = ControlStats()
    c.add(L[:50_000])
    c.add(np.append(L[50_000:], np.nan))
    assert c.n == 100_000
    assert c.median_log10p() == pytest.approx(float(np.median(L)), abs=2e-4)
    assert c.lambda_gc() == pytest.approx(1.0, abs=0.02)
    assert c.n_sig == int((L > 5).sum())


def test_add_levels_requires_sign_and_uses_tested_count() -> None:
    df = pd.DataFrame(
        {
            "status": ["tested", "tested", "tested", "variant_absent", "not_tested"],
            "beta_ref": [0.5, 0.5, -0.5, 0.5, 0.5],
            "beta": [0.4, -0.4, -0.1, np.nan, np.nan],
            "log10p": [6.0, 10.0, 1.4, np.nan, np.nan],
        }
    )
    out = add_levels(df, LEVELS)
    assert out["sign_concordant"].tolist()[:3] == [True, False, True]
    assert out["sign_concordant"].isna().tolist() == [False, False, False, True, True]
    assert out["rep_sign"].tolist() == [True, False, True, False, False]
    assert out["rep_strict"].tolist() == [
        True,
        False,
        False,
        False,
        False,
    ]  # p 1e-10 but discordant
    assert out["rep_nominal"].tolist() == [True, False, True, False, False]
    # Bonferroni over 3 tested pairs: p < 0.0167 -> log10p > 1.778; 1.4 fails
    assert out["rep_bonferroni"].tolist() == [True, False, False, False, False]


def test_pair_assoc_stats_raw_rint_and_missing_rint(tmp_path: Path, qpgen_fixtures: Path) -> None:
    df = read_pair_assoc(qpgen_fixtures / "pair_assoc.tsv.gz")
    raw = pair_assoc_stats(df, False)
    rint = pair_assoc_stats(df, True)
    assert len(raw) == len(df) == 9
    assert raw["log10p"].tolist() == df["LOG10P"].tolist()
    assert rint["beta"].tolist() == df["BETA_RINT"].tolist()
    with pytest.raises(StepError, match="skip-rint"):
        pair_assoc_stats(df.drop(columns=["BETA_RINT"]), True)


def test_pair_assoc_short_rows_recomputed(qpgen_fixtures: Path) -> None:
    full = pair_assoc_stats(read_pair_assoc(qpgen_fixtures / "pair_assoc.tsv.gz"), False)
    short = read_pair_assoc(qpgen_fixtures / "pair_assoc.short.tsv.gz")
    assert short.attrs["recomputed_log10p"]
    np.testing.assert_allclose(pair_assoc_stats(short, False)["log10p"], full["log10p"], atol=2e-3)


def test_rect_stats_known_cells_and_controls(qpgen_fixtures: Path) -> None:
    rect = qpgen_fixtures / "rect_assoc.tsv.gz"
    header = gzip.decompress(rect.read_bytes()).decode().splitlines()
    traits = [c[5:] for c in header[0].split("\t") if c.startswith("BETA.")]
    rows = [line.split("\t") for line in header[1:]]
    ids = [r[2] for r in rows]
    pairs = pd.DataFrame(
        {"trait": [traits[0], traits[1], "missing_trait"], "variant": [ids[0], ids[1], ids[0]]}
    )
    chrom = {traits[0]: rows[0][0], traits[1]: "99"}  # other traits without locus
    res = rect_stats(rect, pairs, chrom, chunksize=2)
    assert len(res.known) == 2
    k0 = res.known.set_index("trait").loc[traits[0]]
    col = header[0].split("\t").index(f"LOG10P.{traits[0]}")
    assert k0["log10p"] == pytest.approx(float(rows[0][col]))
    assert res.n_variants == len(ids) and res.n_traits == len(traits)
    assert res.controls.n == len(ids) * len(traits) - 2
    # cross-chromosome controls: trait 0 against variants off its chromosome, trait 1 ("99")
    # against every variant except its known cell; traits without a locus never count
    off0 = sum(1 for r in rows if r[0] != rows[0][0])
    assert res.controls_cross.n == off0 + len(ids) - 1


def test_read_mapped_pairs_types(tmp_path: Path) -> None:
    params = make_run(tmp_path)
    assert params["engine"] == "pair"
    df = read_mapped_pairs(tmp_path / L.dsref_pairs("d", "r"))
    assert df["present"].tolist() == [True] * 5 + [False]
    assert df["pval_ref"].iloc[0] == pytest.approx(1e-20)
    assert df["log10p_ref"].iloc[0] == pytest.approx(20.0)
    assert df["z_ref"].tolist() == pytest.approx([10.0, -8.0, 6.0, 4.0, 2.0, 6.0])  # beta / se


def test_read_mapped_pairs_z_from_pval_without_se(tmp_path: Path) -> None:
    f = tmp_path / "pairs.tsv"
    f.write_text(
        "trait\tvariant\tbeta\tse\tpval\tref_trait\tpresent\n"
        "A\t1:1:A:G\t-0.2\tNA\t0.05\ta\t1\nB\t1:2:A:G\t0.2\tNA\tNA\tb\t1\n"
    )
    df = read_mapped_pairs(f)
    assert df["z_ref"].iloc[0] == pytest.approx(-stats.norm.isf(0.025))
    assert math.isnan(df["z_ref"].iloc[1]) and math.isnan(df["log10p_ref"].iloc[1])


def test_read_mapped_pairs_log10p_below_double_range_or_from_z(tmp_path: Path) -> None:
    f = tmp_path / "pairs.tsv"
    f.write_text(
        "trait\tvariant\tbeta\tse\tpval\tref_trait\tpresent\n"
        "A\t1:1:A:G\t0.5\t0.01\t1e-400\ta\t1\n"  # underflows to 0 as a float
        "B\t1:2:A:G\t0.5\t0.01\t2.5E-310\tb\t1\n"
        "C\t1:3:A:G\t-0.2\t0.1\tNA\tc\t1\n"  # no p-value: from z_ref = -2
        "D\t1:4:A:G\t0.2\t0.1\t0\td\t1\n"
    )
    df = read_mapped_pairs(f)
    lp = df["log10p_ref"].tolist()
    assert lp[0] == pytest.approx(400.0) and lp[1] == pytest.approx(310 - math.log10(2.5))
    assert lp[2] == pytest.approx(-math.log10(2 * stats.norm.sf(2.0)))
    assert lp[3] == pytest.approx(-math.log10(2 * stats.norm.sf(2.0)))
    assert df["z_ref"].iloc[0] == pytest.approx(50.0)


def test_z_from_log10p_matches_normal_and_survives_underflow() -> None:
    from anqchor.steps.cis import z_from_log10p

    L = np.array([-math.log10(0.05), 0.0, 400.0, np.nan])
    z = z_from_log10p(L, np.array([-1.0, 1.0, 1.0, 1.0]))
    assert z[0] == pytest.approx(-1.959964, abs=1e-6) and z[1] == 0.0
    assert np.isfinite(z[2]) and z[2] > 40 and math.isnan(z[3])


def test_wald_z_falls_back_to_the_p_value() -> None:
    from anqchor.steps.cis import wald_z

    # beta / se where se is usable; else the signed normal quantile of the p-value
    z = wald_z(
        np.array([0.6, -0.3, -0.3, 0.2]),
        np.array([0.1, 0.0, np.nan, np.nan]),
        np.array([3.0, -math.log10(0.05), -math.log10(0.05), np.nan]),
    )
    assert z[0] == pytest.approx(6.0)
    assert z[1] == pytest.approx(-1.959964, abs=1e-6) and z[2] == pytest.approx(z[1])
    assert math.isnan(z[3])
    # a strong pair: beta / se keeps its size where the p-value's normal quantile shrinks
    from anqchor.steps.cis import z_from_log10p

    t, dof = 59.0, 660
    lp = -(stats.t.logsf(t, dof) + math.log(2)) / math.log(10)
    assert (
        z_from_log10p(np.array([lp]), np.array([1.0]))[0]
        < 40
        < wald_z(np.array([5.9]), np.array([0.1]), np.array([lp]))[0]
    )


def test_recommend_tie_breaks() -> None:
    def row(label: str, strict: int, bonf: int, n: int = 10) -> dict[str, Any]:
        return {
            "label": label,
            "run": label,
            "n_tested": n,
            "n_rep_strict": strict,
            "n_rep_bonferroni": bonf,
            "n_rep_fdr": 0,
            "n_rep_nominal": 0,
            "n_rep_sign": 0,
        }

    lv = ["sign", "nominal", "fdr", "bonferroni", "strict"]
    r = recommend([row("a", 5, 6), row("b", 7, 1)], lv, LEVELS)
    assert r["label"] == "b" and r["run"] == "b" and "2 more than a" in r["reason"]
    r = recommend([row("a", 5, 6), row("b", 5, 8)], lv, LEVELS)
    assert r["label"] == "b" and "tied with a" in r["reason"] and "bonferroni" in r["reason"]
    r = recommend([row("a", 5, 6), row("b", 5, 6)], lv, LEVELS)
    assert r["label"] == "a" and "every level" in r["reason"]
    assert recommend([row("a", 0, 0, n=0)], lv, LEVELS)["label"] is None


# ------------------------------------------------------------------ full step
def test_summarize_outputs(tmp_path: Path) -> None:
    params = make_run(tmp_path)
    run_summarize(tmp_path, params)
    for n in SUMMARY_FILES:
        assert (tmp_path / L.cis_summary("d", n)).exists()
    long = pd.read_csv(tmp_path / L.cis_summary("d", "pairs.long.tsv.gz"), sep="\t")
    assert len(long) == 6 * 2  # 6 mapped pairs x 2 strategies (one transform per run)
    assert list(long.columns[:18]) == [
        "run",
        "trait",
        "ref_trait",
        "variant",
        "beta_ref",
        "se_ref",
        "z_ref",
        "pval_ref",
        "log10p_ref",
        "status",
        "a1freq",
        "n",
        "beta",
        "se",
        "tstat",
        "z",
        "log10p",
        "sign_concordant",
    ]
    assert [c for c in long.columns if c.startswith("rep_")] == [
        "rep_sign",
        "rep_nominal",
        "rep_fdr",
        "rep_bonferroni",
        "rep_strict",
    ]
    none = long[long.run == "none"].set_index("trait")
    assert none.loc["T4", "status"] == "not_tested" and pd.isna(none.loc["T4", "z"])
    assert none.loc["T6", "status"] == "variant_absent"
    assert none.loc["T5", "sign_concordant"] == 0 and none.loc["T5", "rep_nominal"] == 0
    assert pd.isna(none.loc["T6", "sign_concordant"])
    # z: beta / se (the Wald statistic, as z_ref), aligned with beta (and so with z_ref)
    assert none.loc["T2", "z"] == pytest.approx(-0.30 / 0.05)
    assert none.loc["T5", "z"] < 0 < none.loc["T5", "z_ref"]
    assert none.loc["T1", "log10p_ref"] == pytest.approx(20.0)

    st = pd.read_csv(tmp_path / L.cis_summary("d", "strategies.tsv"), sep="\t").set_index("label")
    assert list(st.index) == ["full", "none"] and "transform" not in st.columns
    full = st.loc["full"]
    assert full.rint == 0
    assert (full.n_pairs, full.n_tested, full.n_variant_absent, full.n_not_tested) == (6, 5, 1, 0)
    assert (full.n_rep_sign, full.n_rep_nominal, full.n_rep_strict) == (5, 5, 3)
    assert full.n_cov == 10
    # Bonferroni 0.05/5 = 0.01 -> log10p > 2: 30, 15, 8, 4.5 pass, 1.5 fails
    assert full.n_rep_bonferroni == 4
    assert full.frac_rep_strict == pytest.approx(3 / 5)
    assert full.median_log10p == pytest.approx(8.0)
    assert full.rel_neff == pytest.approx(1.0)  # the baseline
    none_row = st.loc["none"]
    assert (none_row.n_tested, none_row.n_not_tested, none_row.n_rep_sign) == (4, 1, 3)
    assert none_row.frac_sign_discordant == pytest.approx(0.25)
    assert 0 < none_row.rel_neff < 1  # weaker statistics on the same pairs

    s = read_json(tmp_path / L.cis_summary("d", "summary.json"))
    assert (
        s["schema_version"] == SCHEMA_VERSION == 3 and s["engine"] == "pair" and s["rint"] is False
    )
    assert s["recommended"]["label"] == "full" and s["baseline"] == "full"
    assert [r["label"] for r in s["runs"]] == ["full", "none"]
    assert s["n_pairs_unmapped"] == 1 and s["n_variant_absent"] == 1
    assert s["n_absent_variants_indel"] == 1
    assert s["unresolved_fraction"] == pytest.approx(2 / 6)
    assert s["ranking"]["strict"][0] == "full"
    assert set(s["per_run"]) == {"full", "none"}
    assert any("not map" in w for w in s["warnings"])
    assert any("no (finite) qpgen result" in w for w in s["warnings"])
    assert s["level_order"] == ["sign", "nominal", "fdr", "bonferroni", "strict"]

    traits = pd.read_csv(tmp_path / L.cis_summary("d", "traits.tsv"), sep="\t").set_index("trait")
    assert traits.loc["T1", "best_log10p:full"] == pytest.approx(30.0)
    assert traits.loc["T1", "delta_log10p:none"] == pytest.approx(12.0 - 30.0)
    assert traits.loc["T1", "level:full"] == "strict"
    assert traits.loc["T5", "level:none"] == "none"
    assert pd.isna(traits.loc["T6", "level:full"]) and pd.isna(traits.loc["T6", "best_run"])
    assert traits.loc["T1", "best_run"] == "full"
    assert "delta_log10p:full" not in traits.columns  # baseline has no delta

    curve = pd.read_csv(tmp_path / L.cis_summary("d", "power_curve.tsv"), sep="\t")
    c = curve[curve.run == "full"]
    assert c.threshold_log10p.iloc[0] == 0 and c.n_replicated.iloc[0] == 5
    assert (np.diff(c.n_replicated) <= 0).all()
    assert c.threshold_log10p.max() == pytest.approx(30.0)
    # p < 10^-t: the pair at exactly -log10 p = 8 does not count at t = 8
    assert int(c[c.threshold_log10p == 8.0].n_replicated.iloc[0]) == 2


def test_summarize_fails_when_too_much_is_unresolved(tmp_path: Path) -> None:
    params = make_run(tmp_path, weak=[("T1", "1:100:A:G", 0.4, 12.0)])
    params["max_unresolved_fraction"] = 0.3
    with pytest.raises(StepError, match="max_unresolved_fraction"):
        run_summarize(tmp_path, params)


def test_summarize_rint_and_disabled_levels(tmp_path: Path) -> None:
    params = make_run(tmp_path, rint=True)
    params["levels"] = {"nominal": None, "strict": 1e-5, "bonferroni": False, "fdr": None}
    run_summarize(tmp_path, params)
    st = pd.read_csv(tmp_path / L.cis_summary("d", "strategies.tsv"), sep="\t")
    assert st["label"].tolist() == ["full", "none"] and st["rint"].tolist() == [1, 1]
    assert "n_rep_bonferroni" not in st.columns and "n_rep_strict" in st.columns
    s = read_json(tmp_path / L.cis_summary("d", "summary.json"))
    assert s["rint"] is True and all(r["rint"] for r in s["runs"])
    assert s["recommended"]["label"] == "full" and s["baseline"] == "full"
    long = pd.read_csv(tmp_path / L.cis_summary("d", "pairs.long.tsv.gz"), sep="\t")
    # RINT columns were 0.9 x raw in the fake output
    t1 = long[(long.run == "full") & (long.trait == "T1")]
    assert t1["log10p"].iloc[0] == pytest.approx(27.0)


def test_power_curve_empty_strategy() -> None:
    long = pd.DataFrame(
        {
            "run": ["a"],
            "status": ["variant_absent"],
            "log10p": [np.nan],
            "sign_concordant": pd.array([pd.NA], dtype="boolean"),
        }
    )
    c = power_curve(long, ["a"])
    assert len(c) == 1 and math.isnan(c["frac_replicated"].iloc[0])


def test_power_curve_counts_p_below_the_threshold() -> None:
    # p < 10^-t, as the levels: a pair exactly at -log10 p = 5 does not count at t = 5
    long = pd.DataFrame(
        {
            "run": ["a"] * 3,
            "status": ["tested"] * 3,
            "log10p": [5.0, 5.0001, 0.5],
            "sign_concordant": [True, True, True],
        }
    )
    c = power_curve(long, ["a"], step=0.5)
    n = dict(zip(c["threshold_log10p"], c["n_replicated"], strict=True))
    assert n[0.0] == 3 and n[0.5] == 2 and n[5.0] == 1


def test_trait_table_multiple_pairs_per_trait() -> None:
    df = pd.DataFrame(
        {
            "run": ["s"] * 3,
            "trait": ["A", "A", "B"],
            "ref_trait": ["a", "a", "b"],
            "variant": ["v1", "v2", "v3"],
            "status": ["tested", "tested", "tested"],
            "beta_ref": [1.0, 1.0, 1.0],
            "beta": [1.0, -1.0, 1.0],
            "log10p": [2.0, 9.0, 0.5],
        }
    )
    df = add_levels(df, LEVELS)
    t = trait_table(df, ["s"], active_levels(LEVELS)).set_index("trait")
    assert t.loc["A", "n_pairs"] == 2
    assert t.loc["A", "best_log10p:s"] == pytest.approx(9.0)  # best p regardless of sign
    assert t.loc["A", "level:s"] == "bonferroni"  # concordant pair: p 0.01 < 0.05/3
    assert t.loc["B", "level:s"] == "sign"


def test_pair_assoc_short_rows_use_covariate_df(qpgen_fixtures: Path) -> None:
    """Old builds without raw LOG10P: the recomputation uses N - 2 - n_cov df (and the
    TSTAT rescaled from the N - 2 residual variance), not N - 2."""
    from anqchor.qpgen.outputs import t_log10p

    base = read_pair_assoc(qpgen_fixtures / "pair_assoc.short.tsv.gz")
    k = 5
    adj = read_pair_assoc(qpgen_fixtures / "pair_assoc.short.tsv.gz", n_cov=k)
    n = base["N"].astype(float)
    dof = n - 2 - k
    expected = t_log10p(base["TSTAT"].astype(float) * np.sqrt(dof / (n - 2)), dof)
    np.testing.assert_allclose(adj["LOG10P"], expected)
    assert (adj["LOG10P"] < base["LOG10P"]).all()


def test_levels_ranked_by_threshold_when_bonferroni_is_stricter() -> None:
    """With n_tested > 5000, 0.05/n_tested < 1e-5: Bonferroni is the strongest level."""
    from anqchor.steps.cis import level_cutoffs, pair_level, stringency_order

    n = 10_000
    lp = np.full(n, 0.1)
    lp[0], lp[1] = 5.2, 9.0  # p = 6.3e-6 (strict only) and 1e-9 (both)
    df = pd.DataFrame(
        {
            "run": "s",
            "trait": [f"T{i}" for i in range(n)],
            "ref_trait": "x",
            "variant": [f"v{i}" for i in range(n)],
            "status": "tested",
            "beta_ref": 1.0,
            "beta": 1.0,
            "log10p": lp,
        }
    )
    df = add_levels(df, LEVELS)
    levels = active_levels(LEVELS)
    order = stringency_order(level_cutoffs(df, levels, LEVELS))
    assert order[-1] == "bonferroni" and order.index("strict") < order.index("bonferroni")
    lv = pair_level(df, levels, order)
    assert lv.iloc[0] == "strict" and lv.iloc[1] == "bonferroni"
    t = trait_table(df, ["s"], levels, LEVELS).set_index("trait")
    assert t.loc["T0", "level:s"] == "strict" and t.loc["T1", "level:s"] == "bonferroni"
    rows = [
        {
            "run": s,
            "label": s,
            "n_tested": n,
            "n_rep_sign": 9,
            "n_rep_nominal": 9,
            "n_rep_fdr": 9,
            "n_rep_bonferroni": b,
            "n_rep_strict": st,
        }
        for s, b, st in (("a", 1, 9), ("b", 2, 3))
    ]
    assert recommend(rows, levels, LEVELS, order)["label"] == "b"  # keyed on Bonferroni
    assert recommend(rows, levels, LEVELS)["label"] == "a"  # default order: strict first


def test_long_rows_na_statistics_are_not_tested() -> None:
    """A qpgen row without a finite estimate must not count in n_tested (as in trans)."""
    from anqchor.steps.cis import long_rows, strategy_row

    pairs = pd.DataFrame(
        {
            "trait": ["A", "B"],
            "ref_trait": ["a", "b"],
            "variant": ["v1", "v2"],
            "beta_ref": [1.0, 1.0],
            "se_ref": [0.1, 0.1],
            "z_ref": [10.0, 10.0],
            "pval_ref": [1e-9, 1e-9],
            "log10p_ref": [9.0, 9.0],
            "present": [True, True],
        }
    )
    st = pd.DataFrame(
        {
            "trait": ["A", "B"],
            "variant": ["v1", "v2"],
            "a1freq": [0.3, 0.0],
            "n": [100.0, 100.0],
            "beta": [1.0, np.nan],
            "se": [0.1, np.nan],
            "tstat": [10.0, np.nan],
            "log10p": [9.0, np.nan],
        }
    )
    df = long_rows(pairs, st, "s", LEVELS)
    assert df["status"].tolist() == ["tested", "not_tested"]
    row = strategy_row(df, "s", False, active_levels(LEVELS), None)
    assert row["n_tested"] == 1 and row["n_not_tested"] == 1 and row["frac_rep_strict"] == 1.0
