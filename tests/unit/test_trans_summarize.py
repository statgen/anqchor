"""steps.trans: statistics helpers and ``summarize_trans`` on a small synthetic rectangle.

The rectangle (window 1 Mb, bins 1 Mb) has six variants and five genes:

=========  =================  ===========================================
gene       locus              note
=========  =================  ===========================================
gA         chr1:1,000
gB         chr1:5,000,000
gC         chr2:1,000
gD         chr2:6,800,000     in cis with bin chr2:5-6 Mb (variant v6)
gE         none
=========  =================  ===========================================

Variants v1 chr1:1,100 and v2 chr1:900,000 (bin chr1:0-1 Mb), v3 chr1:6,000,000,
v4 chr2:2,000, v5 chr2:5,500,000 and v6 chr2:5,900,000 (bin chr2:5-6 Mb). Reference pairs:
gA-v1 (cis), gC-v1 (trans, other chromosome), gD-v4 (trans, same chromosome 6.8 Mb away,
opposite sign) and gA with an absent variant (cis).
"""

from __future__ import annotations

import gzip
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from anqchor.pipelines.layout import Layout as L
from anqchor.qpgen.outputs import iter_rect_assoc, read_rect_header
from anqchor.stats import passes
from anqchor.steps import StepContext
from anqchor.steps import trans as T

LEVELS = {"nominal": 1e-5, "strict": 1e-8, "bonferroni": True}


# ------------------------------------------------------------------ helpers
def test_lambda_gc_matches_chi_square_definition() -> None:
    assert T.lambda_gc(-math.log10(0.5)) == pytest.approx(1.0, abs=1e-9)
    p = 0.2
    want = stats.chi2.isf(p, 1) / stats.chi2.ppf(0.5, 1)
    assert T.lambda_gc(-math.log10(p)) == pytest.approx(want, rel=1e-9)
    assert math.isfinite(T.lambda_gc(400.0))  # far below the double range of p
    assert math.isnan(T.lambda_gc(float("nan")))


def test_qq_points_are_thinned_keep_the_tail_and_have_a_band() -> None:
    lp = np.sort(-np.log10(np.random.default_rng(2).uniform(size=50_000)))[::-1]
    q = T.qq_points(lp, max_points=500, head=100)
    assert len(q) <= 501
    assert q["rank"].iloc[0] == 1 and q["rank"].iloc[-1] == len(lp)
    assert (q["rank"].iloc[:100] == np.arange(1, 101)).all()
    assert q["expected_log10p"].is_monotonic_decreasing
    assert (q["band_lo_log10p"] <= q["expected_log10p"] + 1e-9).all()
    assert (q["band_hi_log10p"] >= q["expected_log10p"] - 0.31).all()
    # the band of the smallest p of n uniforms: Beta(1, n) quantiles
    n = len(lp)
    assert q["band_hi_log10p"].iloc[0] == pytest.approx(
        -math.log10(stats.beta.ppf(0.025, 1, n)), rel=1e-9
    )
    assert q["band_lo_log10p"].iloc[0] == pytest.approx(
        -math.log10(stats.beta.ppf(0.975, 1, n)), rel=1e-9
    )
    assert len(T.qq_points(lp[:10])) == 10
    assert T.qq_points(np.zeros(0)).empty


def test_grid_thresholds_merge_levels_and_add_bonferroni() -> None:
    th = T.grid_thresholds({"nominal": 1e-5, "strict": 1e-8, "bonferroni": True}, 100_000)
    alphas = [t.alpha for t in th]
    assert alphas == sorted(alphas, reverse=True)
    keys = [t.key for t in th]
    assert keys == ["p1e-4", "p1e-5", "p1e-6", "bonferroni", "p5e-8", "p1e-8", "p1e-10"]
    by = {t.key: t for t in th}
    assert by["p1e-5"].roles == ["nominal"] and by["p1e-8"].roles == ["strict"]
    assert by["bonferroni"].alpha == pytest.approx(5e-7) and by["bonferroni"].label.startswith(
        "Bonferroni"
    )
    # a level outside the fixed list is added; no Bonferroni without trans cells or when off
    th = T.grid_thresholds({"nominal": 1e-3, "strict": None, "bonferroni": False}, 10)
    assert (
        th[0].key == "p0.001"
        and th[0].roles == ["nominal"]
        and "bonferroni" not in [t.key for t in th]
    )
    assert T.grid_k(4) == [1, 2, 3, 4, 5, 10, 20] and T.grid_k(5) == [1, 2, 3, 5, 10, 20]


def test_pair_class() -> None:
    assert T.pair_class("1", 1100, "1", 1000, 1_000_000) == ("cis", "cis")
    assert T.pair_class("1", 1_001_000, "1", 1000, 1_000_000) == ("cis", "cis")  # exactly 1 Mb
    assert T.pair_class("1", 1_001_001, "1", 1000, 1_000_000) == ("trans", "same_chrom_distal")
    assert T.pair_class("1", 1100, "2", 1000, 1_000_000) == ("trans", "cross_chrom")
    assert T.pair_class("1", 1100, None, None, 1_000_000) == ("no_locus", "no_locus")


# --------------------------------------------------------- synthetic run dir
TRAITS = ["gA", "gB", "gC", "gD", "gE"]
LOCI = {
    "gA": ("1", 1000),
    "gB": ("1", 5_000_000),
    "gC": ("2", 1000),
    "gD": ("2", 6_800_000),
    "gE": (None, None),
}
VARIANTS = [
    ("1", 1100, "A", "G"),
    ("1", 900_000, "C", "T"),
    ("1", 6_000_000, "G", "A"),
    ("2", 2000, "G", "A"),
    ("2", 5_500_000, "T", "C"),
    ("2", 5_900_000, "A", "C"),
]
# -log10 p per (variant, gene); NaN = not tested
LP = np.array(
    [
        [12.0, 5.5, 8.5, 0.3, 0.2],  # v1: gA cis (ref), gB distal, gC cross (ref trans), gD cross
        [3.0, 0.4, 7.0, 6.5, 0.1],  # v2
        [4.5, 2.0, 0.6, 0.2, 0.3],  # v3: gA distal, gB cis (exactly 1 Mb)
        [0.4, 0.2, 15.0, 3.0, 0.5],  # v4: gC cis, gD distal (ref trans, opposite sign)
        [5.2, 0.1, 0.8, 9.0, np.nan],  # v5: gD distal but in cis with the bin (via v6)
        [0.25, 0.35, 0.9, 11.0, 0.2],  # v6: gD cis
    ]
)
BETA = np.where(np.isnan(LP), np.nan, 0.3)
BETA[3, 3] = -0.3  # gD-v4: opposite sign to the reference
EXPECTED_CELLS = {"cis": 5, "cross_chrom": 12, "same_chrom_distal": 7, "no_locus": 5}


def write_rect(path: Path, lp: np.ndarray, beta: np.ndarray) -> None:
    head = [
        "#CHROM",
        "GENPOS",
        "ID",
        "ALLELE0",
        "ALLELE1",
        "A1FREQ",
        "N",
        "N_RR",
        "N_RA",
        "N_AA",
        "TEST",
    ]
    for t in TRAITS:
        head += [f"BETA.{t}", f"SE.{t}", f"TSTAT.{t}", f"LOG10P.{t}"]
    lines = ["\t".join(head)]
    for i, (c, pos, a0, a1) in enumerate(VARIANTS):
        row = [
            c,
            str(pos),
            f"{c}:{pos}:{a0}:{a1}",
            a0,
            a1,
            "0.3",
            "100",
            "50",
            "40",
            "10",
            "Linear",
        ]
        for j in range(len(TRAITS)):
            if np.isnan(lp[i, j]):
                row += ["NA"] * 4
            else:
                row += [f"{beta[i, j]:.6g}", "0.1", f"{beta[i, j] / 0.1:.6g}", f"{lp[i, j]:.6g}"]
        lines.append("\t".join(row))
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


RECT = L.rect_assoc("trans", "ds", "default", "raw")


def make_run(tmp: Path, lp: np.ndarray | None = None) -> Path:
    run = tmp / "run"
    files = {
        L.geno_json: json.dumps(
            {
                "chrom_style": "plain",
                "chromosomes": ["1", "2"],
                "entries": [{"chrom": "1"}, {"chrom": "2"}],
            }
        ),
        L.ds_traits("ds"): "trait\tchrom\tstart\tend\n"
        + "".join(f"{t}\t{c or 'NA'}\t{s or 'NA'}\t{s or 'NA'}\n" for t, (c, s) in LOCI.items()),
        L.trans_variants("ds"): "".join(f"{c}:{p}:{a}:{b}\n" for c, p, a, b in VARIANTS),
        L.trans_traits("ds"): "\n".join(TRAITS) + "\n",
        L.trans_lists_json("ds"): json.dumps(
            {
                "traits": {
                    "trait_chroms": "all",
                    "n_kept": 5,
                    "n_candidates": 5,
                    "n_dropped_no_locus": 0,
                    "allowed_chroms": None,
                    "source": "all",
                },
                "variants": {"n_kept": 6},
            }
        ),
        L.strategy_json("trans", "ds", "default"): json.dumps(
            {"n_cov": 6, "cov_source": "dataset"}
        ),
    }
    for p, text in files.items():
        (run / p).parent.mkdir(parents=True, exist_ok=True)
        (run / p).write_text(text)
    pairs = (
        "trait\tvariant\tbeta\tse\tpval\tref_trait\tpresent\n"
        "gA\t1:1100:A:G\t0.4\t0.1\t1e-8\tGA\t1\n"
        "gC\t1:1100:A:G\t0.4\t0.1\t1e-9\tGC\t1\n"
        "gD\t2:2000:G:A\t0.8\t0.1\t1e-20\tGD\t1\n"
        "gA\t1:999:A:C\t0.1\t0.1\t1e-3\tGA\t0\n"
    )
    (run / L.dsref_pairs("ds", "ref")).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(run / L.dsref_pairs("ds", "ref"), "wt") as fh:
        fh.write(pairs)
    write_rect(run / RECT, LP if lp is None else lp, BETA)
    return run


def params(**kw: Any) -> dict[str, Any]:
    p = {
        "dataset": "ds",
        "reference": "ref",
        "run": "myrun",
        "strategy": "default",
        "rect": str(RECT),
        "rint": False,
        "variants": "reference",
        "traits": "all",
        "trait_chroms": "all",
        "known_pairs": True,
        "cis_window_bp": 1_000_000,
        "min_maf": 0.01,
        "min_mac": 1.0,
        "levels": dict(LEVELS),
        "hotspot_bin_bp": 1_000_000,
        "hotspot_min_traits": 2,
        "max_sig_rows": 1000,
        "schema_version": T.SCHEMA_VERSION,
    }
    p.update(kw)
    return p


def summarize(run: Path, **kw: Any) -> dict[str, Any]:
    ctx = StepContext(
        name="trans.ds.summarize", run_dir=run, params=params(**kw), inputs=[], outputs=[]
    )
    T.summarize_trans(ctx)

    def tsv(name: str) -> pd.DataFrame:
        return pd.read_csv(
            run / L.trans_summary("ds", name), sep="\t", dtype={"chrom": str, "trait_chrom": str}
        )

    return {
        "summary": json.loads((run / L.trans_summary("ds", "summary.json")).read_text()),
        "pairs": tsv("reference_pairs.tsv"),
        "bins": tsv("hotspots.tsv"),
        "genes": tsv("hotspot_genes.tsv"),
        "qq": tsv("qq.tsv"),
        "variants": tsv("variants.tsv"),
        "sig": tsv("pairs.sig.tsv.gz"),
    }


def test_cell_classes_and_run_facts(tmp_path: Path) -> None:
    sm = summarize(make_run(tmp_path))["summary"]
    assert sm["schema_version"] == 2 and sm["run"] == "myrun" and sm["rint"] is False
    assert sm["cells"] == EXPECTED_CELLS
    assert sm["n_cells"] == 29 and sm["n_cells_trans"] == 19
    assert (sm["n_variants_tested"], sm["n_traits_tested"], sm["n_samples"], sm["n_cov"]) == (
        6,
        5,
        100,
        6,
    )
    assert set(T.DEFINITIONS) <= set(sm["definitions"])


def test_replication_by_class(tmp_path: Path) -> None:
    out = summarize(make_run(tmp_path))
    rep = out["summary"]["replication"]
    assert rep["levels"] == {"nominal": 0.05, "strict": 1e-5}
    assert {
        k: rep["cis"][k]
        for k in ("n_pairs", "n_tested", "n_rep_sign", "n_rep_nominal", "n_rep_strict")
    } == {"n_pairs": 2, "n_tested": 1, "n_rep_sign": 1, "n_rep_nominal": 1, "n_rep_strict": 1}
    t = rep["trans"]
    assert (
        t["n_pairs"],
        t["n_tested"],
        t["n_rep_sign"],
        t["n_rep_nominal"],
        t["n_rep_strict"],
    ) == (2, 2, 1, 1, 1)
    assert t["n_pairs_cross_chrom"] == 1 and t["n_pairs_same_chrom_distal"] == 1
    k = out["pairs"].set_index(["trait", "variant"])
    assert k.loc[("gC", "1:1100:A:G"), "pair_class"] == "trans"
    assert k.loc[("gC", "1:1100:A:G"), "trans_kind"] == "cross_chrom"
    assert k.loc[("gD", "2:2000:G:A"), "trans_kind"] == "same_chrom_distal"
    assert k.loc[("gD", "2:2000:G:A"), "rep_sign"] == 0  # opposite sign
    assert k.loc[("gA", "1:999:A:C"), "status"] == "variant_absent"
    assert k.loc[("gA", "1:999:A:C"), "pair_class"] == "cis"  # classified from the variant ID
    assert pd.isna(k.loc[("gA", "1:999:A:C"), "sign_concordant"])
    assert k.loc[("gA", "1:1100:A:G"), "z"] == pytest.approx(3.0)
    assert (out["pairs"]["run"] == "myrun").all()


def test_hotspot_bins_count_distinct_genes_and_skip_cis_genes(tmp_path: Path) -> None:
    out = summarize(make_run(tmp_path))
    b = out["bins"].set_index(["chrom", "bin_start"])
    assert list(b.index) == [("1", 0), ("1", 6_000_000), ("2", 0), ("2", 5_000_000)]
    assert b["n_variants"].tolist() == [2, 1, 1, 2]
    # chr1:0-1 Mb: gB 5.5, gC 8.5 (twice, counted once), gD 6.5; gA is in cis with the bin
    assert b.loc[("1", 0), "n_genes_p1e-5"] == 3 and b.loc[("1", 0), "n_genes_p1e-6"] == 2
    assert b.loc[("1", 0), "n_genes_p1e-8"] == 1 and b.loc[("1", 0), "top_trait"] == "gC"
    assert b.loc[("1", 0), "top_variant"] == "1:1100:A:G" and b.loc[("1", 0), "top_log10p"] == 8.5
    # chr2:5-6 Mb: gD (9 from v5) is in cis with v6, so only gA counts
    assert b.loc[("2", 5_000_000), "n_genes_p1e-5"] == 1
    assert b.loc[("2", 5_000_000), "top_trait"] == "gA"
    # chr2:0-1 Mb: gD at 3 passes only Bonferroni (0.05 / 19 trans cells: -log10 = 2.58)
    assert b.loc[("2", 0), "n_genes_bonferroni"] == 1 and b.loc[("2", 0), "n_genes_p1e-4"] == 0
    g = out["genes"]
    assert not ((g["trait"] == "gD") & (g["bin_start"] == 5_000_000)).any()
    assert set(g[g["chrom"] == "1"]["trait"]) == {"gA", "gB", "gC", "gD"}  # gA from v3 (distal)
    h = out["summary"]["hotspots"]
    keys = [t["key"] for t in h["thresholds"]]
    assert keys[0] == "bonferroni" and h["k"] == [1, 2, 3, 5, 10, 20]
    grid = dict(zip(keys, h["counts"], strict=True))
    assert grid["p1e-5"] == [2, 1, 1, 0, 0, 0]
    assert grid["bonferroni"] == [4, 1, 1, 0, 0, 0]
    assert h["default"] == {"threshold": "p1e-5", "k": 2, "n_hotspots": 1}
    assert h["bonferroni_p"] == pytest.approx(0.05 / 19, rel=1e-3)


def test_bin_width_changes_the_bins(tmp_path: Path) -> None:
    out = summarize(make_run(tmp_path), hotspot_bin_bp=10_000_000)
    b = out["bins"]
    assert b["n_variants"].tolist() == [3, 3]
    assert out["summary"]["hotspots"]["bin_bp"] == 10_000_000


def test_calibration_leaves_out_reference_trans_regions(tmp_path: Path) -> None:
    out = summarize(make_run(tmp_path))
    inf = out["summary"]["inflation"]
    # reference trans pairs gC-v1 and gD-v4: variants within 1 Mb (v1, v2, v4) and genes gC, gD
    assert inf["n_reference_trans_pairs"] == 2
    assert inf["n_variants_excluded"] == 3 and inf["n_genes_excluded"] == 2
    assert inf["n_cells"] == 4 and inf["n_cells_excluded"] == 8  # v5/v6 x gA/gB remain
    assert inf["median_log10p"] == pytest.approx(0.3)
    assert inf["lambda_gc"] == pytest.approx(T.lambda_gc(0.3), rel=1e-3)
    assert [(c["observed"], c["expected"]) for c in inf["counts"]] == [
        (1, pytest.approx(0.004)),
        (1, pytest.approx(4e-5)),
    ]
    q = out["qq"]
    assert q["n_cells"].iloc[0] == 4 and q["observed_log10p"].tolist() == pytest.approx(
        [5.2, 0.35, 0.25, 0.1]
    )
    v = out["variants"].set_index("variant")
    assert v["excluded_from_calibration"].tolist() == [1, 1, 0, 1, 0, 0]


def test_without_reference_pairs_nothing_is_left_out(tmp_path: Path) -> None:
    out = summarize(make_run(tmp_path), known_pairs=False)
    inf = out["summary"]["inflation"]
    assert (
        inf["n_cells"] == 12
        and inf["n_cells_excluded"] == 0
        and inf["n_reference_trans_pairs"] == 0
    )
    rep = out["summary"]["replication"]
    assert rep["cis"]["n_pairs"] == 0 and rep["trans"]["n_pairs"] == 0
    assert len(out["pairs"]) == 0 and "rep_strict" in out["pairs"].columns


def test_variants_and_significant_cells(tmp_path: Path) -> None:
    out = summarize(make_run(tmp_path))
    v = out["variants"].set_index("variant")
    assert v["n_cis"].tolist() == [1, 1, 1, 1, 0, 1]
    assert v["n_cross_chrom"].tolist() == [2, 2, 2, 2, 2, 2]
    assert v["n_sig_trans"].tolist() == [2, 2, 0, 0, 2, 0]  # p < 1e-5
    assert v.loc["2:5500000:T:C", "top_trans_trait"] == "gD"
    assert v.loc["1:1100:A:G", "max_log10p_cis"] == 12.0
    sig = out["sig"]
    assert len(sig) == 9 and out["summary"]["n_sig_cells"] == 9
    assert sig["log10p"].tolist() == sorted(sig["log10p"].tolist(), reverse=True)
    assert set(sig["cell_class"]) == {"cis", "cross_chrom", "same_chrom_distal"}
    assert (
        sig[(sig["trait"] == "gC") & (sig["variant"] == "1:1100:A:G")]["reference_pair"].iloc[0]
        == 1
    )


def test_sig_rows_cap_is_reported(tmp_path: Path) -> None:
    out = summarize(make_run(tmp_path), max_sig_rows=2)
    assert out["sig"]["log10p"].tolist() == [15.0, 12.0]
    assert out["summary"]["n_sig_cells"] == 9 and out["summary"]["sig_rows_written"] == 2
    assert any("capped" in w for w in out["summary"]["warnings"])


def test_inflation_warning(tmp_path: Path) -> None:
    lp = np.where(np.isnan(LP), np.nan, LP * 3 + 1)
    out = summarize(make_run(tmp_path, lp))
    assert any("λGC" in w and "inflated" in w for w in out["summary"]["warnings"])


def test_streaming_is_chunk_size_invariant(tmp_path: Path) -> None:
    """Every chunk size retains the same cells, including the Bonferroni cells kept by
    the running bound (gD-v4 at -log10 p 3 is below every fixed threshold)."""
    run = make_run(tmp_path)
    path = run / RECT
    _, names = read_rect_header(path)
    loci = T.read_loci(run / L.ds_traits("ds"))
    pairs = T.read_known(run / L.dsref_pairs("ds", "ref"))
    excl = T.build_exclusion(pairs, loci, "plain", 1_000_000)
    res = []
    for cs in (1, 2, 100):
        res.append(
            T.summarize_rect(
                iter_rect_assoc(path, ("BETA", "SE", "LOG10P"), chunksize=cs),
                trait_names=names,
                loci=loci,
                known={"1:1100:A:G": ["gA", "gC"], "2:2000:G:A": ["gD"]},
                exclusion=excl,
                cis_window_bp=1_000_000,
                keep_log10p=4.0,
                nominal_log10p=5.0,
                bonferroni=True,
            )
        )
    thresholds = T.grid_thresholds(LEVELS, 19)
    grids = []
    for r in res:
        assert r.n_cells == res[0].n_cells
        np.testing.assert_allclose(r.cross_desc, res[0].cross_desc)
        assert r.known_cells == res[0].known_cells
        bins = T.bins_of(r.meta, 1_000_000, ["1", "2"])
        grids.append(T.hotspot_tables(r, bins, thresholds, [1, 2], 1_000_000)[2])
        kept = {(int(a), int(b)) for a, b in zip(r.kept["row"], r.kept["col"], strict=True)}
        assert (3, 3) in kept  # gD-v4, -log10 p 3 > final Bonferroni 2.58
    assert grids[0] == grids[1] == grids[2]
    assert res[0].cross_desc.dtype == np.float32  # the documented memory bound


def test_levels_match_cis_at_thresholds() -> None:
    """p exactly at a level (LOG10P printed as 5 for p = 1e-5) does not replicate (p < alpha)."""
    lpa = np.array([5.0, -math.log10(0.05), 5.000001])
    assert passes(lpa, -math.log10(T.REPLICATION_LEVELS["strict"])).tolist() == [False, False, True]
    assert passes(lpa, -math.log10(T.REPLICATION_LEVELS["nominal"])).tolist() == [True, False, True]


def test_summary_json_keeps_small_thresholds() -> None:
    """The Bonferroni p (~1e-7 on real data) must not be rounded to 0 in summary.json."""
    assert T._round(0.05 / 473_928) == pytest.approx(1.055e-07, rel=1e-3)
    assert T._round(1.01234567) == 1.0123
    assert T._round(0.0) == 0.0
    assert T._round(float("nan")) is None
    assert T._round(3) == 3
    assert T._round({"a": [np.float32(0.5), np.int64(2)]}) == {"a": [0.5, 2]}
