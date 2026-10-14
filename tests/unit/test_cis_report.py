"""report_cis: page order and content, chart data, single run vs comparison."""

from __future__ import annotations

import itertools
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from test_cis_summarize import make_run, run_summarize

from anqchor.pipelines.layout import Layout as L
from anqchor.report.cis import (
    DASHES,
    TITLE,
    neff_interval,
    power_counts,
    report_cis,
    run_colors,
    run_dash,
)
from anqchor.report.theme import LIGHT, series_color
from anqchor.steps import StepContext


def report_params(**extra: Any) -> dict[str, Any]:
    return {
        "project": "p",
        "datasets": ["d"],
        "reference": "r",
        "engine": "pair",
        "plotlyjs": "cdn",
        **extra,
    }


def render(run: Path, single: bool = False, **extra: Any) -> str:
    """Summarise the fixture run (two runs ``full`` and ``none``, or ``full`` alone) and
    render its report."""
    params = make_run(run)
    if single:
        params["strategies"] = params["strategies"][:1]
    run_summarize(run, params)
    report_cis(
        StepContext(
            name="report.cis",
            run_dir=run,
            params=report_params(**extra),
            inputs=[],
            outputs=[L.report("cis")],
        )
    )
    return (run / L.report("cis")).read_text(encoding="utf-8")


def page_data(html: str) -> dict[str, Any]:
    m = re.search(r"window\.ANQ_DATA = (.*?);\n", html)
    assert m
    data: dict[str, Any] = json.loads(m.group(1))
    return data


def visible(html: str) -> str:
    """The page without its scripts and styles (what a reader sees, plus markup)."""
    return re.sub(r"<(script|style)\b.*?</\1>", "", html, flags=re.S)


def test_comparison_page_order_and_words(tmp_path: Path) -> None:
    html = render(tmp_path)
    page = visible(html)
    assert f"<title>Comparing 2 runs: {TITLE}</title>" in html
    order = [page.find(f'id="{s}"') for s in ("power", "strength", "hits", "agreement")]
    order += [page.find('id="pairs"'), page.find('id="details"')]
    assert -1 not in order and order == sorted(order)  # power curve first, details last
    assert "fig-neff" in page and "Effective sample size relative to full" in page
    assert 'data-k="zb"' in page  # z against the baseline run
    text = re.sub(r'href="[^"]*"', "", page)  # file links may name summary files
    for word in ("strateg", "Strateg", "Dataset", "dataset", "default", "Raw vs RINT"):
        assert word not in text, word
    assert 'id="rint"' not in page and 'href="index.html"' not in page
    assert "Baseline full" in page


def test_hits_rows_are_runs_with_three_levels(tmp_path: Path) -> None:
    page = visible(render(tmp_path))
    hits = page[page.find('id="hits"') : page.find('id="agreement"')]
    heads = re.findall(r'<span class="hh"[^>]*>([^<]+)</span>', hits)
    assert heads == ["Same sign", "p &lt; 0.05", "p &lt; 1e-5"]
    assert "FDR" not in hits and "Bonferroni" not in hits
    rows = re.findall(r'<span class="hn">.*?</span>([^<]+)</span>', hits)
    assert rows == ["full", "none"]
    full = re.findall(r"<b>(\d+)</b>", hits.split("none</span>")[0])
    assert full == ["5", "5", "3"]  # sign, p < 0.05, p < 1e-5 of the fixture's 'full'
    # FDR and Bonferroni stay in the run details
    details = page[page.find('id="details"') :]
    assert "Bonferroni" in details and "FDR" in details


def test_chart_data(tmp_path: Path) -> None:
    data = page_data(render(tmp_path))
    (P,) = data["panels"]
    assert [r["label"] for r in P["runs"]] == ["full", "none"] and P["baseline"] == "full"
    assert P["runs"][0]["color"] == series_color(0) and P["runs"][0]["dash"] == "solid"
    assert [c["label"] for c in P["cutoffs"]] == ["p < 0.05", "p < 1e-5"]
    grid = P["curve"]["grid"]
    assert grid[0] == 0 and grid == sorted(grid)
    full = P["curve"]["counts"]["full"]
    assert full[0] == 5 and all(a >= b for a, b in itertools.pairwise(full))
    i5 = grid.index(5.0)
    assert full[i5] == 3  # same count as the strict level
    # effective sample size: the baseline is 1 with no interval; the weaker run is below
    assert P["neff"][0]["v"] == pytest.approx(1.0) and P["neff"][0]["lo"] is None
    assert 0 < P["neff"][1]["v"] < 1 and P["neff"][1]["lo"] <= P["neff"][1]["hi"]
    assert P["nCommon"] == 4  # tested in both runs
    # distributions over the pairs tested in every run, the reference study on the same
    assert P["dist"]["z"][0]["n"] == 4 and P["dist"]["lp"][1]["n"] == 4
    assert P["ref"]["z"]["n"] == 4
    # agreement: one value per drawn pair and run; opposite sign coded 0
    assert len(P["pairs"]["trait"]) == len(P["vals"]["none"]["z"]) == 6
    assert 0 in P["vals"]["none"]["sign"]
    assert P["selected"] == "none"  # a comparison opens on a run other than the baseline
    # the "All pairs" table: one row per pair, per run -log10 p, sign and why not tested
    T = P["table"]
    assert T["n"] == len(T["trait"]) == 6 and set(T["lp"]) == {"full", "none"}
    assert 0 in T["sign"]["none"] and "variant not in genotypes" in T["why"]["full"]
    assert P["file"] == "cis_pairs.tsv"


def test_single_run_page(tmp_path: Path) -> None:
    html = render(tmp_path, single=True)
    page = visible(html)
    assert f"<title>{TITLE}</title>" in html
    assert "fig-neff" not in page and 'data-k="zb"' not in page
    assert 'class="conc-grid n3"' in page and "Run full" in page
    assert 'data-seg="conc"' not in page  # no run selector
    (P,) = page_data(html)["panels"]
    assert P["neff"] == [] and not P["comparison"]


def test_run_details_words_and_links(tmp_path: Path) -> None:
    page = visible(render(tmp_path, covariates={"full": "all", "none": "none"}))
    details = page[page.find('id="details"') :]
    assert "10 (all)" in details and "<td>none</td>" in details
    assert "Effect r vs reference" in details and "<dt>Test</dt>" in details
    # only the levels that appear nowhere else are defined; no raw summary columns
    assert "<dt>FDR</dt>" in details and "<dt>Bonferroni</dt>" in details
    assert "<dt>Nominal</dt>" not in details and "<dt>Strict</dt>" not in details
    assert "rel_neff" not in details and "mean_z2_common" not in details
    assert 'href="../cis/d/strategies.tsv"' in details


def test_power_counts_include_level_cutoffs() -> None:
    tested = pd.DataFrame(
        {
            "run": ["a"] * 4,
            "log10p": [0.5, 2.0, 6.0, 30.0],
            "sign_concordant": [1, 1, 0, 1],
        }
    )
    c = power_counts(tested, ["a"], [-np.log10(0.05), 5.0])
    g = c["grid"]
    assert 5.0 in g and any(abs(x - 1.301) < 1e-3 for x in g)
    n = dict(zip(g, c["counts"]["a"], strict=True))
    assert n[0.0] == 3 and n[5.0] == 1  # 6.0 has the opposite sign
    # p < 10^-t as the levels: a pair exactly at the threshold does not count
    tied = pd.DataFrame({"run": ["a"], "log10p": [5.0], "sign_concordant": [1]})
    n = dict(
        zip(
            *(lambda c: (c["grid"], c["counts"]["a"]))(power_counts(tied, ["a"], [5.0])),
            strict=True,
        )
    )
    assert n[5.0] == 0


def test_neff_interval() -> None:
    rng = np.random.default_rng(3)
    z0 = rng.normal(0, 4, 400)
    zmat = np.column_stack([z0, z0, 2 * z0])
    ci = neff_interval(zmat, reps=200)
    assert ci[0] == pytest.approx((1.0, 1.0)) and ci[1] == pytest.approx((1.0, 1.0))
    lo, hi = ci[2]
    assert lo < 4.2 and hi > 4.0 and lo > 3.5  # 4 (mean z^2 - 1) / (mean z^2 - 1) ~ 4.3
    assert neff_interval(np.empty((1, 2)))[0] != neff_interval(np.empty((1, 2)))[0]  # NaN


def test_run_colors_fold_extra_with_dashes(tmp_path: Path) -> None:
    c = run_colors([f"s{i}" for i in range(20)])
    assert c["s0"] == series_color(0) and c["s7"] == series_color(7)
    assert c["s8"] == LIGHT["muted"] and c["s9"] == LIGHT["ink2"]
    assert DASHES[0] != "solid" and run_dash(0) == "solid"
    # past the eighth run: 10 distinct (grey, dash) styles
    styles = [(c[f"s{i}"], run_dash(i)) for i in range(8, 18)]
    assert len(set(styles)) == 10


def test_rect_report_has_local_controls(tmp_path: Path) -> None:
    run_summarize(tmp_path, make_run(tmp_path))
    st_path = tmp_path / L.cis_summary("d", "strategies.tsv")
    st = pd.read_csv(st_path, sep="\t")
    st["n_control_cells"] = 100
    st["lambda_gc_control"] = 1.1
    st["n_control_p_lt_1e5"] = 0
    st.to_csv(st_path, sep="\t", index=False)
    report_cis(
        StepContext(
            name="report.cis",
            run_dir=tmp_path,
            params=report_params(engine="rect"),
            inputs=[],
            outputs=[],
        )
    )
    page = visible((tmp_path / L.report("cis")).read_text(encoding="utf-8"))
    assert "Local controls" in page and "tbl-controls" in page


def test_multi_dataset_ids_are_suffixed(tmp_path: Path) -> None:
    run_summarize(tmp_path, make_run(tmp_path))
    shutil.copytree(tmp_path / L.cis_dir("d"), tmp_path / L.cis_dir("e"))
    report_cis(
        StepContext(
            name="report.cis",
            run_dir=tmp_path,
            params=report_params(datasets=["d", "e"]),
            inputs=[],
            outputs=[],
        )
    )
    html = (tmp_path / L.report("cis")).read_text(encoding="utf-8")
    assert 'id="power-d"' in html and 'id="power-e"' in html
    ids = re.findall(r'<div class="plotly-fig" id="([^"]+)"', html)
    assert len(ids) == len(set(ids)) and len(ids) == 2 * 7  # power dist neff + 4 agreement


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_inline_scripts_parse(tmp_path: Path) -> None:
    html = render(tmp_path)
    for i, js in enumerate(re.findall(r"<script>(.*?)</script>", html, re.S)):
        f = tmp_path / f"s{i}.js"
        f.write_text(js, encoding="utf-8")
        proc = subprocess.run(
            ["node", "--check", str(f)], capture_output=True, text=True, check=False
        )
        assert proc.returncode == 0, proc.stderr


def test_index_link_only_when_run_has_index(tmp_path: Path) -> None:
    assert 'href="index.html"' not in render(tmp_path)  # single-analysis runs: no index
    report_cis(
        StepContext(
            name="report.cis",
            run_dir=tmp_path,
            params=report_params(has_index=True),
            inputs=[],
            outputs=[],
        )
    )
    html = (tmp_path / L.report("cis")).read_text(encoding="utf-8")
    assert '<a href="index.html">Index</a>' in html
