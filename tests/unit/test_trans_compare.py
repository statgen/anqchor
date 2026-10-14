"""steps.trans_compare and the comparison page of report.trans, on synthetic trans runs."""

from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from test_trans_report import page_data
from test_trans_summarize import LP, TRAITS, make_run, params
from typer.testing import CliRunner

from anqchor.cli.main import app
from anqchor.pipelines.layout import Layout as L
from anqchor.report import trans as R
from anqchor.steps import StepContext, StepError
from anqchor.steps import trans_compare as TC
from anqchor.steps.compare import comparison_kind, default_labels
from anqchor.steps.trans import summarize_trans

runner = CliRunner()


def resolved_yaml(run: Path) -> str:
    return (
        "version: 1\nproject: demo\n"
        f"out_dir: {run}\n"
        "genotypes: {pgen_list: /data/geno.list.tsv}\n"
        "references: {ref: {kind: cis, pairs: /data/pairs.tsv.gz}}\n"
        "datasets: {ds: {pheno: /data/pheno.tsv.gz}}\n"
        "pipelines: {trans: {enabled: true}}\n"
    )


def trans_run(
    base: Path, name: str, lp: np.ndarray | None = None, ref: str = "a", **kw: Any
) -> Path:
    """A finished synthetic trans run directory ``base/name/run`` labelled ``name``."""
    run = make_run(base / name, lp)
    summarize_trans(
        StepContext(
            name="trans.ds.summarize",
            run_dir=run,
            params=params(**{"run": name, **kw}),
            inputs=[],
            outputs=[],
        )
    )
    (run / L.resolved_config).write_text(resolved_yaml(run))
    f = run / L.ref_pairs("ref")
    f.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(f, "wt") as fh:
        fh.write(f"trait\tvariant\tbeta\tse\tpval\n{ref}\n")
    return run


def test_compare_writes_tables_settings_and_manifest(tmp_path: Path) -> None:
    a, b = trans_run(tmp_path, "one"), trans_run(tmp_path, "two", rint=True)
    out = tmp_path / "cmp"
    s = TC.compare_runs([a, b], out, labels=["nocov", "peer"], argv=["anqchor", "trans", "compare"])
    assert (
        s["kind"] == "comparison" and s["pipeline"] == "trans" and s["labels"] == ["nocov", "peer"]
    )
    for name in TC.TABLES:
        assert (out / L.compare_table(name)).exists(), name
    assert comparison_kind(out) == "trans"
    m = json.loads((out / L.compare_manifest).read_text())
    assert m["kind"] == "trans" and m["labels"] == ["nocov", "peer"] and m["argv"][1] == "trans"
    assert (out / L.compare_settings("peer")).read_text() == (b / L.resolved_config).read_text()
    # identical rectangles and references: no comparability warning
    c = s["comparability"]
    assert c["same_reference"] and c["same_rectangle"] and c["same_cells"]
    assert c["rint"] == {"nocov": False, "peer": True}
    assert not [w for w in s["warnings"] if "different" in w]
    runs = pd.read_csv(out / "runs.tsv", sep="\t")
    assert list(runs["label"]) == ["nocov", "peer"] and list(runs["rint"]) == [0, 1]
    r0 = runs.iloc[0]
    # reference pairs: gA-v1 cis tested (the absent-variant gA pair is not), gC-v1 and gD-v4 trans
    assert (r0["cis_n_pairs"], r0["cis_n_tested"], r0["trans_n_pairs"], r0["trans_n_tested"]) == (
        2,
        1,
        2,
        2,
    )
    assert (
        r0["trans_n_common"] == 2 and r0["trans_n_rep_sign"] == r0["trans_n_rep_sign_common"] == 1
    )
    assert r0["n_cells"] == 29 and r0["n_cells_cross_chrom"] == 12
    assert {"lambda_gc", "observed_p1e-3", "expected_p1e-5", "n_hotspots", "bonferroni_p"} <= set(
        runs
    )
    grid = pd.read_csv(out / "hotspot_grid.tsv", sep="\t")
    per = s["per_run"]["nocov"]["hotspots"]
    assert len(grid) == 2 * len(per["thresholds"]) * len(per["k"])
    row = grid[(grid["run"] == "nocov") & (grid["threshold"] == "p1e-5") & (grid["k"] == 2)]
    assert int(row["n_hotspots"].iloc[0]) == per["default"]["n_hotspots"] == 1
    for name in ("hotspots.tsv", "hotspot_genes.tsv", "qq.tsv", "reference_pairs.tsv"):
        t = pd.read_csv(out / name, sep="\t")
        assert t.columns[0] == "run" and set(t["run"]) == {"nocov", "peer"}, name
    rp = pd.read_csv(out / "reference_pairs.tsv", sep="\t")
    assert rp["tested_in_all"].sum() == 2 * 3


def test_labels(tmp_path: Path) -> None:
    a, b = trans_run(tmp_path, "x"), trans_run(tmp_path, "y")
    s = TC.compare_runs([a, b], tmp_path / "c1")
    assert s["labels"] == default_labels([a, b]) == ["x_run", "y_run"]  # both named 'run'
    with pytest.raises(StepError, match="1 --label value"):
        TC.compare_runs([a, b], tmp_path / "c2", labels=["one"])
    with pytest.raises(StepError, match="duplicate run label"):
        TC.compare_runs([a, b], tmp_path / "c3", labels=["s", "s"])


def test_different_cells_and_reference_are_warned(tmp_path: Path) -> None:
    lp = LP.copy()
    lp[0, TRAITS.index("gC")] = np.nan  # the trans reference pair gC-v1 is not tested in run b
    a, b = trans_run(tmp_path, "a"), trans_run(tmp_path, "b", lp=lp, ref="other")
    s = TC.compare_runs([a, b], tmp_path / "cmp")
    w = " | ".join(s["warnings"])
    assert "different references" in w
    assert "the runs tested different variant–gene tests (a_run 29, b_run 28)" in w
    assert "different trans reference pairs: 1 tested in every run, up to 2" in w
    runs = pd.read_csv(tmp_path / "cmp" / "runs.tsv", sep="\t")
    assert list(runs["trans_n_tested"]) == [2, 1] and list(runs["trans_n_common"]) == [1, 1]
    assert not s["comparability"]["same_reference"] and not s["comparability"]["same_cells"]


def test_not_a_trans_run_and_old_schema(tmp_path: Path) -> None:
    with pytest.raises(StepError, match="not an anqchor run directory"):
        TC.compare_runs([tmp_path], tmp_path / "c")
    a = trans_run(tmp_path, "a")
    f = a / L.trans_summary("ds", "summary.json")
    f.write_text(json.dumps({**json.loads(f.read_text()), "schema_version": 1}))
    with pytest.raises(StepError, match="re-run 'anqchor trans run"):
        TC.compare_runs([a], tmp_path / "c")
    f.unlink()
    with pytest.raises(StepError, match="is not a trans run directory"):
        TC.compare_runs([a], tmp_path / "c")


def render(tmp_path: Path) -> tuple[Path, str]:
    a, b = trans_run(tmp_path, "nocov"), trans_run(tmp_path, "peer60", rint=True)
    out = tmp_path / "cmp"
    TC.compare_runs([a, b], out)
    ctx = StepContext(
        name="report.trans",
        run_dir=out,
        inputs=[],
        outputs=[],
        params={"comparison": True, "plotlyjs": "cdn"},
    )
    R.report_trans(ctx)
    return out, (out / L.report("trans")).read_text()


def test_comparison_page(tmp_path: Path) -> None:
    _, html = render(tmp_path)
    assert "<title>Comparing 2 runs: trans-QTL evaluation with anQChor</title>" in html
    assert (
        "2 runs · 6 variants × 5 genes · 29 variant–gene tests, 12 across chromosomes · 100 samples"
        in html
    )
    order = [
        html.index(f'id="{s}-cmp"') for s in ("replication", "hotspots", "inflation", "details")
    ]
    assert order == sorted(order)
    visible = re.sub(r"<(script|style).*?</(script|style)>", "", html, flags=re.S)
    for word in ("strategy", "Strategy", "dataset", "default"):
        assert word not in visible, word
    assert "No trans pairs in the reference" not in html
    d = page_data(html)
    assert d["panels"] == []
    c = d["comparison"]
    assert [r["label"] for r in c["runs"]] == ["nocov_run", "peer60_run"]
    assert c["default"] == {"threshold": "p1e-5", "k": 2}
    assert [t["key"] for t in c["thresholds"]][:2] == ["bonferroni", "p1e-4"]
    assert c["runs"][0]["hot"]["bins"]["n"]["p1e-5"] == [3, 0, 0, 1]
    assert c["runs"][1]["qq"]["n"] == 4 and c["runs"][0]["color"] != c["runs"][1]["color"]
    assert 'href="../settings/nocov_run.anqchor.yaml"' in html


def test_trans_report_cli_renders_a_comparison_and_cis_report_refuses_it(tmp_path: Path) -> None:
    out, _ = render(tmp_path)
    (out / L.report("trans")).unlink()
    r = runner.invoke(app, ["trans", "report", "--out", str(out), "--plotlyjs", "cdn"])
    assert r.exit_code == 0, r.output
    assert (out / L.report("trans")).exists()
    r = runner.invoke(app, ["cis", "report", "--out", str(out)])
    assert r.exit_code == 2 and "trans comparison" in r.output


def test_trans_compare_cli(tmp_path: Path) -> None:
    a, b = trans_run(tmp_path, "a"), trans_run(tmp_path, "b")
    out = tmp_path / "cmp"
    r = runner.invoke(
        app,
        [
            "trans",
            "compare",
            str(a),
            str(b),
            "--out",
            str(out),
            "--label",
            "x",
            "--label",
            "y",
            "--no-report",
        ],
    )
    assert r.exit_code == 0, r.output
    assert "compared 2 run(s): x, y" in r.output
    assert not (out / L.report("trans")).exists()


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_comparison_scripts_parse(tmp_path: Path) -> None:
    _, html = render(tmp_path)
    for i, js in enumerate(re.findall(r"<script>(.*?)</script>", html, re.S)):
        f = tmp_path / f"s{i}.js"
        f.write_text(js)
        proc = subprocess.run(
            ["node", "--check", str(f)], capture_output=True, text=True, check=False
        )
        assert proc.returncode == 0, proc.stderr


def test_missing_table_is_a_clean_error(tmp_path: Path) -> None:
    a = trans_run(tmp_path, "a")
    (a / L.trans_summary("ds", "hotspot_genes.tsv")).unlink()
    with pytest.raises(StepError, match=r"missing .*hotspot_genes\.tsv"):
        TC.compare_runs([a], tmp_path / "c")


def test_out_dir_guards_and_stale_files(tmp_path: Path) -> None:
    a, b = trans_run(tmp_path, "a"), trans_run(tmp_path, "b")
    with pytest.raises(StepError, match="is a run directory"):
        TC.compare_runs([a, b], a)
    with pytest.raises(StepError, match="given more than once"):
        TC.compare_runs([a, a], tmp_path / "c")
    cis = tmp_path / "ciscmp"
    cis.mkdir()
    (cis / L.compare_manifest).write_text(json.dumps({"kind": "cis"}))
    with pytest.raises(StepError, match="holds a cis comparison"):
        TC.compare_runs([a, b], cis)
    out = tmp_path / "c"
    TC.compare_runs([a, b], out, labels=["x", "y"])
    report = out / L.report("trans")
    report.parent.mkdir(parents=True)
    report.write_text("old")
    TC.compare_runs([a], out, labels=["x"])
    assert not (out / L.compare_settings("y")).exists() and (out / L.compare_settings("x")).exists()
    assert not report.exists()


def test_shared_hotspot_count_and_opening_warning(tmp_path: Path) -> None:
    a, b = trans_run(tmp_path, "a"), trans_run(tmp_path, "b", hotspot_min_traits=3)
    s = TC.compare_runs([a, b], tmp_path / "c")
    runs = pd.read_csv(tmp_path / "c" / "runs.tsv", sep="\t")
    assert "n_hotspots_shared" in runs and runs["n_hotspots_shared"].notna().all()
    assert any("different thresholds or minimum genes" in w for w in s["warnings"])


def test_default_labels_keep_symlink_names(tmp_path: Path) -> None:
    a = trans_run(tmp_path, "a")
    link = tmp_path / "linked"
    link.symlink_to(a)
    assert default_labels([link]) == ["linked"]


def test_reference_pairs_table_matches_by_reference_trait_and_variant() -> None:
    pairs = pd.DataFrame(
        {
            "run": ["a", "b", "a"],
            "trait": ["g1", "g1_renamed", "g2"],
            "ref_trait": ["G1", "G1", "G2"],
            "variant": ["1:1:A:C", "1:1:A:C", "2:5:G:T"],
            "pair_class": ["cis", "trans", "cis"],
            "log10p": [3.14159, 8.0, 12.0],
            "sign_concordant": [1, 0, 1],
        }
    )
    wide = R._pairs_wide(pairs, ["a", "b"])
    assert len(wide) == 2 and list(wide["Gene"]) == ["g2", "g1"]
    row = wide[wide["Gene"] == "g1"].iloc[0]
    assert row["Class"] == "cis/trans" and row["a"] == 3.1 and row["b"] == -8.0
