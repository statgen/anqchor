"""report.trans: page built from synthetic trans summaries (no qpgen needed)."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest
from test_trans_summarize import make_run, params

from anqchor.pipelines.layout import Layout as L
from anqchor.report import trans as R
from anqchor.steps import StepContext, StepError
from anqchor.steps.trans import summarize_trans

REPORT_PARAMS = {
    "project": "demo",
    "datasets": ["ds"],
    "reference": "ref",
    "reference_file": "pairs.tsv.gz",
    "covariates": "sex, PC1-5",
    "plotlyjs": "cdn",
}


def render(tmp_path: Path, **kw: object) -> tuple[Path, str]:
    run = make_run(tmp_path)
    summarize_trans(
        StepContext(
            name="trans.ds.summarize", run_dir=run, params=params(**kw), inputs=[], outputs=[]
        )
    )
    ctx = StepContext(
        name="report.trans", run_dir=run, inputs=[], outputs=[], params=dict(REPORT_PARAMS)
    )
    R.report_trans(ctx)
    return run, (run / L.report("trans")).read_text()


def page_data(html: str) -> dict:
    m = re.search(r"window\.ANQ_DATA = (\{.*?\});\n", html, re.S)
    assert m
    return json.loads(m.group(1))


def test_title_context_and_sections_in_order(tmp_path: Path) -> None:
    _, html = render(tmp_path)
    assert "<title>trans-QTL evaluation with anQChor</title>" in html
    assert "6 variants × 5 genes" in html and "29 variant–gene tests, 12 across chromosomes" in html
    assert "100 samples" in html and "Run myrun" in html and "Covariates 6 (sex, PC1-5)" in html
    assert "pairs.tsv.gz" in html
    order = [html.index(f'id="{s}"') for s in ("replication", "hotspots", "inflation", "details")]
    assert order == sorted(order)
    for word in ("strategy", "Strategy", "dataset"):
        assert word not in re.sub(r"<script.*?</script>", "", html, flags=re.S), word


def test_replication_rows_cis_and_trans(tmp_path: Path) -> None:
    _, html = render(tmp_path)
    assert "No trans pairs in the reference" not in html
    assert "1 other chromosome · 1 same, &gt; 1 Mb" in html
    _, html = render(tmp_path / "b", known_pairs=False)
    assert "No cis pairs in the reference" in html and "No trans pairs in the reference" in html
    assert "None left out: the reference has no trans pairs." in html


def test_page_data_for_hotspots_and_qq(tmp_path: Path) -> None:
    _, html = render(tmp_path)
    (p,) = page_data(html)["panels"]
    h = p["hot"]
    assert [t["key"] for t in h["thresholds"]][:2] == ["bonferroni", "p1e-4"]
    assert h["default"] == {"threshold": "p1e-5", "k": 2, "n_hotspots": 1}
    assert h["bins"]["chrom"] == ["1", "1", "2", "2"] and h["bins"]["n"]["p1e-5"] == [3, 0, 0, 1]
    # bins along the genome: chromosome blocks in order, x increasing
    assert h["bins"]["x"] == sorted(h["bins"]["x"]) and [c["name"] for c in h["chroms"]] == [
        "1",
        "2",
    ]
    # genes per bin, strongest first, short names
    assert h["bins"]["genes"][0][0][:2] == ["gC", 8.5]
    assert p["qq"]["n"] == 4 and len(p["qq"]["exp"]) == 4
    assert "Left out: 8 tests of 3 variants within 1 Mb of the 2 reference trans pairs" in html


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_inline_scripts_parse(tmp_path: Path) -> None:
    _, html = render(tmp_path)
    scripts = re.findall(r"<script>(.*?)</script>", html, re.S)
    assert scripts
    for i, js in enumerate(scripts):
        f = tmp_path / f"s{i}.js"
        f.write_text(js)
        proc = subprocess.run(
            ["node", "--check", str(f)], capture_output=True, text=True, check=False
        )
        assert proc.returncode == 0, proc.stderr


def test_helpers() -> None:
    assert R.lambda_status(1.02) == "good" and R.lambda_status(1.08) == "warning"
    assert R.lambda_status(1.6) == "serious" and R.lambda_status(0.8) == "serious"
    assert R.lambda_status(None) == ""
    assert R.short_trait("ENSG1.2:APP") == "APP" and R.short_trait("ENSG1") == "ENSG1"
    assert (
        R.fmt_bp(1_000_000) == "1 Mb"
        and R.fmt_bp(500_000) == "500 kb"
        and R.fmt_bp(2_500_000) == "2.5 Mb"
    )


def test_old_summaries_ask_for_a_rerun(tmp_path: Path) -> None:
    run, _ = render(tmp_path)
    f = run / L.trans_summary("ds", "summary.json")
    f.write_text(json.dumps({"schema_version": 1}))
    with pytest.raises(StepError, match="re-run 'anqchor trans run'"):
        R.load_dataset(run, "ds")


def test_index_link_only_when_run_has_index(tmp_path: Path) -> None:
    run, html = render(tmp_path)
    assert 'href="index.html"' not in html  # single-analysis runs write no index page
    ctx = StepContext(
        name="report.trans",
        run_dir=run,
        inputs=[],
        outputs=[],
        params={**REPORT_PARAMS, "has_index": True},
    )
    R.report_trans(ctx)
    assert '<a href="index.html">Index</a>' in (run / L.report("trans")).read_text()
