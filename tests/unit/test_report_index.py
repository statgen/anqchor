"""report_index: scanning report/, reading summaries defensively, and the page."""

from __future__ import annotations

import json
from pathlib import Path

from anqchor.pipelines.layout import Layout as L
from anqchor.report.index import (
    build_index,
    facts,
    first,
    fmt_value,
    read_summary,
    report_index,
    scan_reports,
    summary_datasets,
)
from anqchor.report.render import PLOTLY_CDN, render
from anqchor.steps import StepContext
from anqchor.util.log import get_logger

MATCH = {
    "verdict": "728 of 729 mapped samples confirmed; 1 mismatch (NA11919 best matches HG00237)",
    "verdict_status": "critical",
    "group_counts": {"PASS": 725, "WARN": 3, "FAIL": 1, "NO_TEST": 0},
    "n_mapped": 729,
    "n_pheno": 729,
    "n_rescue_candidates": 0,
    "n_related": 46,
    "lambda_mode": "independent",
    "lambda_used": None,
    "warnings": [],
}
CIS = {
    "recommended": {
        "strategy": "peer60",
        "transform": "raw",
        "label": "peer60/raw",
        "key_level": "strict",
        "reason": "peer60/raw replicates the most reference pairs at the strict level.",
    },
    "per_strategy": {
        "peer60/raw": {"n_rep_strict": 455, "n_tested": 614},
        "default/raw": {"n_rep_strict": 400, "n_tested": 614},
    },
    "n_pairs_mapped": 660,
    "n_ref_pairs": 660,
    "strategies": [{"name": "peer60"}, {"name": "default"}],
    "warnings": ["w"],
}
TRANS = {
    "recommended": {
        "strategy": "peer60",
        "transform": "raw",
        "lambda_gc_trans": 1.11,
        "known_replicated": 549,
    },
    "strategies": [
        {
            "strategy": "default",
            "transform": "raw",
            "lambda_gc_trans": 1.4,
            "n_hotspots": 3,
            "n_sig_trans_bonferroni": 9,
        },
        {
            "strategy": "peer60",
            "transform": "raw",
            "lambda_gc_trans": 1.11,
            "n_hotspots": 0,
            "n_sig_trans_bonferroni": 1,
        },
    ],
    "n_variants": 588,
    "n_traits": 806,
    "n_known_pairs": 660,
}


def _run(tmp_path: Path) -> Path:
    run = tmp_path / "run"
    for p, ds, data in (("match", "rna", MATCH), ("cis", "rna", CIS), ("trans", "rna", TRANS)):
        d = run / p / ds
        d.mkdir(parents=True)
        (d / "summary.json").write_text(json.dumps(data))
    (run / "trans" / "broken").mkdir()
    (run / "trans" / "broken" / "summary.json").write_text("{not json")
    (run / "report").mkdir()
    for n in ("trans", "match", "custom"):
        (run / "report" / f"{n}.html").write_text("<html></html>")
    return run


def test_first_and_fmt_value() -> None:
    d = {"a": None, "b": {"c": 2}, "e": 0}
    assert first(d, "a", "b.c") == 2
    assert first(d, "x", "b.x") is None
    assert first(d, "e") == 0
    assert fmt_value(True) == "yes" and fmt_value(1234567) == "1,234,567"
    assert fmt_value(1.1142) == "1.11" and fmt_value(float("nan")) == "NA"
    assert fmt_value(2e-8) == "2.00e-08" and fmt_value(250.4) == "250"
    assert fmt_value({"PASS": 1, "FAIL": 2}) == "PASS 1, FAIL 2"
    assert fmt_value(list(range(8))).endswith("…")


def test_facts_per_pipeline() -> None:
    m = dict(facts("match", MATCH))
    assert m["groups"] == "PASS 725, WARN 3, FAIL 1, NO_TEST 0"
    assert m["mapped samples"] == "729 of 729 (100%)" and m["λ"] == "none (independent traits)"
    assert m["possible relatives"] == "46" and "warnings" not in m
    c = dict(facts("cis", CIS))
    assert c["recommended"] == "peer60/raw"
    assert c["replicated (strict)"] == "455 of 614 (74%)"
    assert c["strategies"] == "2" and c["warnings"] == "1"
    # cis summary schema 2: per_run and runs
    c2 = dict(
        facts(
            "cis",
            {
                "recommended": {"run": "peer60", "label": "peer60", "key_level": "strict"},
                "per_run": {"peer60": {"n_rep_strict": 421, "n_tested": 479}},
                "runs": [{"label": "peer60", "rint": False}],
            },
        )
    )
    assert c2["recommended"] == "peer60" and c2["replicated (strict)"] == "421 of 479 (88%)"
    assert c2["runs"] == "1"
    t = dict(facts("trans", TRANS))
    assert t["recommended"] == "peer60/raw"
    assert t["λ_GC trans"] == "1.11" and t["hotspots"] == "0"  # the recommended run, not the first
    assert t["Bonferroni trans cells"] == "1"
    assert t["variants × traits"] == "588 × 806"
    assert t["known pairs replicated"] == "549 of 660 (83%)"


def test_facts_tolerate_missing_and_other_keys() -> None:
    assert facts("cis", {}) == []
    assert facts("trans", {}) == []
    assert facts("match", {}) == []
    # alternative top-level keys
    assert dict(facts("trans", {"lambda_gc": 1.02, "n_hotspots": 2})) == {
        "λ_GC trans": "1.02",
        "hotspots": "2",
    }
    assert dict(facts("cis", {"recommended_strategy": "peer15"}))["recommended"] == "peer15"
    # recommendation without a label
    assert (
        dict(facts("cis", {"recommended": {"strategy": "s", "transform": "rint"}}))["recommended"]
        == "s/rint"
    )
    # unknown pipeline: only generic facts
    assert facts("other", {"warnings": ["a", "b"]}) == [("warnings", "2")]


def test_read_summary_problems(tmp_path: Path) -> None:
    assert read_summary(tmp_path / "nope.json") == (None, "summary not written yet")
    bad = tmp_path / "bad.json"
    bad.write_text("{x")
    s, why = read_summary(bad)
    assert s is None and "unreadable" in why
    arr = tmp_path / "arr.json"
    arr.write_text("[1]")
    assert read_summary(arr) == (None, "summary is not a JSON object")


def test_scan_reports_and_datasets(tmp_path: Path) -> None:
    run = _run(tmp_path)
    assert scan_reports(run / "report") == ["match", "trans", "custom"]
    assert scan_reports(tmp_path / "missing") == []
    assert summary_datasets(run, "trans", ["rna"]) == ["rna", "broken"]
    assert summary_datasets(run, "cis", ["other"]) == ["other", "rna"]


def test_build_index_page(tmp_path: Path) -> None:
    run = _run(tmp_path)
    page = build_index(run, "proj", ["match", "cis"], {"match": ["rna"], "cis": ["rna"]})
    html = render(page, "cdn")
    ids = [s.id for s in page.sections]
    assert ids == ["reports", "pipe-match", "pipe-cis", "pipe-trans"]  # trans found in report/
    assert 'href="match.html#ds-rna"' in html and 'href="custom.html"' in html
    assert "Requested but not rendered yet: cis." in html
    assert "NA11919 best matches HG00237" in html and "status-critical" in html
    assert "peer60/raw replicates the most reference pairs" in html  # cis reason as text
    assert "unreadable summary" in html  # trans/broken
    # cis has a summary but no report: no link on its card
    cis = next(s for s in page.sections if s.id == "pipe-cis")
    assert "<a " not in cis.body


def test_report_index_step(tmp_path: Path) -> None:
    run = _run(tmp_path)
    ctx = StepContext(
        name="report.index",
        run_dir=run,
        params={
            "project": "proj",
            "pipelines": ["match"],
            "datasets": {"match": ["rna"]},
            "plotlyjs": "embed",
        },
        inputs=[],
        outputs=[L.report("index")],
        log=get_logger("test.index"),
    )
    report_index(ctx)
    html = (run / L.report("index")).read_text(encoding="utf-8")
    assert html.startswith("<!doctype html>")
    assert PLOTLY_CDN not in html
    assert "Plotly.newPlot" not in html.split("</head>")[0]  # plotly.js is not embedded either
    assert len(html) < 200_000
    assert "proj: anqchor reports" in html


def test_trans_facts_of_one_run_summary() -> None:
    """trans summaries since schema 2 describe one run (replication by class, hotspots,
    cross-chromosome calibration)."""
    t = dict(
        facts(
            "trans",
            {
                "schema_version": 2,
                "run": "peer60",
                "inflation": {"lambda_gc": 1.002},
                "hotspots": {"default": {"threshold": "p1e-5", "k": 5, "n_hotspots": 0}},
                "replication": {
                    "cis": {"n_pairs": 367, "n_tested": 367, "n_rep_strict": 251},
                    "trans": {"n_pairs": 0, "n_tested": 0, "n_rep_strict": 0},
                },
                "n_variants_tested": 346,
                "n_traits_tested": 367,
            },
        )
    )
    assert t["run"] == "peer60" and t["λ_GC cross-chromosome"] == "1.00"
    assert t["hotspots (≥ 5 genes)"] == "0"
    assert t["cis pairs replicated (p < 1e-5)"] == "251 of 367 (68%)"
    assert "trans pairs replicated (p < 1e-5)" not in t
    assert t["variants × genes"] == "346 × 367"
