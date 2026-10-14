"""report_match: page content, embedded data, wording, and the browser-side rules."""

from __future__ import annotations

import inspect
import json
import re
import shutil
import subprocess
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

import pytest
from test_match_helpers import summarized_run

from anqchor.pipelines.layout import Layout as L
from anqchor.report import render as render_mod
from anqchor.report.match import (
    MATCH_RULES_JS,
    ROW_FIELDS,
    TITLE,
    _ds_labels,
    build_page,
    category_of,
    display_flags,
    load_dataset,
    report_match,
    run_categories,
    scoring_detail,
    scoring_text,
    tile_counts,
)
from anqchor.report.render import render
from anqchor.steps import StepContext, StepError
from anqchor.steps.match import (
    CLASSIFY_TEST_VECTORS,
    classify,
    close_candidates,
    group_of,
    is_duplicate,
    is_related,
)
from anqchor.util.log import get_logger

needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")


def _ctx(run: Path, datasets: list[str], **kw: Any) -> StepContext:
    p = {
        "datasets": datasets,
        "omics": dict.fromkeys(datasets, "rna"),
        "project": "proj",
        "z_threshold": 1.96,
        "z_diff": 2.0,
        "plotlyjs": "cdn",
    }
    p.update(kw)
    return StepContext(
        name="report.match",
        run_dir=run,
        params=p,
        inputs=[],
        outputs=[],
        log=get_logger("test.report"),
    )


def _anq_data(html: str) -> dict[str, Any]:
    m = re.search(r"window\.ANQ_DATA = (\{.*?\});\n", html, re.S)
    assert m
    return json.loads(m.group(1))


class _Text(HTMLParser):
    """Visible text: text nodes outside <script>, <style> and <title>."""

    def __init__(self) -> None:
        super().__init__()
        self.skip = 0
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: Any) -> None:
        if tag in ("script", "style"):
            self.skip += 1

    def handle_endtag(self, tag: str) -> None:
        if tag in ("script", "style"):
            self.skip -= 1

    def handle_data(self, data: str) -> None:
        if not self.skip:
            self.parts.append(data)


def visible_text(html: str, drop_codes_table: bool = True) -> str:
    if drop_codes_table:  # the one place where qpgen codes are listed on purpose
        html = re.sub(r'<table class="codes">.*?</table>', "", html, flags=re.S)
        html = re.sub(
            r"<details><summary>qpgen (command|log excerpt)</summary>.*?</details>",
            "",
            html,
            flags=re.S,
        )
    p = _Text()
    p.feed(html)
    return " ".join(p.parts)


@pytest.fixture
def two_dataset_run(tmp_path: Path) -> Path:
    summarized_run(tmp_path, "rna")
    summarized_run(tmp_path, "prot.v2")
    return tmp_path


def test_report_step_writes_page(two_dataset_run: Path) -> None:
    run = two_dataset_run
    report_match(_ctx(run, ["rna", "prot.v2"]))
    html = (run / L.report("match")).read_text()
    assert TITLE == "Sample identity matching with anQChor"
    assert f"<title>{TITLE}</title>" in html and f"<h1>{TITLE}</h1>" in html
    assert '<div class="brand">' not in html
    # one panel per dataset; the second starts hidden; the switch is labelled with the dataset
    # names because both datasets have the omics name "rna"
    assert 'id="panel-rna"' in html and 'id="panel-prot-v2" data-ds="prot-v2" hidden' in html
    assert 'data-ds-set="rna"' in html and ">prot.v2</button>" in html
    assert 'id="anq-zt-r"' in html and 'value="1.96"' in html and 'id="anq-copy"' in html
    data = _anq_data(html)
    assert data["fields"] == ROW_FIELDS and data["config"] == {"z_threshold": 1.96, "z_diff": 2.0}
    assert [d["name"] for d in data["datasets"]] == ["rna", "prot.v2"]
    assert len(data["datasets"][0]["rows"]) == 16
    assert data["datasets"][1]["fig"] == "fig-match-prot-v2"
    for sid in ("check-rna", "chart-rna", "dups-rna", "all-rna", "fig-weights-rna"):
        assert f'id="{sid}"' in html
    # no section nav, no index link (a single-analysis run writes no report/index.html)
    assert 'class="toc"' not in html and 'href="index.html"' not in html


def test_index_link_only_when_run_has_index(tmp_path: Path) -> None:
    summarized_run(tmp_path, "rna")
    report_match(_ctx(tmp_path, ["rna"], has_index=True))
    assert '<a href="index.html">Index</a>' in (tmp_path / L.report("match")).read_text()


def test_dashboard_tiles_at_run_thresholds(tmp_path: Path) -> None:
    summarized_run(tmp_path, "rna")
    report_match(_ctx(tmp_path, ["rna"]))
    html = (tmp_path / L.report("match")).read_text()
    tiles = dict(re.findall(r'data-tile="(\w+)".*?<span class="t-value">([\d,]+)</span>', html))
    # S1-S6 strict, S10 lenient, S7-S9 another sample, S11 own ID too weak, S13 (not scored
    # by qpgen) and the sample without observed traits not scored, S12 matches several
    # (two genotype IDs about equally, neither its own: counted on the two-IDs tile);
    # U1-U3 have no assigned ID (U1: new ID found)
    assert tiles == {
        "examined": "16",
        "assigned": "13",
        "strict": "6",
        "lenient": "7",
        "different": "3",
        "twoids": "1",
        "several": "1",
        "none": "1",
        "nodata": "1",
        "unassigned": "3",
    }
    assert ">1 new ID found · 1 no clear match · 1 not scored<" in html and ">strict + 1<" in html
    assert ">S11<" in html and "3 without" in html and "not even lenient" not in html
    assert (
        'data-tile="examined"' in html
    )  # requested by the user, even though the context line repeats it
    assert not re.search(r'data-tile="twoids"[^>]*hidden>', html)
    # the "no covariates" warning reaches the page
    assert "1 sample not scored: no covariates" in html
    # one dataset: no switch and no dataset name in the visible text
    assert '<div class="seg"' not in html
    assert "rna" not in visible_text(html).replace(
        "rna: ", ""
    )  # only the omics prefix of the context line


def test_banned_words_absent_from_visible_text(tmp_path: Path) -> None:
    summarized_run(tmp_path, "dataset")
    report_match(_ctx(tmp_path, ["dataset"], omics={"dataset": "omics"}))
    html = (tmp_path / L.report("match")).read_text()
    text = visible_text(html)
    for word in (
        "PASS",
        "WARN",
        "FAIL",
        "NO_TEST",
        "SELF_BEST",
        "SINGLE_NEW_BEST",
        "rescue",
        "verdict",
        "mapped",
        "z_threshold",
        "z_diff",
        "dataset",
        "omics:",
        "Possible relatives",
    ):
        assert not re.search(rf"\b{re.escape(word)}\b", text), word
    assert "Dataset dataset" not in html
    assert "Duplicates" in text and "Samples to check" in text and "Run details" in text
    assert "Possible duplicates" not in text and "relatives" not in text and "Margin" not in html
    # the codes are listed in the run details
    full = visible_text(html, drop_codes_table=False)
    assert "SELF_BEST" in full and "MULTI_NEW_BEST" in full


def test_scoring_chip_and_hint(tmp_path: Path) -> None:
    summarized_run(tmp_path, "rna")
    d = load_dataset(_ctx(tmp_path, ["rna"]), "rna", "rna")
    html = render(build_page("p", 1.96, 2.0, [d]), "cdn")
    assert "Scoring: independent traits" in html and "try Mahalanobis scoring" in html
    # no link: the docs site may not be published, and a local report should not depend on it
    assert "<code>--mahalanobis</code>" in html and "identity-matching/#mahalanobis" not in html
    d["summary"].update(mahalanobis=True, lambda_mode="auto", lambda_used=0.5795)
    html = render(build_page("p", 1.96, 2.0, [d]), "cdn")
    assert ">Scoring: Mahalanobis<" in html and "try Mahalanobis scoring" not in html
    # λ only in the run details
    assert "<dd>Mahalanobis (λ 0.58, auto)</dd>" in html
    assert scoring_detail({"mahalanobis": True, "lambda_mode": "fixed", "lambda_used": 0.5}) == (
        "Mahalanobis (λ 0.50, fixed)"
    )
    assert scoring_text({"mahalanobis": False}) == "Scoring: independent traits"


def test_disagreement_line_is_computed_by_the_page(tmp_path: Path) -> None:
    summarized_run(tmp_path, "rna")
    d = load_dataset(_ctx(tmp_path, ["rna"]), "rna", "rna")
    html = render(build_page("p", 1.96, 2.0, [d]), "cdn")
    # an empty, hidden line the page script fills when its statuses differ from the saved ones
    assert '<p class="m-alert js-disagree" role="note" hidden></p>' in html
    assert "the page's result differs from the saved result" in html


def test_report_missing_summary_is_step_error(tmp_path: Path) -> None:
    with pytest.raises(StepError, match="cannot read the match summaries of dataset nope"):
        report_match(_ctx(tmp_path, ["nope"]))


def test_default_table_page_size_is_500() -> None:
    assert inspect.signature(render_mod.table).parameters["page_size"].default == 500
    assert '[[100, "100"], [500, "500"], [0, "All"]]' in render_mod._static("match_page.js")
    assert "ds.pageSize = window.innerWidth < 600 ? 100 : 500" in render_mod._static(
        "match_page.js"
    )


def test_tile_partition_python(two_dataset_run: Path) -> None:
    d = load_dataset(_ctx(two_dataset_run, ["rna"]), "rna", "rna")
    cats, assigned = run_categories(d["samples"])
    t = tile_counts(cats, assigned)
    assert t["assigned"] == _assigned_sum(t)
    assert t["examined"] == t["assigned"] + t["unassigned"]
    assert (
        t["unassigned"] == t["newid"] + t["taken"] + t["un_twoids"] + t["unclear"] + t["un_nodata"]
    )


def _assigned_sum(t: dict[str, int]) -> int:
    return t["lenient"] + t["different"] + t["twoids"] + t["several"] + t["none"] + t["nodata"]


@pytest.mark.parametrize(("zt", "zd"), [(1.96, 2.0), (0.5, 0.3), (3.0, 1.0), (6.0, 4.0)])
def test_tile_partition_on_classifier_vectors(zt: float, zd: float) -> None:
    cats = []
    assigned = []
    for v in CLASSIFY_TEST_VECTORS:
        st = classify(
            mapped=v["mapped"],
            n_traits=v["n_traits"],
            rank_self=v["rank_self"],
            z_self=v["z_self"],
            z_top=v["z_top"],
            z_threshold=zt,
            z_diff=zd,
        )
        cats.append(category_of(st, v["mapped"]))
        assigned.append(v["mapped"])
    t = tile_counts(cats, assigned)
    assert t["assigned"] == _assigned_sum(t) == sum(assigned)


def _run_node(tmp_path: Path, payload: dict[str, Any], body: str) -> Any:
    (tmp_path / "rules.js").write_text(MATCH_RULES_JS)
    (tmp_path / "in.json").write_text(json.dumps(payload))
    (tmp_path / "run.js").write_text(
        "const M = require('./rules.js'); const P = require('./in.json');\n"
        + body
        + "\nprocess.stdout.write(JSON.stringify(out));\n"
    )
    proc = subprocess.run(
        ["node", "run.js"], cwd=tmp_path, capture_output=True, text=True, check=False
    )
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


@needs_node
@pytest.mark.parametrize(("zt", "zd"), [(1.96, 2.0), (0.5, 0.3), (3.0, 1.0), (6.0, 4.0)])
def test_js_tile_counts_match_python_on_vectors(tmp_path: Path, zt: float, zd: float) -> None:
    vecs = CLASSIFY_TEST_VECTORS
    out = _run_node(
        tmp_path,
        {"v": vecs, "zt": zt, "zd": zd},
        (
            "const st = P.v.map(v => M.classify(v.mapped, v.n_traits, v.rank_self, v.z_self, v.z_top, P.zt, P.zd));"
            "const cats = st.map((s, i) => M.category(s, P.v[i].mapped));"
            "const out = {cats: cats, tiles: M.tileCounts(cats, P.v.map(v => v.mapped))};"
        ),
    )
    py_cats = [
        category_of(
            classify(
                mapped=v["mapped"],
                n_traits=v["n_traits"],
                rank_self=v["rank_self"],
                z_self=v["z_self"],
                z_top=v["z_top"],
                z_threshold=zt,
                z_diff=zd,
            ),
            v["mapped"],
        )
        for v in vecs
    ]
    assert out["cats"] == py_cats
    t = out["tiles"]
    assert t == tile_counts(py_cats, [v["mapped"] for v in vecs])
    assert t["assigned"] == _assigned_sum(t)


@needs_node
def test_rendered_scripts_parse(two_dataset_run: Path, tmp_path: Path) -> None:
    datasets = [load_dataset(_ctx(two_dataset_run, ["rna"]), "rna", "rna")]
    html = render(build_page("p", 1.96, 2.0, datasets), "cdn")
    scripts = re.findall(r"<script>(.*?)</script>", html, re.S)
    assert scripts
    for i, js in enumerate(scripts):
        f = tmp_path / f"s{i}.js"
        f.write_text(js)
        proc = subprocess.run(
            ["node", "--check", str(f)], capture_output=True, text=True, check=False
        )
        assert proc.returncode == 0, proc.stderr


@needs_node
@pytest.mark.parametrize(("zt", "zd"), [(1.96, 2.0), (3.0, 1.0), (0.5, 0.3), (6.5, 4.0)])
def test_js_evaluate_matches_python_on_embedded_rows(
    two_dataset_run: Path, tmp_path: Path, zt: float, zd: float
) -> None:
    """The page's per-sample recomputation equals classify()/group_of() for every embedded row."""
    datasets = [load_dataset(_ctx(two_dataset_run, ["rna"]), "rna", "rna")]
    page = build_page("p", 1.96, 2.0, datasets)
    entry = page.data["datasets"][0]
    objs = [dict(zip(page.data["fields"], r, strict=True)) for r in entry["rows"]]
    out = _run_node(
        tmp_path,
        {"rows": objs, "zt": zt, "zd": zd},
        "const out = M.evaluateAll(P.rows, P.zt, P.zd);",
    )
    for r, e in zip(objs, out, strict=True):
        mapped = bool(r["geno_id_assigned"])
        if r["qpgen_status"] == "NOT_SCORED":
            assert e["status"] == "NOT_SCORED" and e["group"] == "NO_TEST"
            continue
        py = classify(
            mapped=mapped,
            n_traits=r["n_traits"],
            rank_self=r["rank_self"],
            z_self=r["z_self"],
            z_top=[r[f"z_{k}"] for k in range(1, 6)],
            z_threshold=zt,
            z_diff=zd,
        )
        assert e["status"] == py, r["pheno_id"]
        assert e["group"] == group_of(py, mapped)
        zs = [r[f"z_{k}"] for k in range(1, 6)]
        close = close_candidates([r[f"id_{k}"] or None for k in range(1, 6)], zs, zt)
        dup = is_duplicate(zs, zt, zd)
        assert e["dup"] == dup, r["pheno_id"]
        assert e["related"] == (not dup and is_related(len(close))), r["pheno_id"]
        _, two_ids, taken = display_flags({**r, "proposed_geno_id": e["proposed"]}, zt, zd)
        assert e["cat"] == category_of(py, mapped, two_ids, taken), r["pheno_id"]
    if (zt, zd) == (1.96, 2.0):
        # at the run's thresholds the page reproduces qpgen's statuses, the proposed IDs and
        # the categories of the saved results
        samples = datasets[0]["samples"]
        assert [x["status"] for x in out] == list(samples["qpgen_status"])
        assert [x["proposed"] for x in out] == list(samples["proposed_geno_id"].fillna(""))
        assert [x["cat"] for x in out] == run_categories(samples, zt, zd)[0]
        assert [r["pheno_id"] for r, x in zip(objs, out, strict=True) if x["dup"]] == ["S12"]


def test_dataset_switch_labels_never_show_default_names() -> None:
    def ds(name: str, omics: str) -> dict[str, Any]:
        return {"name": name, "omics": omics}

    assert _ds_labels([ds("a", "rna"), ds("b", "protein")]) == ["rna", "protein"]
    assert _ds_labels([ds("a", "rna"), ds("b", "rna")]) == ["a", "b"]
    assert _ds_labels([ds("a", "omics"), ds("b", "rna")]) == ["a", "b"]
    assert _ds_labels([ds("dataset", "omics"), ds("x", "omics")]) == ["Dataset 1", "Dataset 2"]


def test_categories_for_two_ids_and_used_ids() -> None:
    assert category_of("MULTI_NEW_BEST", True, two_ids=True) == "twoids"
    assert category_of("MULTI_NEW_BEST", False, two_ids=True) == "twoids"
    assert category_of("SINGLE_NEW_BEST", True, two_ids=True) == "twoids"
    assert category_of("SELF_BEST", True, two_ids=True) == "strict"  # own ID first: still strict
    assert category_of("MULTI_NEW_BEST", True) == "several"
    assert category_of("SINGLE_NEW_BEST", False, taken=True) == "taken"
    assert category_of("SINGLE_NEW_BEST", False) == "newid"
    assert category_of("NO_SELF", False) == "noid"
    zs = {"z_1": 7.0, "z_2": 6.83, "z_3": 2.56}
    r = {
        "qpgen_status": "MULTI_NEW_BEST",
        "geno_id_assigned": "g",
        "id_1": "g_twin",
        "id_2": "g",
        **zs,
    }
    assert display_flags(r, 1.96, 2.0) == (True, True, False)
    # matches two genotypes about equally, neither its own: a mismatch, shown as "fits several"
    assert display_flags({**r, "id_2": "h"}, 1.96, 2.0) == (True, False, False)
    u = {
        "qpgen_status": "SINGLE_NEW_BEST",
        "geno_id_assigned": "",
        "proposed_geno_id": "",
        "z_1": 5.0,
        "z_2": 2.9,
        "z_3": 2.8,
    }
    assert display_flags(u, 1.96, 2.0) == (False, False, True)
    assert display_flags({**u, "proposed_geno_id": "g"}, 1.96, 2.0) == (False, False, False)


def test_table_and_duplicates_markup(tmp_path: Path) -> None:
    summarized_run(tmp_path, "rna")
    report_match(_ctx(tmp_path, ["rna"]))
    html = (tmp_path / L.report("match")).read_text()
    js = render_mod._static("match_page.js")
    # page scroll with a sticky header and fixed column widths: no inner scroll box
    assert '<table class="mtable"><colgroup></colgroup>' in html and "mtable-wrap" not in html
    css = render_mod._static("report.css")
    assert "table-layout: fixed" in css and ".mtable-wrap" not in css
    # plain column names, no info icons in the table header
    for label in ("Own score", "Own rank", "Best other", "Its score", "Lead", "Suggested ID"):
        assert f'label: "{label}"' in js
    assert "help:" not in js and "Margin" not in js
    # short headers on one line; numeric columns with fixed pixel widths
    assert "white-space: nowrap; line-height: 1.25; height: 44px" in css
    assert "w: 124, pw: 64" in js
    # no info icons anywhere on the page (sliders have a plain hint line)
    assert 'class="tip"' not in html and ".tip {" not in css
    assert "How far the best genotype must beat the next one" in html and ">used 2.00<" in html
    # the duplicates section: one rule (ANQ_MATCH.isDuplicate), compact empty state, no relatives
    assert "No sample matches two genotype IDs." in js and "relatives" not in js


@needs_node
def test_js_two_genotype_ids_get_no_relabel_and_are_not_mismatches(tmp_path: Path) -> None:
    """Strong signal (12.4 vs 10, next 3): the page shows "Matches two genotype IDs", with
    no suggested ID and no mismatch in the verdict, for assigned and unassigned samples."""
    base = {"qpgen_status": "", "n_traits": 50, "cor_self": None}
    tops = {"z_1": 12.4, "z_2": 10.0, "z_3": 3.0, "z_4": 2.9, "z_5": 2.8}
    rows = [
        {
            **base,
            **tops,
            "pheno_id": "A",
            "geno_id_assigned": "G2",
            "rank_self": 2,
            "z_self": 10.0,
            "id_1": "G1",
            "id_2": "G2",
            "id_3": "X",
            "id_4": "Y",
            "id_5": "Z",
        },
        {
            **base,
            **tops,
            "pheno_id": "_U",
            "geno_id_assigned": "",
            "rank_self": None,
            "z_self": None,
            "id_1": "H1",
            "id_2": "H2",
            "id_3": "X",
            "id_4": "Y",
            "id_5": "Z",
        },
    ]
    out = _run_node(
        tmp_path,
        {"rows": rows},
        "const e = M.evaluateAll(P.rows, 1.96, 2.0); const out = {e: e, v: M.verdict(P.rows, e)};",
    )
    assert [x["status"] for x in out["e"]] == ["SINGLE_NEW_BEST", "SINGLE_NEW_BEST"]
    assert [x["cat"] for x in out["e"]] == ["twoids", "twoids"]
    assert [x["proposed"] for x in out["e"]] == ["", ""]
    assert "mismatch" not in out["v"]["text"] and out["v"]["status"] == "warning"
    assert "0 rescue candidates" in out["v"]["text"]
    t = tile_counts([x["cat"] for x in out["e"]], [True, False], [True, True])
    assert (t["twoids"], t["un_twoids"], t["dups"]) == (1, 1, 2)
