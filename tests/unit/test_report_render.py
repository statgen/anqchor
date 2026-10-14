"""Report page rendering: structure, embedded JSON, escaping, plotly.js modes."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from anqchor.report.render import (
    PLOTLY_CDN,
    Page,
    Section,
    _json,
    callout,
    code_block,
    definition_list,
    details,
    figure,
    paragraph,
    render,
    table,
    tiles,
    write_page,
)
from anqchor.report.theme import CATEGORICAL_LIGHT, base_layout, css_variables, series_color

SCRIPT_RE = re.compile(
    r'<script type="application/json" class="(fig|table)-data" data-target="([^"]+)">(.*?)</script>',
    re.S,
)


def embedded(html: str) -> dict[str, dict]:
    """Every embedded figure/table JSON payload, parsed like the browser does."""
    return {m.group(2): json.loads(m.group(3)) for m in SCRIPT_RE.finditer(html)}


def sample_page() -> Page:
    df = pd.DataFrame(
        {
            "sample": ["A", "B", "C"],
            "z": [1.5, float("nan"), np.inf],
            "n": np.array([1, 2, 3], dtype=np.int64),
            "group": pd.array(["PASS", None, "FAIL"], dtype="string"),
        }
    )
    fig = go.Figure(
        go.Scatter(x=np.arange(3), y=np.array([1.0, np.nan, 2.0])), layout=base_layout()
    )
    body = (
        tiles([("Samples", "729", "mapped", "good"), ("Fail", "1", "")])
        + callout("One sample matches another genotype.", "critical", title="Mismatch")
        + figure(fig, "fig-z", caption="Z scores <b>raw</b>")
        + table(df, "tbl-samples", caption="Per-sample status", status_columns=["group"])
    )
    return Page(
        title="Identity <matching>",
        subtitle="MAGE & friends",
        sections=[
            Section("overview", "Overview", body, intro="Lead <text>"),
            Section(
                "params", "Parameters", definition_list([("lambda", 0.5)]) + code_block("a < b")
            ),
        ],
        nav=[("Index", "index.html"), ("cis", "cis.html")],
        data={"thresholds": {"z": 1.96}, "evil": "</script><script>alert(1)</script>"},
        extra_js="ANQ.onReady(function () { window.__ok = true; });",
    )


def test_page_structure() -> None:
    html = render(sample_page(), plotlyjs="cdn")
    assert html.startswith("<!doctype html>")
    assert '<section id="overview">' in html and '<section id="params">' in html
    assert "<h2>Overview</h2>" in html and "<h2>Parameters</h2>" in html
    assert '<a href="#overview">Overview</a>' in html  # table of contents
    assert '<a href="index.html">Index</a>' in html
    # page text is escaped
    assert "<title>Identity &lt;matching&gt;</title>" in html
    assert "MAGE &amp; friends" in html and "Lead &lt;text&gt;" in html
    assert "Z scores &lt;b&gt;raw&lt;/b&gt;" in html
    assert "a &lt; b" in html
    assert 'class="tiles"' in html and "status-good" in html and "callout-critical" in html


def test_css_and_js_not_html_escaped() -> None:
    """<style>/<script> are raw text: entities there would break the CSS and JS."""
    html = render(sample_page(), plotlyjs="cdn")
    assert '"use strict";' in html
    assert ':root:not([data-theme="light"])' in html
    assert "&#34;" not in html and "&quot;use strict" not in html
    assert "ANQ.onReady(function () { window.__ok = true; });" in html


def test_embedded_json_valid_and_clean() -> None:
    html = render(sample_page(), plotlyjs="cdn")
    data = embedded(html)
    assert set(data) == {"fig-z", "tbl-samples"}
    tbl = data["tbl-samples"]
    assert tbl["columns"] == ["sample", "z", "n", "group"]
    assert tbl["rows"] == [["A", 1.5, 1, "PASS"], ["B", None, 2, None], ["C", None, 3, "FAIL"]]
    assert tbl["status"] == ["group"] and tbl["filename"] == "tbl-samples.tsv"
    fig = data["fig-z"]
    assert fig["layout"]["height"] == 420
    assert fig["data"][0]["type"] == "scatter"
    # page data: '</script>' cannot terminate the element
    m = re.search(r"window\.ANQ_DATA = (.*?);\n", html)
    assert m is not None
    assert "</script>" not in m.group(1)
    assert json.loads(m.group(1))["evil"] == "</script><script>alert(1)</script>"
    assert html.count("<script>alert(1)</script>") == 0


def test_json_escaping_and_nan() -> None:
    text = _json(
        {
            "s": "</script><!-- & >",
            "nan": float("nan"),
            "inf": -np.inf,
            "np": np.float32(np.nan),
            "na": pd.NA,
            "nat": pd.NaT,
            "arr": np.array([1.0, np.nan]),
            "ser": pd.Series([1, 2]),
            "ts": pd.Timestamp("2026-10-11"),
            "int": np.int64(3),
            "bool": np.bool_(True),
            "path": Path("/a/b"),
            "tuple": (1, 2),
            "sep": " ",
        }
    )
    assert "<" not in text and ">" not in text and "&" not in text and " " not in text
    back = json.loads(text)
    assert back["s"] == "</script><!-- & >"
    assert back["nan"] is None and back["inf"] is None and back["np"] is None
    assert back["na"] is None and back["nat"] is None
    assert back["arr"] == [1.0, None] and back["ser"] == [1, 2]
    assert back["ts"].startswith("2026-10-11")
    assert back["int"] == 3 and back["bool"] is True and back["path"] == "/a/b"
    assert back["tuple"] == [1, 2] and back["sep"] == " "


def test_table_truncation_caption() -> None:
    df = pd.DataFrame({"x": range(30)})
    out = table(df, "t", caption="Values.", max_rows=10)
    assert "Showing the first 10 of 30 rows" in out
    spec = embedded(out)["t"]
    assert len(spec["rows"]) == 10
    full = table(df, "t2", caption="Values.")
    assert "Showing the first" not in full
    assert len(embedded(full)["t2"]["rows"]) == 30
    only_note = table(df, "t3", max_rows=5)
    assert "Showing the first 5 of 30 rows" in only_note


def test_figure_from_dict_does_not_mutate() -> None:
    layout = {"title": {"text": "t"}}
    spec = {
        "data": [{"type": "bar", "x": np.array(["a", "b"]), "y": np.array([1, 2])}],
        "layout": layout,
    }
    out = figure(spec, "f", height=300)
    assert layout == {"title": {"text": "t"}}
    parsed = embedded(out)["f"]
    assert parsed["data"][0]["x"] == ["a", "b"] and parsed["layout"]["height"] == 300
    assert 'style="min-height:300px"' in out


def test_plotlyjs_modes(tmp_path: Path) -> None:
    page = Page(title="t", sections=[Section("a", "A", paragraph("x"))])
    cdn = render(page, plotlyjs="cdn")
    assert f'<script src="{PLOTLY_CDN}"' in cdn
    embed = render(page, plotlyjs="embed")
    assert PLOTLY_CDN not in embed
    assert len(embed) > len(cdn) + 1_000_000  # plotly.js inlined
    with pytest.raises(ValueError, match="embed' or 'cdn"):
        render(page, plotlyjs="local")
    out = write_page(tmp_path / "report" / "x.html", page, plotlyjs="cdn")
    assert out.read_text(encoding="utf-8").startswith("<!doctype html>")
    assert not [p for p in out.parent.iterdir() if p.name.endswith(".tmp")]


def test_single_section_has_no_toc() -> None:
    html = render(
        Page(title="t", sections=[Section("only", "Only", details("more", "<p>x</p>"))]),
        plotlyjs="cdn",
    )
    assert 'class="toc"' not in html
    assert "<details><summary>more</summary><p>x</p></details>" in html


def test_theme_helpers() -> None:
    css = css_variables()
    assert (
        "--surface:" in css
        and "--series-1: #2a78d6;" in css
        and "prefers-color-scheme: dark" in css
    )
    assert series_color(0) == CATEGORICAL_LIGHT[0]
    with pytest.raises(ValueError):
        series_color(8)
    lay = base_layout(xaxis={"title": {"text": "x"}}, height=200)
    assert lay["xaxis"]["title"] == {"text": "x"} and lay["xaxis"]["showgrid"] is True
    assert lay["height"] == 200


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_inline_scripts_parse_as_javascript(tmp_path: Path) -> None:
    """Every inline <script> of a rendered page is syntactically valid JavaScript."""
    html = render(sample_page(), plotlyjs="cdn")
    scripts = re.findall(r"<script>(.*?)</script>", html, re.S)
    assert len(scripts) == 3  # theme/data, runtime, page script
    for i, js in enumerate(scripts):
        f = tmp_path / f"s{i}.js"
        f.write_text(js, encoding="utf-8")
        proc = subprocess.run(
            ["node", "--check", str(f)], capture_output=True, text=True, check=False
        )
        assert proc.returncode == 0, proc.stderr


def test_tile_status_has_a_word_for_screen_readers() -> None:
    html = tiles([("Unresolved", "7%", "note", "warning"), ("Plain", "1", "n")])
    assert 'class="sr-only">status warning: </span>' in html
    assert html.count("status-dot") == 1


def test_offline_cdn_page_explains_missing_charts() -> None:
    """With plotlyjs: cdn and no network, each chart area says why it is empty."""
    from anqchor.report import render as R

    js = (Path(R.__file__).parent / "static" / "report.js").read_text()
    assert "--plotlyjs embed" in js and "plotly-missing" in js
