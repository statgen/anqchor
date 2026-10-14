"""Assemble self-contained HTML report pages.

A page is a list of :class:`Section` objects whose bodies are built with the helpers
here (:func:`figure`, :func:`table`, :func:`tiles`, :func:`callout`). Data for figures
and tables is embedded as JSON; ``static/report.js`` renders it, provides sortable and
searchable tables with TSV download, and re-themes Plotly figures for dark mode.
Pages need no network access when plotly.js is embedded (the default).
"""

from __future__ import annotations

import datetime as _dt
import html
import json
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Any

import jinja2
import numpy as np
import pandas as pd
from markupsafe import Markup

from anqchor import __version__
from anqchor.report.theme import CATEGORICAL_DARK, CATEGORICAL_LIGHT, DARK, LIGHT, css_variables
from anqchor.util.files import atomic_write

PLOTLY_CDN = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"


@dataclass
class Section:
    """One titled block of a page. ``body`` is trusted HTML built with this module."""

    id: str
    title: str
    body: str
    intro: str = ""
    """Optional plain-text lead paragraph (escaped)."""


@dataclass
class Page:
    """A report page: title, sections, navigation links and page-specific script/data."""

    title: str
    subtitle: str = ""
    sections: list[Section] = field(default_factory=list)
    nav: list[tuple[str, str]] = field(default_factory=list)
    """(label, href) links to sibling reports."""
    extra_js: str = ""
    """Page-specific JavaScript (runs after the shared script)."""
    data: dict[str, Any] = field(default_factory=dict)
    """Extra JSON exposed to ``extra_js`` as ``window.ANQ_DATA``."""
    header_html: str = ""
    """Trusted HTML placed under the title (context line, dataset switch)."""
    hide_brand: bool = False
    """Omit the small "anqchor" line above the title."""
    body_class: str = ""
    """CSS class of ``<body>`` (page-specific styling)."""
    toc: bool = True
    """Show the section links under the header (pages with more than one section)."""


def index_nav(params: dict[str, Any]) -> list[tuple[str, str]]:
    """The navigation link to ``report/index.html``, only when the run writes that page.

    :func:`anqchor.pipelines.build_graph` sets the ``has_index`` step parameter on every
    report step of a run with an index page; runs of a single analysis have none, and a
    link to it would be broken."""
    return [("Index", "index.html")] if params.get("has_index") else []


def esc(text: object) -> str:
    """HTML-escape ``text`` (``None`` → empty string)."""
    return html.escape("" if text is None else str(text), quote=True)


_SCRIPT_ESCAPES = {
    "<": "\\u003c",
    ">": "\\u003e",
    "&": "\\u0026",
    "\u2028": "\\u2028",
    "\u2029": "\\u2029",
}


def _clean(o: Any) -> Any:
    """Convert ``o`` into plain JSON types: non-finite floats, ``pd.NA`` and ``NaT``
    become ``None``; numpy arrays, pandas Series/Index, sets and tuples become lists."""
    if o is None or isinstance(o, (bool, str, int)):
        return o
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set, frozenset)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.ndarray, pd.Series, pd.Index)):
        return [_clean(v) for v in o.tolist()]
    if o is pd.NA or o is pd.NaT:
        return None
    if isinstance(o, (pd.Timestamp, _dt.datetime, _dt.date)):
        return o.isoformat()
    if isinstance(o, np.generic):  # numpy scalar
        return _clean(o.item())
    if isinstance(o, Path):
        return str(o)
    return o


def _json(obj: Any) -> str:
    """JSON safe to embed in a ``<script>`` element.

    NaN/inf become ``null``; ``<``, ``>``, ``&`` and the JavaScript line separators are
    written as ``\\uXXXX`` escapes, so neither ``</script>`` nor ``<!--`` can end or
    alter the element and the text stays valid JSON for ``JSON.parse``."""
    text = json.dumps(_clean(obj), allow_nan=False, separators=(",", ":"), ensure_ascii=False)
    for ch, rep in _SCRIPT_ESCAPES.items():
        text = text.replace(ch, rep)
    return text


def figure(fig: Any, fig_id: str, caption: str = "", height: int = 420) -> str:
    """Embed a Plotly figure. ``fig`` is a ``plotly.graph_objects.Figure`` or a dict
    with ``data`` and ``layout``; the caption is plain text."""
    spec = fig.to_plotly_json() if hasattr(fig, "to_plotly_json") else dict(fig)
    layout = dict(spec.get("layout") or {})  # copy: never mutate the caller's figure
    layout.setdefault("height", height)
    payload = _json({"data": spec.get("data", []), "layout": layout})
    cap = f"<figcaption>{esc(caption)}</figcaption>" if caption else ""
    return (
        f'<figure class="viz"><div class="plotly-fig" id="{esc(fig_id)}" style="min-height:{int(height)}px"></div>'
        f'<script type="application/json" class="fig-data" data-target="{esc(fig_id)}">{payload}</script>{cap}</figure>'
    )


def table(
    df: pd.DataFrame,
    table_id: str,
    caption: str = "",
    *,
    max_rows: int = 20000,
    page_size: int = 500,
    digits: int = 4,
    searchable: bool = True,
    column_help: dict[str, str] | None = None,
    status_columns: Sequence[str] = (),
) -> str:
    """A sortable, searchable table with a TSV download of the embedded rows.

    Rows beyond ``max_rows`` are not embedded; the caption says so (no silent caps).
    ``status_columns`` are rendered as labelled status chips (PASS/WARN/FAIL/NO_TEST).
    """
    n = len(df)
    shown = df.head(max_rows)
    cols = [str(c) for c in shown.columns]
    rows = [list(rec) for rec in shown.itertuples(index=False, name=None)]  # cleaned by _json
    spec = {
        "columns": cols,
        "rows": rows,
        "pageSize": page_size,
        "digits": digits,
        "searchable": searchable,
        "help": column_help or {},
        "status": list(status_columns),
        "filename": f"{table_id}.tsv",
    }
    note = ""
    if n > max_rows:
        note = f" Showing the first {max_rows:,} of {n:,} rows; the full table is in the run directory."
    cap = f'<p class="table-caption">{esc(caption)}{esc(note)}</p>' if (caption or note) else ""
    return (
        f'{cap}<div class="anq-table" id="{esc(table_id)}"></div>'
        f'<script type="application/json" class="table-data" data-target="{esc(table_id)}">{_json(spec)}</script>'
    )


def tiles(items: Iterable[tuple[str, str, str]] | Iterable[tuple[str, str, str, str]]) -> str:
    """Stat tiles: (label, value, note[, status]) where status is good/warning/serious/critical."""
    out = ['<div class="tiles">']
    for it in items:
        label, value, note, *rest = it
        status = rest[0] if rest else ""
        dot = (
            (
                f'<span class="status-dot status-{esc(status)}" aria-hidden="true"></span>'
                f'<span class="sr-only">status {esc(status)}: </span>'
            )
            if status
            else ""
        )
        out.append(
            f'<div class="tile"><div class="tile-label">{dot}{esc(label)}</div>'
            f'<div class="tile-value">{esc(value)}</div><div class="tile-note">{esc(note)}</div></div>'
        )
    out.append("</div>")
    return "".join(out)


STATUS_ICON = {"good": "✓", "warning": "!", "serious": "!", "critical": "✕", "info": "i"}


def callout(text: str, status: str = "info", title: str = "") -> str:
    """A verdict or note box. ``status``: good, warning, serious, critical or info."""
    icon = STATUS_ICON.get(status, "i")
    t = f"<strong>{esc(title)}</strong> " if title else ""
    return f'<div class="callout callout-{esc(status)}" role="note"><span class="callout-icon" aria-hidden="true">{icon}</span><div>{t}{esc(text)}</div></div>'


def paragraph(text: str) -> str:
    """A ``<p>`` with escaped text."""
    return f"<p>{esc(text)}</p>"


def definition_list(items: Iterable[tuple[str, object]]) -> str:
    """A key/value list (parameters, provenance); keys and values are escaped."""
    rows = "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in items)
    return f'<dl class="params">{rows}</dl>'


def code_block(text: str) -> str:
    """Preformatted, escaped text."""
    return f'<pre class="code"><code>{esc(text)}</code></pre>'


def details(summary: str, body_html: str) -> str:
    """A collapsible block; ``body_html`` is trusted HTML."""
    return f"<details><summary>{esc(summary)}</summary>{body_html}</details>"


def _static(name: str) -> str:
    # one joinpath argument at a time: Python 3.10's Traversable.joinpath takes only one
    return (
        resources.files("anqchor.report")
        .joinpath("static")
        .joinpath(name)
        .read_text(encoding="utf-8")
    )


def _plotly_tag(mode: str) -> str:
    if mode not in ("embed", "cdn"):
        raise ValueError(f"plotlyjs must be 'embed' or 'cdn', got {mode!r}")
    if mode == "cdn":
        return f'<script src="{PLOTLY_CDN}" charset="utf-8"></script>'
    from plotly.offline import get_plotlyjs

    return "<script>" + str(get_plotlyjs()) + "</script>"


def render(page: Page, plotlyjs: str = "embed") -> str:
    """Render a :class:`Page` to a complete HTML document."""
    env = jinja2.Environment(
        loader=jinja2.PackageLoader("anqchor.report", "templates"),
        autoescape=jinja2.select_autoescape(["html"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    tpl = env.get_template("base.html")
    theme = {
        "light": LIGHT,
        "dark": DARK,
        "seriesMap": dict(zip(CATEGORICAL_LIGHT, CATEGORICAL_DARK, strict=True)),
    }
    # The CSS and JS are trusted package files placed inside <style>/<script>, where
    # HTML entities are not decoded: they must not be autoescaped.
    return tpl.render(
        page=page,
        css=Markup(css_variables() + _static("report.css")),
        js=Markup(_static("report.js")),
        theme_json=Markup(_json(theme)),
        data_json=Markup(_json(page.data)),
        extra_js=Markup(page.extra_js),
        header_html=Markup(page.header_html),
        plotly_tag=Markup(_plotly_tag(plotlyjs)),
        sections=[
            {"id": s.id, "title": s.title, "intro": s.intro, "body": Markup(s.body)}
            for s in page.sections
        ],
        version=__version__,
        generated=_dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z"),
    )


def write_page(path: str | Path, page: Page, plotlyjs: str = "embed") -> Path:
    """Render ``page`` and write it atomically to ``path``."""
    path = Path(path)
    with atomic_write(path) as fh:
        fh.write(render(page, plotlyjs))
    return path
