"""The report index (``report/index.html``, step kind ``report_index``).

The page links every report present in ``report/`` (the directory is scanned, so
reports from earlier invocations with other pipelines are listed too) and shows one
tile per pipeline and dataset taken from that dataset's ``summary.json``:

* match: the verdict and the group counts;
* cis: the recommended strategy, its replication at the deciding level and the reason;
* trans: the run's cross-chromosome λ_GC, hotspots and replicated reference pairs
  (older summaries: the recommended run with its trans λ_GC and Bonferroni cells).

The cis and trans summaries are read defensively (``.get`` with several candidate
keys): whatever exists is shown, and a missing or unreadable summary is noted
rather than failing the step.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from anqchor.pipelines.layout import Layout
from anqchor.report.render import PLOTLY_CDN, Page, Section, callout, esc, render
from anqchor.steps import StepContext, register
from anqchor.util.files import atomic_write

PIPELINE_ORDER = ("match", "cis", "trans")
PIPELINE_TITLES = {
    "match": "Sample identity (match)",
    "cis": "cis-xQTL power by strategy (cis)",
    "trans": "trans-QTL evaluation (trans)",
}
SUMMARY_DIR = {"match": Layout.match_dir, "cis": Layout.cis_dir, "trans": Layout.trans_dir}
STATUS_WORDS = {
    "good": "good",
    "pass": "good",
    "ok": "good",
    "warning": "warning",
    "warn": "warning",
    "serious": "serious",
    "critical": "critical",
    "fail": "critical",
    "error": "critical",
}

INDEX_CSS = """
<style>
.idx-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 12px; }
.idx-card { border: 1px solid var(--grid); border-radius: 8px; padding: 12px 14px; background: var(--surface); }
.idx-card h3 { margin: 0 0 6px; font-size: 15px; display: flex; align-items: center; gap: 6px; }
.idx-card p { margin: 4px 0; color: var(--ink2); font-size: 14px; }
.idx-card dl { display: grid; grid-template-columns: max-content 1fr; gap: 2px 12px; margin: 6px 0 0; font-size: 13px; }
.idx-card dt { color: var(--muted); }
.idx-card dd { margin: 0; font-variant-numeric: tabular-nums; }
.idx-links { display: flex; flex-wrap: wrap; gap: 8px 16px; margin: 0 0 4px; font-size: 15px; }
</style>
"""


def first(d: dict[str, Any], *keys: str) -> Any:
    """The first present, non-null value among ``keys`` (dotted keys reach into dicts)."""
    for k in keys:
        cur: Any = d
        for part in k.split("."):
            cur = cur.get(part) if isinstance(cur, dict) else None
            if cur is None:
                break
        if cur is not None:
            return cur
    return None


def fmt_value(v: Any) -> str:
    """Short display text for a summary value."""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if not math.isfinite(v):
            return "NA"
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e6):
            return f"{v:.2e}"
        return f"{v:.3g}" if abs(v) < 100 else f"{v:,.0f}"
    if isinstance(v, int):
        return f"{v:,}"
    if isinstance(v, dict):
        return ", ".join(f"{k} {fmt_value(x)}" for k, x in list(v.items())[:6])
    if isinstance(v, list):
        return ", ".join(fmt_value(x) for x in v[:6]) + (" …" if len(v) > 6 else "")
    return str(v)


def _status(summary: dict[str, Any]) -> str:
    s = first(summary, "verdict_status", "status")
    return STATUS_WORDS.get(str(s).lower(), "") if s is not None else ""


def _label(rec: Any) -> str | None:
    """``strategy/transform`` of a recommendation dict (or the value itself)."""
    if isinstance(rec, dict):
        lab = first(rec, "label", "name")
        if lab is None and rec.get("strategy") is not None:
            lab = rec["strategy"] + (f"/{rec['transform']}" if rec.get("transform") else "")
        return None if lab is None else str(lab)
    return None if rec is None else str(rec)


def _run_of(summary: dict[str, Any], rec: Any) -> dict[str, Any]:
    """The per-run record of the recommended run (cis ``per_run``, formerly ``per_strategy``,
    keyed by label; trans ``strategies`` list), or an empty dict."""
    lab = _label(rec)
    per = summary.get("per_run", summary.get("per_strategy"))
    if isinstance(per, dict) and isinstance(per.get(lab), dict):
        found: dict[str, Any] = per[lab]
        return found
    runs = summary.get("strategies")
    if isinstance(runs, list) and isinstance(rec, dict):
        for r in runs:
            if (
                isinstance(r, dict)
                and r.get("strategy") == rec.get("strategy")
                and r.get("transform", rec.get("transform")) == rec.get("transform")
            ):
                return r
    return {}


def _of(n: Any, total: Any) -> str | None:
    """``"n of total (p%)"`` when both are numbers, else what exists."""
    if isinstance(n, (int, float)) and isinstance(total, (int, float)) and total:
        return f"{fmt_value(n)} of {fmt_value(total)} ({100.0 * n / total:.0f}%)"
    return None if n is None else fmt_value(n)


def facts(pipeline: str, summary: dict[str, Any]) -> list[tuple[str, str]]:
    """The key facts of one dataset summary for its tile (only those present).

    Keys are read defensively: the cis and trans summaries evolve independently, so
    every value is optional and alternative key names are tried."""
    out: list[tuple[str, Any]] = []
    if pipeline == "match":
        gc = summary.get("group_counts") or {}
        if gc:
            out.append(
                (
                    "groups",
                    ", ".join(
                        f"{g} {gc[g]}" for g in ("PASS", "WARN", "FAIL", "NO_TEST") if g in gc
                    ),
                )
            )
        out += [
            ("mapped samples", _of(summary.get("n_mapped"), summary.get("n_pheno"))),
            ("rescue candidates", summary.get("n_rescue_candidates")),
            ("possible relatives", summary.get("n_related")),
            (
                "λ",
                "none (independent traits)"
                if summary.get("lambda_mode") == "independent"
                else summary.get("lambda_used"),
            ),
        ]
    elif pipeline == "cis":
        rec = first(summary, "recommended", "recommended_strategy", "best_strategy")
        run = _run_of(summary, rec)
        level = first(rec, "key_level", "level") if isinstance(rec, dict) else None
        out.append(("recommended", _label(rec)))
        if level is not None:
            rep = _of(run.get(f"n_rep_{level}"), run.get("n_tested"))
            out.append((f"replicated ({level})", rep))
        out.append(
            (
                "reference pairs mapped",
                _of(summary.get("n_pairs_mapped"), summary.get("n_ref_pairs")),
            )
        )
        if isinstance(summary.get("runs"), list):
            out.append(("runs", len(summary["runs"])))
        else:
            strategies = summary.get("strategies")
            out.append(
                ("strategies", len(strategies) if isinstance(strategies, list) else strategies)
            )
    elif pipeline == "trans" and isinstance(summary.get("inflation"), dict):
        # schema 2: one run; replication by class, hotspots, cross-chromosome calibration
        out.append(("run", summary.get("run")))
        lam = summary["inflation"].get("lambda_gc")
        out.append(
            ("λ_GC cross-chromosome", f"{lam:.2f}" if isinstance(lam, (int, float)) else lam)
        )
        hot: dict[str, Any] = summary.get("hotspots") or {}
        dflt: dict[str, Any] = hot.get("default") or {}
        if dflt.get("n_hotspots") is not None:
            out.append((f"hotspots (≥ {dflt.get('k')} genes)", dflt.get("n_hotspots")))
        replication: dict[str, Any] = summary.get("replication") or {}
        for cls in ("cis", "trans"):
            rc: dict[str, Any] = replication.get(cls) or {}
            if rc.get("n_pairs"):
                out.append(
                    (
                        f"{cls} pairs replicated (p < 1e-5)",
                        _of(rc.get("n_rep_strict"), rc.get("n_tested")),
                    )
                )
        nv, nt = summary.get("n_variants_tested"), summary.get("n_traits_tested")
        if nv is not None and nt is not None:
            out.append(("variants × genes", f"{fmt_value(nv)} × {fmt_value(nt)}"))
    elif pipeline == "trans":
        rec = first(summary, "recommended")
        run = _run_of(summary, rec)
        out.append(("recommended", _label(rec)))
        out.append(
            (
                "λ_GC trans",
                first(run, "lambda_gc_trans", "lambda_gc")
                if run
                else first(summary, "lambda_gc_trans", "lambda_gc", "recommended.lambda_gc_trans"),
            )
        )
        if not run.get("lambda_gc_trans") and run.get("lambda_gc_same_chrom") is not None:
            out.append(("λ_GC same-chrom", run.get("lambda_gc_same_chrom")))
        out.append(
            (
                "hotspots",
                first(run, "n_hotspots") if run else first(summary, "n_hotspots", "hotspots"),
            )
        )
        out.append(("Bonferroni trans cells", first(run, "n_sig_trans_bonferroni")))
        nv, nt = summary.get("n_variants"), summary.get("n_traits")
        if nv is not None and nt is not None:
            out.append(("variants × traits", f"{fmt_value(nv)} × {fmt_value(nt)}"))
        if isinstance(rec, dict) and rec.get("known_replicated") is not None:
            # Denominator: the known pairs actually tested in the recommended run (pairs
            # with an absent variant cannot replicate), as on the trans report itself.
            known = run.get("known") if run else None
            n_tested = known.get("n_tested") if isinstance(known, dict) else None
            out.append(
                (
                    "known pairs replicated",
                    _of(rec.get("known_replicated"), n_tested or summary.get("n_known_pairs")),
                )
            )
    warnings = summary.get("warnings")
    if isinstance(warnings, list) and warnings:
        out.append(("warnings", len(warnings)))
    return [(k, fmt_value(v)) for k, v in out if v is not None and v != ""]


def read_summary(path: Path) -> tuple[dict[str, Any] | None, str]:
    """``(summary, problem)``: the parsed JSON, or ``None`` and why."""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None, "summary not written yet"
    except (OSError, ValueError) as e:
        return None, f"unreadable summary ({e.__class__.__name__})"
    if not isinstance(data, dict):
        return None, "summary is not a JSON object"
    return data, ""


def _card(
    pipeline: str, ds: str, summary: dict[str, Any] | None, problem: str, href: str | None
) -> str:
    status = _status(summary) if summary else ""
    dot = (
        f'<span class="status-dot status-{esc(status)}" aria-hidden="true"></span>'
        if status
        else ""
    )
    title = f"{esc(ds)}"
    if href:
        anchor = f"#ds-{re.sub(r'[^A-Za-z0-9_-]', '-', ds)}" if pipeline == "match" else ""
        title = f'<a href="{esc(href)}{esc(anchor)}">{title}</a>'
    parts = [f'<div class="idx-card"><h3>{dot}{title}</h3>']
    if summary is None:
        parts.append(f"<p>{esc(problem)}</p></div>")
        return "".join(parts)
    verdict = first(summary, "verdict", "headline", "message", "recommended.reason")
    if verdict:
        label = f"{status.upper()}: " if status in ("critical",) else ""
        parts.append(f"<p>{esc(label)}{esc(verdict)}</p>")
    f = facts(pipeline, summary)
    if f:
        parts.append(
            "<dl>" + "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in f) + "</dl>"
        )
    parts.append("</div>")
    return "".join(parts)


def scan_reports(report_dir: Path) -> list[str]:
    """Pipeline report names present in ``report_dir`` (``index`` excluded), known
    pipelines first."""
    if not report_dir.is_dir():
        return []
    names = sorted(p.stem for p in report_dir.glob("*.html") if p.stem != "index")
    return [n for n in PIPELINE_ORDER if n in names] + [n for n in names if n not in PIPELINE_ORDER]


def summary_datasets(run_dir: Path, pipeline: str, declared: Sequence[str]) -> list[str]:
    """Datasets with a summary for ``pipeline``: the declared ones first, then any other
    ``<pipeline>/<ds>/summary.json`` found in the run directory."""
    out = list(declared)
    base = run_dir / pipeline
    if base.is_dir():
        for p in sorted(base.glob("*/summary.json")):
            if p.parent.name not in out:
                out.append(p.parent.name)
    return out


def build_index(
    run_dir: Path, project: str, pipelines: Sequence[str], datasets: dict[str, list[str]]
) -> Page:
    """The index page for the run in ``run_dir``."""
    reports = scan_reports(run_dir / "report")
    shown = [p for p in PIPELINE_ORDER if p in pipelines or p in reports]
    shown += [p for p in reports if p not in shown]
    sections: list[Section] = []
    links = "".join(
        f'<a href="{esc(n)}.html">{esc(PIPELINE_TITLES.get(n, n))}</a>' for n in reports
    )
    if links:
        body = INDEX_CSS + f'<div class="idx-links">{links}</div>'
    else:
        body = INDEX_CSS + callout("No pipeline report has been written yet.", "info")
    missing = [p for p in pipelines if p not in reports]
    if missing:
        body += callout(f"Requested but not rendered yet: {', '.join(missing)}.", "warning")
    sections.append(Section(id="reports", title="Reports", body=body))
    for p in shown:
        if p not in SUMMARY_DIR:
            continue
        href = f"{p}.html" if p in reports else None
        cards = []
        for ds in summary_datasets(run_dir, p, datasets.get(p, [])):
            summary, problem = read_summary(run_dir / SUMMARY_DIR[p](ds) / "summary.json")
            cards.append(_card(p, ds, summary, problem, href))
        if not cards:
            cards.append('<div class="idx-card"><p>No dataset summaries found.</p></div>')
        sections.append(
            Section(
                id=f"pipe-{p}",
                title=PIPELINE_TITLES.get(p, p),
                body='<div class="idx-grid">' + "".join(cards) + "</div>",
            )
        )
    return Page(
        title=f"{project}: anqchor reports",
        subtitle="Genetically anchored quality control of the omics datasets in this run.",
        sections=sections,
    )


@register("report_index")
def report_index(ctx: StepContext) -> None:
    """Render ``report/index.html``."""
    p = ctx.params
    page = build_index(
        ctx.run_dir,
        p.get("project", "anqchor"),
        list(p.get("pipelines") or []),
        dict(p.get("datasets") or {}),
    )
    write_index(ctx.path(Layout.report("index")), page)


def write_index(path: Path, page: Page) -> Path:
    """Write the index page without plotly.js (it has no figures): rendered in CDN mode
    and the CDN script tag removed, so the page stays small and makes no network request."""
    html = render(page, "cdn").replace(f'<script src="{PLOTLY_CDN}" charset="utf-8"></script>', "")
    with atomic_write(path) as fh:
        fh.write(html)
    return path
