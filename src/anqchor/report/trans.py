"""``report.trans`` (kind ``report_trans``): the trans report (``report/trans.html``).

Built only from the summaries of ``summarize_trans`` (``trans/<dataset>/…``), never from
qpgen output, so ``anqchor trans report`` re-renders it in about a second. One page shows
one run. From top to bottom:

1. header: title, a context line (cells tested, cross-chromosome cells, samples, reference
   file, covariates, RINT) and the summary warnings;
2. replicated reference pairs: cis and trans pairs side by side (tested, same sign,
   p < 0.05, p < 1e-5), with a "no trans pairs in the reference" state;
3. hotspots: a p-value threshold and a minimum gene count chosen with buttons, the number
   of hotspots, the genes with a trans signal per bin along the genome, the grid of hotspot
   counts over thresholds × minimum genes (click a cell to select it) and the hotspot bins;
4. inflation: λ_GC and observed vs expected counts of the cross-chromosome cells, the QQ
   plot with its 95% band, and what was left out;
5. run details (collapsed): settings, cell counts per class, the reference pairs table.

A comparison directory (``anqchor trans compare``, params ``comparison: true``) gets one page
for every run, built from the comparison tables only (DESIGN 11.2-R): the replicated pairs
per run (cis and trans blocks, runs as rows), the hotspot count of every run under one
threshold and minimum, the bins, grid and table of a selected run, λ_GC per run, one QQ
line per run, observed vs expected per run, and the compared runs (collapsed).

The page script (``static/trans_page.js``) draws the charts from ``ANQ_DATA.panels`` or
``ANQ_DATA.comparison``.
"""

from __future__ import annotations

import math
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from anqchor.pipelines.layout import Layout
from anqchor.report.cis import (
    _cov_text,
    _details,
    _empty_fig,
    _r4,
    _seg,
    _series_var,
    chart_layout,
    fmt_alpha,
    fmt_int,
    fmt_pct,
    run_colors,
    run_dash,
)
from anqchor.report.render import Page, Section, definition_list, esc, index_nav, table, write_page
from anqchor.report.theme import LIGHT, STATUS
from anqchor.steps import StepContext, StepError, register
from anqchor.util.files import read_json

TITLE = "trans-QTL evaluation with anQChor"
"""Title of a single-run page (a comparison: ``Comparing N runs: trans-QTL evaluation …``)."""

TRANS_PAGE_JS = (
    resources.files("anqchor.report")
    .joinpath("static")
    .joinpath("trans_page.js")
    .read_text("utf-8")
)
"""The page script (``static/trans_page.js``)."""

MAX_GENES_PER_BIN = 50
"""Genes listed per bin in tooltips and the hotspot table (the strongest; counts are exact)."""

LAMBDA_GOOD = 1.05
LAMBDA_WARN = 1.1

CLASS_COLOR = {"cis": "var(--series-1)", "trans": "var(--series-2)"}
"""cis is blue and trans orange everywhere on the page."""


# ------------------------------------------------------------------- helpers
def short_trait(trait: str) -> str:
    """Gene name of ``ENSG…:NAME`` IDs, the ID otherwise."""
    s = str(trait)
    return s.split(":", 1)[1] if ":" in s and s.split(":", 1)[1] else s


def fmt_bp(bp: int) -> str:
    """``1000000`` → ``1 Mb``, ``500000`` → ``500 kb``."""
    if bp % 1_000_000 == 0:
        return f"{bp // 1_000_000} Mb"
    if bp >= 1_000_000:
        return f"{bp / 1e6:g} Mb"
    if bp % 1000 == 0:
        return f"{bp // 1000} kb"
    return f"{bp:,} bp"


def fmt_mb(bp: float) -> str:
    v = bp / 1e6
    return f"{v:.0f}" if abs(v - round(v)) < 1e-9 else f"{v:g}"


def lambda_status(lam: Any) -> str:
    """good up to 1.05, warning up to 1.1, serious above (deflation mirrored)."""
    try:
        f = float(lam)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(f):
        return ""
    d = max(f, 1 / f) if f > 0 else math.inf
    return "good" if d <= LAMBDA_GOOD else "warning" if d <= LAMBDA_WARN else "serious"


def _read_tsv(path: Path, **kw: Any) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def load_dataset(run: Path, ds: str) -> dict[str, Any]:
    """The summary files of one dataset."""

    def p(name: str) -> Path:
        return run / Layout.trans_summary(ds, name)

    str_cols = {
        "variant": str,
        "trait": str,
        "chrom": str,
        "trait_chrom": str,
        "top_trait": str,
        "top_variant": str,
        "variant_tested": str,
        "ref_trait": str,
        "run": str,
    }
    try:
        summary = read_json(p("summary.json"))
        if int(summary.get("schema_version", 0)) < 2:
            raise StepError(
                f"the trans summaries of dataset {ds} are from an older anqchor (schema "
                f"{summary.get('schema_version')}); re-run 'anqchor trans run' on this directory"
            )
        return {
            "summary": summary,
            "pairs": _read_tsv(p("reference_pairs.tsv"), dtype=str_cols),
            "bins": _read_tsv(p("hotspots.tsv"), dtype=str_cols),
            "genes": _read_tsv(p("hotspot_genes.tsv"), dtype=str_cols),
            "qq": _read_tsv(p("qq.tsv")),
        }
    except (OSError, ValueError) as e:
        raise StepError(f"cannot read the trans summaries of dataset {ds}: {e}") from e


# ------------------------------------------------------------- 1. replication
def replication_html(summary: dict[str, Any], window_bp: int) -> str:
    """Rows = cis and trans reference pairs; columns = same sign, p < 0.05, p < 1e-5 with
    bars on the scale of the class's tested pairs."""
    rep = summary.get("replication") or {}
    lv = rep.get("levels") or {"nominal": 0.05, "strict": 1e-5}
    cols = [("sign", "Same sign", "Effect in the reference's direction")]
    for name in ("nominal", "strict"):
        if name in lv:
            lab = f"p < {fmt_alpha(float(lv[name]))}"
            cols.append((name, lab, f"Same sign and {lab}"))
    head = (
        '<div class="hits-row hits-head"><span></span>'
        + "".join(f'<span class="hh" title="{esc(tip)}">{esc(t)}</span>' for _, t, tip in cols)
        + "</div>"
    )
    rows = []
    for cls, name in (("cis", "cis"), ("trans", "trans")):
        r = rep.get(cls) or {}
        n_pairs, n_t = int(r.get("n_pairs") or 0), int(r.get("n_tested") or 0)
        sub = f"{n_t:,} tested" if n_t == n_pairs else f"{n_t:,} of {n_pairs:,} tested"
        sub2 = ""
        if cls == "trans" and n_pairs:
            parts = []
            for key, what in (
                ("n_pairs_cross_chrom", "other chromosome"),
                ("n_pairs_same_chrom_distal", f"same, > {fmt_bp(window_bp)}"),
            ):
                if r.get(key):
                    parts.append(f"{int(r[key]):,} {what}")
            sub2 = " · ".join(parts)
        label = (
            f'<span class="hn"><span class="key" aria-hidden="true"></span>'
            f'<span class="hn-t">{esc(name)}<span class="hn-sub">{esc(sub if n_pairs else "")}</span>'
            + (f'<span class="hn-sub">{esc(sub2)}</span>' if sub2 else "")
            + "</span></span>"
        )
        if not n_pairs:
            rows.append(
                f'<div class="hits-row is-empty" style="--c:{CLASS_COLOR[cls]}">{label}'
                f'<span class="hits-none">No {name} pairs in the reference</span></div>'
            )
            continue
        cells = []
        for k, t, _ in cols:
            n = int(r.get(f"n_rep_{k}") or 0)
            w = 100.0 * n / n_t if n_t else 0.0
            cells.append(
                f'<span class="hc" title="{esc(name)}: {n:,} of {n_t:,} tested pairs, {esc(t)}">'
                f'<span class="hv"><b>{n:,}</b><span class="hp">{fmt_pct(n, n_t)}</span></span>'
                f'<span class="track"><span class="fill" style="width:{w:.1f}%"></span></span></span>'
            )
        rows.append(
            f'<div class="hits-row" style="--c:{CLASS_COLOR[cls]}">{label}{"".join(cells)}</div>'
        )
    return (
        f'<div class="hits rep-hits" style="--ncol:{len(cols)}">' + head + "".join(rows) + "</div>"
    )


# ---------------------------------------------------------------- 2. hotspots
def hotspot_data(
    summary: dict[str, Any], bins: pd.DataFrame, genes: pd.DataFrame
) -> dict[str, Any]:
    """Bins along the genome with their gene counts per threshold, the strongest genes per
    bin and the precomputed grid."""
    h = summary.get("hotspots") or {}
    ths = h.get("thresholds") or []
    bin_bp = int(h.get("bin_bp") or 1_000_000)
    chroms: list[dict[str, Any]] = []
    x: list[float] = []
    offset = 0.0
    gap = 0.0
    if len(bins):
        span = float(bins.groupby("chrom", sort=False)["bin_end"].max().sum())
        gap = max(bin_bp, 0.01 * span)
        for c, sub in bins.groupby("chrom", sort=False):
            end = float(sub["bin_end"].max())
            chroms.append({"name": str(c), "start": offset, "end": offset + end})
            x += (offset + (sub["bin_start"].to_numpy(dtype=float) + bin_bp / 2)).tolist()
            offset += end + gap
    lists: list[list[list[Any]]] = [[] for _ in range(len(bins))]
    if len(genes) and len(bins):
        key = {
            (str(c), int(s)): i
            for i, (c, s) in enumerate(zip(bins["chrom"], bins["bin_start"], strict=True))
        }
        g = genes.sort_values("log10p", ascending=False, kind="stable")
        for c, s, t, lp, v in zip(
            g["chrom"], g["bin_start"], g["trait"], g["log10p"], g["variant"], strict=True
        ):
            i = key.get((str(c), int(s)))
            if i is not None and len(lists[i]) < MAX_GENES_PER_BIN:
                lists[i].append([short_trait(t), float(f"{float(lp):.4g}"), str(v), str(t)])
    return {
        "binBp": bin_bp,
        "thresholds": [
            {
                "key": t["key"],
                "label": t["label"],
                "lp": t["log10p"],
                "short": "Bonferroni"
                if t["key"] == "bonferroni"
                else t["label"].replace("p < ", ""),
                "roles": t.get("roles") or [],
            }
            for t in ths
        ],
        "k": h.get("k") or [],
        "counts": h.get("counts") or [],
        "default": h.get("default") or {},
        "bins": {
            "chrom": bins["chrom"].astype(str).tolist() if len(bins) else [],
            "start": bins["bin_start"].astype(int).tolist() if len(bins) else [],
            "nVar": bins["n_variants"].astype(int).tolist() if len(bins) else [],
            "x": x,
            "n": {
                t["key"]: bins[f"n_genes_{t['key']}"].astype(int).tolist() if len(bins) else []
                for t in ths
            },
            "genes": lists,
        },
        "chroms": chroms,
    }


def hotspots_html(sfx: str, hd: dict[str, Any], summary: dict[str, Any]) -> str:
    ths, ks = hd["thresholds"], hd["k"]
    d = hd["default"]
    t_sel = d.get("threshold")
    k_sel = d.get("k")
    t_btn = "".join(
        f'<button type="button" data-v="{esc(t["key"])}" aria-pressed="{"true" if t["key"] == t_sel else "false"}"'
        f' title="{esc(t["label"])}">{esc(t["short"])}</button>'
        for t in ths
    )
    k_btn = "".join(
        f'<button type="button" data-v="{k}" aria-pressed="{"true" if k == k_sel else "false"}">{k}</button>'
        for k in ks
    )
    bin_txt = fmt_bp(int(hd["binBp"]))
    n_bins = int((summary.get("hotspots") or {}).get("n_bins") or 0)
    return (
        f'<section class="m-sec" id="hotspots{sfx}"><div class="m-sec-head"><h2>Hotspots</h2></div>'
        '<div class="hot-ctl">'
        f'<div class="ctl"><span class="ctl-l">Trans signal</span><span class="ctl-pre">p &lt;</span>'
        f'<div class="seg seg-sm" role="group" aria-label="p-value threshold" data-hot-t="{esc(sfx)}">{t_btn}</div></div>'
        f'<div class="ctl"><span class="ctl-l">Genes per bin</span><span class="ctl-pre">≥</span>'
        f'<div class="seg seg-sm" role="group" aria-label="Minimum genes" data-hot-k="{esc(sfx)}">{k_btn}</div></div>'
        "</div>"
        f'<div class="hot-hero" aria-live="polite"><b class="num" data-hot-n="{esc(sfx)}"></b>'
        f'<span data-hot-say="{esc(sfx)}"></span></div>'
        f'<div class="chart-card">{_empty_fig(f"fig-man{sfx}", 360)}</div>'
        f'<p class="m-foot">Each bar: one {esc(bin_txt)} bin of the {n_bins:,} with a tested variant; height: '
        "genes with a trans signal from any variant in the bin, genes within the cis window of the bin not "
        "counted. Dotted: the minimum.</p>"
        '<div class="hot-grid-wrap">'
        '<h3 class="m-h3">Hotspots by threshold</h3>'
        f'<div class="hot-grid-scroll"><table class="hot-grid" data-hot-grid="{esc(sfx)}"></table></div>'
        '<p class="m-foot">Rows: p-value threshold; columns: minimum genes per bin. Click a cell to show it.</p>'
        "</div>"
        f'<h3 class="m-h3">Hotspot bins <span class="cnt" data-hot-cnt="{esc(sfx)}"></span></h3>'
        f'<div class="hot-tbl-wrap"><table class="mtable hot-tbl" data-hot-tbl="{esc(sfx)}"><colgroup></colgroup>'
        "<thead></thead><tbody></tbody></table></div>"
        f'<div class="tfoot" data-hot-foot="{esc(sfx)}"></div>'
        "</section>"
    )


# ------------------------------------------------------------- 3. inflation
def inflation_html(sfx: str, summary: dict[str, Any], window_bp: int) -> str:
    inf = summary.get("inflation") or {}
    lam = inf.get("lambda_gc")
    st = lambda_status(lam)
    lam_txt = (
        f"{float(lam):.2f}" if isinstance(lam, (int, float)) and math.isfinite(float(lam)) else "–"
    )
    n = int(inf.get("n_cells") or 0)
    verdict = {
        "good": "calibrated",
        "warning": "slightly off",
        "serious": "inflated" if (isinstance(lam, (int, float)) and float(lam) > 1) else "deflated",
    }.get(st, "")
    tiles = [
        f'<div class="itile"><div class="it-v num">{esc(lam_txt)}'
        + (f'<span class="it-of"> {esc(verdict)}</span>' if verdict else "")
        + '</div><div class="it-l"><span>λ<sub>GC</sub> · 1 when calibrated</span></div></div>'
    ]
    for c in inf.get("counts") or []:
        exp = float(c["expected"])
        exp_txt = f"{exp:,.0f}" if exp >= 10 else f"{exp:.1f}" if exp >= 0.1 else f"{exp:.2g}"
        tiles.append(
            f'<div class="itile"><div class="it-v num">{int(c["observed"]):,}<span class="it-of"> vs {esc(exp_txt)}</span>'
            f'</div><div class="it-l">{esc(c["label"])} · observed vs expected</div></div>'
        )
    n_ex = int(inf.get("n_cells_excluded") or 0)
    n_pairs = int(inf.get("n_reference_trans_pairs") or 0)
    if n_pairs:
        ex = (
            f"Left out: {n_ex:,} tests of {int(inf.get('n_variants_excluded') or 0):,} variants within "
            f"{fmt_bp(window_bp)} of the {n_pairs:,} reference trans pairs, or of their "
            f"{int(inf.get('n_genes_excluded') or 0):,} genes."
        )
    else:
        ex = "None left out: the reference has no trans pairs."
    if n == 0:
        body = '<p class="hits-none">No tests across chromosomes: every gene is on its variant\'s chromosome.</p>'
    else:
        body = (
            '<div class="infl-grid">'
            f'<div class="chart-card qq-card">{_empty_fig(f"fig-qq{sfx}", 460)}</div>'
            f'<div class="itiles">{"".join(tiles)}</div></div>'
        )
    return (
        f'<section class="m-sec" id="inflation{sfx}"><div class="m-sec-head"><h2>Inflation</h2></div>'
        f'<p class="m-lead">{n:,} tests with the variant and the gene on different chromosomes. '
        "Almost all should be null.</p>"
        + body
        + f'<p class="m-foot">{esc(ex)} Grey band: 95% range under no association.</p></section>'
    )


def qq_data(qq: pd.DataFrame) -> dict[str, Any]:
    if qq.empty:
        return {}
    return {
        "exp": _r4(qq["expected_log10p"]),
        "obs": _r4(qq["observed_log10p"]),
        "lo": _r4(qq["band_lo_log10p"]),
        "hi": _r4(qq["band_hi_log10p"]),
        "n": int(qq["n_cells"].iloc[0]),
    }


# ------------------------------------------------------------- run details
def run_details_html(sfx: str, d: dict[str, Any], params: dict[str, Any]) -> str:
    s = d["summary"]
    pr = s.get("params") or {}
    window = int(pr.get("cis_window_bp") or 1_000_000)
    cells = s.get("cells") or {}
    hot = s.get("hotspots") or {}
    lv = pr.get("levels") or {}
    items: list[tuple[str, object]] = [
        ("Run", s.get("run", "")),
        ("Covariates", _cov_text(s.get("n_cov"), params.get("covariates") or "")),
        ("Phenotypes", "RINT after covariate adjustment" if s.get("rint") else "as given"),
        ("Samples", fmt_int(s.get("n_samples"))),
        (
            "Variants × genes",
            f"{fmt_int(s.get('n_variants_tested'))} × {fmt_int(s.get('n_traits_tested'))} "
            f"(listed {fmt_int(s.get('n_variants'))} × {fmt_int(s.get('n_traits'))})",
        ),
        (
            "cis cells",
            f"{fmt_int(cells.get('cis'))} · variant within {fmt_bp(window)} of the gene start",
        ),
        (
            "trans cells",
            f"{fmt_int(cells.get('cross_chrom'))} other chromosome · "
            f"{fmt_int(cells.get('same_chrom_distal'))} same chromosome, farther",
        ),
    ]
    if cells.get("no_locus"):
        items.append(("Gene without a locus", fmt_int(cells.get("no_locus"))))
    b = hot.get("bonferroni_p")
    if b:
        items.append(("Bonferroni", f"p < {_sci(float(b), 1)} = 0.05 / trans tests"))
    items.append(
        (
            "Hotspot bins",
            f"{fmt_bp(int(hot.get('bin_bp') or 1_000_000))}, {fmt_int(hot.get('n_bins'))} "
            "with a tested variant",
        )
    )
    items.append(
        (
            "Significant cells file",
            f"p < {fmt_alpha(float(lv.get('nominal', 1e-5)))}: "
            f"{fmt_int(s.get('n_sig_cells'))} (pairs.sig.tsv.gz)",
        )
    )
    items.append(("Minimum MAF", str(pr.get("min_maf"))))
    items.append(("Test", "qpgentools rect-assoc"))
    parts = [definition_list(items)]
    pairs = d["pairs"]
    if len(pairs):
        names = {
            "trait": "Gene",
            "variant": "Variant",
            "pair_class": "Class",
            "status": "Status",
            "distance": "Distance",
            "log10p": "−log10 p",
            "z": "z",
            "beta_ref": "Reference effect",
            "beta": "Effect",
            "sign_concordant": "Same sign",
        }
        cols = [c for c in names if c in pairs.columns]
        t = pairs[cols].rename(columns=names)
        parts.append(f"<h3>Reference pairs ({len(pairs):,})</h3>")
        parts.append(table(t, f"tbl-refpairs{sfx}", page_size=500, digits=3))
    rel = Layout.trans_dir(s.get("dataset", "dataset")).as_posix()
    links = " · ".join(
        f'<a href="../{rel}/{n}" download>{n}</a>'
        for n in (
            "summary.json",
            "reference_pairs.tsv",
            "hotspots.tsv",
            "hotspot_genes.tsv",
            "qq.tsv",
            "variants.tsv",
            "pairs.sig.tsv.gz",
        )
    )
    parts.append(
        f'<p class="rd-note">Tables: {links} · <a href="../anqchor.resolved.yaml">settings</a></p>'
    )
    return _details("Run details", "".join(parts), f"details{sfx}")


# ---------------------------------------------------------------------- page
def context_html(d: dict[str, Any], params: dict[str, Any]) -> str:
    s = d["summary"]
    cells = s.get("cells") or {}
    bits = [
        f"{fmt_int(s.get('n_variants_tested'))} variants × {fmt_int(s.get('n_traits_tested'))} genes",
        f"{fmt_int(s.get('n_cells'))} variant–gene tests, {fmt_int(cells.get('cross_chrom'))} across chromosomes",
    ]
    if s.get("n_samples"):
        bits.append(f"{fmt_int(s.get('n_samples'))} samples")
    pills = [
        f'<span class="pill">Run {esc(s.get("run", ""))}{" · RINT" if s.get("rint") else ""}</span>'
    ]
    pills.append(
        f'<span class="pill">Covariates {esc(_cov_text(s.get("n_cov"), params.get("covariates") or ""))}'
        "</span>"
    )
    if params.get("reference_file"):
        pills.append(
            f'<span class="pill" title="reference pairs">{esc(params["reference_file"])}</span>'
        )
    alerts = "".join(
        f'<p class="m-alert" role="note">⚠ {esc(w)}</p>' for w in s.get("warnings") or []
    )
    return (
        f'<div class="ctx"><span class="num">{esc(" · ".join(bits))}</span>{"".join(pills)}</div>'
        + alerts
    )


def panel(
    ds: str, sfx: str, d: dict[str, Any], params: dict[str, Any]
) -> tuple[str, dict[str, Any]]:
    s = d["summary"]
    window = int((s.get("params") or {}).get("cis_window_bp") or 1_000_000)
    hd = hotspot_data(s, d["bins"], d["genes"])
    out = [
        f'<section class="m-sec" id="replication{sfx}"><div class="m-sec-head"><h2>Replicated reference pairs</h2>'
        "</div>"
        + replication_html(s, window)
        + f'<p class="m-foot">cis: same chromosome, variant within {esc(fmt_bp(window))} of the gene start; '
        "trans: all other pairs. Every column needs the same sign as the reference.</p></section>",
        hotspots_html(sfx, hd, s),
        inflation_html(sfx, s, window),
        run_details_html(sfx, d, params),
    ]
    entry = {
        "sfx": sfx,
        "hot": hd,
        "qq": qq_data(d["qq"]),
        "lambda": (s.get("inflation") or {}).get("lambda_gc"),
        "colors": {
            "trans": "#eb6834",
            "cis": "#2a78d6",
            "muted": LIGHT["muted"],
            "axis": LIGHT["axis"],
            "ink": LIGHT["ink"],
            "ink2": LIGHT["ink2"],
            "grid": LIGHT["grid"],
            "surface": LIGHT["surface"],
            "critical": STATUS["critical"],
        },
        "layout": chart_layout(),
    }
    return "".join(out), entry


def build_page(datasets: list[tuple[str, dict[str, Any]]], params: dict[str, Any]) -> Page:
    multi = len(datasets) > 1
    bodies, entries, header = [], [], []
    for ds, d in datasets:
        sfx = f"-{ds}" if multi else ""
        html, entry = panel(ds, sfx, d, params)
        header.append(context_html(d, params))
        bodies.append(
            f'<div class="ds-block" id="ds{sfx}"><h2 class="ds-title">{esc(ds)}</h2>{html}</div>'
            if multi
            else html
        )
        entries.append(entry)
    return Page(
        title=TITLE,
        sections=[Section(id="trans", title="", body="".join(bodies))],
        nav=index_nav(params),
        extra_js=TRANS_PAGE_JS,
        data={"panels": entries},
        header_html="".join(header),
        hide_brand=True,
        body_class="page-match page-cis page-trans",
        toc=False,
    )


# ---------------------------------------------------------------- comparison
def load_comparison(cdir: Path) -> dict[str, Any]:
    """The tables of a ``trans compare`` directory."""
    T = Layout.compare_table
    try:
        summary = read_json(cdir / T("summary.json"))
        if summary.get("kind") != "comparison" or summary.get("pipeline") != "trans":
            raise StepError(f"{cdir} is not a trans comparison directory")
        if int(summary.get("schema_version", 0)) < 2:
            raise StepError(
                f"the comparison in {cdir} is from an older anqchor; run 'anqchor trans compare' again"
            )
        return {
            "summary": summary,
            "runs": _read_tsv(cdir / T("runs.tsv"), dtype={"label": str}),
            "bins": _read_tsv(cdir / T("hotspots.tsv"), dtype={"run": str, "chrom": str}),
            "genes": _read_tsv(
                cdir / T("hotspot_genes.tsv"),
                dtype={"run": str, "chrom": str, "trait": str, "variant": str},
            ),
            "qq": _read_tsv(cdir / T("qq.tsv"), dtype={"run": str}),
            "pairs": _read_tsv(
                cdir / T("reference_pairs.tsv"),
                dtype={"run": str, "trait": str, "variant": str, "ref_trait": str},
            ),
        }
    except (OSError, ValueError) as e:
        raise StepError(f"cannot read the trans comparison in {cdir}: {e}") from e


def _of_run(df: pd.DataFrame, label: str) -> pd.DataFrame:
    if df.empty or "run" not in df:
        return pd.DataFrame()
    return df[df["run"].astype(str) == label].drop(columns=["run"]).reset_index(drop=True)


def cmp_replication_html(runs: pd.DataFrame, labels: list[str]) -> str:
    """Two blocks, cis and trans; in each, runs as rows and same sign, p < 0.05, p < 1e-5 as
    columns with bars on the scale of the tested pairs. When the runs tested different pairs
    of a class, that block counts only the pairs tested in every run (``*_common``)."""
    cols = [
        ("sign", "Same sign", "Effect in the reference's direction"),
        ("nominal", "p < 0.05", "Same sign and p < 0.05"),
        ("strict", "p < 1e-5", "Same sign and p < 1e-5"),
    ]
    rows_by = {str(r["label"]): r for _, r in runs.iterrows()}
    blocks = []
    for cls in ("cis", "trans"):
        n_pairs = [int(rows_by[lab].get(f"{cls}_n_pairs") or 0) for lab in labels]
        n_tested = [int(rows_by[lab].get(f"{cls}_n_tested") or 0) for lab in labels]
        if not any(n_pairs):
            blocks.append(
                f'<div class="rep-blk"><h3 class="m-h3">{cls}</h3>'
                f'<div class="hits"><p class="hits-none rep-none">No {cls} pairs in the reference</p></div></div>'
            )
            continue
        n_common = [int(rows_by[lab].get(f"{cls}_n_common") or 0) for lab in labels]
        hi = max(n_tested)
        use_common = len(labels) > 1 and 0 < min(n_common) < hi
        sfx_col = "_common" if use_common else ""
        lo = min(n_tested)
        own = f"{hi:,} tested" if lo == hi else f"{lo:,}–{hi:,} tested"
        tested = f"{min(n_common):,} tested in every run" if use_common else own
        head = (
            '<div class="hits-row hits-head"><span></span>'
            + "".join(f'<span class="hh" title="{esc(tip)}">{esc(t)}</span>' for _, t, tip in cols)
            + "</div>"
        )
        body = []
        for i, lab in enumerate(labels):
            r = rows_by[lab]
            n_t = int(r.get(f"{cls}_n_common" if use_common else f"{cls}_n_tested") or 0)
            cells = []
            for k, t, _ in cols:
                n = int(r.get(f"{cls}_n_rep_{k}{sfx_col}") or 0)
                w = 100.0 * n / n_t if n_t else 0.0
                cells.append(
                    f'<span class="hc" title="{esc(lab)}: {n:,} of {n_t:,} tested {cls} pairs, {esc(t)}">'
                    f'<span class="hv"><b>{n:,}</b><span class="hp">{fmt_pct(n, n_t)}</span></span>'
                    f'<span class="track"><span class="fill" style="width:{w:.1f}%"></span></span></span>'
                )
            body.append(
                f'<div class="hits-row" style="--c:{_series_var(i)}">'
                f'<span class="hn" title="{esc(lab)}"><span class="key" aria-hidden="true"></span>'
                f'<span class="hn-l">{esc(lab)}</span></span>' + "".join(cells) + "</div>"
            )
        blocks.append(
            f'<div class="rep-blk"><h3 class="m-h3">{cls} <span class="cnt">{esc(tested)}</span></h3>'
            f'<div class="hits" style="--ncol:{len(cols)}">' + head + "".join(body) + "</div></div>"
        )
    return '<div class="rep-cmp">' + "".join(blocks) + "</div>"


def cmp_hotspots_html(
    sfx: str,
    thresholds: list[dict[str, Any]],
    ks: list[int],
    default: dict[str, Any],
    labels: list[str],
    bin_txt: str,
) -> str:
    t_sel, k_sel = default.get("threshold"), default.get("k")
    t_btn = "".join(
        f'<button type="button" data-v="{esc(t["key"])}" aria-pressed="{"true" if t["key"] == t_sel else "false"}"'
        f' title="{esc(t["label"])}">{esc(t["short"])}</button>'
        for t in thresholds
    )
    k_btn = "".join(
        f'<button type="button" data-v="{k}" aria-pressed="{"true" if k == k_sel else "false"}">{k}</button>'
        for k in ks
    )
    first = labels[0] if labels else ""
    if len(labels) > 5 or sum(len(lab) for lab in labels) > 40:
        pick = (
            f'<label class="run-pick">Run <select data-pick="hot{esc(sfx)}">'
            + "".join(f'<option value="{i}">{esc(lab)}</option>' for i, lab in enumerate(labels))
            + "</select></label>"
        )
    else:
        pick = _seg(f"hot{sfx}", [(str(i), lab) for i, lab in enumerate(labels)], "Run", "0")
    return (
        f'<section class="m-sec" id="hotspots{sfx}"><div class="m-sec-head"><h2>Hotspots</h2></div>'
        '<div class="hot-ctl">'
        f'<div class="ctl"><span class="ctl-l">Trans signal</span><span class="ctl-pre">p &lt;</span>'
        f'<div class="seg seg-sm" role="group" aria-label="p-value threshold" data-hot-t="{esc(sfx)}">{t_btn}</div></div>'
        f'<div class="ctl"><span class="ctl-l">Genes per bin</span><span class="ctl-pre">≥</span>'
        f'<div class="seg seg-sm" role="group" aria-label="Minimum genes" data-hot-k="{esc(sfx)}">{k_btn}</div></div>'
        "</div>"
        f'<p class="hot-say" data-hot-say="{esc(sfx)}" aria-live="polite"></p>'
        f'<div class="hits hot-runs" data-hot-runs="{esc(sfx)}"></div>'
        '<p class="m-foot">Genes near the bin (cis window) are not counted. Bar: share of bins with a '
        "tested variant. Click a run to see its bins.</p>"
        '<div class="m-sec-head m-h3-head"><h3 class="m-h3">Bins along the genome '
        f'<span class="cnt" data-hot-bin="{esc(sfx)}">{esc(bin_txt)}</span></h3>{pick}</div>'
        f'<div class="chart-card">{_empty_fig(f"fig-man{sfx}", 360)}</div>'
        '<p class="m-foot">Each bar: one bin; height: genes with a trans signal. Dotted: the minimum.</p>'
        '<div class="hot-grid-wrap">'
        f'<h3 class="m-h3">Hotspots by threshold <span class="cnt" data-hot-run="{esc(sfx)}">{esc(first)}</span></h3>'
        f'<div class="hot-grid-scroll"><table class="hot-grid" data-hot-grid="{esc(sfx)}"></table></div>'
        '<p class="m-foot">Rows: p-value threshold; columns: minimum genes per bin. Click a cell to apply it to every run.</p>'
        "</div>"
        f'<h3 class="m-h3">Hotspot bins <span class="cnt" data-hot-cnt="{esc(sfx)}"></span> '
        f'<span class="cnt" data-hot-run="{esc(sfx)}">{esc(first)}</span></h3>'
        f'<div class="hot-tbl-wrap"><table class="mtable hot-tbl" data-hot-tbl="{esc(sfx)}"><colgroup></colgroup>'
        "<thead></thead><tbody></tbody></table></div>"
        f'<div class="tfoot" data-hot-foot="{esc(sfx)}"></div>'
        "</section>"
    )


def cmp_inflation_html(
    sfx: str, runs: pd.DataFrame, labels: list[str], per_run: dict[str, Any]
) -> str:
    rows_by = {str(r["label"]): r for _, r in runs.iterrows()}
    n_cells = [int(rows_by[lab].get("n_inflation_cells") or 0) for lab in labels]
    lo, hi = (min(n_cells), max(n_cells)) if n_cells else (0, 0)
    n_txt = f"{hi:,}" if lo == hi else f"{lo:,}–{hi:,}"
    # observed vs expected: one column per p-value of the runs' inflation counts
    alphas: list[tuple[float, str]] = []
    for lab in labels:
        for c in (per_run.get(lab, {}).get("inflation") or {}).get("counts") or []:
            a = float(c["alpha"])
            if a not in [x for x, _ in alphas]:
                alphas.append((a, str(c["label"])))
    alphas.sort(key=lambda x: -x[0])
    head = (
        '<div class="hits-row hits-head"><span></span>'
        + "".join(f'<span class="hh">p &lt; {esc(_sci(a))}</span>' for a, _ in alphas)
        + "</div>"
    )
    body = []
    for i, lab in enumerate(labels):
        counts = {
            float(c["alpha"]): c
            for c in (per_run.get(lab, {}).get("inflation") or {}).get("counts") or []
        }
        cells = []
        for a, t in alphas:
            c = counts.get(a)
            if c is None:
                cells.append('<span class="hc"><span class="hv"><b>–</b></span></span>')
                continue
            obs, exp = int(c["observed"]), float(c["expected"])
            exp_txt = f"{exp:,.0f}" if exp >= 100 else f"{exp:.1f}" if exp >= 0.1 else f"{exp:.2g}"
            fold = (
                obs / exp
                if exp > 0 and (exp >= MIN_EXPECTED or obs >= 2 * MIN_EXPECTED)
                else math.nan
            )
            if not math.isfinite(fold):
                fold_txt = ""
            elif obs == 0:
                fold_txt = "×0"
            else:
                fold_txt = f"×{fold:,.0f}" if fold >= 10 else f"×{fold:.1f}"
            hot = math.isfinite(fold) and fold > 2 and obs >= 5
            cells.append(
                f'<span class="hc" title="{esc(lab)}: {obs:,} tests at {esc(t)}; {esc(exp_txt)} expected">'
                f'<span class="hv"><b>{obs:,}</b><span class="hp">vs {esc(exp_txt)}</span>'
                + (
                    f'<span class="fold{" bad" if hot else ""}">{esc(fold_txt)}</span>'
                    if fold_txt
                    else ""
                )
                + "</span></span>"
            )
        body.append(
            f'<div class="hits-row" style="--c:{_series_var(i)}">'
            f'<span class="hn" title="{esc(lab)}"><span class="key" aria-hidden="true"></span>'
            f'<span class="hn-l">{esc(lab)}</span></span>' + "".join(cells) + "</div>"
        )
    oe = (
        f'<div class="hits oe-hits" style="--ncol:{max(1, len(alphas))}">'
        + head
        + "".join(body)
        + "</div>"
        if alphas
        else ""
    )
    excl = [int(rows_by[lab].get("n_inflation_cells_excluded") or 0) for lab in labels]
    ex = (
        "None left out: the reference has no trans pairs."
        if not any(excl)
        else "Left out: the tests near the reference trans pairs ("
        + ", ".join(f"{lab} {n:,}" for lab, n in zip(labels, excl, strict=True))
        + ")."
    )
    n_runs = len(labels)
    if hi == 0:
        body_html = '<p class="hits-none">No tests across chromosomes: every gene is on its variant\'s chromosome.</p>'
    else:
        body_html = (
            '<h3 class="m-h3">λ<sub>GC</sub></h3>'
            f'<div class="chart-card">{_empty_fig(f"fig-lambda{sfx}", 72 + 48 * n_runs)}</div>'
            '<p class="m-foot">1 = calibrated (dotted). Above 1.1 is flagged.</p>'
            '<div class="m-sec-head m-h3-head"><h3 class="m-h3">QQ plot</h3>'
            + _seg(f"qq{sfx}", [("fit", "Zoom"), ("full", "Full")], "QQ plot range", "fit")
            + "</div>"
            f'<div class="chart-card qq-card">{_empty_fig(f"fig-qq{sfx}", 480)}</div>'
            '<p class="m-foot">Observed against expected −log10 p, one line per run. Grey band: 95% range with '
            "no association. Click a run in the legend to hide it.</p>"
            '<h3 class="m-h3">Observed vs expected</h3>'
            + oe
            + f'<p class="m-foot">Ratio shown when at least {MIN_EXPECTED} are expected or {2 * MIN_EXPECTED} observed.</p>'
        )
    return (
        f'<section class="m-sec" id="inflation{sfx}"><div class="m-sec-head"><h2>Inflation</h2></div>'
        f'<p class="m-lead">{n_txt} tests with the variant and the gene on different chromosomes. '
        "Almost all should be null.</p>" + body_html + f'<p class="m-foot">{esc(ex)}</p></section>'
    )


def cmp_details_html(sfx: str, c: dict[str, Any]) -> str:
    s, runs = c["summary"], c["runs"]
    per_run = s.get("per_run") or {}
    rows = []
    for i, (_, r) in enumerate(runs.iterrows()):
        lab = str(r["label"])
        pr = per_run.get(lab) or {}
        b = _num_or_none(r.get("bonferroni_p"))
        rows.append(
            {
                "Run": lab,
                "Covariates": _cov_text(
                    r.get("n_cov"),
                    ""
                    if str(r.get("cov_include")) in ("nan", "None")
                    else str(r.get("cov_include") or ""),
                ),
                "Rank-normalized": "yes" if str(r.get("rint")) in ("1", "True", "true") else "no",
                "Samples": fmt_int(r.get("n_samples")),
                "Variants × genes": f"{fmt_int(r.get('n_variants_tested'))} × {fmt_int(r.get('n_traits_tested'))}",
                "Tests": fmt_int(r.get("n_cells")),
                "Across chromosomes": fmt_int(r.get("n_cells_cross_chrom")),
                "Bonferroni": f"p < {_sci(b, 1)}" if b else "",
                "Bin": fmt_bp(int((pr.get("hotspots") or {}).get("bin_bp") or 1_000_000)),
                "Settings": Layout.compare_settings(lab).as_posix(),
                "_i": i,
            }
        )
    text_cols = (
        "Run",
        "Covariates",
        "Rank-normalized",
        "Bin",
        "Settings",
        "Bonferroni",
        "Variants × genes",
    )
    show = [k for k in rows[0] if k != "_i"] if rows else []
    head = "".join(
        f"<th>{esc(k)}</th>" if k in text_cols else f'<th class="num">{esc(k)}</th>' for k in show
    )
    body = ""
    for rec in rows:
        cells = []
        for k in show:
            v = str(rec[k])
            if k == "Run":
                cells.append(
                    f'<td class="rn" style="--c:{_series_var(int(rec["_i"]))}" title="{esc(v)}">'
                    f'<span class="key" aria-hidden="true"></span>{esc(v)}</td>'
                )
            elif k == "Settings":
                cells.append(f'<td><a href="../{esc(v)}" title="{esc(Path(v).name)}">yaml</a></td>')
            else:
                cells.append(
                    f"<td>{esc(v)}</td>" if k in text_cols else f'<td class="num">{esc(v)}</td>'
                )
        body += "<tr>" + "".join(cells) + "</tr>"
    parts = [
        '<div class="rd-scroll"><table class="codes runs-t"><thead><tr>'
        + head
        + "</tr></thead><tbody>"
        + body
        + "</tbody></table></div>"
    ]
    first: dict[str, Any] = next(iter(per_run.values()), {}) if per_run else {}
    window = int((first.get("params") or {}).get("cis_window_bp") or 1_000_000)
    parts.append(
        definition_list(
            [
                ("cis", f"same chromosome, variant within {fmt_bp(window)} of the gene start"),
                ("trans", "every other pair: another chromosome, or the same chromosome farther"),
            ]
        )
    )
    pairs = c["pairs"]
    labels = [str(x) for x in runs["label"]]
    if len(pairs):
        wide = _pairs_wide(pairs, labels)
        parts.append(f"<h3>Reference pairs ({len(wide):,})</h3>")
        parts.append(
            '<p class="rd-note">−log10 p per run; negative: opposite sign to the reference.</p>'
        )
        parts.append(table(wide, f"tbl-refpairs{sfx}", page_size=500, digits=6))
    links = " · ".join(
        f'<a href="../{n}" download>{n}</a>'
        for n in (
            "runs.tsv",
            "hotspot_grid.tsv",
            "hotspots.tsv",
            "hotspot_genes.tsv",
            "qq.tsv",
            "reference_pairs.tsv",
            "summary.json",
        )
    )
    parts.append(f'<p class="rd-note">Tables: {links}</p>')
    return _details("Run details", "".join(parts), f"details{sfx}")


MIN_EXPECTED = 5
"""Observed / expected ratios are shown only when at least this many tests are expected (or
twice as many observed): a ratio over an expected count below 1 is noise."""


def _sci(a: float, digits: int = 0) -> str:
    """``0.001`` → ``1e-3``; ``4.139e-07`` with ``digits=1`` → ``4.1e-7``."""
    m, x = f"{a:.{digits}e}".split("e")
    return f"{m}e{int(x)}"


def _pairs_wide(pairs: pd.DataFrame, labels: list[str]) -> pd.DataFrame:
    """Reference pairs × runs: signed −log10 p per run (negative: opposite sign), matched by
    ``(ref_trait, variant)``; sorted by the first run's value, largest first."""
    key = ["ref_trait", "variant"]
    p = pairs.copy()
    if "ref_trait" not in p:
        p["ref_trait"] = p["trait"]
    lp = pd.to_numeric(p["log10p"], errors="coerce")
    opp = p["sign_concordant"].astype(str).isin(["0", "0.0", "False"])
    p["v"] = np.where(opp, -lp, lp).round(1)
    first = p.drop_duplicates(key)[[*key, "trait"]]
    cls = (
        p.groupby(key, sort=False)["pair_class"]
        .agg(lambda x: "/".join(sorted({str(v) for v in x})))
        .reset_index()
    )
    wide = first.merge(cls, on=key, how="left")
    for lab in labels:
        sub = p[p["run"] == lab].drop_duplicates(key)[[*key, "v"]].rename(columns={"v": lab})
        wide = wide.merge(sub, on=key, how="left")
    if labels:
        wide = wide.sort_values(labels[0], ascending=False, na_position="last", kind="stable")
    wide["trait"] = wide["trait"].map(short_trait)
    return wide.drop(columns=["ref_trait"]).rename(
        columns={"trait": "Gene", "variant": "Variant", "pair_class": "Class"}
    )[["Gene", "Variant", "Class", *labels]]


def _num_or_none(v: Any) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def cmp_context_html(c: dict[str, Any]) -> str:
    s, runs = c["summary"], c["runs"]

    def span(col: str) -> str:
        v = [int(x) for x in pd.to_numeric(runs[col], errors="coerce").dropna()]
        if not v:
            return "NA"
        return f"{min(v):,}" if min(v) == max(v) else f"{min(v):,}–{max(v):,}"

    bits = [
        f"{len(runs)} runs",
        f"{span('n_variants_tested')} variants × {span('n_traits_tested')} genes",
        f"{span('n_cells')} variant–gene tests, {span('n_cells_cross_chrom')} across chromosomes",
        f"{span('n_samples')} samples",
    ]
    refs = sorted(
        {Path(str(x)).name for x in runs.get("reference_pairs", []) if str(x) not in ("", "nan")}
    )
    pills = [f'<span class="pill" title="reference pairs">{esc(r)}</span>' for r in refs[:3]]
    alerts = "".join(
        f'<p class="m-alert" role="note">⚠ {esc(w)}</p>' for w in s.get("warnings") or []
    )
    return (
        f'<div class="ctx"><span class="num">{esc(" · ".join(bits))}</span>{"".join(pills)}</div>'
        + alerts
    )


def comparison_panel(c: dict[str, Any], params: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    """The comparison page's sections and its ``ANQ_DATA.comparison`` entry."""
    s, runs = c["summary"], c["runs"]
    labels = [str(x) for x in runs["label"]]
    per_run = s.get("per_run") or {}
    colors = run_colors(labels)
    sfx = "-cmp"
    # one set of controls for every run: the union of the runs' thresholds and minimums
    ths = [
        {
            "key": t["key"],
            "label": "Bonferroni" if t["key"] == "bonferroni" else t["label"],
            "short": "Bonferroni"
            if t["key"] == "bonferroni"
            else str(t["label"]).replace("p < ", ""),
            "lp": t["log10p"],
        }
        for t in s.get("thresholds") or []
    ]
    ks = sorted(
        {int(k) for pr in per_run.values() for k in (pr.get("hotspots") or {}).get("k") or []}
    )
    d0: dict[str, Any] = (
        ((per_run.get(labels[0]) or {}).get("hotspots") or {}).get("default") or {}
        if labels
        else {}
    )
    default: dict[str, Any] = {"threshold": d0.get("threshold"), "k": d0.get("k")}
    rlist = []
    for i, lab in enumerate(labels):
        pr = per_run.get(lab) or {}
        hd = hotspot_data(pr, _of_run(c["bins"], lab), _of_run(c["genes"], lab))
        rlist.append(
            {
                "label": lab,
                "color": colors[lab],
                "css": _series_var(i),
                "dash": run_dash(i),
                "hot": hd,
                "nBins": int((pr.get("hotspots") or {}).get("n_bins") or 0),
                "qq": qq_data(_of_run(c["qq"], lab)),
                "lambda": (pr.get("inflation") or {}).get("lambda_gc"),
            }
        )
    bin_bps = {
        int((pr.get("hotspots") or {}).get("bin_bp") or 1_000_000) for pr in per_run.values()
    }
    bin_txt = f"{fmt_bp(next(iter(bin_bps)))} bins" if len(bin_bps) == 1 else ""
    entry = {
        "sfx": sfx,
        "runs": rlist,
        "sameBin": len(bin_bps) <= 1,
        "thresholds": ths,
        "k": ks,
        "default": default,
        "colors": {
            "muted": LIGHT["muted"],
            "axis": LIGHT["axis"],
            "ink": LIGHT["ink"],
            "ink2": LIGHT["ink2"],
            "grid": LIGHT["grid"],
            "surface": LIGHT["surface"],
            "critical": STATUS["critical"],
        },
        "layout": chart_layout(),
    }
    out = [
        f'<section class="m-sec" id="replication{sfx}"><div class="m-sec-head"><h2>Replicated reference pairs</h2>'
        "</div>"
        + cmp_replication_html(runs, labels)
        + '<p class="m-foot">cis: same chromosome, variant near the gene start; trans: all other pairs. '
        "Every column needs the same sign as the reference; % of the run's tested pairs.</p></section>",
        cmp_hotspots_html(sfx, ths, ks, default, labels, bin_txt),
        cmp_inflation_html(sfx, runs, labels, per_run),
        cmp_details_html(sfx, c),
    ]
    return "".join(out), entry


def build_comparison_page(c: dict[str, Any], params: dict[str, Any]) -> Page:
    html, entry = comparison_panel(c, params)
    n = len(entry["runs"])
    return Page(
        title=f"Comparing {n} runs: {TITLE}" if n > 1 else TITLE,
        sections=[Section(id="trans", title="", body=html)],
        nav=[],
        extra_js=TRANS_PAGE_JS,
        data={"panels": [], "comparison": entry},
        header_html=cmp_context_html(c),
        hide_brand=True,
        body_class="page-match page-cis page-trans",
        toc=False,
    )


@register("report_trans")
def report_trans(ctx: StepContext) -> None:
    """Render ``report/trans.html`` from a run's trans summaries, or (``comparison``) from
    the tables of a ``trans compare`` directory."""
    p = ctx.params
    if p.get("comparison"):
        page = build_comparison_page(load_comparison(ctx.run_dir), p)
    else:
        datasets = [(ds, load_dataset(ctx.run_dir, ds)) for ds in p["datasets"]]
        page = build_page(datasets, p)
    write_page(ctx.path(Layout.report("trans")), page, p.get("plotlyjs", "embed"))


__all__ = [
    "TITLE",
    "build_comparison_page",
    "build_page",
    "hotspot_data",
    "lambda_status",
    "load_comparison",
    "load_dataset",
    "report_trans",
    "short_trait",
]
