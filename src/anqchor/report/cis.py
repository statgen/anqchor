"""``report.cis`` (kind ``report_cis``): the cis replication report (``report/cis.html``).

Built only from the summaries of ``summarize_cis`` (``cis/<dataset>/…``) or of
``cis compare``, never from qpgen output. One page shows one run or a comparison of
runs (the first run is the baseline). From top to bottom:

1. header: title, a context line (pairs, samples, reference file) and qpgen/summary
   warnings;
2. power curve: reference pairs replicated (same sign and ``p < 10^-t``, as the levels)
   against the threshold ``t`` on a log1p axis with natural tick labels, one line per
   run, the nominal and strict thresholds marked;
3. signal strength: per run the distribution of ``|z|`` or ``-log10 p`` (the reference
   study as a hollow row) over the pairs tested in every run and, for a comparison, the
   effective sample size relative to the baseline (``rel_neff`` from ``strategies.tsv``,
   with a 95 % paired bootstrap interval over the same pairs) on a log axis around 1;
4. significant hits: per run (rows) the pairs with the same sign, ``p < 0.05`` and
   ``p < 1e-5`` (the configured nominal and strict alphas), as in-row bars;
5. agreement with the reference: for one run at a time (selector), effect size, z-score
   and ``-log10 p`` against the reference study, and z against the baseline run;
6. all pairs (sortable, searchable, 100 rows a page, header and first column pinned) and
   the run details (collapsed).

The page script (``static/cis_page.js``) draws the charts from compact per-pair arrays
in ``ANQ_DATA.panels``. Runs keep one color everywhere (fixed categorical order).
"""

from __future__ import annotations

import math
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from anqchor.pipelines.layout import Layout
from anqchor.report.render import (
    Page,
    Section,
    definition_list,
    esc,
    figure,
    index_nav,
    table,
    write_page,
)
from anqchor.report.theme import CATEGORICAL_LIGHT, LIGHT, STATUS, base_layout, series_color
from anqchor.steps import StepContext, StepError, register
from anqchor.util.files import read_json

TITLE = "cis-QTL replication with anQChor"
"""Title of a single-run page (a comparison: ``Comparing N runs: cis-QTL replication …``)."""

MAX_SCATTER_POINTS = 20000
"""Upper bound on pairs drawn in the agreement charts (pairs with an opposite sign in any
run are kept first, then the most significant in the reference; the page says when
pairs were not drawn)."""

MAX_TABLE_ROWS = 20000
"""Rows embedded in the "All pairs" table (the page says when there are more; the full
table is ``pairs.long.tsv.gz``)."""

CURVE_POINTS = 400
"""Thresholds per power curve (evenly spaced on the log1p axis, plus the level cutoffs)."""

BOOTSTRAP_REPS = 1000
"""Paired bootstrap replicates of the effective sample size ratio."""

CIS_PAGE_JS = (
    resources.files("anqchor.report").joinpath("static").joinpath("cis_page.js").read_text("utf-8")
)
"""The page script (``static/cis_page.js``)."""


# ------------------------------------------------------------------- helpers
DASHES = ["dash", "dot", "dashdot", "longdash", "longdashdot"]
"""Line dashes of the runs past the eighth (two greys): color is never the only cue."""

EXTRA_GREYS = ("muted", "ink2")
"""Theme colors alternated by the runs past the eighth (with :data:`DASHES`: 10 styles)."""


def run_colors(runs: list[str]) -> dict[str, str]:
    """Run → light-mode color in order; past the eighth, two alternating greys (dashed)."""
    n = len(CATEGORICAL_LIGHT)
    return {
        s: series_color(i) if i < n else LIGHT[EXTRA_GREYS[(i - n) % 2]] for i, s in enumerate(runs)
    }


def run_dash(i: int) -> str:
    """Line dash of run ``i``: solid for the first eight, then one dash per grey pair."""
    n = len(CATEGORICAL_LIGHT)
    return "solid" if i < n else DASHES[((i - n) // 2) % len(DASHES)]


def _series_var(i: int) -> str:
    """CSS color of run ``i`` (follows the light/dark theme)."""
    n = len(CATEGORICAL_LIGHT)
    return f"var(--series-{i + 1})" if i < n else f"var(--{EXTRA_GREYS[(i - n) % 2]})"


def fmt_int(v: Any) -> str:
    try:
        return f"{int(v):,}"
    except (TypeError, ValueError):
        return "NA"


def fmt_pct(n: float, d: float) -> str:
    return f"{100.0 * n / d:.0f}%" if d else ""


def fmt_alpha(a: float) -> str:
    """``1e-05`` → ``1e-5``, ``0.05`` → ``0.05``."""
    if a >= 1e-3:
        return f"{a:g}"
    m, e = f"{a:.0e}".split("e")
    return f"{m}e{int(e)}"


def _num(v: Any) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _r4(values: Any) -> list[float | None]:
    """Floats rounded to 4 significant digits (``None`` for missing): compact JSON."""
    out: list[float | None] = []
    for v in np.asarray(values, dtype=float).tolist():
        out.append(float(f"{v:.4g}") if math.isfinite(v) else None)
    return out


def load_dataset(run: Path, ds: str, summary_dir: str | None = None) -> dict[str, Any]:
    """The five summary files of one dataset (from ``summary_dir`` when given: a
    comparison directory keeps them at its top level)."""

    def p(name: str) -> Path:
        return (
            (run / summary_dir / name)
            if summary_dir is not None
            else run / Layout.cis_summary(ds, name)
        )

    href = (
        f"../{summary_dir}/strategies.tsv".replace("/./", "/")
        if summary_dir is not None
        else f"../{Layout.cis_summary(ds, 'strategies.tsv').as_posix()}"
    )
    try:
        return {
            "strategies_href": href,
            "summary": read_json(p("summary.json")),
            "strategies": pd.read_csv(p("strategies.tsv"), sep="\t"),
            "traits": pd.read_csv(p("traits.tsv"), sep="\t"),
            "curve": pd.read_csv(p("power_curve.tsv"), sep="\t"),
            "long": pd.read_csv(
                p("pairs.long.tsv.gz"),
                sep="\t",
                dtype={"trait": str, "ref_trait": str, "variant": str},
            ),
        }
    except (OSError, ValueError) as e:
        raise StepError(f"cannot read the cis summaries of dataset {ds}: {e}") from e


def level_cutoffs(summary: dict[str, Any]) -> dict[str, tuple[float, str]]:
    """``nominal`` / ``strict`` → (``-log10 alpha``, ``p < alpha`` label), when active."""
    out = {}
    for name in ("nominal", "strict"):
        a = _num((summary.get("levels") or {}).get(name, {}).get("alpha"))
        if a is not None and 0 < a < 1:
            out[name] = (-math.log10(a), f"p < {fmt_alpha(a)}")
    return out


# ------------------------------------------------------------------ data prep
def power_counts(tested: pd.DataFrame, runs: list[str], cutoffs: list[float]) -> dict[str, Any]:
    """Replicated pairs (same sign, ``-log10 p > t``: ``p < 10^-t`` as the levels) per run
    on a log1p-even grid."""
    lp = tested["log10p"].to_numpy(dtype=float)
    top = float(np.nanmax(lp)) if np.isfinite(lp).any() else 1.0
    top = max(top, 1.0, *cutoffs)
    grid = np.expm1(np.linspace(0.0, math.log1p(top), CURVE_POINTS))
    grid = np.unique(np.concatenate([grid, np.asarray(cutoffs, dtype=float)]))
    counts = {}
    for r in runs:
        sub = tested[tested["run"] == r]
        L = sub["log10p"].to_numpy(dtype=float)
        ok = (sub["sign_concordant"] == 1).to_numpy() & np.isfinite(L)
        Ls = np.sort(L[ok])
        counts[r] = (len(Ls) - np.searchsorted(Ls, grid, side="right")).astype(int).tolist()
    return {"grid": _r4(grid), "counts": counts}


def box_stats(values: np.ndarray) -> dict[str, Any] | None:
    """5th, 25th, 50th, 75th and 95th percentiles, the mean and the mean square."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if not v.size:
        return None
    q = np.percentile(v, [5, 25, 50, 75, 95])
    return {
        "q05": q[0],
        "q25": q[1],
        "med": q[2],
        "q75": q[3],
        "q95": q[4],
        "mean": float(v.mean()),
        "msq": float(np.mean(v * v)),
        "n": int(v.size),
    }


def neff_interval(
    zmat: np.ndarray, reps: int = BOOTSTRAP_REPS, seed: int = 1
) -> list[tuple[float, float]]:
    """95 % paired bootstrap interval of ``(mean z^2 - 1) / (baseline mean z^2 - 1)`` for
    every column of ``zmat`` (pairs × runs; column 0 is the baseline). Pairs are resampled
    jointly, so the correlation between runs on the same samples is kept."""
    n, k = zmat.shape
    if n < 2 or k < 1:
        return [(float("nan"), float("nan"))] * k
    sq = zmat * zmat
    rng = np.random.default_rng(seed)
    chunk = max(1, min(reps, 2_000_000 // max(n, 1)))
    ratios = []
    done = 0
    while done < reps:
        b = min(chunk, reps - done)
        w = rng.multinomial(n, np.full(n, 1.0 / n), size=b).astype(float)
        m = (w @ sq) / n  # b × k mean z^2
        base = m[:, :1] - 1.0
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios.append(np.where(base > 0, (m - 1.0) / base, np.nan))
        done += b
    r = np.vstack(ratios)
    out = []
    for j in range(k):
        col = r[:, j]
        col = col[np.isfinite(col)]
        out.append(
            (float(np.percentile(col, 2.5)), float(np.percentile(col, 97.5)))
            if col.size
            else (float("nan"), float("nan"))
        )
    return out


def wide_pairs(long: pd.DataFrame, runs: list[str]) -> pd.DataFrame:
    """One row per mapped (reference trait, variant) with the reference columns and, per
    run, ``beta``, ``z``, ``log10p``, ``sign`` and ``status`` (suffix ``@<run>``)."""
    key = ["ref_trait", "variant"]
    base = long.drop_duplicates(key)[
        ["trait", "ref_trait", "variant", "beta_ref", "z_ref", "log10p_ref"]
    ].reset_index(drop=True)
    for r in runs:
        sub = long[long["run"] == r].drop_duplicates(key)
        sub = sub[[*key, "beta", "z", "log10p", "sign_concordant", "status"]].rename(
            columns={
                "beta": f"beta@{r}",
                "z": f"z@{r}",
                "log10p": f"log10p@{r}",
                "sign_concordant": f"sign@{r}",
                "status": f"status@{r}",
            }
        )
        base = base.merge(sub, on=key, how="left")
    return base


def common_z(wide: pd.DataFrame, runs: list[str]) -> np.ndarray:
    """z of the pairs tested (finite z) in every run: pairs × runs."""
    if not runs:
        return np.empty((0, 0))
    z = wide[[f"z@{r}" for r in runs]].to_numpy(dtype=float)
    out: np.ndarray = z[np.isfinite(z).all(axis=1)]
    return out


def scatter_rows(wide: pd.DataFrame, runs: list[str]) -> tuple[pd.DataFrame, int]:
    """The pairs drawn in the agreement charts (at most :data:`MAX_SCATTER_POINTS`)."""
    if len(wide) <= MAX_SCATTER_POINTS:
        return wide, 0
    disc = np.zeros(len(wide), dtype=bool)
    for r in runs:
        disc |= (wide[f"sign@{r}"] == 0).to_numpy()
    order = wide.assign(_d=disc, _l=wide["log10p_ref"].fillna(-1.0)).sort_values(
        ["_d", "_l"], ascending=[False, False]
    )
    keep = order.head(MAX_SCATTER_POINTS).drop(columns=["_d", "_l"])
    return keep, len(wide) - len(keep)


# ---------------------------------------------------------------------- figures
def chart_layout(**overrides: Any) -> dict[str, Any]:
    """Base layout of the page's charts: larger tick and title fonts."""
    tick = {"tickfont": {"color": LIGHT["ink2"], "size": 14}}
    title_font = {"color": LIGHT["ink2"], "size": 15}
    lay = base_layout(
        margin={"l": 64, "r": 24, "t": 16, "b": 56},
        xaxis={**tick, "title": {"font": title_font}},
        yaxis={**tick, "title": {"font": title_font}},
        hoverlabel={"font": {"size": 14}},
        legend={"font": {"color": LIGHT["ink2"], "size": 14}},
    )
    for k, v in overrides.items():
        if isinstance(v, dict) and isinstance(lay.get(k), dict):
            lay[k] = {**lay[k], **v}
        else:
            lay[k] = v
    return lay


def _empty_fig(fig_id: str, height: int) -> str:
    return figure(
        {"data": [], "layout": chart_layout(height=height)}, fig_id, height=height
    ).replace('<figure class="viz">', '<figure class="viz" style="margin:0">')


def _seg(group: str, options: list[tuple[str, str]], label: str, current: str) -> str:
    """A segmented control (page script: ``data-seg``)."""
    return (
        f'<div class="seg seg-sm" role="group" aria-label="{esc(label)}" data-seg="{esc(group)}">'
        + "".join(
            f'<button type="button" data-v="{esc(v)}" aria-pressed="{"true" if v == current else "false"}">'
            f"{esc(t)}</button>"
            for v, t in options
        )
        + "</div>"
    )


def hits_html(strat: pd.DataFrame, summary: dict[str, Any], colors: dict[str, int]) -> str:
    """Rows = runs; columns = same sign, nominal and strict; in-row bars on a common
    scale (the run's tested pairs)."""
    cols: list[tuple[str, str, str]] = [
        ("sign", "Same sign", "Effect in the reference's direction")
    ]
    cut = level_cutoffs(summary)
    for name in ("nominal", "strict"):
        if name in cut and f"n_rep_{name}" in strat.columns:
            cols.append((name, cut[name][1], f"Same sign and {cut[name][1]}"))
    head = (
        '<div class="hits-row hits-head"><span></span>'
        + "".join(f'<span class="hh" title="{esc(tip)}">{esc(t)}</span>' for _, t, tip in cols)
        + "</div>"
    )
    rows = []
    for _, r in strat.iterrows():
        lab = str(r["label"])
        n_t = int(r["n_tested"]) if _num(r["n_tested"]) is not None else 0
        cells = []
        for k, t, _ in cols:
            n = int(r[f"n_rep_{k}"])
            w = 100.0 * n / n_t if n_t else 0.0
            cells.append(
                f'<span class="hc" title="{esc(lab)}: {n:,} of {n_t:,} tested pairs, {esc(t)}">'
                f'<span class="hv"><b>{n:,}</b>'
                f'<span class="hp">{fmt_pct(n, n_t)}</span></span>'
                f'<span class="track"><span class="fill" style="width:{w:.1f}%"></span></span></span>'
            )
        rows.append(
            f'<div class="hits-row" style="--c:{_series_var(colors[str(r["run"])])}">'
            f'<span class="hn"><span class="key" aria-hidden="true"></span>{esc(lab)}</span>'
            + "".join(cells)
            + "</div>"
        )
    return f'<div class="hits" style="--ncol:{len(cols)}">' + head + "".join(rows) + "</div>"


def pairs_table(wide: pd.DataFrame, runs: list[str]) -> dict[str, Any]:
    """The "All pairs" table (drawn by the page script): trait, variant, reference
    ``-log10 p`` and per run ``-log10 p``, the sign code (1 same, 0 opposite, ``None``
    not tested) and the reason a pair was not tested."""
    why = {"variant_absent": "variant not in genotypes", "not_tested": "no qpgen result"}
    rows = wide.head(MAX_TABLE_ROWS)
    return {
        "trait": rows["trait"].astype(str).tolist(),
        "variant": rows["variant"].astype(str).tolist(),
        "ref": _r4(rows["log10p_ref"]),
        "lp": {r: _r4(rows[f"log10p@{r}"]) for r in runs},
        "sign": {
            r: [None if not math.isfinite(x) else int(x) for x in rows[f"sign@{r}"].astype(float)]
            for r in runs
        },
        "why": {
            r: [
                "" if str(x) == "tested" else why.get(str(x), "not tested")
                for x in rows[f"status@{r}"]
            ]
            for r in runs
        },
        "n": len(wide),
    }


# ------------------------------------------------------------------ run details
def _details(title: str, body: str, id_: str) -> str:
    return (
        f'<details class="m-details" id="{esc(id_)}"><summary><h2>{esc(title)}</h2></summary>'
        f'<div class="rd-body">{body}</div></details>'
    )


def _cov_text(n_cov: Any, names: str) -> str:
    """``none``, ``66 (all)`` or ``6 (sex, PC1-5)``: one style for runs and comparisons."""
    n = fmt_int(n_cov)
    names = ", ".join(x.strip() for x in str(names or "").split(",") if x.strip())
    if names in ("", "NA", "none") or n == "0":
        return "none" if n in ("0", "NA") else n
    return f"{n} ({names})"


def run_details_html(
    sfx: str,
    d: dict[str, Any],
    runs_df: pd.DataFrame | None,
    params: dict[str, Any],
    long: pd.DataFrame,
) -> str:
    """Per-run settings and counts, the FDR and Bonferroni rules, local controls (rect
    engine), untested pairs and a link to every per-run number (``strategies.tsv``)."""
    summary, strat = d["summary"], d["strategies"]
    comparison = runs_df is not None
    parts: list[str] = []
    # ---- runs
    rows = []
    cov_text: dict[str, str] = dict(params.get("covariates") or {})
    rd = runs_df.set_index("label") if runs_df is not None and "label" in runs_df else None
    for i, (_, r) in enumerate(strat.iterrows()):
        lab = str(r["label"])
        rint = bool(int(r["rint"])) if "rint" in r and pd.notna(r["rint"]) else False
        names = cov_text.get(lab, "")
        samples = ""
        settings = ""
        if rd is not None and lab in rd.index:
            x = rd.loc[lab]
            names = str(x.get("cov_include", "") or "")
            samples = fmt_int(x.get("n_samples_analysable", ""))
            settings = f"../settings/{lab}.anqchor.yaml"
        elif not comparison:
            settings = "../anqchor.resolved.yaml"
        rows.append(
            {
                "Run": lab,
                "Covariates": _cov_text(r["n_cov"], names),
                "RINT": "yes" if rint else "no",
                "Samples": samples if samples != "NA" else "",
                "Tested pairs": fmt_int(r["n_tested"]),
                "FDR": fmt_int(r.get("n_rep_fdr")) if "n_rep_fdr" in r else "",
                "Bonferroni": fmt_int(r.get("n_rep_bonferroni")) if "n_rep_bonferroni" in r else "",
                "Mean z²": f"{float(r['mean_z2']):.1f}"
                if _num(r.get("mean_z2")) is not None
                else "",
                "Effect r vs reference": f"{float(r['beta_corr']):.3f}"
                if _num(r.get("beta_corr")) is not None
                else "",
                "Settings": settings,
                "_i": i,
            }
        )
    rt = pd.DataFrame(rows)
    show = [c for c in rt.columns if c != "_i" and (rt[c] != "").any()]
    text_cols = ("Run", "Covariates", "RINT", "Settings")
    head = "".join(
        f"<th>{esc(c)}</th>" if c in text_cols else f'<th class="num">{esc(c)}</th>' for c in show
    )
    body = ""
    for rec in rt.to_dict("records"):
        cells = []
        for c in show:
            v = str(rec[c])
            if c == "Run":
                cells.append(
                    f'<td class="rn" style="--c:{_series_var(int(rec["_i"]))}" title="{esc(v)}">'
                    f'<span class="key" aria-hidden="true"></span>{esc(v)}</td>'
                )
            elif c == "Settings" and v:
                cells.append(f'<td><a href="{esc(v)}" title="{esc(Path(v).name)}">yaml</a></td>')
            else:
                cells.append(
                    f"<td>{esc(v)}</td>" if c in text_cols else f'<td class="num">{esc(v)}</td>'
                )
        body += "<tr>" + "".join(cells) + "</tr>"
    parts.append(
        '<div class="rd-scroll"><table class="codes runs-t"><thead><tr>'
        + head
        + "</tr></thead><tbody>"
        + body
        + "</tbody></table></div>"
    )
    # ---- the two levels shown only here, and the shared settings
    lv = summary.get("levels") or {}
    items: list[tuple[str, object]] = []
    if comparison:
        items.append(("Baseline", str(summary.get("baseline", ""))))
    a = _num((lv.get("fdr") or {}).get("alpha"))
    if "fdr" in lv and a:
        items.append(("FDR", f"same sign and Benjamini-Hochberg q ≤ {fmt_alpha(a)}"))
    a = _num((lv.get("bonferroni") or {}).get("alpha"))
    if "bonferroni" in lv and a:
        items.append(("Bonferroni", f"same sign and p < {fmt_alpha(a)} / tested pairs"))
    items.append(
        (
            "Reference pairs",
            f"{fmt_int(summary.get('n_pairs_mapped'))} of {fmt_int(summary.get('n_ref_pairs'))} map onto the phenotypes",
        )
    )
    items.append(
        ("Test", f"qpgentools {'pair-assoc' if summary.get('engine') == 'pair' else 'rect-assoc'}")
    )
    if summary.get("recomputed_log10p"):
        items.append(("−log10 p recomputed", ", ".join(summary["recomputed_log10p"])))
    parts.append(definition_list(items))
    # ---- local controls (rect engine)
    if "n_control_cells" in strat.columns:
        names_c = {
            "label": "Run",
            "n_control_cells": "Cells",
            "lambda_gc_control": "λGC",
            "n_control_p_lt_1e5": "p < 1e-5",
            "n_control_cross_chrom": "Other-chromosome cells",
            "lambda_gc_control_cross_chrom": "λGC, other chromosomes",
            "n_control_cross_chrom_p_lt_1e5": "p < 1e-5, other chromosomes",
        }
        ccols = [c for c in names_c if c in strat.columns]
        parts.append("<h3>Local controls</h3>")
        parts.append(
            '<p class="rd-note">Other cells of the variant × trait rectangle: in LD with nearby '
            "QTLs, so not a null; the other-chromosome subset is closer to one.</p>"
        )
        parts.append(
            table(
                strat[ccols].rename(columns=names_c),
                f"tbl-controls{sfx}",
                searchable=False,
                digits=3,
            )
        )
    # ---- untested pairs
    un = long[long["status"] != "tested"]
    if len(un):
        g = (
            un.groupby(["trait", "variant", "status"], sort=False)["run"]
            .apply(lambda x: ", ".join(map(str, x)))
            .reset_index()
            .rename(
                columns={"trait": "Trait", "variant": "Variant", "status": "Reason", "run": "Runs"}
            )
        )
        g["Reason"] = (
            g["Reason"]
            .map({"variant_absent": "variant not in genotypes", "not_tested": "no qpgen result"})
            .fillna(g["Reason"])
        )
        parts.append(f"<h3>Pairs not tested ({len(g):,})</h3>")
        parts.append(table(g, f"cis_untested{sfx}", page_size=100))
    # ---- every per-run number: the summary table itself
    href = d.get("strategies_href")
    if href:
        parts.append(
            f'<p class="rd-note"><a href="{esc(href)}" download>All per-run numbers (TSV)</a></p>'
        )
    return _details("Run details", "".join(parts), f"details{sfx}")


# ---------------------------------------------------------------------- panel
def panel(
    ds: str,
    sfx: str,
    d: dict[str, Any],
    params: dict[str, Any],
    runs_df: pd.DataFrame | None,
) -> tuple[str, dict[str, Any]]:
    """One dataset's sections and its ``ANQ_DATA.panels`` entry."""
    summary, strat, long = d["summary"], d["strategies"], d["long"]
    runs = [str(r) for r in dict.fromkeys(strat["run"])]
    labels = {str(r): str(lab) for r, lab in zip(strat["run"], strat["label"], strict=True)}
    cidx = {r: i for i, r in enumerate(runs)}
    colors = run_colors(runs)
    comparison = len(runs) > 1
    tested = long[long["status"] == "tested"]
    cut = level_cutoffs(summary)
    baseline = runs[0] if runs else ""

    # ---- data
    curve = power_counts(tested, runs, [c[0] for c in cut.values()])
    wide = wide_pairs(long, runs)
    # Distributions over one set of pairs: those tested in every run (the pairs of the
    # effective sample size); the reference study over the same pairs.
    zw = wide[[f"z@{r}" for r in runs]].to_numpy(dtype=float)
    on = np.isfinite(zw).all(axis=1) if runs else np.zeros(len(wide), dtype=bool)
    wc = wide[on]
    dist: dict[str, Any] = {"z": [], "lp": []}
    for r in runs:
        dist["z"].append(box_stats(np.abs(wc[f"z@{r}"].to_numpy(dtype=float))))
        dist["lp"].append(box_stats(wc[f"log10p@{r}"].to_numpy(dtype=float)))
    ref_z = box_stats(np.abs(wc["z_ref"].to_numpy(dtype=float)))
    ref_lp = box_stats(wc["log10p_ref"].to_numpy(dtype=float))
    zc = common_z(wide, runs)
    neff: list[dict[str, Any]] = []
    if comparison:
        ci = neff_interval(zc)
        for i, r in enumerate(runs):
            row = strat[strat["run"] == r].iloc[0]
            neff.append(
                {
                    "v": _num(row.get("rel_neff")),
                    "lo": _num(ci[i][0]) if i else None,
                    "hi": _num(ci[i][1]) if i else None,
                }
            )
    pts, dropped = scatter_rows(wide, runs)
    vals = {
        r: {
            "beta": _r4(pts[f"beta@{r}"]),
            "z": _r4(pts[f"z@{r}"]),
            "lp": _r4(pts[f"log10p@{r}"]),
            "sign": [
                None if not math.isfinite(s) else int(s) for s in pts[f"sign@{r}"].astype(float)
            ],
        }
        for r in runs
    }
    # open on the first run after the baseline (the page explains no other choice)
    first = runs[1] if comparison else baseline
    entry = {
        "sfx": sfx,
        "comparison": comparison,
        "runs": [
            {
                "run": r,
                "label": labels[r],
                "color": colors[r],
                "dash": run_dash(i),
                "nTested": int(strat[strat["run"] == r]["n_tested"].iloc[0]),
            }
            for i, r in enumerate(runs)
        ],
        "baseline": baseline,
        "cutoffs": [{"t": c[0], "label": c[1]} for c in cut.values()],
        "curve": curve,
        "dist": dist,
        "ref": {"z": ref_z, "lp": ref_lp},
        "neff": neff,
        "nCommon": int(zc.shape[0]),
        "pairs": {
            "trait": pts["trait"].astype(str).tolist(),
            "variant": pts["variant"].astype(str).tolist(),
            "beta": _r4(pts["beta_ref"]),
            "z": _r4(pts["z_ref"]),
            "lp": _r4(pts["log10p_ref"]),
        },
        "vals": vals,
        "dropped": dropped,
        "selected": first,
        "table": pairs_table(wide, runs),
        "file": f"cis_pairs{sfx}.tsv",
        "colors": {
            "critical": STATUS["critical"],
            "muted": LIGHT["muted"],
            "axis": LIGHT["axis"],
            "ink": LIGHT["ink"],
            "ink2": LIGHT["ink2"],
            "surface": LIGHT["surface"],
        },
        "layout": chart_layout(),
    }

    # ---- sections
    n_runs = len(runs)
    bl = labels.get(baseline, baseline)
    n_on = int(on.sum())
    over = (
        f"over the {n_on:,} pairs tested in every run"
        if comparison
        else f"over the {n_on:,} tested pairs"
    )
    out: list[str] = []
    out.append(
        f'<section class="m-sec" id="power{sfx}"><div class="m-sec-head"><h2>Replicated pairs by p-value threshold</h2>'
        + _seg(f"power{sfx}", [("n", "Count"), ("f", "%")], "Show", "n")
        + "</div>"
        f'<div class="chart-card">{_empty_fig(f"fig-power{sfx}", 440)}</div>'
        '<p class="m-foot">Same sign as the reference and p below the threshold.</p></section>'
    )
    measure = _seg(f"dist{sfx}", [("z", "|z|"), ("lp", "−log10 p")], "Measure", "z")
    sig = (
        f'<section class="m-sec" id="strength{sfx}"><div class="m-sec-head"><h2>Signal strength</h2>'
        + ("" if comparison else measure)
        + "</div>"
    )
    if comparison:
        sig += (
            f'<h3 class="m-h3">Effective sample size relative to {esc(bl)}</h3>'
            f'<div class="chart-card">{_empty_fig(f"fig-neff{sfx}", 72 + 48 * n_runs)}</div>'
            f'<p class="m-foot">Ratio of (mean z² − 1) to {esc(bl)}\'s, {over}; '
            "strong pairs weigh most. Whiskers: 95% bootstrap interval.</p>"
            '<div class="m-sec-head m-h3-head"><h3 class="m-h3">Distribution</h3>'
            + measure
            + "</div>"
        )
    sig += (
        f'<div class="chart-card">{_empty_fig(f"fig-dist{sfx}", 120 + 52 * (n_runs + 1))}</div>'
        f'<p class="m-foot">{over[0].upper() + over[1:]}. Box: middle half; whiskers: 5th to 95th '
        "percentile; line: median. Hollow: reference study.</p>"
    )
    out.append(sig + "</section>")
    out.append(
        f'<section class="m-sec" id="hits{sfx}"><div class="m-sec-head"><h2>Significant hits</h2></div>'
        + hits_html(strat, summary, cidx)
        + '<p class="m-foot">Every column requires the same sign as the reference; '
        "% of tested pairs.</p></section>"
    )
    if comparison and n_runs > 5:
        sel = (
            f'<label class="run-pick">Run <select data-pick="conc{sfx}">'
            + "".join(
                f'<option value="{esc(r)}"{" selected" if r == first else ""}>{esc(labels[r])}</option>'
                for r in runs
            )
            + "</select></label>"
        )
    elif comparison:
        sel = _seg(f"conc{sfx}", [(r, labels[r]) for r in runs], "Run", first)
    else:
        sel = ""
    cards = [
        ("beta", "Effect size", "x: reference study · dashed: fit through 0"),
        ("z", "Z-score", "x: reference study · dashed: fit through 0"),
        ("lp", "−log10 p", "x: reference study."),
    ]
    if comparison:
        cards.append(
            (
                "zb",
                f"Z-score vs {bl}",
                f"x: {bl} · dashed: y = x",
            )
        )
    grid = "".join(
        f'<div class="chart-card conc-card" data-k="{k}"><div class="cc-head"><h3>{esc(t)}</h3>'
        f'<span class="cc-stat"></span></div>{_empty_fig(f"fig-conc-{k}{sfx}", 360)}'
        f'<p class="cc-foot">{esc(foot)}</p></div>'
        for k, t, foot in cards
    )
    note = (
        f'<p class="m-foot">{dropped:,} pairs not drawn (only the {MAX_SCATTER_POINTS:,} most significant in the reference '
        "and every pair with an opposite sign are).</p>"
        if dropped
        else ""
    )
    out.append(
        f'<section class="m-sec" id="agreement{sfx}"><div class="m-sec-head"><h2>Agreement with the reference</h2>'
        + sel
        + f'</div><div class="conc-grid{"" if comparison else " n3"}">{grid}</div>'
        '<p class="m-foot">y: the selected run. Open diamonds: opposite sign to the reference '
        f"(last plot: the two runs disagree in sign).</p>{note}</section>"
        if comparison
        else f'<section class="m-sec" id="agreement{sfx}"><div class="m-sec-head"><h2>Agreement with the reference</h2>'
        f'</div><div class="conc-grid n3">{grid}</div>'
        '<p class="m-foot">y: this run. Open diamonds: opposite sign to the reference.</p>'
        f"{note}</section>"
    )
    n_rows = len(wide)
    more = (
        f'<p class="rd-note">The first {MAX_TABLE_ROWS:,} of {n_rows:,} pairs; all are in '
        "pairs.long.tsv.gz.</p>"
        if n_rows > MAX_TABLE_ROWS
        else ""
    )
    out.append(
        f'<details class="m-details" id="pairs{sfx}" open><summary><h2>All pairs <span class="cnt">({n_rows:,})</span>'
        '</h2></summary><div class="m-details-body">'
        + (
            '<p class="rd-note">−log10 p per run; red with −: opposite sign to the reference; '
            "– : not tested.</p>"
            if comparison
            else '<p class="rd-note">Red with −: opposite sign to the reference; – : not tested.</p>'
        )
        + more
        + f'<div class="tbar" data-ptbl="{esc(sfx)}"></div>'
        f'<div class="ptbl-wrap" id="ptbl{sfx}"><table class="mtable ptbl"><colgroup></colgroup>'
        "<thead></thead><tbody></tbody></table></div>"
        f'<div class="tfoot" data-ptfoot="{esc(sfx)}"></div>'
        "</div></details>"
    )
    out.append(run_details_html(sfx, d, runs_df, params, long))
    return "".join(out), entry


def context_html(
    d: dict[str, Any], params: dict[str, Any], runs_df: pd.DataFrame | None, ds_label: str = ""
) -> str:
    """The context line: pairs tested, samples, reference file; run and baseline pills."""
    summary, strat, long = d["summary"], d["strategies"], d["long"]
    tested = long[long["status"] == "tested"]
    n_map = int(summary.get("n_pairs_mapped") or 0)
    nt = [int(x) for x in strat["n_tested"]] if len(strat) else [0]
    ref = params.get("reference_file") or ""
    if runs_df is not None and "reference_pairs" in runs_df and len(runs_df):
        refs = {Path(str(x)).name for x in runs_df["reference_pairs"] if str(x) and str(x) != "NA"}
        if len(refs) == 1:
            ref = refs.pop()
    src = f" from {ref}" if ref else " from the reference"
    if min(nt) == max(nt):
        pairs = (
            f"{nt[0]:,} pairs{src} tested"
            if nt[0] == n_map
            else f"{nt[0]:,} of {n_map:,} pairs{src} tested"
        )
    else:
        pairs = f"{min(nt):,}–{max(nt):,} of {n_map:,} pairs{src} tested"
    n = (
        pd.to_numeric(tested.get("n"), errors="coerce")
        if "n" in tested
        else pd.Series([], dtype=float)
    )
    n = n[np.isfinite(n)] if len(n) else n
    samples = f"{int(np.median(n)):,} samples" if len(n) else ""
    bits = [x for x in (pairs, samples) if x]
    labels = strat["label"].astype(str).tolist()
    pills = []
    if ds_label:
        pills.append(f'<span class="pill">{esc(ds_label)}</span>')
    if len(labels) == 1:
        rint = (
            bool(int(strat["rint"].iloc[0]))
            if "rint" in strat and pd.notna(strat["rint"].iloc[0])
            else False
        )
        pills.append(f'<span class="pill">Run {esc(labels[0])}{" · RINT" if rint else ""}</span>')
    elif labels:
        pills.append(f'<span class="pill">Baseline {esc(labels[0])}</span>')
    alerts = "".join(
        f'<p class="m-alert" role="note">⚠ {esc(w)}</p>' for w in summary.get("warnings") or []
    )
    return (
        f'<div class="ctx"><span class="num">{esc(" · ".join(bits))}</span>{"".join(pills)}</div>'
        + alerts
    )


def build_page(
    datasets: list[tuple[str, dict[str, Any]]],
    params: dict[str, Any],
    runs_df: pd.DataFrame | None = None,
) -> Page:
    """The cis report: one panel per dataset (normally one)."""
    multi = len(datasets) > 1
    bodies: list[str] = []
    entries: list[dict[str, Any]] = []
    header: list[str] = []
    n_runs = 0
    for ds, d in datasets:
        sfx = f"-{ds}" if multi else ""
        html, entry = panel(ds, sfx, d, params, runs_df)
        n_runs = max(n_runs, len(entry["runs"]))
        if multi:
            header.append(context_html(d, params, runs_df, ds_label=ds))
            bodies.append(
                f'<div class="ds-block" id="ds{sfx}"><h2 class="ds-title">{esc(ds)}</h2>{html}</div>'
            )
        else:
            header.append(context_html(d, params, runs_df))
            bodies.append(html)
        entries.append(entry)
    title = f"Comparing {n_runs} runs: {TITLE}" if n_runs > 1 else TITLE
    return Page(
        title=title,
        sections=[Section(id="cis", title="", body="".join(bodies))],
        nav=index_nav(params),
        extra_js=CIS_PAGE_JS,
        data={"panels": entries},
        header_html="".join(header),
        hide_brand=True,
        body_class="page-match page-cis",
        toc=False,
    )


@register("report_cis")
def report_cis(ctx: StepContext) -> None:
    """Render ``report/cis.html`` from every dataset's cis summaries."""
    p = ctx.params
    summary_dirs: dict[str, str] = dict(p.get("summary_dirs") or {})
    datasets = [(ds, load_dataset(ctx.run_dir, ds, summary_dirs.get(ds))) for ds in p["datasets"]]
    runs_df = None
    if p.get("comparison"):
        f = ctx.run_dir / summary_dirs.get(datasets[0][0], ".") / "runs.tsv"
        try:
            runs_df = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        except (OSError, ValueError) as e:
            raise StepError(f"cannot read {f}: {e}") from e
    page = build_page(datasets, p, runs_df)
    write_page(ctx.path(Layout.report("cis")), page, p.get("plotlyjs", "embed"))
