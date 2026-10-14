"""The identity-matching report (``report/match.html``, step kind ``report_match``).

The page answers one question: does each omics sample come from the individual whose ID
it carries? From top to bottom (DESIGN.md section 6.3):

* header: title, a context line (samples, genotypes, traits) and the scoring mode;
* dashboard: large tiles with the counts (examined, with an assigned ID, strict match,
  lenient match, matches a different sample, matches two genotype IDs and fits several
  when present, matches nothing, and no assigned ID when present);
* two sliders (lenient cutoff, clear lead) that reclassify every sample in the browser;
* one card per sample to check (one summary card for many samples without an ID), a
  scatter chart of every scored sample's best against its runner-up match score, the
  duplicates (samples that match two genotype IDs about equally), the table of all samples
  and the run details.

"Match score" is qpgen's Z-score of a genotype's polygenic scores against the sample's
measurements; the page never says "z".

Statuses are shown with plain labels (:data:`CATEGORIES`); qpgen's status codes appear
only in tooltips, the TSV download and the code table in the run details.

Interactive thresholds (decision 9): the page embeds every sample's Z-scores and
recomputes statuses in the browser with :data:`MATCH_RULES_JS`, a line-by-line port of
:func:`anqchor.steps.match.classify` and its helpers. The TSV outputs keep the run's
thresholds.
"""

from __future__ import annotations

import math
import re
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from anqchor.pipelines.layout import Layout
from anqchor.report.render import Page, Section, code_block, esc, figure, index_nav, write_page
from anqchor.report.theme import (
    CATEGORICAL_LIGHT,
    LIGHT,
    STATUS,
    WARNING_LINE,
    base_layout,
    series_color,
)
from anqchor.steps import StepContext, StepError, register
from anqchor.steps.match import DUPLICATE_RATIO, NOT_SCORED, RELATED_GAP, is_duplicate
from anqchor.util.files import read_json

TITLE = "Sample identity matching with anQChor"
"""Page title and ``<h1>`` (display name anQChor; the package and CLI stay ``anqchor``)."""

ROW_FIELDS = [
    "pheno_id",
    "geno_id_assigned",
    "qpgen_status",
    "n_traits",
    "rank_self",
    "z_self",
    "cor_self",
    *[f"id_{k}" for k in range(1, 6)],
    *[f"z_{k}" for k in range(1, 6)],
]
"""Field order of the per-sample rows embedded in the page (``ANQ_DATA.fields``)."""

CATEGORIES: list[dict[str, str]] = [
    {"key": "different", "label": "Matches another sample", "icon": "✕", "role": "critical"},
    {"key": "taken", "label": "ID already used", "icon": "!", "role": "critical"},
    {"key": "nomatch", "label": "Own ID too weak", "icon": "?", "role": "serious"},
    {"key": "several", "label": "Matches several", "icon": "…", "role": "serious"},
    {"key": "twoids", "label": "Matches two genotype IDs", "icon": "2", "role": "dup"},
    {"key": "newid", "label": "New ID found", "icon": "+", "role": "info"},
    {"key": "lenient", "label": "Lenient match", "icon": "≈", "role": "warning"},
    {"key": "strict", "label": "Strict match", "icon": "✓", "role": "good"},
    {"key": "noid", "label": "No clear match", "icon": "–", "role": "muted"},
    {"key": "nodata", "label": "Not scored", "icon": "–", "role": "muted"},
]
"""Display categories in severity order: the label, icon and color role of each result.
The same labels appear on tiles, cards, chips, the chart legend and the table.
``newid``, ``taken`` and ``noid`` occur only for samples without an assigned ID."""

CATEGORY_CODES = [
    ("Strict match", "SELF_BEST"),
    ("Lenient match", "SELF_LENIENT"),
    ("Matches another sample / New ID found / ID already used", "SINGLE_NEW_BEST"),
    ("Matches several", "MULTI_NEW_BEST"),
    ("Matches two genotype IDs", "SINGLE_NEW_BEST or MULTI_NEW_BEST"),
    ("Own ID too weak", "UNCLEAR"),
    ("No clear match", "NO_SELF"),
    ("Not scored", f"NO_OBS_TRAITS, {NOT_SCORED}"),
]
"""Result names and the qpgen codes in ``samples.tsv`` (run details)."""

ROLE_COLORS = {
    **STATUS,
    "info": CATEGORICAL_LIGHT[0],
    "dup": CATEGORICAL_LIGHT[6],
    "muted": LIGHT["muted"],
}
"""Color of each category role; ``dup`` (violet, series 7) is not a status color: two
genotype IDs of one person are a data issue, not a failed sample."""


def category_of(status: str, assigned: bool, two_ids: bool = False, taken: bool = False) -> str:
    """Display category (:data:`CATEGORIES` key) of a qpgen status. Mirrors
    ``ANQ_MATCH.category``.

    ``two_ids``: the sample matches two genotype IDs about equally
    (:func:`anqchor.steps.match.is_duplicate`) and, when it has an assigned ID, that ID is
    one of the two; it then shows as ``twoids`` unless its own ID ranks first.
    ``taken``: a sample without an assigned ID whose clearly best genotype is the confirmed
    match of another sample (no ID is suggested)."""
    if status == "SELF_BEST":
        return "strict"
    if two_ids and status in ("SINGLE_NEW_BEST", "MULTI_NEW_BEST"):
        return "twoids"
    if status == "SELF_LENIENT":
        return "lenient"
    if status == "SINGLE_NEW_BEST":
        if assigned:
            return "different"
        return "taken" if taken else "newid"
    if status == "MULTI_NEW_BEST":
        return "several"
    if status == "UNCLEAR":
        return "nomatch"
    if status == "NO_SELF":
        return "noid"
    return "nodata"


TILE_KEYS = (
    "examined",
    "assigned",
    "strict",
    "lenient_only",
    "different",
    "twoids",
    "several",
    "none",
    "nodata",
    "unassigned",
    "newid",
    "taken",
    "un_twoids",
    "un_nodata",
    "unclear",
    "dups",
)


def tile_counts(
    categories: list[str], assigned: list[bool], dups: list[bool] | None = None
) -> dict[str, int]:
    """Dashboard counts. Mirrors ``ANQ_MATCH.tileCounts``.

    Every sample with an assigned ID is in exactly one of ``lenient`` (strict and lenient
    matches), ``different``, ``twoids``, ``several``, ``none`` (own ID too weak) and
    ``nodata`` (not scored): ``assigned == lenient + different + twoids + several + none +
    nodata``. qpgen's MULTI_NEW_BEST overrides SELF_LENIENT, so ``several`` is its own
    tile. Samples without an assigned ID count in ``unassigned`` and one of ``newid``,
    ``taken``, ``un_twoids``, ``un_nodata`` and ``unclear`` (no clear match). ``dups``:
    samples that match two genotype IDs about equally (the flags ``dups``, assigned or
    not, including strict matches whose own ID ties with another), the number on the
    "Matches two genotype IDs" tile.
    """
    t = dict.fromkeys(TILE_KEYS, 0)
    t["examined"] = len(categories)
    for c, a in zip(categories, assigned, strict=True):
        if a:
            t["assigned"] += 1
            key = {"strict": "strict", "lenient": "lenient_only", "nomatch": "none"}.get(c, c)
            if key not in ("different", "twoids", "several", "nodata", "strict", "lenient_only"):
                key = "none"
            t[key] += 1
        else:
            t["unassigned"] += 1
            key = {"twoids": "un_twoids", "nodata": "un_nodata"}.get(c, c)
            t[key if key in ("newid", "taken", "un_twoids", "un_nodata") else "unclear"] += 1
    t["lenient"] = t["strict"] + t["lenient_only"]
    t["dups"] = sum(bool(d) for d in dups) if dups is not None else t["twoids"] + t["un_twoids"]
    return t


MATCH_RULES_JS = r"""
/* Port of anqchor.steps.match: classify, group_of, is_rescue, close_candidates,
   is_related, is_duplicate, proposed genotype and verdict; and of anqchor.report.match:
   category_of, tile_counts, display_flags. Keep line-by-line in sync. */
var ANQ_MATCH = (function () {
  "use strict";
  var MISSING_Z = -9999;
  var RELATED_GAP = __RELATED_GAP__;
  var DUPLICATE_RATIO = __DUPLICATE_RATIO__;
  var STATUS_GROUP = { SELF_BEST: "PASS", SELF_LENIENT: "WARN", UNCLEAR: "WARN",
    SINGLE_NEW_BEST: "FAIL", MULTI_NEW_BEST: "FAIL", NO_SELF: "NO_TEST", NO_OBS_TRAITS: "NO_TEST" };
  function isMissing(v) { return v === null || v === undefined || typeof v !== "number" || isNaN(v); }
  function zval(v) { return isMissing(v) ? MISSING_Z : v; }
  function classify(mapped, nTraits, rankSelf, zSelf, zTop, zThreshold, zDiff) {
    if (isMissing(nTraits) || nTraits === 0) return "NO_OBS_TRAITS";
    var z = [];
    for (var k = 0; k < 5; k++) z.push(zval(zTop[k]));
    if (mapped && !isMissing(rankSelf) && Math.trunc(rankSelf) === 1) return "SELF_BEST";
    if (z[0] > z[1] + zDiff) return "SINGLE_NEW_BEST";
    for (var j = 2; j < 5; j++) {
      if (z[j - 1] > z[j] + zDiff) return "MULTI_NEW_BEST";
    }
    if (mapped && zval(zSelf) > zThreshold) return "SELF_LENIENT";
    return mapped ? "UNCLEAR" : "NO_SELF";
  }
  function groupOf(status, mapped) {
    if (!mapped || status === "NOT_SCORED") return "NO_TEST";
    return STATUS_GROUP[status] || "NO_TEST";
  }
  function isRescue(status, mapped) { return !mapped && status === "SINGLE_NEW_BEST"; }
  function closeCandidates(ids, zs, zThreshold, gap) {
    if (gap === undefined) gap = RELATED_GAP;
    var pairs = [];
    for (var i = 0; i < ids.length; i++) {
      if (ids[i] && zval(zs[i]) > MISSING_Z) pairs.push([ids[i], zval(zs[i])]);
    }
    var bestK = 0, bestDrop = -Infinity;
    for (var k = 2; k < Math.min(pairs.length, 5); k++) {
      var drop = pairs[k - 1][1] - pairs[k][1];
      if (drop > bestDrop) { bestK = k; bestDrop = drop; }
    }
    if (bestK === 0 || bestDrop < gap || pairs[bestK - 1][1] < zThreshold) return [];
    return pairs.slice(0, bestK);
  }
  function isRelated(nClose) { return nClose >= 2; }
  /* steps.match.is_duplicate: matches two genotype IDs about equally. */
  function isDuplicate(zTop, zThreshold, zDiff, ratio) {
    if (ratio === undefined) ratio = DUPLICATE_RATIO;
    var z1 = zval(zTop[0]), z2 = zval(zTop[1]), z3 = zval(zTop[2]);
    if (z1 <= 0 || z2 <= MISSING_Z) return false;
    return z2 > z3 + zDiff && z2 > zThreshold && z2 >= ratio * z1;
  }
  /* Display category of a status (anqchor.report.match.category_of). */
  function category(status, mapped, twoIds, taken) {
    if (status === "SELF_BEST") return "strict";
    if (twoIds && (status === "SINGLE_NEW_BEST" || status === "MULTI_NEW_BEST")) return "twoids";
    if (status === "SELF_LENIENT") return "lenient";
    if (status === "SINGLE_NEW_BEST") return mapped ? "different" : (taken ? "taken" : "newid");
    if (status === "MULTI_NEW_BEST") return "several";
    if (status === "UNCLEAR") return "nomatch";
    if (status === "NO_SELF") return "noid";
    return "nodata";
  }
  /* Dashboard counts (anqchor.report.match.tile_counts);
     assigned = lenient + different + twoids + several + none + nodata. */
  function tileCounts(cats, mapped, dups) {
    var t = { examined: cats.length, assigned: 0, strict: 0, lenient_only: 0, different: 0, twoids: 0,
              several: 0, none: 0, nodata: 0, unassigned: 0, newid: 0, taken: 0, un_twoids: 0, un_nodata: 0,
              unclear: 0, dups: 0 };
    for (var i = 0; i < cats.length; i++) {
      var c = cats[i];
      if (mapped[i]) {
        t.assigned++;
        if (c === "strict") t.strict++;
        else if (c === "lenient") t.lenient_only++;
        else if (c === "different" || c === "twoids" || c === "several" || c === "nodata") t[c]++;
        else t.none++;
      } else {
        t.unassigned++;
        if (c === "newid" || c === "taken") t[c]++;
        else if (c === "twoids") t.un_twoids++;
        else if (c === "nodata") t.un_nodata++;
        else t.unclear++;
      }
    }
    t.lenient = t.strict + t.lenient_only;
    if (dups) { for (var j = 0; j < dups.length; j++) if (dups[j]) t.dups++; }
    else t.dups = t.twoids + t.un_twoids;
    return t;
  }
  /* One sample. r is an object with the ROW_FIELDS keys. */
  function evaluate(r, zThreshold, zDiff) {
    var mapped = !!r.geno_id_assigned;
    var ids = [r.id_1, r.id_2, r.id_3, r.id_4, r.id_5];
    var zs = [r.z_1, r.z_2, r.z_3, r.z_4, r.z_5];
    var status = r.qpgen_status === "NOT_SCORED" ? "NOT_SCORED"
      : classify(mapped, r.n_traits, r.rank_self, r.z_self, zs, zThreshold, zDiff);
    var group = groupOf(status, mapped);
    var gap = (isMissing(r.z_1) || isMissing(r.z_2)) ? null : r.z_1 - r.z_2;
    var rescue = isRescue(status, mapped);
    var proposed = ((group === "FAIL" || rescue) && gap !== null && gap >= zDiff && r.id_1) ? r.id_1 : "";
    var close = closeCandidates(ids, zs, zThreshold);
    var dup = status !== "NOT_SCORED" && isDuplicate(zs, zThreshold, zDiff);
    /* steps.match.two_genotype_ids: one person under two genotype IDs, no relabel */
    var twoIds = dup && (!mapped || r.geno_id_assigned === r.id_1 || r.geno_id_assigned === r.id_2);
    if (twoIds) proposed = "";
    return { status: status, group: group, mapped: mapped, gap: gap, rescue: rescue, twoIds: twoIds,
             proposed: proposed, close: close, related: !dup && isRelated(close.length), dup: dup };
  }
  /* Every row; then no relabel where the best genotype is the confirmed (PASS) match of
     another sample (steps.match.proposal_conflict); then the display category
     (report.match.display_flags and category_of). */
  function evaluateAll(rows, zThreshold, zDiff) {
    var evals = rows.map(function (r) { return evaluate(r, zThreshold, zDiff); });
    var byGeno = {};
    rows.forEach(function (r, i) {
      if (r.geno_id_assigned) (byGeno[r.geno_id_assigned] = byGeno[r.geno_id_assigned] || []).push(i);
    });
    evals.forEach(function (e, i) {
      if (e.proposed && (byGeno[e.proposed] || []).some(function (j) { return j !== i && evals[j].group === "PASS"; })) {
        e.proposed = "";
        e.conflict = true;
      }
    });
    evals.forEach(function (e, i) {
      e.taken = !e.mapped && e.status === "SINGLE_NEW_BEST" && !!e.conflict;
      e.cat = category(e.status, e.mapped, e.twoIds, e.taken);
    });
    return evals;
  }
  function failPhrase(phenoId, assigned, best) {
    var label = (assigned && assigned !== phenoId) ? " (labelled " + assigned + ")" : "";
    return phenoId + label + " best matches " + best;
  }
  /* summary.json's verdict sentence (steps.match.make_verdict); not shown on the page. */
  function verdict(rows, evals, maxListed) {
    maxListed = maxListed || 3;
    var nMapped = 0, selfBest = 0, lenient = 0, unclear = 0, noObs = 0, notScored = 0, rescue = 0, fails = [],
        dups = [];
    for (var i = 0; i < rows.length; i++) {
      var e = evals[i];
      if (e.status === "NO_OBS_TRAITS") noObs++;
      if (e.status === "NOT_SCORED") notScored++;
      if (e.dup) dups.push(rows[i].pheno_id + ": " + rows[i].id_1 + ", " + rows[i].id_2);
      if (e.mapped) {
        nMapped++;
        if (e.status === "SELF_BEST") selfBest++;
        if (e.status === "SELF_LENIENT") lenient++;
        if (e.status === "UNCLEAR") unclear++;
        if (e.group === "FAIL" && !e.twoIds) fails.push(failPhrase(rows[i].pheno_id, rows[i].geno_id_assigned, rows[i].id_1));
      } else if (e.status === "SINGLE_NEW_BEST" && e.proposed !== "") {
        rescue++;
      }
    }
    var parts = [];
    if (nMapped) {
      var lead = (selfBest + lenient) + " of " + nMapped + " mapped samples confirmed";
      if (lenient) lead += " (" + lenient + " only leniently)";
      parts.push(lead);
    } else {
      parts.push("no phenotyped sample is mapped to a genotype");
    }
    if (fails.length) {
      var more = fails.length > maxListed ? "; and " + (fails.length - maxListed) + " more" : "";
      parts.push(fails.length + " mismatch" + (fails.length !== 1 ? "es" : "") + " (" +
                 fails.slice(0, maxListed).join("; ") + more + ")");
    }
    if (unclear) parts.push(unclear + " unclear");
    if (noObs) parts.push(noObs + " without observed traits");
    if (notScored) parts.push(notScored + " not scored (no covariates)");
    if (dups.length) {
      var moreD = dups.length > maxListed ? "; and " + (dups.length - maxListed) + " more" : "";
      parts.push(dups.length + " sample" + (dups.length !== 1 ? "s" : "") + " matching two genotype IDs about equally (" +
                 dups.slice(0, maxListed).join("; ") + moreD + ")");
    }
    var nUnmapped = rows.length - nMapped;
    if (nUnmapped) {
      parts.push(nUnmapped + " unmapped sample" + (nUnmapped !== 1 ? "s" : "") + ", " + rescue +
                 " rescue candidate" + (rescue !== 1 ? "s" : ""));
    }
    var status = fails.length ? "critical" : ((unclear || !nMapped || dups.length) ? "warning" : "good");
    return { text: parts.join("; "), status: status };
  }
  return { MISSING_Z: MISSING_Z, RELATED_GAP: RELATED_GAP, DUPLICATE_RATIO: DUPLICATE_RATIO, classify: classify,
           groupOf: groupOf, isRescue: isRescue, closeCandidates: closeCandidates, isRelated: isRelated,
           isDuplicate: isDuplicate, category: category, tileCounts: tileCounts,
           evaluate: evaluate, evaluateAll: evaluateAll, verdict: verdict };
})();
if (typeof module !== "undefined" && module.exports) { module.exports = ANQ_MATCH; }
""".replace("__RELATED_GAP__", repr(RELATED_GAP)).replace(
    "__DUPLICATE_RATIO__", repr(DUPLICATE_RATIO)
)

MATCH_PAGE_JS = (
    resources.files("anqchor.report")
    .joinpath("static")
    .joinpath("match_page.js")
    .read_text(encoding="utf-8")
)
"""The page script (``static/match_page.js``): dashboard, sliders, cards, chart,
duplicates and the sample table, recomputed with :data:`MATCH_RULES_JS`."""


# --------------------------------------------------------------------------- helpers
def slug(name: str) -> str:
    """An HTML-id-safe version of a dataset name."""
    return re.sub(r"[^A-Za-z0-9_-]", "-", name)


def _fmt(v: Any, digits: int = 2) -> str:
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "NA"
    if isinstance(v, float):
        return f"{v:.{digits}f}"
    return str(v)


def _num(v: Any) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _str(v: Any) -> str:
    if v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v)):
        return ""
    return str(v)


def _yes(v: Any) -> str:
    return "yes" if v else "no"


def read_tsv(path: Path) -> pd.DataFrame:
    """A summary TSV with string ID columns and ``NA`` as missing."""
    return pd.read_csv(
        path,
        sep="\t",
        na_values=["NA"],
        keep_default_na=False,
        dtype={
            c: str
            for c in (
                "pheno_id",
                "geno_id_assigned",
                "best_id",
                "second_id",
                "proposed_geno_id",
                "note",
                "qpgen_status",
                "group",
                *[f"id_{k}" for k in range(1, 6)],
                "partner_pheno_id",
                "partner_geno_id_assigned",
                "partner_best_id",
                "trait",
                "ref_trait",
                "candidate_geno_id",
                "candidate_assigned_to",
            )
        },
    )


def embed_rows(samples: pd.DataFrame) -> list[list[Any]]:
    """Per-sample rows in :data:`ROW_FIELDS` order for the page script."""
    out = []
    for r in samples.to_dict("records"):
        row: list[Any] = []
        for f in ROW_FIELDS:
            v = r.get(f)
            if f in ("n_traits", "rank_self"):
                n = _num(v)
                row.append(int(n) if n is not None else None)
            elif f.startswith("z_") or f == "cor_self":
                n = _num(v)
                row.append(round(n, 6) if n is not None else None)
            else:
                row.append(_str(v))
        out.append(row)
    return out


def display_flags(r: dict[str, Any], z_threshold: float, z_diff: float) -> tuple[bool, bool, bool]:
    """``(dup, two_ids, taken)`` of one ``samples.tsv`` row for :func:`category_of`.
    Mirrors ``ANQ_MATCH.evaluateAll``. ``taken``: no assigned ID, one clearly best
    genotype, and no proposed ID (the genotype is another sample's confirmed match)."""
    status = _str(r.get("qpgen_status"))
    assigned = _str(r.get("geno_id_assigned"))
    dup = status != NOT_SCORED and is_duplicate(
        [_num(r.get(f"z_{k}")) for k in (1, 2, 3)], z_threshold, z_diff
    )
    two_ids = dup and (not assigned or assigned in (_str(r.get("id_1")), _str(r.get("id_2"))))
    taken = not assigned and status == "SINGLE_NEW_BEST" and not _str(r.get("proposed_geno_id"))
    return dup, two_ids, taken


def run_categories(
    samples: pd.DataFrame, z_threshold: float = 1.96, z_diff: float = 2.0
) -> tuple[list[str], list[bool]]:
    """Display categories of the run's statuses (``samples.tsv`` ``qpgen_status``) at the
    run's thresholds."""
    cats: list[str] = []
    assigned: list[bool] = []
    for r in samples.to_dict("records"):
        a = bool(_str(r.get("geno_id_assigned")))
        _, two_ids, taken = display_flags(r, z_threshold, z_diff)
        cats.append(category_of(_str(r.get("qpgen_status")), a, two_ids, taken))
        assigned.append(a)
    return cats, assigned


def qpgen_version(summary: dict[str, Any]) -> str:
    """qpgen version from the log excerpt (``qpgentools 0.2.0 (git 5b811e2)``), or ``unknown``."""
    for line in summary.get("qpgen_log_excerpt") or []:
        m = re.search(r"qpgentools\s+(\S+(?:\s+\(git [0-9a-f]+\))?)", str(line))
        if m:
            return m.group(1)
    return "unknown"


def scoring_text(summary: dict[str, Any]) -> str:
    """The header chip: ``Scoring: independent traits`` or ``Scoring: Mahalanobis``."""
    return "Scoring: Mahalanobis" if summary.get("mahalanobis") else "Scoring: independent traits"


def scoring_detail(summary: dict[str, Any]) -> str:
    """Run details: ``Independent traits`` or ``Mahalanobis (λ 0.58, auto)``."""
    if not summary.get("mahalanobis"):
        return "Independent traits"
    lam = _num(summary.get("lambda_used"))
    mode = {"auto": "auto", "fixed": "fixed", "fallback": "fallback"}.get(
        str(summary.get("lambda_mode")), ""
    )
    inner = ", ".join(x for x in (f"λ {lam:.2f}" if lam is not None else "", mode) if x)
    return "Mahalanobis" + (f" ({inner})" if inner else "")


def weight_figure(weights: pd.DataFrame) -> dict[str, Any]:
    """Histogram of trait weights (computed bins, one bar trace)."""
    w = (
        pd.to_numeric(weights["weight"], errors="coerce").dropna().to_numpy(dtype=float)
        if "weight" in weights.columns
        else np.array([], dtype=float)
    )
    layout = base_layout(
        height=300,
        margin={"t": 16},
        xaxis={"title": {"text": "Trait weight"}},
        yaxis={"title": {"text": "Traits"}},
        showlegend=False,
    )
    if not len(w):
        return {"data": [], "layout": layout}
    lo, hi = float(w.min()), float(w.max())
    if hi <= lo:
        hi = lo + 1e-6
    edges = np.linspace(lo, hi, 41)
    counts, _ = np.histogram(w, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2
    return {
        "data": [
            {
                "type": "bar",
                "x": centers.tolist(),
                "y": counts.tolist(),
                "width": float(edges[1] - edges[0]) * 0.8,
                "marker": {"color": series_color(0)},
                "customdata": np.stack([edges[:-1], edges[1:]], axis=1).tolist(),
                "hovertemplate": "weight %{customdata[0]:.3f} to %{customdata[1]:.3f}: %{y} traits<extra></extra>",
            }
        ],
        "layout": layout,
    }


def chart_layout() -> dict[str, Any]:
    """Base layout of the match chart; the page script sets ranges, lines and labels."""
    tick = {"tickfont": {"color": LIGHT["ink2"], "size": 14}}
    title_font = {"color": LIGHT["ink2"], "size": 15}
    return base_layout(
        height=460,
        margin={"l": 60, "r": 24, "t": 12, "b": 56},
        showlegend=False,
        xaxis={**tick, "title": {"text": "Best match score", "font": title_font}},
        yaxis={**tick, "title": {"text": "Runner-up score", "font": title_font}},
        hoverlabel={"font": {"size": 14}},
    )


# ------------------------------------------------------------------------------ page
_TILES = [
    # key, label, icon, color role, tooltip; the result tiles use the CATEGORIES labels
    ("examined", "Samples examined", "", "", "Omics samples in this report."),
    ("assigned", "With an assigned ID", "", "", "Samples labelled with a genotype ID."),
    (
        "strict",
        "Strict match to own ID",
        "✓",
        "good",
        "The own ID's genotype fits the sample better than any other.",
    ),
    (
        "lenient",
        "Lenient match to own ID",
        "≈",
        "warning",
        "Strict matches plus samples whose own ID is not the top match but scores above the "
        "lenient cutoff.",
    ),
    (
        "different",
        "Matches another sample",
        "✕",
        "critical",
        "Another sample's genotype fits clearly better than the own ID.",
    ),
    (
        "twoids",
        "Matches two genotype IDs",
        "2",
        "dup",
        "Two genotypes fit about equally well and clearly better than the rest: probably one "
        "person under two genotype IDs.",
    ),
    (
        "several",
        "Matches several",
        "…",
        "serious",
        "A few genotypes fit about equally well, clearly better than the rest.",
    ),
    (
        "none",
        "Own ID too weak",
        "?",
        "serious",
        "The own ID is not the top match and scores at or below the lenient cutoff, and no "
        "other genotype leads clearly.",
    ),
    ("nodata", "Not scored", "–", "muted", "Samples qpgen did not score."),
    ("unassigned", "No assigned ID", "+", "info", "Samples without a genotype ID."),
]
TILE_OPTIONAL = ("twoids", "several", "nodata")
"""Tiles shown only when their count (at the run's or the sliders' values) is above 0."""
TILE_VALUE = {"twoids": "dups"}
"""Tiles whose number is another :func:`tile_counts` key."""


def _pct(n: int, d: int) -> str:
    return f"{100.0 * n / d:.1f}%" if d else ""


def _short_list(ids: list[str]) -> str:
    return ", ".join(ids) if len(ids) <= 3 else f"{ids[0]}, {ids[1]} and {len(ids) - 2} more"


def tile_subs(
    t: dict[str, int],
    different_ids: list[str],
    twoids_ids: list[str] | None = None,
    none_ids: list[str] | None = None,
) -> dict[str, str]:
    """The small line under each tile's number. Mirrors ``drawTiles`` in the page script."""
    return {
        "examined": "",
        "assigned": f"{t['examined'] - t['assigned']:,} without"
        if t["examined"] > t["assigned"]
        else "",
        "strict": _pct(t["strict"], t["assigned"]),
        "lenient": f"strict + {t['lenient_only']:,}" if t["lenient_only"] else "",
        "different": _short_list(different_ids) if t["different"] else "",
        "several": "",
        "twoids": _short_list(twoids_ids or []) if t["dups"] else "",
        "none": _short_list(none_ids or []) if t["none"] else "",
        "nodata": "",
        "unassigned": unassigned_sub(t),
    }


def unassigned_sub(t: dict[str, int]) -> str:
    """The line under the "No assigned ID" tile, in the category labels. Mirrors
    ``unassignedSub`` in the page script."""
    parts = []
    for key, label in (
        ("newid", "new ID found"),
        ("taken", "ID already used"),
        ("un_twoids", "two genotype IDs"),
        ("unclear", "no clear match"),
        ("un_nodata", "not scored"),
    ):
        if t[key]:
            parts.append(f"{t[key]:,} {label}")
    return " · ".join(parts)


def tile_html(
    ds: str,
    t: dict[str, int],
    different_ids: list[str],
    twoids_ids: list[str] | None = None,
    none_ids: list[str] | None = None,
) -> str:
    """The dashboard tiles with the run's counts (the page script updates them)."""
    subs = tile_subs(t, different_ids, twoids_ids, none_ids)
    out = []
    for key, label, icon, role, tip in _TILES:
        if key == "unassigned" and not t["unassigned"]:
            continue
        v = t[TILE_VALUE.get(key, key)]
        hidden = key in TILE_OPTIONAL and not v
        on = role and (v or role == "info")
        cls = f"dash-tile has-role r-{role}" if on else "dash-tile"
        out.append(
            f'<button type="button" class="{cls}" data-tile="{key}" aria-pressed="false" '
            f'title="{esc(tip)}"{" hidden" if hidden else ""}>'
            f'<span class="t-label"><span class="ic" aria-hidden="true">{esc(icon)}</span>{esc(label)}</span>'
            f'<span class="t-value">{v:,}</span><span class="t-run" hidden></span>'
            f'<span class="t-sub">{esc(subs[key])}</span></button>'
        )
    return (
        f'<div class="dash" role="group" aria-label="Summary" data-ds="{esc(ds)}">'
        + "".join(out)
        + "</div>"
    )


MATCH_SCORE_HELP = (
    "Match score: how well a genotype predicts a sample's measurements. "
    "Unrelated genotypes score about 2-4 by chance; most true matches score 5-9."
)
"""The one explanation of "match score" (under the slider panel)."""


def thresholds_html(z_threshold: float, z_diff: float) -> str:
    """The slider panel (moved into the visible dataset's panel by the page script)."""

    def row(key: str, label: str, hint: str, value: float, vmax: float, step: float) -> str:
        frac = min(max(value / vmax, 0.0), 1.0) if vmax > 0 else 0.0
        return (
            f'<div class="th-row"><label for="anq-{key}"><span class="th-l">{esc(label)}</span>'
            f'<span class="th-hint">{esc(hint)}</span></label>'
            f'<div class="th-track"><input type="range" id="anq-{key}-r" min="0" max="{vmax:g}" step="{step:g}" '
            f'value="{value:g}" aria-label="{esc(label)}" aria-valuetext="{value:.2f}">'
            f'<span class="th-mark" style="left:calc(8px + (100% - 16px) * {frac:.4f})" title="run value"></span></div>'
            f'<input type="number" id="anq-{key}" min="0" step="{step:g}" value="{value:.2f}" '
            f'inputmode="decimal" aria-label="{esc(label)}">'
            f'<span class="th-run">used {value:.2f}</span></div>'
        )

    zt_max = max(6.0, math.ceil(z_threshold * 1.5))
    zd_max = max(5.0, math.ceil(z_diff * 1.5))
    return (
        '<div class="th-band" id="anq-th">'
        + row(
            "zt",
            "Lenient cutoff",
            "Own ID score that still counts when another genotype ranks higher",
            z_threshold,
            zt_max,
            0.01,
        )
        + row(
            "zd",
            "Clear lead",
            "How far the best genotype must beat the next one",
            z_diff,
            zd_max,
            0.05,
        )
        + f'<p class="th-help">{esc(MATCH_SCORE_HELP)}</p>'
        + '<div class="th-foot"><span class="th-note" id="anq-th-note" hidden>'
        '<button type="button" class="linkbtn" id="anq-copy" title="YAML for a run with these values">'
        "Copy as config</button> "
        '<span id="anq-copy-msg" role="status" aria-live="polite"></span></span>'
        '<button type="button" class="btn" id="anq-reset" disabled>Reset</button></div></div>'
    )


def context_html(ds: str, omics: str, summary: dict[str, Any], n_pheno: int, hidden: bool) -> str:
    """The context line: counts, scoring chip and (for independent scoring) a short hint."""
    prefix = f"{omics}: " if omics and omics != "omics" else ""
    counts = (
        f"{prefix}{n_pheno:,} omics samples · {int(summary.get('n_prs_samples') or 0):,} genotypes · "
        f"{int(summary.get('n_traits_used') or 0):,} traits"
    )
    hint = (
        ""
        if summary.get("mahalanobis")
        else '<span class="hint">Tip: try Mahalanobis scoring (<code>--mahalanobis</code>)</span>'
    )
    return (
        f'<div class="ctx" data-ctx="{esc(ds)}"{" hidden" if hidden else ""}><span class="num">{esc(counts)}</span>'
        f'<span class="pill">{esc(scoring_text(summary))}</span>{hint}</div>'
    )


def alerts_html(summary: dict[str, Any]) -> str:
    """qpgen's warnings; the page script adds a line when its own statuses differ from
    the saved ones."""
    return (
        "".join(
            f'<p class="m-alert" role="note">⚠ {esc(w)}</p>' for w in summary.get("warnings") or []
        )
        + '<p class="m-alert js-disagree" role="note" hidden></p>'
    )


def _short_argv(argv: list[Any]) -> str:
    """The qpgen command with the binary shown by its file name."""
    if not argv:
        return ""
    return " ".join([Path(str(argv[0])).name, *(str(a) for a in argv[1:])])


def run_details_html(s: str, summary: dict[str, Any], weights: pd.DataFrame, n_pheno: int) -> str:
    """The collapsed run details: settings, trait weights, result codes, qpgen command and log."""
    dc = summary.get("dataset_counts") or {}
    n_cov = dc.get("n_cov_columns")
    cov = _yes(summary.get("use_cov")) + (
        f" ({n_cov} columns)" if summary.get("use_cov") and n_cov else ""
    )
    items = [
        ("Scoring", scoring_detail(summary)),
        ("Covariates", cov),
        ("Measurements rank-normalized", _yes(summary.get("rint"))),
        (
            "Lenient cutoff / clear lead",
            f"{float(summary.get('z_threshold', 0)):.2f} / {float(summary.get('z_diff', 0)):.2f}",
        ),
        ("Traits used", f"{int(summary.get('n_traits_used') or 0):,}"),
        ("Genotypes", f"{int(summary.get('n_prs_samples') or 0):,}"),
        ("Samples scored", f"{int(summary.get('n_scored') or 0):,} of {n_pheno:,}"),
        ("qpgen version", qpgen_version(summary)),
        (
            "Method",
            "Polygenic scores of known xQTLs, matched to the measured traits (qpgentools match-prs-pheno).",
        ),
    ]
    dl = (
        '<dl class="params">'
        + "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in items)
        + "</dl>"
    )
    n_le0 = int(summary.get("n_traits_weight_le0") or 0)
    wnote = f'<p class="rd-note">{n_le0:,} traits ignored (weight ≤ 0)</p>' if n_le0 else ""
    codes = (
        '<table class="codes"><thead><tr><th>Result</th><th>Code in samples.tsv</th></tr></thead><tbody>'
        + "".join(
            f"<tr><td>{esc(a)}</td><td><code>{esc(b)}</code></td></tr>" for a, b in CATEGORY_CODES
        )
        + "</tbody></table>"
    )
    argv = summary.get("qpgen_argv") or []
    excerpt = summary.get("qpgen_log_excerpt") or []
    nested = ""
    if argv:
        nested += (
            f"<details><summary>qpgen command</summary>{code_block(_short_argv(argv))}</details>"
        )
    if excerpt:
        nested += f"<details><summary>qpgen log excerpt</summary>{code_block(chr(10).join(excerpt))}</details>"
    return (
        '<details class="m-details run-details"><summary><h2>Run details</h2></summary><div class="rd-body">'
        + dl
        + "<h3>Trait weights</h3>"
        + '<div class="chart-card">'
        + figure(weight_figure(weights), f"fig-weights-{s}", height=300).replace(
            '<figure class="viz">', '<figure class="viz" style="margin:0">'
        )
        + "</div>"
        + wnote
        + "<h3>Result names in the files</h3>"
        + codes
        + f"<div>{nested}</div></div></details>"
    )


def dataset_panel(
    d: dict[str, Any], hidden: bool, z_threshold: float = 1.96, z_diff: float = 2.0
) -> tuple[str, dict[str, Any]]:
    """One dataset's panel (dashboard to run details) and its ``ANQ_DATA.datasets`` entry."""
    s = slug(d["name"])
    samples, summary = d["samples"], d["summary"]
    cats, assigned = run_categories(samples, z_threshold, z_diff)
    dups = [display_flags(r, z_threshold, z_diff)[0] for r in samples.to_dict("records")]
    t = tile_counts(cats, assigned, dups)
    diff_ids = [
        f"{p} → {b}"
        for p, b, c in zip(
            samples["pheno_id"].map(_str), samples["id_1"].map(_str), cats, strict=True
        )
        if c == "different"
    ]
    twoids_ids = [_str(p) for p, d in zip(samples["pheno_id"], dups, strict=True) if d]
    none_ids = [
        _str(p)
        for p, c, a in zip(samples["pheno_id"], cats, assigned, strict=True)
        if c == "nomatch" and a
    ]
    n = len(samples)
    # the "Assigned ID" column is shown only when some assigned ID differs from the sample ID
    same_ids = all(
        not _str(g) or _str(g) == _str(p)
        for p, g in zip(samples["pheno_id"], samples["geno_id_assigned"], strict=True)
    )
    body = [
        f'<div class="ds-panel" id="panel-{esc(s)}" data-ds="{esc(s)}"{" hidden" if hidden else ""}>',
        alerts_html(summary),
        tile_html(s, t, diff_ids, twoids_ids, none_ids),
        '<div class="th-slot"></div>',
        f'<section class="m-sec" id="check-{esc(s)}" aria-labelledby="check-h-{esc(s)}">'
        f'<div class="m-sec-head"><h2 id="check-h-{esc(s)}">Samples to check <span class="cnt"></span></h2></div>'
        '<div class="check-grid"></div><div class="more-link"></div></section>',
        f'<section class="m-sec" id="chart-{esc(s)}"><div class="m-sec-head"><h2>How clearly each sample '
        "matches one genotype</h2></div>"
        '<div class="chart-card"><div class="legend-row" role="group" aria-label="Results shown"></div>'
        + figure({"data": [], "layout": chart_layout()}, f"fig-match-{s}", height=460).replace(
            '<figure class="viz">', '<figure class="viz" style="margin:0">'
        )
        + "</div></section>",
        f'<section class="m-sec" id="dups-{esc(s)}"><div class="m-sec-head"><h2>Duplicates '
        '<span class="cnt"></span></h2></div><div class="dup-list"></div></section>',
        f'<details class="m-details all-samples" id="all-{esc(s)}" open><summary><h2>All samples '
        f'<span class="cnt">({n:,})</span></h2></summary><div class="m-details-body"><div class="tbar"></div>'
        '<table class="mtable"><colgroup></colgroup><thead></thead><tbody></tbody></table>'
        '<div class="tfoot"></div></div></details>',
        run_details_html(s, summary, d["weights"], n),
        "</div>",
    ]
    entry = {
        "name": d["name"],
        "slug": s,
        "omics": d["omics"],
        "rows": embed_rows(samples),
        "fig": f"fig-match-{s}",
        "sameIds": same_ids,
    }
    return "".join(body), entry


_DEFAULT_NAMES = ("omics", "dataset", "")


def _ds_labels(datasets: list[dict[str, Any]]) -> list[str]:
    """Dataset switch labels: omics names, else dataset names, else ``Dataset 1``, ...;
    never an internal default name."""
    for key in ("omics", "name"):
        names = [str(d[key]) for d in datasets]
        if len(set(names)) == len(names) and not any(n in _DEFAULT_NAMES for n in names):
            return names
    return [f"Dataset {i + 1}" for i in range(len(datasets))]


def build_page(
    project: str,
    z_threshold: float,
    z_diff: float,
    datasets: list[dict[str, Any]],
    nav: list[tuple[str, str]] | None = None,
) -> Page:
    """The match report page.

    ``datasets`` items: ``name``, ``omics``, ``summary`` (summary.json), ``samples``,
    ``swaps`` and ``weights`` (DataFrames read from the summary TSVs). ``project`` is
    kept for the signature; the page shows no internal names."""
    del project
    panels: list[str] = []
    entries: list[dict[str, Any]] = []
    header: list[str] = []
    for i, d in enumerate(datasets):
        html, entry = dataset_panel(d, hidden=i > 0, z_threshold=z_threshold, z_diff=z_diff)
        panels.append(html)
        entries.append(entry)
        header.append(
            context_html(
                entry["slug"], str(d["omics"]), d["summary"], len(d["samples"]), hidden=i > 0
            )
        )
    if len(datasets) > 1:
        header.append(
            '<div class="seg" role="group" aria-label="Dataset">'
            + "".join(
                f'<button type="button" data-ds-set="{esc(e["slug"])}" aria-pressed="{"true" if i == 0 else "false"}">'
                f"{esc(lbl)}</button>"
                for i, (e, lbl) in enumerate(zip(entries, _ds_labels(datasets), strict=True))
            )
            + "</div>"
        )
    body = thresholds_html(z_threshold, z_diff) + "".join(panels)
    data = {
        "config": {"z_threshold": z_threshold, "z_diff": z_diff},
        "fields": ROW_FIELDS,
        "categories": CATEGORIES,
        "colors": {
            "role": ROLE_COLORS,
            "warnLine": WARNING_LINE,
            "surface": LIGHT["surface"],
            "axis": LIGHT["axis"],
            "muted": LIGHT["muted"],
            "ink2": LIGHT["ink2"],
            "ink": LIGHT["ink"],
        },
        "chartLayout": chart_layout(),
        "datasets": entries,
        "notScored": NOT_SCORED,
    }
    return Page(
        title=TITLE,
        sections=[Section(id="match", title="", body=body)],
        nav=list(nav or []),
        extra_js=MATCH_RULES_JS + "\n" + MATCH_PAGE_JS,
        data=data,
        header_html="".join(header),
        hide_brand=True,
        body_class="page-match",
        toc=False,
    )


def load_dataset(run: StepContext, ds: str, omics: str) -> dict[str, Any]:
    """Read one dataset's summaries for :func:`build_page`."""
    L = Layout
    try:
        summary = read_json(run.path(L.match_summary(ds, "summary.json")))
        samples = read_tsv(run.path(L.match_summary(ds, "samples.tsv")))
        swaps = read_tsv(run.path(L.match_summary(ds, "swaps.tsv")))
        weights = read_tsv(run.path(L.match_summary(ds, "weights.tsv")))
    except (OSError, ValueError) as e:
        raise StepError(f"cannot read the match summaries of dataset {ds}: {e}") from e
    return {
        "name": ds,
        "omics": summary.get("omics") or omics,
        "summary": summary,
        "samples": samples,
        "swaps": swaps,
        "weights": weights,
    }


@register("report_match")
def report_match(ctx: StepContext) -> None:
    """Render ``report/match.html`` from every dataset's summaries."""
    p = ctx.params
    omics = p.get("omics") or {}
    datasets = [load_dataset(ctx, ds, omics.get(ds, "omics")) for ds in p["datasets"]]
    page = build_page(
        p.get("project", "anqchor"),
        float(p["z_threshold"]),
        float(p["z_diff"]),
        datasets,
        nav=index_nav(ctx.params),
    )
    write_page(ctx.path(Layout.report("match")), page, p.get("plotlyjs", "embed"))
