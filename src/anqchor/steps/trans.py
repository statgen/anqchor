"""``trans.*`` step implementations: rectangle lists and the trans summarizer.

``trans_lists`` writes the variant × trait rectangle tested by ``rect-assoc``.

``summarize_trans`` streams the run's ``rect-assoc`` matrix once and classifies every
(variant, trait) cell with the trait loci of ``prepare/datasets/<ds>/traits.tsv``:

* ``cis``: same chromosome and ``|pos - start| <= cis_window_bp`` (start = gene TSS);
* ``trans``: every other cell whose trait has a locus, split into ``cross_chrom``
  (different chromosomes) and ``same_chrom_distal`` (same chromosome beyond the window);
* ``no_locus``: the trait has no locus (counted only).

It writes three things for the report (DESIGN 11.2-Q):

1. **Reference-pair replication by class.** Every reference pair is classified the same
   way (cis or trans) and counted as tested / same sign / same sign and p < 0.05 / same
   sign and p < 1e-5 (the cis report's levels).
2. **Hotspots.** Trans cells only. Variants are grouped into bins of ``hotspot_bin_bp`` by
   position; a bin's count at a p-value threshold is the number of distinct genes with a
   trans cell below the threshold from any variant of the bin, leaving out genes in cis
   with the bin (gene start within ``cis_window_bp`` of the bin's tested variants). A
   hotspot is a bin with at least k genes. Counts are precomputed for a grid of thresholds
   (:data:`HOTSPOT_ALPHAS`, the configured levels and Bonferroni over the trans cells) × k
   (:data:`HOTSPOT_K` and ``hotspot_min_traits``).
3. **Calibration.** The cross-chromosome cells, leaving out (when the reference has trans
   pairs) every cell of a variant within ``cis_window_bp`` of a reference trans pair's
   variant and every cell of that pair's gene: λ_GC (median χ² / 0.4549), observed vs
   expected counts at p < 1e-3 and 1e-5, and thinned QQ quantiles with a 95% band.

Memory bound
------------
Besides one read chunk (``CHUNK_ROWS`` rows × traits × 3 statistics as float64): the
``-log10 p`` of every cross-chromosome cell used for calibration as float32 (4 bytes per
cell, plus a sorted copy at the end), O(variants + traits) accumulators, the reference
pairs, and the cells with ``p`` below the loosest grid threshold (or below Bonferroni over
the trans cells seen so far, when that is looser), with beta and SE.
"""

from __future__ import annotations

import math
import time
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import special, stats

from anqchor.io.genotypes import normalize_chrom
from anqchor.io.pairs import parse_cpra, read_pairs_table, read_variant_list, sort_variants
from anqchor.io.tables import examples, read_lines, write_lines
from anqchor.pipelines.layout import Layout
from anqchor.qpgen.outputs import QpgenOutputError, RectChunk, iter_rect_assoc, read_rect_header
from anqchor.stats import passes
from anqchor.steps import StepContext, StepError, register
from anqchor.steps.prepare import input_errors
from anqchor.util.files import atomic_write, read_header, read_json, write_json

SCHEMA_VERSION = 2
"""2: one run per directory; reference pairs by cis/trans class, hotspot bins and grid,
cross-chromosome calibration with exclusion (DESIGN 11.2-Q)."""

CELL_CLASSES = ("cis", "cross_chrom", "same_chrom_distal", "no_locus")
"""Disjoint cell classes; ``trans`` = ``cross_chrom`` + ``same_chrom_distal``."""
CIS, CROSS, DISTAL, NO_LOCUS = range(4)
NOT_TESTED = -1

CHI2_MEDIAN_1DF = 0.454936423119572
"""Median of the chi-square distribution with one degree of freedom."""
BONFERRONI_ALPHA = 0.05
LAMBDA_WARN = 1.1
CHUNK_ROWS = 16
QQ_HEAD = 250
QQ_MAX_POINTS = 2000
QQ_BAND = 0.95

REPLICATION_LEVELS = {"nominal": 0.05, "strict": 1e-5}
"""Reference-pair replication levels (same sign and p below each), as in the cis report."""
HOTSPOT_ALPHAS = (1e-4, 1e-5, 1e-6, 5e-8, 1e-8, 1e-10)
"""p-value thresholds of the hotspot grid (plus the configured levels and Bonferroni)."""
HOTSPOT_K = (1, 2, 3, 5, 10, 20)
"""Minimum genes per bin of the hotspot grid (plus ``hotspot_min_traits``)."""
INFLATION_ALPHAS = (1e-3, 1e-5)
"""p-value thresholds of the observed vs expected cross-chromosome counts."""

DEFINITIONS = {
    "cis": "same chromosome and variant within cis_window_bp of the gene start",
    "trans": "every other pair with a gene locus: another chromosome, or the same chromosome beyond "
    "cis_window_bp",
    "replication": "same sign as the reference, and p below each level",
    "hotspot": "a bin of hotspot_bin_bp with at least k distinct genes that have a trans signal "
    "(p below the threshold) from any variant of the bin; genes in cis with the bin (start within "
    "cis_window_bp of its tested variants) are not counted",
    "calibration": "cross-chromosome cells; when the reference has trans pairs, the cells of variants "
    "within cis_window_bp of a reference trans variant and of the reference trans genes are left out",
    "lambda_gc": "median chi-square (1 df) of the cells / 0.4549",
    "bonferroni": "p < 0.05 / trans cells tested",
}


def _log(ctx: StepContext, msg: str, *args: object) -> None:
    if ctx.log is not None:
        ctx.log.info(msg, *args)


def _warn(ctx: StepContext, msg: str, *args: object) -> None:
    if ctx.log is not None:
        ctx.log.warning(msg, *args)


# ======================================================================= lists
def _genotype_chroms(geno: dict[str, Any]) -> list[str]:
    """Chromosomes in pgen-list order."""
    chroms = [str(e["chrom"]) for e in geno.get("entries") or [] if "chrom" in e]
    return list(dict.fromkeys(chroms or [str(c) for c in geno.get("chromosomes") or []]))


@dataclass
class TraitSelection:
    """Outcome of the ``trait_chroms`` restriction."""

    kept: list[str]
    allowed: list[str] | None
    """Allowed chromosomes (genotype style); ``None`` when ``trait_chroms`` is ``all``."""
    n_candidates: int
    n_off_chrom: int
    n_no_locus: int
    examples_off: list[str] = field(default_factory=list)
    examples_no_locus: list[str] = field(default_factory=list)

    def as_json(self) -> dict[str, Any]:
        return {
            "allowed_chroms": self.allowed,
            "n_candidates": self.n_candidates,
            "n_kept": len(self.kept),
            "n_dropped_off_chrom": self.n_off_chrom,
            "n_dropped_no_locus": self.n_no_locus,
            "examples_dropped_off_chrom": self.examples_off,
            "examples_dropped_no_locus": self.examples_no_locus,
        }


def trait_chrom_map(loci: pd.DataFrame, style: str) -> dict[str, str]:
    """Trait → chromosome (genotype style) for the traits with a locus."""
    out: dict[str, str] = {}
    for t, c in zip(loci["trait"].astype(str), loci["chrom"], strict=True):
        if isinstance(c, str) and c and c != "NA":
            out[t] = normalize_chrom(c, style)  # type: ignore[arg-type]
    return out


def no_trait_kept_message(ds: str, tc: Any, mode: str, n_traits: int, sel: TraitSelection) -> str:
    """Error text when ``trait_chroms`` keeps no trait (examples only where there are some)."""
    hint = (
        "no trait has a locus: give the dataset a trait_annot, or set "
        "pipelines.trans.trait_chroms: all"
        if sel.n_no_locus == n_traits
        else "set pipelines.trans.trait_chroms to a list of chromosomes, or 'all'"
    )
    off = f"{sel.n_off_chrom} traits on other chromosomes"
    if sel.examples_off:
        off += f" (e.g. {examples(sel.examples_off, 3)})"
    return (
        f"trans lists for {ds}: trait_chroms {tc!r} keeps none of the {n_traits} traits "
        f"(traits: {mode}); allowed chromosomes {examples(sel.allowed or [], 5)}, "
        f"{off}, {sel.n_no_locus} without a locus; {hint}"
    )


def select_trait_chroms(
    traits: list[str],
    trait_chrom: dict[str, str],
    trait_chroms: str | list[str],
    geno_chroms: list[str],
    style: str,
) -> TraitSelection:
    """Apply ``trait_chroms`` to a candidate trait list (order kept).

    ``genotyped`` keeps traits on the genotype set's chromosomes, a list keeps traits on the
    listed chromosomes (``chr`` or plain spelling), ``all`` keeps everything. Except with
    ``all``, traits without a locus are dropped."""
    if trait_chroms == "all":
        return TraitSelection(
            kept=list(traits), allowed=None, n_candidates=len(traits), n_off_chrom=0, n_no_locus=0
        )
    if trait_chroms == "genotyped":
        allowed = list(geno_chroms)
    elif isinstance(trait_chroms, str):
        raise StepError(f"trait_chroms must be 'genotyped', 'all' or a list, got {trait_chroms!r}")
    else:
        allowed = list(dict.fromkeys(normalize_chrom(str(c), style) for c in trait_chroms))  # type: ignore[arg-type]
    ok = set(allowed)
    kept: list[str] = []
    off: list[str] = []
    none: list[str] = []
    for t in traits:
        c = trait_chrom.get(t)
        if c is None:
            none.append(t)
        elif c in ok:
            kept.append(t)
        else:
            off.append(t)
    return TraitSelection(
        kept=kept,
        allowed=allowed,
        n_candidates=len(traits),
        n_off_chrom=len(off),
        n_no_locus=len(none),
        examples_off=off[:5],
        examples_no_locus=none[:5],
    )


@register("trans_lists")
@input_errors
def trans_lists(ctx: StepContext) -> None:
    """Write the rectangle's variant list (sorted, deduplicated, genotype chromosome style),
    its trait list (dataset trait IDs, restricted by ``trait_chroms``) and their counts."""
    p = ctx.params
    ds, ref = p["dataset"], p["reference"]
    geno = read_json(ctx.path(Layout.geno_json))
    style = geno["chrom_style"]
    chroms = _genotype_chroms(geno)
    order = {c: i for i, c in enumerate(chroms)}
    vsrc = Layout.ref_variants_present(ref) if p["variants"] == "reference" else Path(p["variants"])
    found, n_invalid, n_renamed = read_variant_list(ctx.path(vsrc), style)
    off = [v for v in found.values() if v.chrom not in order]
    keep = [v for v in found.values() if v.chrom in order]
    if not keep:
        raise StepError(
            f"trans lists for {ds}: none of the {len(found)} variants in {vsrc} is on a genotyped "
            f"chromosome ({examples(chroms, 5)}); variants look like "
            f"{examples([v.cpra for v in found.values()], 3)}"
        )
    variants = sort_variants(keep, order)
    mode = p["traits"]
    if mode == "reference":
        traits = list(dict.fromkeys(read_lines(ctx.path(Layout.dsref_traits(ds, ref)))))
        missing: list[str] = []
    else:
        ds_traits = read_header(ctx.path(Layout.ds_pheno(ds)))[2:]
        if mode == "all":
            traits, missing = ds_traits, []
        else:
            wanted = list(dict.fromkeys(line.split()[0] for line in read_lines(ctx.path(mode))))
            have = set(ds_traits)
            want_set = set(wanted)
            traits = [t for t in ds_traits if t in want_set]
            missing = [t for t in wanted if t not in have]
    if not traits:
        raise StepError(
            f"trans lists for {ds}: no trait to test (traits: {mode}); "
            + (
                f"none of the listed traits ({examples(missing, 3)}) is in the dataset"
                if missing
                else "the reference maps no trait onto the dataset"
            )
        )
    tc = p.get("trait_chroms", "all")
    loci = read_loci(ctx.path(Layout.ds_traits(ds)))
    sel = select_trait_chroms(traits, trait_chrom_map(loci, style), tc, chroms, style)
    if not sel.kept:
        raise StepError(no_trait_kept_message(ds, tc, mode, len(traits), sel))
    write_lines(ctx.path(Layout.trans_variants(ds)), [v.cpra for v in variants])
    write_lines(ctx.path(Layout.trans_traits(ds)), sel.kept)
    write_json(
        ctx.path(Layout.trans_lists_json(ds)),
        {
            "schema_version": SCHEMA_VERSION,
            "dataset": ds,
            "reference": ref,
            "variants": {
                "source": p["variants"],
                "n_listed": len(found) + n_invalid,
                "n_invalid": n_invalid,
                "n_renamed_chrom": n_renamed,
                "n_off_genome": len(off),
                "n_kept": len(variants),
            },
            "traits": {
                "source": mode,
                "trait_chroms": tc,
                "n_missing_from_dataset": len(missing),
                **sel.as_json(),
            },
        },
    )
    if off:
        _warn(
            ctx,
            "trans lists for %s: %d variant(s) on chromosomes without genotypes dropped, e.g. %s",
            ds,
            len(off),
            examples([v.cpra for v in off], 3),
        )
    if missing:
        _warn(
            ctx,
            "trans lists for %s: %d listed trait(s) are not in the dataset, e.g. %s",
            ds,
            len(missing),
            examples(missing, 3),
        )
    if sel.n_no_locus:
        _warn(
            ctx,
            "trans lists for %s: %d trait(s) without a locus dropped (trait_chroms %s), e.g. %s",
            ds,
            sel.n_no_locus,
            tc,
            examples(sel.examples_no_locus, 3),
        )
    _log(
        ctx,
        "trans lists for %s: %d variants (%d invalid, %d renamed, %d off-genome); "
        "%d of %d traits kept by trait_chroms %s (%d on other chromosomes, %d without locus)",
        ds,
        len(variants),
        n_invalid,
        n_renamed,
        len(off),
        len(sel.kept),
        len(traits),
        tc,
        sel.n_off_chrom,
        sel.n_no_locus,
    )


# ================================================================== statistics
def lambda_gc(median_log10p: float | None) -> float:
    """λ_GC from the median ``-log10 p``: the 1-df chi-square quantile of the median p
    divided by 0.4549 (computed in log space, so tiny p-values stay finite)."""
    if median_log10p is None or not math.isfinite(median_log10p):
        return float("nan")
    log_half_p = -median_log10p * math.log(10.0) - math.log(2.0)
    z = float(special.ndtri_exp(log_half_p))
    return z * z / CHI2_MEDIAN_1DF


def median_sorted(desc: np.ndarray) -> float:
    """Median of an array sorted in either order (NaN when empty)."""
    n = len(desc)
    if n == 0:
        return float("nan")
    if n % 2:
        return float(desc[n // 2])
    return float((float(desc[n // 2 - 1]) + float(desc[n // 2])) / 2.0)


def qq_points(
    desc: np.ndarray, max_points: int = QQ_MAX_POINTS, head: int = QQ_HEAD
) -> pd.DataFrame:
    """Expected vs observed ``-log10 p`` for cells sorted by decreasing ``-log10 p``, with
    the 95% band of the null.

    Rank ``r`` (1 = smallest p) of ``n`` uniform p-values follows Beta(r, n − r + 1); the
    band is its 2.5% and 97.5% quantiles on the ``-log10`` scale (``band_lo`` ≤ ``band_hi``).
    The ``head`` strongest cells are kept individually (the informative tail); the rest
    are thinned on a log-spaced rank grid so that at most ``max_points`` remain."""
    n = len(desc)
    cols = ["rank", "expected_log10p", "observed_log10p", "band_lo_log10p", "band_hi_log10p"]
    if n == 0:
        return pd.DataFrame({c: [] for c in cols})
    if n <= max_points:
        ranks = np.arange(1, n + 1)
    else:
        h = min(head, max_points // 2)
        tail = np.unique(np.round(np.geomspace(h + 1, n, max_points - h)).astype(np.int64))
        ranks = np.unique(np.concatenate([np.arange(1, h + 1), tail, [n]]))
    expected = -np.log10((ranks - 0.5) / n)
    observed = np.asarray(desc[ranks - 1], dtype=np.float64)
    a, b = ranks.astype(float), (n - ranks + 1).astype(float)
    tail_p = (1 - QQ_BAND) / 2
    lo = -np.log10(stats.beta.ppf(1 - tail_p, a, b))
    hi = -np.log10(stats.beta.ppf(tail_p, a, b))
    return pd.DataFrame(dict(zip(cols, [ranks, expected, observed, lo, hi], strict=True)))


def fmt_alpha(a: float) -> str:
    """``1e-05`` → ``1e-5``, ``5e-08`` → ``5e-8``, ``0.05`` → ``0.05``."""
    if a >= 1e-3:
        return f"{a:g}"
    m, e = f"{a:.1e}".split("e")
    m = m.rstrip("0").rstrip(".")
    return f"{m}e{int(e)}"


@dataclass
class Threshold:
    """One p-value threshold of the hotspot grid."""

    key: str
    """Column suffix in ``hotspots.tsv`` (``p1e-5``, ``bonferroni``)."""
    label: str
    alpha: float
    roles: list[str] = field(default_factory=list)
    """Configured levels at this threshold (``nominal``, ``strict``, ``bonferroni``)."""

    @property
    def log10p(self) -> float:
        return -math.log10(self.alpha)

    def as_json(self) -> dict[str, Any]:
        return {
            "key": self.key,
            "label": self.label,
            "alpha": self.alpha,
            "log10p": self.log10p,
            "roles": self.roles,
        }


def grid_thresholds(levels: dict[str, Any], n_trans: int) -> list[Threshold]:
    """The hotspot grid's thresholds, loosest first: :data:`HOTSPOT_ALPHAS`, the nominal and
    strict levels and (with ``bonferroni``) 0.05 / trans cells. Equal alphas merge."""
    out: dict[float, Threshold] = {}

    def add(
        alpha: float, role: str | None, key: str | None = None, label: str | None = None
    ) -> None:
        a = float(f"{alpha:.6g}")
        if a not in out:
            out[a] = Threshold(key or f"p{fmt_alpha(a)}", label or f"p < {fmt_alpha(a)}", a)
        if role:
            out[a].roles.append(role)

    for a in HOTSPOT_ALPHAS:
        add(a, None)
    for role in ("nominal", "strict"):
        if levels.get(role) is not None:
            add(float(levels[role]), role)
    if levels.get("bonferroni", True) and n_trans > 0:
        b = BONFERRONI_ALPHA / n_trans
        add(
            b,
            "bonferroni",
            key="bonferroni",
            label=f"Bonferroni (p < {fmt_alpha(float(f'{b:.2g}'))})",
        )
    return sorted(out.values(), key=lambda t: -t.alpha)


def grid_k(hotspot_min_traits: int) -> list[int]:
    return sorted({*HOTSPOT_K, int(hotspot_min_traits)})


# ================================================================= the stream
@dataclass
class TraitInfo:
    """Trait loci aligned to a rectangle's trait order."""

    names: list[str]
    chrom_code: np.ndarray
    """int32 chromosome code, -1 when the trait has no locus."""
    start: np.ndarray
    chrom: list[str | None]


@dataclass
class KnownPair:
    trait: str
    variant: str
    beta_ref: float
    se_ref: float
    pval_ref: float
    ref_trait: str
    present: bool


def _trait_info(names: list[str], loci: pd.DataFrame, codes: dict[str, int]) -> TraitInfo:
    by = loci.drop_duplicates("trait").set_index("trait")
    cc = np.full(len(names), -1, dtype=np.int32)
    st = np.zeros(len(names), dtype=np.int64)
    chroms: list[str | None] = []
    for i, t in enumerate(names):
        if t in by.index:
            row = by.loc[t]
            c, s = row["chrom"], row["start"]
            if isinstance(c, str) and c and c != "NA" and pd.notna(s):
                cc[i] = codes.setdefault(c, len(codes))
                st[i] = int(s)
                chroms.append(c)
                continue
        chroms.append(None)
    return TraitInfo(names=names, chrom_code=cc, start=st, chrom=chroms)


def _row_keys(meta: pd.DataFrame) -> list[tuple[str, str]]:
    """(ID, CHROM:GENPOS:ALLELE0:ALLELE1) of each row, used to find reference pairs."""
    cpra = (
        meta["CHROM"].astype(str)
        + ":"
        + meta["GENPOS"].astype(str)
        + ":"
        + meta["ALLELE0"].astype(str)
        + ":"
        + meta["ALLELE1"].astype(str)
    )
    return list(zip(meta["ID"].astype(str), cpra, strict=True))


def pair_class(
    vchrom: str | None, vpos: int | None, tchrom: str | None, tstart: int | None, window: int
) -> tuple[str, str]:
    """(class, subclass) of a variant-gene pair: ``cis``/``cis``, ``trans``/``cross_chrom``,
    ``trans``/``same_chrom_distal`` or ``no_locus``/``no_locus`` (gene or variant position
    unknown)."""
    if vchrom is None or vpos is None or tchrom is None or tstart is None:
        return "no_locus", "no_locus"
    if vchrom != tchrom:
        return "trans", "cross_chrom"
    if abs(int(vpos) - int(tstart)) <= window:
        return "cis", "cis"
    return "trans", "same_chrom_distal"


@dataclass
class Exclusion:
    """Regions and genes of the reference trans pairs, left out of the calibration."""

    regions: dict[str, np.ndarray]
    """chromosome → sorted positions of reference trans variants."""
    genes: set[str]
    window: int
    n_pairs: int

    def rows(self, chroms: list[str], pos: np.ndarray) -> np.ndarray:
        """Rows (variants) within ``window`` of a reference trans variant."""
        out = np.zeros(len(pos), dtype=bool)
        for i, (c, p) in enumerate(zip(chroms, pos, strict=True)):
            arr = self.regions.get(c)
            if arr is None or not len(arr):
                continue
            j = int(np.searchsorted(arr, p))
            near = [arr[k] for k in (j - 1, j) if 0 <= k < len(arr)]
            out[i] = any(abs(int(p) - int(x)) <= self.window for x in near)
        return out


def build_exclusion(
    pairs: list[KnownPair], loci: pd.DataFrame, style: str, window: int
) -> Exclusion:
    """The reference trans pairs (with a parseable variant and a gene locus) as exclusion
    regions and genes."""
    by = loci.drop_duplicates("trait").set_index("trait")
    regions: dict[str, list[int]] = {}
    genes: set[str] = set()
    n = 0
    for kp in pairs:
        v = parse_cpra(kp.variant)
        if v is None or kp.trait not in by.index:
            continue
        c, s = by.loc[kp.trait, "chrom"], by.loc[kp.trait, "start"]
        if not (isinstance(c, str) and c and c != "NA" and pd.notna(s)):
            continue
        vchrom = normalize_chrom(v.chrom, style)  # type: ignore[arg-type]
        if pair_class(vchrom, v.pos, c, int(s), window)[0] != "trans":
            continue
        regions.setdefault(vchrom, []).append(v.pos)
        genes.add(kp.trait)
        n += 1
    return Exclusion(
        {c: np.array(sorted(p), dtype=np.int64) for c, p in regions.items()}, genes, window, n
    )


@dataclass
class RectResult:
    """Everything the writers need from the rectangle."""

    traits: TraitInfo
    meta: pd.DataFrame
    n_cells: dict[str, int]
    """Tested cells per :data:`CELL_CLASSES` entry."""
    cross_desc: np.ndarray
    """``-log10 p`` of the cross-chromosome cells used for calibration, sorted decreasingly."""
    n_cross_excluded: int
    row_excluded: np.ndarray
    per_row: dict[str, Any]
    kept: dict[str, np.ndarray]
    """Retained cells: ``row``, ``col``, ``cls``, ``lp``, ``beta``, ``se``, ``ref`` (reference pair)."""
    known_cells: dict[tuple[int, int], tuple[float, float, float]]
    seconds: float


def summarize_rect(
    chunks: Iterable[RectChunk],
    *,
    trait_names: list[str],
    loci: pd.DataFrame,
    known: dict[str, list[str]],
    exclusion: Exclusion,
    cis_window_bp: int,
    keep_log10p: float,
    nominal_log10p: float,
    bonferroni: bool,
) -> RectResult:
    """Stream the rectangle once (see the module docstring).

    Cells with ``-log10 p > keep_log10p`` are retained, and also, when ``bonferroni``,
    the trans cells with ``-log10 p > log10(trans cells seen so far / 0.05)``: the final
    Bonferroni threshold is at least as strict, so every Bonferroni-significant cell is
    kept whatever the row order."""
    t0 = time.time()
    codes: dict[str, int] = {}
    tinfo = _trait_info(trait_names, loci, codes)
    tidx = {t: i for i, t in enumerate(trait_names)}
    known_idx = {
        v: np.array(sorted({tidx[t] for t in ts if t in tidx}), dtype=np.int64)
        for v, ts in known.items()
    }
    known_idx = {v: a for v, a in known_idx.items() if len(a)}
    n_traits = len(trait_names)
    has_locus = tinfo.chrom_code >= 0
    col_excl = np.array([t in exclusion.genes for t in trait_names], dtype=bool)
    names_arr = np.array(trait_names, dtype=object)

    n_cells = {c: 0 for c in CELL_CLASSES}
    cross_parts: list[np.ndarray] = []
    n_cross_excl = 0
    metas: list[pd.DataFrame] = []
    row_excl_parts: list[np.ndarray] = []
    row_stats: dict[str, list[Any]] = {
        k: []
        for k in (
            "n_cis",
            "n_cross_chrom",
            "n_same_chrom_distal",
            "n_sig_trans",
            "max_log10p_trans",
            "top_trans_trait",
            "max_log10p_cis",
            "top_cis_trait",
        )
    }
    kept: dict[str, list[np.ndarray]] = {
        k: [] for k in ("row", "col", "cls", "lp", "beta", "se", "ref")
    }
    known_cells: dict[tuple[int, int], tuple[float, float, float]] = {}
    n_trans_seen = 0
    n_rows = 0

    for ch in chunks:
        meta = ch.meta
        r = len(meta)
        if r == 0:
            continue
        L = ch.stats["LOG10P"]
        B = ch.stats["BETA"]
        S = ch.stats["SE"]
        vchrom = meta["CHROM"].astype(str).tolist()
        vc = np.array([codes.setdefault(c, len(codes)) for c in vchrom], dtype=np.int32)
        pos = meta["GENPOS"].to_numpy(dtype=np.int64)
        fin = np.isfinite(L)
        same = tinfo.chrom_code[None, :] == vc[:, None]
        cis = same & (np.abs(pos[:, None] - tinfo.start[None, :]) <= cis_window_bp)
        cls = np.full((r, n_traits), NO_LOCUS, dtype=np.int8)
        cls[has_locus[None, :] & ~same] = CROSS
        cls[same & ~cis] = DISTAL
        cls[cis] = CIS
        cls[~fin] = NOT_TESTED
        ref = np.zeros((r, n_traits), dtype=bool)
        for i, (vid, cpra) in enumerate(_row_keys(meta)):
            idx = known_idx.get(vid)
            if idx is None:
                idx = known_idx.get(cpra)
            if idx is not None:
                ref[i, idx] = True
                for j in idx:
                    known_cells[(n_rows + i, int(j))] = (
                        float(L[i, j]),
                        float(B[i, j]),
                        float(S[i, j]),
                    )
        for k, c in enumerate(CELL_CLASSES):
            n_cells[c] += int(np.count_nonzero(cls == k))
        is_cross = cls == CROSS
        is_trans = is_cross | (cls == DISTAL)
        is_cis = cls == CIS
        rex = exclusion.rows(vchrom, pos)
        row_excl_parts.append(rex)
        use = is_cross & ~rex[:, None] & ~col_excl[None, :]
        cross_parts.append(L[use].astype(np.float32))
        n_cross_excl += int(np.count_nonzero(is_cross)) - int(np.count_nonzero(use))
        # per-variant statistics
        row_stats["n_cis"].append(is_cis.sum(axis=1))
        row_stats["n_cross_chrom"].append(is_cross.sum(axis=1))
        row_stats["n_same_chrom_distal"].append((cls == DISTAL).sum(axis=1))
        row_stats["n_sig_trans"].append((is_trans & passes(L, nominal_log10p)).sum(axis=1))
        for key, mask in (("trans", is_trans), ("cis", is_cis)):
            Lm = np.where(mask, L, -np.inf)
            am = Lm.argmax(axis=1) if n_traits else np.zeros(r, np.int64)
            mx = Lm[np.arange(r), am] if n_traits else np.full(r, -np.inf)
            ok = np.isfinite(mx)
            row_stats[f"max_log10p_{key}"].append(np.where(ok, mx, np.nan))
            row_stats[f"top_{key}_trait"].append(
                np.where(ok, names_arr[am], None) if n_traits else np.full(r, None, dtype=object)
            )
        # retained cells
        n_trans_seen += int(np.count_nonzero(is_trans))
        t_seen = (
            math.log10(n_trans_seen / BONFERRONI_ALPHA) if bonferroni and n_trans_seen else math.inf
        )
        keep = fin & (passes(L, keep_log10p) | (is_trans & passes(L, t_seen)))
        ri, ci = np.nonzero(keep)
        for k2, v in (
            ("row", n_rows + ri),
            ("col", ci.astype(np.int32)),
            ("cls", cls[ri, ci]),
            ("lp", L[ri, ci]),
            ("beta", B[ri, ci]),
            ("se", S[ri, ci]),
            ("ref", ref[ri, ci]),
        ):
            kept[k2].append(v)
        metas.append(meta)
        n_rows += r

    meta_all = (
        pd.concat(metas, ignore_index=True)
        if metas
        else pd.DataFrame(columns=["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N"])
    )
    cross = np.concatenate(cross_parts) if cross_parts else np.zeros(0, np.float32)
    del cross_parts
    cross.sort()
    cross = cross[::-1]
    per_row: dict[str, Any] = {}
    for k3, parts in row_stats.items():
        per_row[k3] = np.concatenate([np.asarray(x) for x in parts]) if parts else np.zeros(0)
    dtypes = {
        "row": np.int64,
        "col": np.int32,
        "cls": np.int8,
        "lp": float,
        "beta": float,
        "se": float,
        "ref": bool,
    }
    kept_arr = {k4: (np.concatenate(v) if v else np.zeros(0, dtypes[k4])) for k4, v in kept.items()}
    return RectResult(
        traits=tinfo,
        meta=meta_all,
        n_cells=n_cells,
        cross_desc=cross,
        n_cross_excluded=n_cross_excl,
        row_excluded=np.concatenate(row_excl_parts) if row_excl_parts else np.zeros(0, bool),
        per_row=per_row,
        kept=kept_arr,
        known_cells=known_cells,
        seconds=time.time() - t0,
    )


# ===================================================================== writers
def _fmt_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Booleans as 0/1 so that every TSV reads back numerically."""
    out = df.copy()
    for c in out.columns:
        if out[c].dtype == bool:
            out[c] = out[c].astype(np.int8)
        elif isinstance(out[c].dtype, pd.BooleanDtype):
            out[c] = out[c].astype("Int8")
    return out


def write_tsv(path: str | Path, df: pd.DataFrame) -> None:
    """Atomically write a TSV (gzip when the name ends in .gz), ``NA`` for missing."""
    with atomic_write(path) as fh:
        _fmt_frame(df).to_csv(fh, sep="\t", index=False, na_rep="NA", float_format="%.6g")


def bins_of(meta: pd.DataFrame, bin_bp: int, chrom_order: list[str]) -> pd.DataFrame:
    """One row per bin with a tested variant: ``bin`` (0-based order along the genome),
    ``chrom``, ``bin_start``, ``bin_end``, ``n_variants``, ``min_pos``, ``max_pos``; and
    the bin of every row (``row_bin``) as the frame's ``attrs``."""
    if not len(meta):
        df = pd.DataFrame(
            columns=["bin", "chrom", "bin_start", "bin_end", "n_variants", "min_pos", "max_pos"]
        )
        df.attrs["row_bin"] = np.zeros(0, np.int64)
        return df
    chrom = meta["CHROM"].astype(str).to_numpy()
    pos = meta["GENPOS"].to_numpy(dtype=np.int64)
    idx = pos // bin_bp
    order = {c: i for i, c in enumerate(chrom_order)}
    rank = np.array([order.get(c, len(order)) for c in chrom])
    keys = pd.DataFrame({"rank": rank, "chrom": chrom, "idx": idx, "pos": pos})
    g = (
        keys.groupby(["rank", "chrom", "idx"], sort=True)["pos"]
        .agg(["size", "min", "max"])
        .reset_index()
    )
    g["bin"] = np.arange(len(g))
    lookup = {(c, int(i)): int(b) for c, i, b in zip(g["chrom"], g["idx"], g["bin"], strict=True)}
    df = pd.DataFrame(
        {
            "bin": g["bin"].to_numpy(),
            "chrom": g["chrom"].to_numpy(),
            "bin_start": (g["idx"] * bin_bp).to_numpy(dtype=np.int64),
            "bin_end": ((g["idx"] + 1) * bin_bp).to_numpy(dtype=np.int64),
            "n_variants": g["size"].to_numpy(dtype=np.int64),
            "min_pos": g["min"].to_numpy(dtype=np.int64),
            "max_pos": g["max"].to_numpy(dtype=np.int64),
        }
    )
    df.attrs["row_bin"] = np.array(
        [lookup[(c, int(i))] for c, i in zip(chrom, idx, strict=True)], dtype=np.int64
    )
    return df


def hotspot_tables(
    res: RectResult, bins: pd.DataFrame, thresholds: list[Threshold], ks: list[int], window: int
) -> tuple[pd.DataFrame, pd.DataFrame, list[list[int]]]:
    """``hotspots.tsv`` (one row per bin), ``hotspot_genes.tsv`` (one row per bin and gene
    with a trans signal below the loosest grid threshold) and the grid counts
    (thresholds × k)."""
    k = res.kept
    t_loose = min(t.log10p for t in thresholds)
    sel = ((k["cls"] == CROSS) | (k["cls"] == DISTAL)) & passes(k["lp"], t_loose)
    rows, cols, lp = k["row"][sel], k["col"][sel], k["lp"][sel]
    beta, se = k["beta"][sel], k["se"][sel]
    row_bin = bins.attrs["row_bin"]
    b = row_bin[rows] if len(rows) else np.zeros(0, np.int64)
    # genes in cis with the bin: same chromosome, start within the window of the bin's variants
    tchrom = np.array([res.traits.chrom[c] or "" for c in cols], dtype=object)
    tstart = res.traits.start[cols] if len(cols) else np.zeros(0, np.int64)
    bchrom = bins["chrom"].to_numpy(dtype=object)[b] if len(b) else np.zeros(0, dtype=object)
    lo = bins["min_pos"].to_numpy()[b] - window if len(b) else np.zeros(0, np.int64)
    hi = bins["max_pos"].to_numpy()[b] + window if len(b) else np.zeros(0, np.int64)
    cis_bin = (tchrom == bchrom) & (tstart >= lo) & (tstart <= hi)
    cells = pd.DataFrame({"bin": b, "col": cols, "row": rows, "lp": lp, "beta": beta, "se": se})[
        ~cis_bin
    ]
    cells = cells.sort_values(["bin", "col", "lp"], ascending=[True, True, False], kind="stable")
    best = cells.drop_duplicates(["bin", "col"]).sort_values(
        ["bin", "lp"], ascending=[True, False], kind="stable"
    )
    ids = res.meta["ID"].astype(str).to_numpy()
    bmeta = bins.set_index("bin")
    genes = pd.DataFrame(
        {
            "chrom": bmeta.loc[best["bin"], "chrom"].to_numpy() if len(best) else [],
            "bin_start": bmeta.loc[best["bin"], "bin_start"].to_numpy() if len(best) else [],
            "bin_end": bmeta.loc[best["bin"], "bin_end"].to_numpy() if len(best) else [],
            "trait": [res.traits.names[c] for c in best["col"]],
            "trait_chrom": [res.traits.chrom[c] for c in best["col"]],
            "trait_start": [int(res.traits.start[c]) for c in best["col"]],
            "variant": ids[best["row"].to_numpy()] if len(best) else [],
            "beta": best["beta"].to_numpy(),
            "se": best["se"].to_numpy(),
            "log10p": best["lp"].to_numpy(),
        }
    )
    out = bins[["chrom", "bin_start", "bin_end", "n_variants"]].copy()
    blp = best["lp"].to_numpy()
    bb = best["bin"].to_numpy()
    for t in thresholds:
        ok = passes(blp, t.log10p)
        out[f"n_genes_{t.key}"] = np.bincount(bb[ok], minlength=len(bins)).astype(np.int64)[
            : len(bins)
        ]
    top = best.drop_duplicates("bin")  # strongest gene per bin (best is sorted by lp within bin)
    tb = top.set_index("bin")
    has = out.index.isin(tb.index)
    out["top_trait"] = [
        res.traits.names[int(tb.loc[i, "col"])] if h else None
        for i, h in zip(out.index, has, strict=True)
    ]
    out["top_variant"] = [
        ids[int(tb.loc[i, "row"])] if h else None for i, h in zip(out.index, has, strict=True)
    ]
    out["top_log10p"] = [
        float(tb.loc[i, "lp"]) if h else np.nan for i, h in zip(out.index, has, strict=True)
    ]
    grid = [[int((out[f"n_genes_{t.key}"] >= kk).sum()) for kk in ks] for t in thresholds]
    return out, genes, grid


def reference_frame(
    pairs: list[KnownPair],
    res: RectResult,
    variant_list: set[str],
    window: int,
    style: str,
    run: str,
) -> pd.DataFrame:
    """``reference_pairs.tsv``: every reference pair with its class and replication.

    ``status``: ``tested``; ``variant_absent`` (not in the genotypes); ``not_in_rectangle``
    (variant or trait outside the lists); ``not_tested`` (no finite cell, e.g. below
    ``min_maf``)."""
    ids = res.meta["ID"].astype(str).tolist()
    row_of: dict[str, int] = {}
    for i, (vid, cpra) in enumerate(_row_keys(res.meta) if len(res.meta) else []):
        row_of.setdefault(vid, i)
        row_of.setdefault(cpra, i)
    tidx = {t: i for i, t in enumerate(res.traits.names)}
    recs = []
    for kp in pairs:
        v = parse_cpra(kp.variant)
        vchrom = normalize_chrom(v.chrom, style) if v is not None else None  # type: ignore[arg-type]
        j = tidx.get(kp.trait)
        tchrom = res.traits.chrom[j] if j is not None else None
        tstart = int(res.traits.start[j]) if tchrom is not None and j is not None else None
        cls, sub = pair_class(vchrom, v.pos if v is not None else None, tchrom, tstart, window)
        rec: dict[str, Any] = {
            "run": run,
            "trait": kp.trait,
            "variant": kp.variant,
            "ref_trait": kp.ref_trait,
            "pair_class": cls,
            "trans_kind": sub if cls == "trans" else None,
            "chrom": vchrom,
            "pos": v.pos if v is not None else None,
            "trait_chrom": tchrom,
            "trait_start": tstart,
            "distance": (v.pos - tstart)
            if v is not None and tstart is not None and vchrom == tchrom
            else None,
            "beta_ref": kp.beta_ref,
            "se_ref": kp.se_ref,
            "pval_ref": kp.pval_ref,
        }
        vals = None
        if not kp.present:
            status = "variant_absent"
        elif kp.variant not in variant_list or j is None:
            status = "not_in_rectangle"
        else:
            r = row_of.get(kp.variant)
            if r is None and v is not None:
                r = row_of.get(v.cpra)
            cell = res.known_cells.get((r, j)) if r is not None else None
            if cell is None or not math.isfinite(cell[0]):
                status = "not_tested"
            else:
                status = "tested"
                vals = (r, cell)
        rec["status"] = status
        if vals is not None:
            r, (lp, b, s) = vals
            m = res.meta.iloc[r]
            rec.update(
                {
                    "variant_tested": ids[r],
                    "a1freq": float(m["A1FREQ"]),
                    "n": int(m["N"]),
                    "beta": b,
                    "se": s,
                    "z": b / s if s and math.isfinite(s) and s > 0 else np.nan,
                    "log10p": lp,
                }
            )
        recs.append(rec)
    cols = [
        "run",
        "trait",
        "variant",
        "ref_trait",
        "pair_class",
        "trans_kind",
        "chrom",
        "pos",
        "trait_chrom",
        "trait_start",
        "distance",
        "beta_ref",
        "se_ref",
        "pval_ref",
        "status",
        "variant_tested",
        "a1freq",
        "n",
        "beta",
        "se",
        "z",
        "log10p",
    ]
    df = pd.DataFrame(recs, columns=cols)
    for c in ("pos", "trait_start", "distance", "n"):
        df[c] = pd.array(pd.to_numeric(df[c], errors="coerce"), dtype="Int64") if len(df) else df[c]
    tested = (df["status"] == "tested").to_numpy()
    br = df["beta_ref"].to_numpy(dtype=float)
    bb = df["beta"].to_numpy(dtype=float)
    lpa = df["log10p"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        sign = tested & np.isfinite(br) & np.isfinite(bb) & (br != 0) & (np.sign(br) == np.sign(bb))
    sc = pd.array(sign, dtype="boolean")
    sc[~tested] = pd.NA
    df["sign_concordant"] = sc
    df["rep_sign"] = sign
    for name, alpha in REPLICATION_LEVELS.items():
        df[f"rep_{name}"] = sign & passes(lpa, -math.log10(alpha))
    return df


def replication_summary(df: pd.DataFrame) -> dict[str, Any]:
    """Counts per pair class (``cis``, ``trans``, ``no_locus``) of ``reference_pairs.tsv``."""
    out: dict[str, Any] = {}
    for cls in ("cis", "trans", "no_locus"):
        sub = df[df["pair_class"] == cls]
        tested = sub[sub["status"] == "tested"]
        d: dict[str, Any] = {
            "n_pairs": len(sub),
            "n_tested": len(tested),
            "n_rep_sign": int(sub["rep_sign"].astype(bool).sum()),
            **{f"n_rep_{n}": int(sub[f"rep_{n}"].astype(bool).sum()) for n in REPLICATION_LEVELS},
            "median_log10p": float(tested["log10p"].median()) if len(tested) else None,
        }
        if cls == "trans":
            for kind in ("cross_chrom", "same_chrom_distal"):
                d[f"n_pairs_{kind}"] = int((sub["trans_kind"] == kind).sum())
        out[cls] = d
    return out


def sig_frame(res: RectResult, t_nominal: float, max_rows: int) -> tuple[pd.DataFrame, int]:
    """``pairs.sig.tsv.gz``: every cell with p < nominal (any class), the strongest
    ``max_rows``; and the number of such cells."""
    k = res.kept
    sel = np.flatnonzero(passes(k["lp"], t_nominal))
    total = len(sel)
    sel = sel[np.argsort(-k["lp"][sel], kind="stable")][:max_rows]
    rows, cols = k["row"][sel], k["col"][sel]
    m = res.meta.iloc[rows] if len(rows) else res.meta.iloc[0:0]
    pos = m["GENPOS"].to_numpy(dtype=np.int64) if len(rows) else np.zeros(0, np.int64)
    tchrom = [res.traits.chrom[c] for c in cols]
    tstart = [int(res.traits.start[c]) if res.traits.chrom[c] is not None else None for c in cols]
    vchrom = m["CHROM"].astype(str).tolist() if len(rows) else []
    cls_names: Any = (
        np.array(CELL_CLASSES, dtype=object)[k["cls"][sel].astype(np.int64)] if len(sel) else []
    )
    df = pd.DataFrame(
        {
            "variant": m["ID"].astype(str).to_numpy() if len(rows) else [],
            "chrom": vchrom,
            "pos": pos,
            "trait": [res.traits.names[c] for c in cols],
            "trait_chrom": tchrom,
            "trait_start": pd.array(tstart, dtype="Int64"),
            "cell_class": cls_names,
            "reference_pair": k["ref"][sel],
            "distance": pd.array(
                [
                    int(p - s) if tc is not None and tc == vc and s is not None else None
                    for p, s, tc, vc in zip(pos, tstart, tchrom, vchrom, strict=True)
                ],
                dtype="Int64",
            ),
            "a1freq": m["A1FREQ"].to_numpy(dtype=float) if len(rows) else [],
            "n": m["N"].to_numpy() if len(rows) else [],
            "beta": k["beta"][sel],
            "se": k["se"][sel],
            "log10p": k["lp"][sel],
        }
    )
    return df, total


def variants_frame(res: RectResult, row_bin_start: np.ndarray) -> pd.DataFrame:
    """``variants.tsv``: one row per tested variant."""
    meta, pr = res.meta, res.per_row
    df = pd.DataFrame(
        {
            "variant": meta["ID"].astype(str).to_numpy(),
            "chrom": meta["CHROM"].astype(str).to_numpy(),
            "pos": meta["GENPOS"].to_numpy(dtype=np.int64),
            "bin_start": row_bin_start,
            "a1freq": meta["A1FREQ"].to_numpy(dtype=float),
            "n": meta["N"].to_numpy(),
        }
    )
    for c in ("n_cis", "n_cross_chrom", "n_same_chrom_distal", "n_sig_trans"):
        df[c] = np.asarray(pr[c], dtype=np.int64)
    df["max_log10p_trans"] = np.asarray(pr["max_log10p_trans"], dtype=float)
    df["top_trans_trait"] = np.asarray(pr["top_trans_trait"], dtype=object)
    df["max_log10p_cis"] = np.asarray(pr["max_log10p_cis"], dtype=float)
    df["top_cis_trait"] = np.asarray(pr["top_cis_trait"], dtype=object)
    df["excluded_from_calibration"] = res.row_excluded
    return df


def inflation_summary(res: RectResult, exclusion: Exclusion) -> dict[str, Any]:
    """λ_GC and observed vs expected counts of the cross-chromosome cells (after exclusion)."""
    arr = res.cross_desc
    n = len(arr)
    med = median_sorted(arr)
    counts = []
    for a in INFLATION_ALPHAS:
        counts.append(
            {
                "alpha": a,
                "label": f"p < {fmt_alpha(a)}",
                "observed": int(np.count_nonzero(passes(arr, -math.log10(a)))),
                "expected": n * a,
            }
        )
    n_ex_rows = int(res.row_excluded.sum()) if len(res.row_excluded) else 0
    tested_genes = set(res.traits.names)
    return {
        "n_cells": n,
        "n_cells_excluded": res.n_cross_excluded,
        "n_reference_trans_pairs": exclusion.n_pairs,
        "n_variants_excluded": n_ex_rows,
        "n_genes_excluded": len(exclusion.genes & tested_genes),
        "exclusion_window_bp": exclusion.window,
        "median_log10p": med,
        "lambda_gc": lambda_gc(med),
        "counts": counts,
    }


# ================================================================ the step
def read_known(path: Path) -> list[KnownPair]:
    """Reference pairs on dataset trait IDs (``Layout.dsref_pairs``)."""
    t = read_pairs_table(path)
    out = []

    def num(x: str) -> float:
        try:
            return float(x)
        except (TypeError, ValueError):
            return float("nan")

    for r in t.itertuples(index=False):
        out.append(
            KnownPair(
                trait=r.trait,
                variant=r.variant,
                beta_ref=num(r.beta),
                se_ref=num(r.se),
                pval_ref=num(r.pval),
                ref_trait=r.ref_trait,
                present=str(r.present) == "1",
            )
        )
    return out


def read_loci(path: Path) -> pd.DataFrame:
    """``trait chrom start end`` of the dataset (strings; NA when unknown)."""
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={"trait": str, "chrom": str},
        keep_default_na=False,
        na_values={"start": ["NA"], "end": ["NA"]},
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    return df


def rect_chunks(path: Path) -> Iterator[RectChunk]:
    """Stream a rect-assoc file with the schema-checked reader (BETA, SE and LOG10P)."""
    return iter_rect_assoc(path, ("BETA", "SE", "LOG10P"), chunksize=CHUNK_ROWS)


def _round(x: Any, nd: int = 4) -> Any:
    """Round floats for summary.json: ``nd`` decimals, but ``nd`` significant digits below 1
    so that thresholds such as the Bonferroni p (~1e-7) keep their precision."""
    if isinstance(x, dict):
        return {k: _round(v, nd) for k, v in x.items()}
    if isinstance(x, list):
        return [_round(v, nd) for v in x]
    if isinstance(x, (np.floating, float)):
        x = float(x)
        if not math.isfinite(x):
            return None
        if x != 0.0 and abs(x) < 1.0:
            return float(f"{x:.{nd}g}")
        return round(x, nd)
    if isinstance(x, np.integer):
        return int(x)
    return x


@register("summarize_trans")
@input_errors
def summarize_trans(ctx: StepContext) -> None:
    """Replication of the reference pairs by class, hotspots and calibration of one run."""
    t0 = time.time()
    p = ctx.params
    ds, ref, run = p["dataset"], p["reference"], str(p.get("run") or p.get("strategy") or "run")
    levels: dict[str, Any] = dict(p["levels"])
    window = int(p["cis_window_bp"])
    bin_bp = int(p["hotspot_bin_bp"])
    hot_min = int(p["hotspot_min_traits"])
    max_sig = int(p["max_sig_rows"])
    geno = read_json(ctx.path(Layout.geno_json))
    style = geno["chrom_style"]
    chrom_order = _genotype_chroms(geno)
    loci = read_loci(ctx.path(Layout.ds_traits(ds)))
    listed_variants = read_lines(ctx.path(Layout.trans_variants(ds)))
    listed_traits = read_lines(ctx.path(Layout.trans_traits(ds)))
    pairs = read_known(ctx.path(Layout.dsref_pairs(ds, ref))) if p.get("known_pairs") else []
    known: dict[str, list[str]] = {}
    for kp in pairs:
        if kp.present:
            known.setdefault(kp.variant, []).append(kp.trait)
            v = parse_cpra(kp.variant)
            if v is not None and v.cpra != kp.variant:
                known.setdefault(v.cpra, []).append(kp.trait)
    exclusion = build_exclusion(pairs, loci, style, window)
    lists = read_json(ctx.path(Layout.trans_lists_json(ds)))
    sj = ctx.path(Layout.strategy_json("trans", ds, p["strategy"]))
    smeta = read_json(sj) if sj.exists() else {}

    t_nominal = -math.log10(float(levels["nominal"]))
    fixed = [-math.log10(a) for a in HOTSPOT_ALPHAS]
    fixed += [
        -math.log10(float(levels[n])) for n in ("nominal", "strict") if levels.get(n) is not None
    ]
    path = ctx.path(p["rect"])
    try:
        _, trait_names = read_rect_header(path)
        res = summarize_rect(
            rect_chunks(path),
            trait_names=trait_names,
            loci=loci,
            known=known,
            exclusion=exclusion,
            cis_window_bp=window,
            keep_log10p=min(fixed),
            nominal_log10p=t_nominal,
            bonferroni=bool(levels.get("bonferroni", True)),
        )
    except QpgenOutputError as e:
        raise StepError(str(e)) from e

    n_trans = res.n_cells["cross_chrom"] + res.n_cells["same_chrom_distal"]
    thresholds = grid_thresholds(levels, n_trans)
    ks = grid_k(hot_min)
    bins = bins_of(res.meta, bin_bp, chrom_order)
    hot, hot_genes, grid = hotspot_tables(res, bins, thresholds, ks, window)
    kdf = reference_frame(pairs, res, set(listed_variants), window, style, run)
    rep = replication_summary(kdf)
    infl = inflation_summary(res, exclusion)
    sig, sig_total = sig_frame(res, t_nominal, max_sig)
    row_bin_start = (
        bins["bin_start"].to_numpy()[bins.attrs["row_bin"]] if len(bins) else np.zeros(0, np.int64)
    )
    vdf = variants_frame(res, row_bin_start)
    default_key = next(t.key for t in thresholds if "nominal" in t.roles)
    n_hot_default = grid[[t.key for t in thresholds].index(default_key)][ks.index(hot_min)]

    warnings: list[str] = []
    lam = infl["lambda_gc"]
    if math.isfinite(lam) and lam > LAMBDA_WARN:
        warnings.append(
            f"p-values inflated (λGC {lam:.2f} across chromosomes, above {LAMBDA_WARN})"
        )
    if math.isfinite(lam) and lam < 1 / LAMBDA_WARN:
        warnings.append(
            f"p-values deflated (λGC {lam:.2f} across chromosomes, below {1 / LAMBDA_WARN:.2f})"
        )
    if len(res.meta) < len(listed_variants):
        warnings.append(
            f"{len(listed_variants) - len(res.meta)} of {len(listed_variants)} listed variants were "
            "not tested (absent or below min_maf/min_mac)"
        )
    if res.n_cells["cross_chrom"] == 0:
        warnings.append(
            "no cross-chromosome cell (every gene is on its variant's chromosome or has no locus): "
            "no calibration"
        )
    n_total = sum(res.n_cells.values())
    if res.n_cells["no_locus"] and n_total and res.n_cells["no_locus"] / n_total > 0.5:
        warnings.append(f"{res.n_cells['no_locus']} of {n_total} cells have no gene locus")
    if sig_total > max_sig:
        warnings.append(
            f"pairs.sig.tsv.gz capped at max_sig_rows = {max_sig} ({sig_total} cells at p < nominal)"
        )
    ds_json = ctx.path(Layout.ds_json(ds))
    if ds_json.exists():
        warnings += [f"dataset {ds}: {w}" for w in read_json(ds_json).get("warnings") or []]
    tsel = lists.get("traits") or {}
    if tsel.get("n_dropped_no_locus"):
        warnings.append(
            f"{tsel['n_dropped_no_locus']} trait(s) without a locus were not tested "
            f"(trait_chroms {tsel.get('trait_chroms')!r}; use 'all' to keep them)"
        )
    n_samples = (
        pd.to_numeric(res.meta["N"], errors="coerce")
        if len(res.meta)
        else pd.Series([], dtype=float)
    )
    n_samples = n_samples[np.isfinite(n_samples)]

    q = qq_points(res.cross_desc)
    q.insert(0, "n_cells", len(res.cross_desc))
    for name, frame in (
        ("reference_pairs.tsv", kdf),
        ("hotspots.tsv", hot),
        ("hotspot_genes.tsv", hot_genes),
        ("qq.tsv", q),
        ("variants.tsv", vdf),
        ("pairs.sig.tsv.gz", sig),
    ):
        write_tsv(ctx.path(Layout.trans_summary(ds, name)), frame)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "dataset": ds,
        "reference": ref,
        "run": run,
        "strategy": p["strategy"],
        "rint": bool(p.get("rint")),
        "n_cov": smeta.get("n_cov"),
        "cov_source": smeta.get("cov_source"),
        "cov_columns": smeta.get("cov_columns"),
        "n_variants": len(listed_variants),
        "n_traits": len(listed_traits),
        "n_variants_tested": len(res.meta),
        "n_traits_tested": len(res.traits.names),
        "n_samples": int(np.median(n_samples)) if len(n_samples) else None,
        "n_cells": n_total,
        "n_cells_trans": n_trans,
        "cells": res.n_cells,
        "variants": p.get("variants"),
        "traits": p.get("traits"),
        "trait_chroms": p.get("trait_chroms"),
        "variant_selection": lists.get("variants"),
        "trait_selection": lists.get("traits"),
        "known_pairs": bool(p.get("known_pairs")),
        "n_reference_pairs": len(pairs),
        "replication": {"levels": REPLICATION_LEVELS, **rep},
        "hotspots": {
            "bin_bp": bin_bp,
            "n_bins": len(bins),
            "thresholds": [t.as_json() for t in thresholds],
            "k": ks,
            "counts": grid,
            "default": {"threshold": default_key, "k": hot_min, "n_hotspots": n_hot_default},
            "bonferroni_p": BONFERRONI_ALPHA / n_trans if n_trans else None,
        },
        "inflation": infl,
        "n_sig_cells": sig_total,
        "sig_rows_written": len(sig),
        "params": {
            "cis_window_bp": window,
            "hotspot_bin_bp": bin_bp,
            "hotspot_min_traits": hot_min,
            "min_maf": p.get("min_maf"),
            "min_mac": p.get("min_mac"),
            "levels": levels,
            "max_sig_rows": max_sig,
            "lambda_warn": LAMBDA_WARN,
        },
        "definitions": DEFINITIONS,
        "warnings": warnings,
        "memory_bound": (
            "float32 -log10 p of the cross-chromosome calibration cells (4 bytes each, plus a sorted "
            "copy), O(variants + traits) accumulators, and the cells below the loosest hotspot "
            "threshold (or Bonferroni over the trans cells seen so far) with beta and SE"
        ),
        "seconds": round(time.time() - t0, 2),
    }
    write_json(ctx.path(Layout.trans_summary(ds, "summary.json")), _round(summary))
    _log(
        ctx,
        "trans %s: %d variants x %d traits, %d cells (%d cis, %d cross-chromosome, %d same-chromosome "
        "distal); lambda_GC %.3f over %d cross-chromosome cells (%d excluded); %d hotspot(s) at p < %s, k >= %d; "
        "%.1fs",
        run,
        len(res.meta),
        len(res.traits.names),
        n_total,
        res.n_cells["cis"],
        res.n_cells["cross_chrom"],
        res.n_cells["same_chrom_distal"],
        lam,
        infl["n_cells"],
        infl["n_cells_excluded"],
        n_hot_default,
        fmt_alpha(float(levels["nominal"])),
        hot_min,
        time.time() - t0,
    )
