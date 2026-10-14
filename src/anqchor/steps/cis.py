"""``cis.<dataset>.summarize`` (kind ``summarize_cis``): replication of known cis pairs.

Reads the reference pairs mapped to the dataset and qpgen's association results of the
run (one phenotype transform per run: as given, or RINT with ``pipelines.cis.rint``) and
writes under ``cis/<dataset>/``. Rows are labelled by ``run``: the run directory name for
a run with one strategy (the usual case), the strategy name otherwise.

``pairs.long.tsv.gz``
    One row per (run, mapped reference pair): ``run trait ref_trait variant beta_ref
    se_ref z_ref pval_ref log10p_ref status a1freq n beta se tstat z log10p
    sign_concordant rep_<level>…``. ``status`` is ``tested``, ``variant_absent`` (the
    variant is not in the genotypes, from ``prepare.reference``) or ``not_tested`` (the
    variant is present but qpgen wrote no result, or one without a finite estimate,
    e.g. no genotype variation among the phenotyped samples). Booleans are written as
    ``1``/``0``; ``sign_concordant`` is ``NA`` for pairs that were not tested. ``beta``
    and ``beta_ref`` refer to the same allele (variants are matched by
    ``CHROM:POS:REF:ALT``), so the signed ``z`` (``beta / se``, else from the p-value) and
    ``z_ref`` (``beta_ref / se_ref``, else from ``pval_ref``) are aligned. ``log10p_ref`` is ``-log10 pval_ref``
    (parsed from the text for p-values below the double range; from ``z_ref`` when the
    reference has no p-value).
``traits.tsv``
    One row per mapped trait: ``n_pairs``, ``n_variant_absent``, then per run
    ``best_log10p:<run>`` and ``level:<run>`` (highest level reached by any of the
    trait's pairs), ``delta_log10p:<run>`` against the first run, and ``best_run``.
``strategies.tsv``
    One row per run with ``rint``, counts, replicated counts and fractions per level,
    ``median_log10p``, ``mean_z2``, ``mean_z2_common``, ``rel_neff``, ``beta_corr``,
    ``frac_sign_discordant`` and ``pi1``. With the ``rect`` engine also local control
    statistics from the non-reference cells.
``power_curve.tsv``
    Fraction of tested pairs that are sign-concordant with ``p < 10^-t`` for
    ``t = 0, 0.5, …``, per run.
``summary.json``
    Counts, levels, ``rint``, per-run key numbers, rankings, the recommended run and
    warnings.

Replication levels (DESIGN decision 11). Sign concordance is required by every level:

* ``sign``: ``sign(beta_ref) == sign(beta)`` (always reported).
* ``nominal`` / ``strict``: sign and ``p < alpha``.
* ``bonferroni``: sign and ``p < 0.05 / n_tested`` (per run).
* ``fdr``: sign and Benjamini-Hochberg ``q <= alpha`` over the tested pairs of the run.

``pi1`` is Storey's ``1 - pi0`` with the single ``lambda = 0.5`` estimator
(``pi0 = #{p > 0.5} / (0.5 m)``): approximate and only indicative for a few hundred pairs.

``rel_neff`` approximates the effective sample size relative to the first (baseline) run:
with ``z = beta / se``, ``E[z^2] = 1 + N beta^2 Var(g) / sigma_e^2`` (the noncentrality), so
removing residual variance ``sigma_e^2`` acts like adding samples. Over the pairs tested in
every run, ``(mean z^2 - 1) / (baseline mean z^2 - 1)`` estimates ``N_eff / N_eff(baseline)``.
It is a ratio of means, so pairs weigh in proportion to their ``z^2`` (the strongest pairs
dominate). The reference pairs are the strongest known signals (selected, partly
in-sample), so it is a diagnostic, not an estimate of the sample size.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import special, stats

from anqchor.pipelines.layout import Layout
from anqchor.qpgen.outputs import QpgenOutputError, iter_rect_assoc, read_pair_assoc
from anqchor.stats import p_level_threshold, passes
from anqchor.steps import StepContext, StepError, register
from anqchor.util.files import atomic_write, read_json, write_json

SCHEMA_VERSION = 3
"""3: ``z = beta / se`` (2 used the normal quantile of the p-value) and ``log10p_ref`` kept
for p-values below the double range or derived from ``z_ref`` when the reference has none."""
LEVEL_ORDER = ("sign", "nominal", "fdr", "bonferroni", "strict")
"""Default order of the levels, weakest to strongest. The order actually used is by
threshold (:func:`stringency_order`): Bonferroni (0.05 / n_tested) is stronger than
``strict`` (1e-5) once more than 5000 pairs are tested."""
RECOMMEND_ORDER = ("strict", "bonferroni", "fdr", "nominal", "sign")
"""Default ranking keys for the recommended strategy (strongest first; see
:func:`stringency_order`)."""
BONFERRONI_ALPHA = 0.05
LONG_COLUMNS = (
    "run",
    "trait",
    "ref_trait",
    "variant",
    "beta_ref",
    "se_ref",
    "z_ref",
    "pval_ref",
    "log10p_ref",
    "status",
    "a1freq",
    "n",
    "beta",
    "se",
    "tstat",
    "z",
    "log10p",
    "sign_concordant",
)
STAT_COLUMNS = ("a1freq", "n", "beta", "se", "tstat", "z", "log10p")
CONTROL_P = 1e-5
"""p-value cutoff counted among the local control cells (rect engine)."""
_HIST_MAX = 20.0
_HIST_BINS = 200_000
_CHI2_MEDIAN = float(stats.chi2.ppf(0.5, 1))


def z_from_log10p(log10p: np.ndarray, sign: np.ndarray) -> np.ndarray:
    """Signed normal quantile of a two-sided p-value given as ``-log10 p``:
    ``sign * Phi^-1(1 - p/2)``, computed in log space so that p-values below the double
    range keep their size. NaN where ``-log10 p`` is not finite."""
    L = np.asarray(log10p, dtype=float)
    out = np.full(L.shape, np.nan)
    ok = np.isfinite(L) & (L >= 0)
    if ok.any():
        log_half_p = -L[ok] * math.log(10.0) - math.log(2.0)
        out[ok] = np.maximum(-special.ndtri_exp(log_half_p), 0.0)
    result: np.ndarray = np.asarray(sign, dtype=float) * out
    return result


def neglog10_text(values: pd.Series) -> np.ndarray:
    """``-log10 p`` of p-values given as text, exact also below the double range: a value
    that parses to 0 (``1e-400``) is read as mantissa and exponent. NaN for missing,
    unparsable, zero or out-of-range values."""
    txt = values.astype(str).str.strip()
    pv = pd.to_numeric(txt, errors="coerce").to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where((pv > 0) & (pv <= 1), -np.log10(pv) + 0.0, np.nan)
    under = np.flatnonzero(pv == 0)
    if under.size:
        parts = txt.iloc[under].str.extract(r"^\+?([0-9]*\.?[0-9]+)[eE]([-+]?[0-9]+)$")
        m = pd.to_numeric(parts[0], errors="coerce").to_numpy(dtype=float)
        e = pd.to_numeric(parts[1], errors="coerce").to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = -(np.log10(m) + e)
        out[under] = np.where(np.isfinite(v) & (m > 0) & (v >= 0), v, np.nan)
    result: np.ndarray = out
    return result


def log10p_from_z(z: np.ndarray) -> np.ndarray:
    """``-log10`` of the two-sided normal p-value of ``z`` (log space: no underflow)."""
    a = np.abs(np.asarray(z, dtype=float))
    out: np.ndarray = np.where(
        np.isfinite(a), -(math.log(2.0) + special.log_ndtr(-a)) / math.log(10.0) + 0.0, np.nan
    )
    return out


def wald_z(beta: np.ndarray, se: np.ndarray, log10p: np.ndarray) -> np.ndarray:
    """Signed z-score ``beta / se`` (the Wald statistic, as ``z_ref``); where ``se`` is
    missing or not positive, the signed normal quantile of the p-value
    (:func:`z_from_log10p`).

    ``beta / se`` rather than the normal quantile of the p-value: the p-value comes from a
    t distribution, so for strong pairs its normal quantile is smaller than ``beta / se``
    (|t| = 59 with 660 degrees of freedom maps to |z| = 35), which would compress ``z^2``
    and bias ``rel_neff`` towards 1. ``(beta / se)^2`` scales with the inverse residual
    variance, i.e. with the effective sample size."""
    b = np.asarray(beta, dtype=float)
    s = np.asarray(se, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(np.isfinite(s) & (s > 0), b / s, np.nan)
    fallback = z_from_log10p(log10p, np.sign(b))
    out: np.ndarray = np.where(np.isfinite(z), z, fallback)
    return out


# ------------------------------------------------------------------------- levels
def active_levels(levels: dict[str, Any]) -> list[str]:
    """Configured levels in :data:`LEVEL_ORDER` (``sign`` always first)."""
    out = ["sign"]
    for name in LEVEL_ORDER[1:]:
        v = levels.get(name)
        if v is None or v is False:
            continue
        out.append(name)
    return out


def level_description(levels: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Machine-readable description of each active level (for ``summary.json``)."""
    d: dict[str, dict[str, Any]] = {"sign": {"rule": "sign(beta_ref) == sign(beta)"}}
    for name in active_levels(levels)[1:]:
        if name in ("nominal", "strict"):
            d[name] = {"rule": "sign and p < alpha", "alpha": float(levels[name])}
        elif name == "bonferroni":
            d[name] = {
                "rule": "sign and p < alpha / n_tested",
                "alpha": BONFERRONI_ALPHA,
                "per": "tested pairs of each run",
            }
        elif name == "fdr":
            d[name] = {
                "rule": "sign and Benjamini-Hochberg q <= alpha",
                "alpha": float(levels[name]),
                "per": "tested pairs of each run",
            }
    return d


def level_text(name: str, levels: dict[str, Any], n_tested: int | None = None) -> str:
    """Short human description, e.g. ``p < 1e-05``."""
    if name == "sign":
        return "sign concordance"
    if name in ("nominal", "strict"):
        return f"p < {float(levels[name]):g}"
    if name == "bonferroni":
        if n_tested:
            return f"p < {BONFERRONI_ALPHA:g}/{n_tested} = {BONFERRONI_ALPHA / n_tested:.2g}"
        return f"p < {BONFERRONI_ALPHA:g}/n_tested"
    return f"BH q ≤ {float(levels[name]):g}"


def bh_neglog10q(log10p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg ``-log10 q`` from ``-log10 p`` (computed in log space, so
    p-values far below the double range keep their order). NaN p-values count as
    ``p = 1`` in the number of tests and get ``q = 1``."""
    L = np.asarray(log10p, dtype=float)
    m = L.size
    if m == 0:
        return np.empty(0)
    L0 = np.where(np.isfinite(L), L, np.where(np.isnan(L), 0.0, L))
    order = np.argsort(-L0, kind="mergesort")  # most significant first
    ranks = np.arange(1, m + 1, dtype=float)
    nlq_sorted = L0[order] - np.log10(m / ranks)
    # q_(i) = min_{j >= i} p_(j) m / j  <=>  -log10 q_(i) = max_{j >= i} (...)
    nlq_sorted = np.maximum.accumulate(nlq_sorted[::-1])[::-1]
    nlq_sorted = np.maximum(nlq_sorted, 0.0)
    out = np.empty(m)
    out[order] = nlq_sorted
    return out


def storey_pi1(log10p: np.ndarray, lam: float = 0.5) -> float:
    """Storey's ``pi1 = 1 - pi0`` with one ``lambda`` (``nan`` without p-values)."""
    L = np.asarray(log10p, dtype=float)
    L = L[np.isfinite(L)]
    if L.size == 0:
        return float("nan")
    p = np.power(10.0, -L)
    pi0 = float(np.mean(p > lam)) / (1.0 - lam)
    return float(1.0 - min(1.0, pi0))


def lambda_gc_from_median_log10p(median_log10p: float) -> float:
    """Genomic control λ from the median ``-log10 p`` (1-df chi-square)."""
    if not math.isfinite(median_log10p):
        return float("nan")
    p = 10.0 ** (-median_log10p)
    return float(stats.chi2.isf(p, 1) / _CHI2_MEDIAN)


def add_levels(df: pd.DataFrame, levels: dict[str, Any]) -> pd.DataFrame:
    """Add ``sign_concordant`` and ``rep_<level>`` columns to one run's rows.

    ``df`` needs ``status``, ``beta_ref``, ``beta`` and ``log10p``. Bonferroni and FDR are
    computed over the rows with ``status == 'tested'``."""
    df = df.copy()
    tested = (df["status"] == "tested").to_numpy()
    br = df["beta_ref"].to_numpy(dtype=float)
    b = df["beta"].to_numpy(dtype=float)
    L = df["log10p"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        sign = tested & np.isfinite(br) & np.isfinite(b) & (br != 0) & (np.sign(br) == np.sign(b))
    sc = pd.array(sign, dtype="boolean")
    sc[~tested] = pd.NA
    df["sign_concordant"] = sc
    n_tested = int(tested.sum())
    Lt = np.where(np.isfinite(L), L, -np.inf)
    for name in active_levels(levels):
        if name == "sign":
            rep = sign
        elif name in ("nominal", "strict"):
            rep = sign & passes(Lt, p_level_threshold(levels[name]))
        elif name == "bonferroni":
            rep = (
                sign & passes(Lt, p_level_threshold(BONFERRONI_ALPHA / n_tested))
                if n_tested
                else sign & False
            )
        else:  # fdr
            nlq = np.full(len(df), 0.0)
            if n_tested:
                nlq[tested] = bh_neglog10q(L[tested])
            rep = sign & passes(nlq, p_level_threshold(levels[name]), "fdr")
        df[f"rep_{name}"] = rep
    return df


def level_cutoffs(
    df: pd.DataFrame, levels: Sequence[str], level_cfg: dict[str, Any]
) -> dict[str, float]:
    """``-log10 p`` cut-off of each active level for one run's rows.

    ``sign`` has none (``-inf``); ``bonferroni`` depends on the number of tested pairs;
    ``fdr`` is the Benjamini-Hochberg cut-off of the tested pairs."""
    tested = df["status"] == "tested"
    n_tested = int(tested.sum())
    out: dict[str, float] = {}
    for name in levels:
        if name == "sign":
            out[name] = -math.inf
        elif name in ("nominal", "strict"):
            out[name] = p_level_threshold(level_cfg[name])
        elif name == "bonferroni":
            out[name] = p_level_threshold(BONFERRONI_ALPHA / n_tested) if n_tested else math.inf
        elif name == "fdr":
            # BH rejects the k smallest p-values, p <= k q / m: the cut-off is -log10(k q / m)
            # (-log10(q / m) when nothing is rejected, the cut-off the best pair had to beat).
            q = float(level_cfg[name])
            if n_tested:
                L = df["log10p"].to_numpy(dtype=float)[tested.to_numpy()]
                k = int(np.count_nonzero(passes(bh_neglog10q(L), p_level_threshold(q), "fdr")))
                out[name] = p_level_threshold(max(k, 1) * q / n_tested)
            else:
                out[name] = math.inf
    return out


def stringency_order(cutoffs: dict[str, float]) -> list[str]:
    """Levels from weakest to strongest by their ``-log10 p`` cut-off (ties: :data:`LEVEL_ORDER`)."""
    return sorted(cutoffs, key=lambda n: (cutoffs[n], LEVEL_ORDER.index(n)))


def pair_level(
    df: pd.DataFrame, levels: Sequence[str], order: Sequence[str] | None = None
) -> pd.Series:
    """Highest level reached by each row (``none`` when tested but not sign-concordant,
    ``NA`` when not tested). ``order`` lists the levels weakest to strongest (default
    :data:`LEVEL_ORDER`; see :func:`stringency_order`)."""
    out = pd.Series(np.where(df["status"] == "tested", "none", None), index=df.index, dtype=object)
    for name in order if order is not None else LEVEL_ORDER:
        if name in levels:
            out[df[f"rep_{name}"].to_numpy(dtype=bool)] = name
    return out


# ---------------------------------------------------------------------- inputs
def read_mapped_pairs(path: Path) -> pd.DataFrame:
    """``prepare.mapping`` pairs → ``trait ref_trait variant beta_ref se_ref z_ref pval_ref
    log10p_ref present``."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    need = ["trait", "variant", "beta", "se", "ref_trait", "present"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise StepError(f"{path}: missing columns {missing} (re-run the prepare.mapping step)")
    if "pval" not in df.columns:
        df["pval"] = "NA"
    out = pd.DataFrame(
        {
            "trait": df["trait"],
            "ref_trait": df["ref_trait"],
            "variant": df["variant"],
            "beta_ref": pd.to_numeric(df["beta"], errors="coerce"),
            "se_ref": pd.to_numeric(df["se"], errors="coerce"),
            "pval_ref": pd.to_numeric(df["pval"].replace({"": "NA"}), errors="coerce"),
            "present": df["present"].astype(str).str.strip() == "1",
        }
    )
    br, se = (out[c].to_numpy(dtype=float) for c in ("beta_ref", "se_ref"))
    log10p_ref = neglog10_text(df["pval"])
    with np.errstate(divide="ignore", invalid="ignore"):
        z_se = np.where(np.isfinite(se) & (se > 0), br / se, np.nan)
    z_ref = np.where(np.isfinite(z_se), z_se, z_from_log10p(log10p_ref, np.sign(br)))
    # no usable reference p-value: the two-sided p of z_ref (so the strongest pairs stay)
    log10p_ref = np.where(np.isfinite(log10p_ref), log10p_ref, log10p_from_z(z_ref))
    out.insert(out.columns.get_loc("pval_ref"), "z_ref", z_ref)
    out.insert(out.columns.get_loc("pval_ref") + 1, "log10p_ref", log10p_ref)
    return out.reset_index(drop=True)


def _empty_stats() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "trait": pd.Series(dtype=str),
            "variant": pd.Series(dtype=str),
            **{c: pd.Series(dtype=float) for c in STAT_COLUMNS},
        }
    )


def pair_assoc_stats(df: pd.DataFrame, rint: bool, wanted: set[str] | None = None) -> pd.DataFrame:
    """Per-(trait, variant) statistics from a ``read_pair_assoc`` frame: the raw columns,
    or the ``*_RINT`` columns with ``rint``.

    The variant key is qpgen's ``ID``; when ``wanted`` is given and an ``ID`` is not in
    it, ``CHROM:GENPOS:ALLELE0:ALLELE1`` is tried instead."""
    suffix = "_RINT" if rint else ""
    cols = [f"{s}{suffix}" for s in ("BETA", "SE", "TSTAT", "LOG10P")]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise StepError(
            f"pair-assoc output has no {'RINT' if rint else 'raw'} columns ({', '.join(missing)}); "
            "was it run with --skip-rint? Re-run the association step"
        )
    if df.empty:
        return _empty_stats()
    variant = df["ID"].astype(str)
    if wanted is not None:
        alt = (
            df["CHROM"].astype(str)
            + ":"
            + df["GENPOS"].astype(str)
            + ":"
            + df["ALLELE0"].astype(str)
            + ":"
            + df["ALLELE1"].astype(str)
        )
        variant = variant.where(variant.isin(wanted), alt)
    out = pd.DataFrame(
        {
            "trait": df["TRAIT"].astype(str),
            "variant": variant,
            "a1freq": df["A1FREQ"].astype(float),
            "n": df["N"].astype(float),
            "beta": df[cols[0]].astype(float),
            "se": df[cols[1]].astype(float),
            "tstat": df[cols[2]].astype(float),
            "log10p": df[cols[3]].astype(float),
        }
    )
    return out.drop_duplicates(["trait", "variant"], keep="first").reset_index(drop=True)


@dataclass
class ControlStats:
    """Streaming summary of ``-log10 p`` over local control cells (histogram median)."""

    n: int = 0
    n_sig: int = 0
    hist: np.ndarray = field(default_factory=lambda: np.zeros(_HIST_BINS, dtype=np.int64))

    def add(self, log10p: np.ndarray) -> None:
        L = log10p[np.isfinite(log10p)]
        if L.size == 0:
            return
        self.n += int(L.size)
        self.n_sig += int((-math.log10(CONTROL_P) < L).sum())
        idx = np.clip((L / _HIST_MAX * _HIST_BINS).astype(np.int64), 0, _HIST_BINS - 1)
        self.hist += np.bincount(idx, minlength=_HIST_BINS)

    def median_log10p(self) -> float:
        if self.n == 0:
            return float("nan")
        cum = np.cumsum(self.hist)
        k = int(np.searchsorted(cum, self.n / 2.0))
        return (k + 0.5) * _HIST_MAX / _HIST_BINS

    def lambda_gc(self) -> float:
        return lambda_gc_from_median_log10p(self.median_log10p())


@dataclass
class RectResult:
    """Known-cell statistics and control summaries from one ``rect-assoc`` file."""

    known: pd.DataFrame
    controls: ControlStats
    controls_cross: ControlStats
    n_variants: int
    n_traits: int


def rect_stats(
    path: Path, pairs: pd.DataFrame, trait_chrom: dict[str, str] | None = None, chunksize: int = 64
) -> RectResult:
    """Extract the reference cells (``pairs`` ``trait``/``variant``) from a ``rect-assoc``
    file and summarise every other cell as a local control.

    ``trait_chrom`` (trait → chromosome) additionally separates cross-chromosome
    controls (variant and trait on different chromosomes)."""
    from anqchor.qpgen.outputs import read_rect_header

    _, traits = read_rect_header(path)
    tidx = {t: i for i, t in enumerate(traits)}
    by_variant: dict[str, list[int]] = {}
    want_trait: dict[str, list[str]] = {}
    for tr, v in zip(pairs["trait"], pairs["variant"], strict=True):
        if tr in tidx:
            by_variant.setdefault(v, []).append(tidx[tr])
            want_trait.setdefault(v, []).append(tr)
    tchrom = np.array([(trait_chrom or {}).get(t, "") for t in traits], dtype=object)
    has_chrom = tchrom != ""
    controls, cross = ControlStats(), ControlStats()
    rows: list[dict[str, Any]] = []
    n_var = 0
    for chunk in iter_rect_assoc(path, ("BETA", "SE", "TSTAT", "LOG10P"), chunksize=chunksize):
        L = chunk.stats["LOG10P"]
        mask = np.ones_like(L, dtype=bool)
        ids = chunk.meta["ID"].astype(str).tolist()
        chroms = chunk.meta["CHROM"].astype(str).tolist()
        n_var += len(ids)
        for r, vid in enumerate(ids):
            cols = by_variant.get(vid)
            if not cols:
                continue
            mask[r, cols] = False
            for tr, c in zip(want_trait[vid], cols, strict=True):
                rows.append(
                    {
                        "trait": tr,
                        "variant": vid,
                        "a1freq": float(chunk.meta.at[r, "A1FREQ"]),
                        "n": float(chunk.meta.at[r, "N"]),
                        "beta": chunk.stats["BETA"][r, c],
                        "se": chunk.stats["SE"][r, c],
                        "tstat": chunk.stats["TSTAT"][r, c],
                        "log10p": L[r, c],
                    }
                )
        controls.add(L[mask])
        if trait_chrom:
            cross_mask = (
                mask
                & has_chrom[None, :]
                & (tchrom[None, :] != np.array(chroms, dtype=object)[:, None])
            )
            cross.add(L[cross_mask])
    known = pd.DataFrame(rows) if rows else _empty_stats()
    known = known.drop_duplicates(["trait", "variant"], keep="first").reset_index(drop=True)
    return RectResult(
        known=known, controls=controls, controls_cross=cross, n_variants=n_var, n_traits=len(traits)
    )


def read_trait_chrom(path: Path) -> dict[str, str]:
    """trait → chromosome from ``prepare/datasets/<ds>/traits.tsv`` (traits without a
    locus are left out)."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "trait" not in df.columns or "chrom" not in df.columns:
        return {}
    df = df[(df["chrom"] != "") & (df["chrom"] != "NA")]
    return dict(zip(df["trait"], df["chrom"], strict=True))


# -------------------------------------------------------------------- summaries
def long_rows(
    pairs: pd.DataFrame, st: pd.DataFrame, run: str, levels: dict[str, Any]
) -> pd.DataFrame:
    """Rows of ``pairs.long`` for one run (with status, ``z`` and levels)."""
    df = pairs.merge(st, on=["trait", "variant"], how="left", indicator=True)
    # A row qpgen wrote without a finite estimate (no genotype variation among the
    # phenotyped samples, failed fit) is 'not_tested', as in the trans pipeline: it
    # cannot replicate, so it must not count in n_tested, Bonferroni or BH.
    has_estimate = (
        (df["_merge"] == "both")
        & np.isfinite(df["log10p"].to_numpy(dtype=float))
        & np.isfinite(df["beta"].to_numpy(dtype=float))
    )
    status = np.where(
        ~df["present"], "variant_absent", np.where(has_estimate, "tested", "not_tested")
    )
    df = df.drop(columns=["_merge", "present"])
    df.insert(0, "run", run)
    df["status"] = status
    df["z"] = wald_z(
        df["beta"].to_numpy(dtype=float),
        df["se"].to_numpy(dtype=float),
        df["log10p"].to_numpy(dtype=float),
    )
    for c in STAT_COLUMNS:
        df.loc[df["status"] != "tested", c] = np.nan
    df = df[[*LONG_COLUMNS[:-1]]]
    return add_levels(df, levels)


def common_keys(long: pd.DataFrame, labels: Sequence[str]) -> set[tuple[str, str]]:
    """``(ref_trait, variant)`` tested in every run of ``labels``."""
    sets = []
    for lab in labels:
        sub = long[(long["run"] == lab) & (long["status"] == "tested")]
        sets.append(set(zip(sub["ref_trait"], sub["variant"], strict=True)))
    return set.intersection(*sets) if sets else set()


def neff_columns(long: pd.DataFrame, labels: Sequence[str]) -> dict[str, dict[str, float]]:
    """Per run: ``mean_z2_common`` (mean z^2 over the pairs tested in every run) and
    ``rel_neff`` = ``(mean_z2_common - 1) / (baseline mean_z2_common - 1)``, the baseline
    being the first run (see the module doc)."""
    keys = common_keys(long, labels)
    means: dict[str, float] = {}
    for lab in labels:
        sub = long[(long["run"] == lab) & (long["status"] == "tested")]
        # one row per (ref_trait, variant), as the report's per-pair tables
        sub = sub.drop_duplicates(["ref_trait", "variant"])
        mask = np.array(
            [k in keys for k in zip(sub["ref_trait"], sub["variant"], strict=True)], dtype=bool
        )
        z = sub["z"].to_numpy(dtype=float)[mask] if len(sub) else np.empty(0)
        z = z[np.isfinite(z)]
        means[lab] = float(np.mean(z * z)) if z.size else float("nan")
    base = means[labels[0]] if labels else float("nan")
    out = {}
    for lab in labels:
        m = means[lab]
        rel = (
            (m - 1.0) / (base - 1.0)
            if math.isfinite(m) and math.isfinite(base) and base > 1.0
            else float("nan")
        )
        out[lab] = {"mean_z2_common": m, "rel_neff": rel}
    return out


def strategy_row(
    df: pd.DataFrame, run: str, rint: bool, levels: Sequence[str], n_cov: int | None
) -> dict[str, Any]:
    """One ``strategies.tsv`` row from a run's ``pairs.long`` rows."""
    tested = df[df["status"] == "tested"]
    n_t = len(tested)
    row: dict[str, Any] = {
        "run": run,
        "label": run,
        "rint": bool(rint),
        "n_cov": n_cov,
        "n_pairs": len(df),
        "n_tested": n_t,
        "n_variant_absent": int((df["status"] == "variant_absent").sum()),
        "n_not_tested": int((df["status"] == "not_tested").sum()),
    }
    for name in levels:
        row[f"n_rep_{name}"] = int(df[f"rep_{name}"].sum())
    for name in levels:
        row[f"frac_rep_{name}"] = row[f"n_rep_{name}"] / n_t if n_t else float("nan")
    L = tested["log10p"].to_numpy(dtype=float)
    row["median_log10p"] = float(np.nanmedian(L)) if np.isfinite(L).any() else float("nan")
    z = tested["z"].to_numpy(dtype=float)
    z = z[np.isfinite(z)]
    row["mean_z2"] = float(np.mean(z * z)) if z.size else float("nan")
    ok = tested[np.isfinite(tested["beta_ref"]) & np.isfinite(tested["beta"])]
    if len(ok) >= 3 and ok["beta"].std() > 0 and ok["beta_ref"].std() > 0:
        row["beta_corr"] = float(np.corrcoef(ok["beta_ref"], ok["beta"])[0, 1])
    else:
        row["beta_corr"] = float("nan")
    row["frac_sign_discordant"] = (1.0 - row["n_rep_sign"] / n_t) if n_t else float("nan")
    row["pi1"] = storey_pi1(L)
    row["bonferroni_log10p"] = (-math.log10(BONFERRONI_ALPHA / n_t)) if n_t else float("nan")
    return row


def power_curve(
    long: pd.DataFrame, labels: Sequence[str], step: float = 0.5, max_points: int = 2000
) -> pd.DataFrame:
    """Fraction of tested pairs that are sign-concordant with ``-log10 p > t`` (``p < 10^-t``,
    as the levels), per run."""
    tested = long[long["status"] == "tested"]
    top = (
        float(np.nanmax(tested["log10p"]))
        if len(tested) and np.isfinite(tested["log10p"]).any()
        else 0.0
    )
    n_steps = min(max_points, math.ceil(top / step) + 1)
    grid = np.arange(n_steps, dtype=float) * step
    out = []
    for lab in labels:
        sub = tested[tested["run"] == lab]
        n_t = len(sub)
        L = sub["log10p"].to_numpy(dtype=float)
        sign = sub["sign_concordant"].to_numpy(dtype=bool)
        Ls = np.sort(np.where(sign & np.isfinite(L), L, -np.inf))
        # count of values > threshold (p < 10^-t, as the levels' p < alpha)
        counts = n_t - np.searchsorted(Ls, grid, side="right")
        for thr, c in zip(grid, counts, strict=True):
            out.append(
                {
                    "run": lab,
                    "threshold_log10p": float(thr),
                    "n_tested": n_t,
                    "n_replicated": int(c),
                    "frac_replicated": (c / n_t) if n_t else float("nan"),
                }
            )
    return pd.DataFrame(
        out, columns=["run", "threshold_log10p", "n_tested", "n_replicated", "frac_replicated"]
    )


def trait_table(
    long: pd.DataFrame,
    labels: Sequence[str],
    levels: Sequence[str],
    level_cfg: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """``traits.tsv``: per-trait best ``-log10 p`` and highest level per run.

    With ``level_cfg`` the levels are ranked by their cut-offs in each run
    (:func:`stringency_order`), otherwise in :data:`LEVEL_ORDER`. ``delta_log10p:<run>``
    is against the first run."""
    first = long[long["run"] == labels[0]]
    base = first.groupby("trait", sort=False).agg(
        ref_trait=("ref_trait", lambda x: ",".join(dict.fromkeys(x))),
        n_pairs=("variant", "size"),
        n_variant_absent=("status", lambda x: int((x == "variant_absent").sum())),
    )
    best_cols: dict[str, str] = {}
    for lab in labels:
        sub = long[long["run"] == lab].copy()
        order = (
            stringency_order(level_cutoffs(sub, levels, level_cfg))
            if level_cfg is not None
            else [lv for lv in LEVEL_ORDER if lv in levels]
        )
        names = ["none", *order]
        rank = {name: i for i, name in enumerate(names)}
        sub["_level"] = pair_level(sub, levels, order)
        sub["_rank"] = sub["_level"].map(
            lambda v, rank=rank: rank.get(v, -1) if isinstance(v, str) else -1
        )
        g = sub.groupby("trait", sort=False)
        best = g["log10p"].max()
        lv = g["_rank"].max().map(lambda r, names=names: names[r] if r >= 0 else None)
        base[f"best_log10p:{lab}"] = best
        base[f"level:{lab}"] = lv
        best_cols[lab] = f"best_log10p:{lab}"
    for lab in labels[1:]:
        base[f"delta_log10p:{lab}"] = base[best_cols[lab]] - base[best_cols[labels[0]]]
    bests = base[list(best_cols.values())]
    has_any = bests.notna().any(axis=1)
    best_lab = pd.Series(None, index=base.index, dtype=object)
    if has_any.any():
        best_lab[has_any] = (
            bests[has_any].idxmax(axis=1).str.replace("best_log10p:", "", regex=False)
        )
    base["best_run"] = best_lab
    return base.reset_index()


def recommend(
    rows: Sequence[dict[str, Any]],
    levels: Sequence[str],
    level_cfg: dict[str, Any],
    order: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Recommended run: most pairs replicated at the most stringent configured level,
    ties broken by the next levels, then by order. ``order`` lists the levels weakest to
    strongest (:func:`stringency_order`); by default ``strict``, ``bonferroni``, ``fdr``,
    ``nominal``, ``sign`` decide in turn. Returns ``{run, label, reason, key_level}``."""
    keys = (
        [lv for lv in reversed(order) if lv in levels]
        if order is not None
        else [lv for lv in RECOMMEND_ORDER if lv in levels]
    )
    candidates = [r for r in rows if r["n_tested"] > 0]
    if not candidates:
        return {
            "run": None,
            "label": None,
            "key_level": None,
            "reason": "no run tested any reference pair",
        }
    ordered = sorted(candidates, key=lambda r: tuple(-r[f"n_rep_{k}"] for k in keys))
    win = ordered[0]
    key0 = keys[0]
    desc = f"{key0} ({level_text(key0, level_cfg, win['n_tested'])})"
    head = (
        f"{win['label']} replicates the most reference pairs at the {desc} level: "
        f"{win[f'n_rep_{key0}']} of {win['n_tested']} tested"
    )
    if len(ordered) == 1:
        reason = head + " (the only run evaluated)."
    else:
        run = ordered[1]
        decided = next((k for k in keys if win[f"n_rep_{k}"] != run[f"n_rep_{k}"]), None)
        if decided == key0:
            reason = (
                head
                + f", {win[f'n_rep_{key0}'] - run[f'n_rep_{key0}']} more than {run['label']} ({run[f'n_rep_{key0}']})."
            )
        elif decided is not None:
            reason = (
                head + f", tied with {run['label']}; it replicates more at {decided} "
                f"({win[f'n_rep_{decided}']} vs {run[f'n_rep_{decided}']})."
            )
        else:
            reason = head + f", tied with {run['label']} at every level; it comes first."
    return {"run": win["run"], "label": win["label"], "key_level": key0, "reason": reason}


def write_tsv(path: Path, df: pd.DataFrame) -> None:
    """Atomic TSV (gzip when the name ends in ``.gz``); booleans as 1/0, NaN as ``NA``."""
    out = df.copy()
    for c in out.columns:
        if out[c].dtype == bool:
            out[c] = out[c].astype(int)
        elif str(out[c].dtype) == "boolean":
            out[c] = out[c].astype("Int64")
    with atomic_write(path) as fh:
        out.to_csv(fh, sep="\t", index=False, na_rep="NA", float_format="%.6g")


def _num(v: Any) -> Any:
    if isinstance(v, (float, np.floating)):
        return None if not math.isfinite(float(v)) else round(float(v), 6)
    if isinstance(v, np.integer):
        return int(v)
    return v


# ------------------------------------------------------------------------- step
@register("summarize_cis")
def summarize_cis(ctx: StepContext) -> None:
    """Replication of the mapped reference pairs by run (see module doc)."""
    p = ctx.params
    ds, ref, engine = p["dataset"], p["reference"], p["engine"]
    rint = bool(p.get("rint", False))
    level_cfg: dict[str, Any] = dict(p.get("levels") or {})
    levels = active_levels(level_cfg)
    pairs = read_mapped_pairs(ctx.path(Layout.dsref_pairs(ds, ref)))
    mapping = read_json(ctx.path(Layout.dsref_json(ds, ref)))
    if pairs.empty:
        raise StepError(f"no reference pair of {ref} maps onto dataset {ds}; nothing to evaluate")
    dup = pairs.duplicated(["trait", "variant"])
    if dup.any():
        pairs = pairs[~dup].reset_index(drop=True)
    wanted = set(pairs["variant"])
    variants = pd.read_csv(
        ctx.path(Layout.ref_variants(ref)), sep="\t", dtype=str, keep_default_na=False
    )
    indel = (
        dict(zip(variants["variant"], variants["is_indel"] == "1", strict=True))
        if "is_indel" in variants
        else {}
    )
    trait_chrom = read_trait_chrom(ctx.path(Layout.ds_traits(ds))) if engine == "rect" else None

    warnings: list[str] = []
    labels: list[str] = []
    longs: list[pd.DataFrame] = []
    rows: list[dict[str, Any]] = []
    run_info: list[dict[str, Any]] = []
    recomputed: list[str] = []
    for strat in p["strategies"]:
        s = strat["name"]
        lab = strat.get("label") or s
        sj_path = ctx.path(strat.get("strategy_json") or Layout.strategy_json("cis", ds, s))
        sj = read_json(sj_path) if sj_path.exists() else {}
        n_cov = sj.get("n_cov")
        f = ctx.path(strat["file"])
        controls: tuple[ControlStats, ControlStats] | None = None
        try:
            if engine == "pair":
                pair_df = read_pair_assoc(f, n_cov=int(n_cov or 0))
                if pair_df.attrs.get("recomputed_log10p") and not rint:
                    recomputed.append(lab)
                st = pair_assoc_stats(pair_df, rint, wanted)
            else:
                rr = rect_stats(f, pairs[pairs["present"]], trait_chrom)
                st = rr.known
                controls = (rr.controls, rr.controls_cross)
        except QpgenOutputError as e:
            raise StepError(str(e)) from e
        df = long_rows(pairs, st, lab, level_cfg)
        longs.append(df)
        labels.append(lab)
        run_info.append({"label": lab, "strategy": s, "rint": rint, "n_cov": n_cov})
        row = strategy_row(df, lab, rint, levels, n_cov)
        if controls is not None:
            c_all, c_cross = controls
            row.update(
                {
                    "n_control_cells": c_all.n,
                    "lambda_gc_control": c_all.lambda_gc(),
                    "n_control_p_lt_1e5": c_all.n_sig,
                    "n_control_cross_chrom": c_cross.n,
                    "lambda_gc_control_cross_chrom": c_cross.lambda_gc(),
                    "n_control_cross_chrom_p_lt_1e5": c_cross.n_sig,
                }
            )
        rows.append(row)
    if len(set(labels)) != len(labels):
        raise StepError(f"duplicate run labels {labels}")
    long = pd.concat(longs, ignore_index=True)
    neff = neff_columns(long, labels)
    for r in rows:
        r.update(neff[r["run"]])
    strategies = pd.DataFrame(rows)

    # ------------------------------------------------------------ resolution check
    n_mapped = len(pairs)
    unresolved = {r["label"]: (r["n_variant_absent"] + r["n_not_tested"]) / n_mapped for r in rows}
    worst_label = max(unresolved, key=lambda k: unresolved[k])
    unresolved_fraction = unresolved[worst_label]
    n_absent = int((~pairs["present"]).sum())
    absent_variants = sorted(set(pairs.loc[~pairs["present"], "variant"]))
    n_absent_indel = sum(1 for v in absent_variants if indel.get(v, False))
    max_unres = float(p.get("max_unresolved_fraction", 0.2))
    if unresolved_fraction > max_unres:
        worst = next(r for r in rows if r["label"] == worst_label)
        raise StepError(
            f"cis {ds}: {unresolved_fraction:.1%} of the {n_mapped} mapped reference pairs could not be "
            f"evaluated for {worst_label} ({worst['n_variant_absent']} variant absent from the genotypes, "
            f"{worst['n_not_tested']} present but without a finite qpgen result); the limit is "
            f"pipelines.cis.max_unresolved_fraction = {max_unres:g}. Absent variants usually mean a "
            f"different genome build, allele order or indel representation (e.g. {', '.join(absent_variants[:3])}); "
            "results missing for present variants point to a qpgen problem (see the association step log)"
        )

    # ---------------------------------------------------------------- warnings
    ds_json = ctx.path(Layout.ds_json(ds))
    if ds_json.exists():
        warnings += list(read_json(ds_json).get("warnings") or [])
    n_unmapped_pairs = int(mapping.get("n_pairs", n_mapped)) - int(
        mapping.get("n_pairs_mapped", n_mapped)
    )
    if n_unmapped_pairs > 0:
        warnings.append(
            f"{n_unmapped_pairs} reference pair(s) of {mapping.get('n_unmapped_traits', '?')} trait(s) did not "
            f"map onto the phenotypes (trait ID rule {mapping.get('rule')!r}) and are excluded"
        )
    if n_absent:
        warnings.append(
            f"{n_absent} mapped pair(s) have a variant absent from the genotypes ({len(absent_variants)} "
            f"variant(s): {n_absent_indel} indel(s), {len(absent_variants) - n_absent_indel} SNV(s)); "
            "allele-flipped or differently normalised IDs are not matched"
        )
    for r in rows:
        if r["n_not_tested"]:
            warnings.append(
                f"{r['label']}: {r['n_not_tested']} pair(s) with a present variant have no (finite) qpgen result "
                "(e.g. no genotype variation among the phenotyped samples; see the association step log)"
            )
        if r["n_tested"] == 0:
            warnings.append(f"{r['label']}: no pair was tested")
    if recomputed:
        warnings.append(
            "pair-assoc rows lack the raw LOG10P field (unpatched qpgen) for "
            f"{', '.join(recomputed)}; raw -log10 p was recomputed from TSTAT with N - 2 - n_cov df "
            "(upgrade qpgen: such builds also do not residualise genotypes on covariates)"
        )
    if dup.any():
        warnings.append(
            f"{int(dup.sum())} duplicated (trait, variant) pair(s) after trait mapping were dropped"
        )

    # ---------------------------------------------------------------- outputs
    # Level stringency by cut-off: Bonferroni beats strict once n_tested > 0.05 / strict.
    ref_label = max(labels, key=lambda lab: next(r["n_tested"] for r in rows if r["run"] == lab))
    order = stringency_order(level_cutoffs(long[long["run"] == ref_label], levels, level_cfg))
    rec = recommend(rows, levels, level_cfg, order)
    ranking = {
        lv: [r["label"] for r in sorted(rows, key=lambda r: -r[f"n_rep_{lv}"])] for lv in levels
    }
    traits = trait_table(long, labels, levels, level_cfg)
    curve = power_curve(long, labels, float(p.get("power_step", 0.5)))
    per_run = {r["label"]: {k: _num(v) for k, v in r.items() if k not in ("label",)} for r in rows}
    summary = {
        "schema_version": SCHEMA_VERSION,
        "dataset": ds,
        "reference": ref,
        "engine": engine,
        "rint": rint,
        "levels": level_description(level_cfg),
        "level_order": order,
        "transform_semantics": p.get("transform_semantics", {}),
        "runs": run_info,
        "baseline": labels[0],
        "n_ref_pairs": int(mapping.get("n_pairs", n_mapped)),
        "n_pairs_mapped": n_mapped,
        "n_pairs_unmapped": n_unmapped_pairs,
        "n_traits_mapped": int(pairs["trait"].nunique()),
        "n_traits_unmapped": mapping.get("n_unmapped_traits"),
        "n_variant_absent": n_absent,
        "n_absent_variants": len(absent_variants),
        "n_absent_variants_indel": n_absent_indel,
        "n_absent_variants_snv": len(absent_variants) - n_absent_indel,
        "absent_variant_examples": absent_variants[:20],
        "unresolved_fraction": unresolved_fraction,
        "unresolved_by_run": unresolved,
        "max_unresolved_fraction": max_unres,
        "recomputed_log10p": recomputed,
        "per_run": per_run,
        "ranking": ranking,
        "recommended": rec,
        "warnings": warnings,
        "sentinels": None,
    }
    write_tsv(ctx.path(Layout.cis_summary(ds, "pairs.long.tsv.gz")), long)
    write_tsv(ctx.path(Layout.cis_summary(ds, "traits.tsv")), traits)
    write_tsv(ctx.path(Layout.cis_summary(ds, "strategies.tsv")), strategies)
    write_tsv(ctx.path(Layout.cis_summary(ds, "power_curve.tsv")), curve)
    write_json(ctx.path(Layout.cis_summary(ds, "summary.json")), summary)
    if ctx.log is not None:
        ctx.log.info("cis %s: %d mapped pairs, recommended %s", ds, n_mapped, rec.get("label"))
