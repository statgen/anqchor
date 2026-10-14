"""Readers for qpgentools output files, with schema checks.

qpgen's dev branch changes; every reader asserts the columns it relies on and raises
:class:`QpgenOutputError` with the file name when they are missing. Readers also
tolerate two defects of older dev builds (``3864e6a`` and earlier):

* ``pair-assoc`` data rows one field short of the header (raw ``LOG10P`` missing):
  ``LOG10P`` is recomputed from ``TSTAT`` with ``N - 2 - n_cov`` degrees of freedom
  (``TSTAT`` rescaled from the old ``N - 2``; see :func:`read_pair_assoc`).
* ``rect-assoc`` ``N`` equal to the allele number: ``N`` is always recomputed as
  ``N_RR + N_RA + N_AA``.
"""

from __future__ import annotations

import re
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats

from anqchor.util.files import open_text, read_header

MATCH_ORDINALS = ("1st", "2nd", "3rd", "4th", "5th")
MATCH_ALL_BASE = (
    "ID.Pheno",
    "N.Traits",
    "MatchStatus",
    "ID.self",
    "Z.self",
    "COR.self",
    "Rank.self",
)
MATCH_STATUSES_ALL = (
    "SELF_BEST",
    "SELF_LENIENT",
    "SINGLE_NEW_BEST",
    "MULTI_NEW_BEST",
    "UNCLEAR",
    "NO_SELF",
    "NO_OBS_TRAITS",
)
MATCH_STATUSES_ASSIGNED = ("BEST_MATCH", "LENIENT_MATCH", "NO_MATCH", "NO_OBS_TRAITS")
PAIR_ASSOC_COLUMNS = (
    "TRAIT",
    "CHROM",
    "GENPOS",
    "ID",
    "ALLELE0",
    "ALLELE1",
    "A1FREQ",
    "N",
    "N_RR",
    "N_RA",
    "N_AA",
    "INFO",
    "BETA",
    "SE",
    "TSTAT",
    "LOG10P",
)
PAIR_ASSOC_RINT_COLUMNS = ("BETA_RINT", "SE_RINT", "TSTAT_RINT", "LOG10P_RINT")
RECT_META_COLUMNS = (
    "#CHROM",
    "GENPOS",
    "ID",
    "ALLELE0",
    "ALLELE1",
    "A1FREQ",
    "N",
    "N_RR",
    "N_RA",
    "N_AA",
    "TEST",
)
RECT_REQUIRED_META = (
    "#CHROM",
    "GENPOS",
    "ID",
    "ALLELE0",
    "ALLELE1",
    "A1FREQ",
    "N_RR",
    "N_RA",
    "N_AA",
)
"""Columns :func:`iter_rect_assoc` reads (``N`` is recomputed from the genotype counts)."""
RECT_STATS = ("BETA", "SE", "TSTAT", "LOG10P")
LAMBDA_RE = re.compile(r"Using automatically tuned lambda = ([0-9.eE+-]+)")


class QpgenOutputError(Exception):
    """A qpgen output file does not have the expected layout."""


def _require(df: pd.DataFrame, cols: Sequence[str], path: Path) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise QpgenOutputError(f"{path}: missing columns {missing}; qpgen output layout changed?")


def _t_log10p_tail(t: np.ndarray, d: np.ndarray) -> np.ndarray:
    """-log10 of the two-sided t p-value from the exact identity
    ``p = I_x(d/2, 1/2) = x^a (1-x)^b / (a B(a, b)) * 2F1(a+b, 1; a+1; x)`` with
    ``x = d / (d + t^2)``, evaluated in log space (finite for any finite ``t``)."""
    a = d / 2.0
    b = 0.5
    x = d / (d + t * t)
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = (
            a * np.log(x)
            + b * np.log1p(-x)
            - np.log(a)
            - special.betaln(a, b)
            + np.log(special.hyp2f1(a + b, 1.0, a + 1.0, x))
        )
    return np.asarray(-logp / np.log(10.0))


def t_log10p(tstat: np.ndarray | pd.Series, df: np.ndarray | pd.Series) -> np.ndarray:
    """Two-sided -log10 p of a t statistic with ``df`` degrees of freedom.

    Computed in log space; where scipy's ``logsf`` underflows (|t| in the tens with
    hundreds of degrees of freedom) an exact tail formula keeps the value finite."""
    t, d = np.broadcast_arrays(np.abs(np.asarray(tstat, dtype=float)), np.asarray(df, dtype=float))
    with np.errstate(invalid="ignore", divide="ignore"):
        logsf = stats.t.logsf(t, d)
    out = np.asarray(-(logsf + np.log(2.0)) / np.log(10.0), dtype=float)
    with np.errstate(invalid="ignore"):
        bad = ~np.isfinite(out) & np.isfinite(t) & (d > 0)
    if bad.any():
        out = out.copy()
        out[bad] = _t_log10p_tail(t[bad], d[bad])
    return out


# --------------------------------------------------------------------------- pair-prs
def read_prs(path: Path) -> pd.DataFrame:
    """PRS matrix (regenie layout) → DataFrame indexed by IID, one column per trait."""
    df = pd.read_csv(path, sep="\t", dtype={"FID": str, "IID": str})
    _require(df, ("FID", "IID"), path)
    return df.drop(columns=["FID"]).set_index("IID")


def read_prs_traits(path: Path) -> list[str]:
    """Trait IDs in a PRS matrix header, without reading the values."""
    header = read_header(path)
    if header[:2] != ["FID", "IID"]:
        raise QpgenOutputError(f"{path}: expected FID IID header, got {header[:2]}")
    return header[2:]


# --------------------------------------------------------------------- match-prs-pheno
def read_match_all(path: Path) -> pd.DataFrame:
    """``<out>.match.all.tsv.gz``: one row per phenotyped sample."""
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={"ID.Pheno": str, "ID.self": str, "MatchStatus": str},
        na_values=["NA", "nan"],
    )
    cols = list(MATCH_ALL_BASE)
    for o in MATCH_ORDINALS:
        cols += [f"ID.{o}", f"Z.{o}", f"COR.{o}"]
    _require(df, cols, path)
    for o in MATCH_ORDINALS:
        df[f"ID.{o}"] = df[f"ID.{o}"].astype("string")
    return df


def read_match_assigned(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={"ID.Pheno": str, "ID.self": str, "MatchStatus": str},
        na_values=["NA", "nan"],
    )
    _require(df, MATCH_ALL_BASE, path)
    return df


def read_match_weights(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"Trait": str})
    _require(df, ("Trait", "Weight", "N.Obs"), path)
    return df


def parse_auto_lambda(log_text: str) -> float | None:
    """λ chosen by ``--auto-lambda``, parsed from the qpgen log (``None`` if absent)."""
    found = LAMBDA_RE.findall(log_text)
    return float(found[-1]) if found else None


# ------------------------------------------------------------------------- pair-assoc
def read_pair_assoc(path: Path, n_cov: int | None = None) -> pd.DataFrame:
    """``pair-assoc`` output with normalised column names (``TRAIT`` without ``#``).

    Handles rows one field short (older builds without raw ``LOG10P``): the raw
    ``LOG10P`` is then recomputed from ``TSTAT`` with the residual degrees of freedom
    ``N - 2 - n_cov`` that current qpgen uses (``n_cov`` = number of covariate columns,
    an upper bound of their rank). Such builds computed ``TSTAT`` with ``N - 2`` df, so
    it is rescaled by ``sqrt((N - 2 - n_cov) / (N - 2))`` first."""
    header = read_header(path)
    if not header:
        raise QpgenOutputError(f"{path}: empty file")
    header = [h.lstrip("#") for h in header]
    _require(pd.DataFrame(columns=header), PAIR_ASSOC_COLUMNS, path)
    rows: list[list[str]] = []
    short = False
    with open_text(path) as fh:
        next(fh)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == len(header):
                rows.append(fields)
            elif len(fields) == len(header) - 1 and "LOG10P" in header:
                k = header.index("LOG10P")
                fields.insert(k, "nan")
                rows.append(fields)
                short = True
            elif line.strip():
                raise QpgenOutputError(
                    f"{path}: row has {len(fields)} fields, header has {len(header)}"
                )
    df = pd.DataFrame(rows, columns=header)
    str_cols = {"TRAIT", "CHROM", "ID", "ALLELE0", "ALLELE1"}
    for c in df.columns:
        if c not in str_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("GENPOS", "N", "N_RR", "N_RA", "N_AA"):
        df[c] = df[c].astype("Int64")
    if short:
        n = df["N"].astype(float)
        n_k = float(n_cov or 0)
        dof = (n - 2 - n_k).where(n - 2 - n_k > 0)
        t = df["TSTAT"].astype(float) * np.sqrt(dof / (n - 2))
        df["LOG10P"] = t_log10p(t, dof)
    df.attrs["recomputed_log10p"] = short
    return df


# ------------------------------------------------------------------------- rect-assoc
@dataclass
class RectChunk:
    """A block of rows from a ``rect-assoc`` matrix."""

    meta: pd.DataFrame
    """Columns ``CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N N_RR N_RA N_AA`` (N recomputed)."""
    stats: dict[str, np.ndarray]
    """Statistic name → float array of shape (rows, traits), traits in header order."""


def read_rect_header(path: Path) -> tuple[list[str], list[str]]:
    """Return (meta columns, trait IDs) of a ``rect-assoc`` file."""
    header = read_header(path)
    meta = [c for c in header if c in RECT_META_COLUMNS]
    if meta[:3] != ["#CHROM", "GENPOS", "ID"]:
        raise QpgenOutputError(f"{path}: unexpected rect-assoc header start {header[:3]}")
    missing = [c for c in RECT_REQUIRED_META if c not in header]
    if missing:
        raise QpgenOutputError(f"{path}: missing columns {missing}; qpgen output layout changed?")
    traits = [c[len("BETA.") :] for c in header if c.startswith("BETA.")]
    for s in RECT_STATS:
        names = [c[len(s) + 1 :] for c in header if c.startswith(s + ".")]
        if names != traits:
            raise QpgenOutputError(f"{path}: {s}.* columns do not match BETA.* columns")
    return meta, traits


def iter_rect_assoc(
    path: Path, stats_wanted: Sequence[str] = ("BETA", "LOG10P"), chunksize: int = 64
) -> Iterator[RectChunk]:
    """Stream a ``rect-assoc`` matrix in row blocks, loading only the requested statistics.

    Memory use is ``chunksize × n_traits × len(stats_wanted)`` doubles."""
    _, traits = read_rect_header(path)
    for s in stats_wanted:
        if s not in RECT_STATS:
            raise ValueError(f"unknown rect-assoc statistic {s}")
    meta_cols = list(RECT_REQUIRED_META)
    usecols = meta_cols + [f"{s}.{t}" for s in stats_wanted for t in traits]
    dtypes = {"#CHROM": str, "ID": str, "ALLELE0": str, "ALLELE1": str}
    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype=dtypes,
        chunksize=chunksize,
        na_values=["NA", "nan", "-nan"],
    )
    for block in reader:
        if block.empty:  # header-only file
            continue
        meta = block[meta_cols].rename(columns={"#CHROM": "CHROM"}).reset_index(drop=True)
        meta["N"] = meta["N_RR"] + meta["N_RA"] + meta["N_AA"]
        out = {}
        for s in stats_wanted:
            out[s] = block[[f"{s}.{t}" for t in traits]].to_numpy(dtype=float)
        yield RectChunk(meta=meta, stats=out)
