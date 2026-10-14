"""Covariate encoding and selection.

Covariates are read once and kept in two forms:

* **raw** (:attr:`CovEncoding.raw`): one column per covariate, categorical covariates
  kept as their text levels (e.g. ``XX``/``XY``). This is what the pipelines write and
  hand to qpgen (``match-prs-pheno``, ``pair-assoc``, ``rect-assoc``), which builds the
  indicator columns itself from ``--cat-cov-names``.
* **encoded** (:attr:`CovEncoding.values`): all numeric. Two levels become one column
  coded 1/2 (sorted levels → 1, 2); k > 2 levels become k - 1 indicator columns
  ``<name>_<level>`` coded 2 (sample has the level) or 1 (otherwise), with the first
  sorted level as the reference. Used to count covariates, to resolve encoded names in
  ``--cov-include``, and by ``anqchor prepare cov --encode``.

Both codings span the same column space as qpgen's L - 1 indicators, so association
results do not depend on which one is used. Constant numeric columns are dropped and
missing values are handled by policy (``error``, ``drop_samples``, ``impute_mean``; for
a categorical covariate ``impute_mean`` uses the most frequent level).
"""

from __future__ import annotations

import os
import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from anqchor.io.pheno import read_matrix, write_regenie
from anqchor.io.tables import InputError, examples

MISSING_POLICIES = ("error", "drop_samples", "impute_mean")
MAX_AUTO_LEVELS = 20
"""A non-numeric covariate not listed in ``cov_categorical`` with more levels than this
(and more than a tenth of the samples) is an error rather than a silent categorical."""
_RANGE_RE = re.compile(r"^(?P<prefix>.*?\D)(?P<a>\d+)-(?P=prefix)?(?P<b>\d+)$")
_NUM_SUFFIX_RE = re.compile(r"^(?P<prefix>.*?\D)(?P<n>\d+)$")


@dataclass
class CovEncoding:
    """Encoded covariates and how they were derived."""

    values: pd.DataFrame
    """samples × encoded columns (float64)."""
    groups: dict[str, list[str]] = field(default_factory=dict)
    """original covariate name → encoded column names (empty when dropped)."""
    categorical: dict[str, dict[str, int]] = field(default_factory=dict)
    """categorical covariate → {level: code} (codes 1..k in sorted level order)."""
    dropped_constant: list[str] = field(default_factory=list)
    dropped_all_missing: list[str] = field(default_factory=list)
    missing_policy: str = "error"
    n_missing_handled: int = 0
    """Samples dropped (``drop_samples``) or values imputed (``impute_mean``)."""
    dropped_samples: list[str] = field(default_factory=list)
    raw: pd.DataFrame | None = None
    """samples × covariates, categorical columns as text levels, numeric columns as
    floats (same samples and missing-value handling as :attr:`values`)."""

    @property
    def categorical_names(self) -> list[str]:
        """Categorical covariates present in :attr:`raw` (qpgen ``--cat-cov-names``)."""
        cols = set(self.raw.columns) if self.raw is not None else set()
        return [n for n in self.categorical if n in cols]

    @property
    def columns(self) -> list[str]:
        return [str(c) for c in self.values.columns]

    def as_json(self) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "columns": self.columns,
            "names": list(self.groups),
            "groups": self.groups,
            "categorical": self.categorical,
            "categorical_names": self.categorical_names,
            "raw_columns": [str(c) for c in self.raw.columns] if self.raw is not None else [],
            "dropped_constant": self.dropped_constant,
            "dropped_all_missing": self.dropped_all_missing,
            "n_samples": int(self.values.shape[0]),
            "missing_policy": self.missing_policy,
            "n_missing_handled": self.n_missing_handled,
            "dropped_samples": self.dropped_samples[:100],
        }


def _sorted_levels(levels: Sequence[str]) -> list[str]:
    """Numeric order when every level is a number, string order otherwise."""
    nums = pd.to_numeric(pd.Series(list(levels), dtype=object), errors="coerce")
    if not nums.isna().any():
        return [lv for _, lv in sorted(zip(nums.tolist(), levels, strict=True))]
    return sorted(levels)


def _safe(level: str) -> str:
    return re.sub(r"\s+", "_", level.strip())


def encode_covariates(
    raw: pd.DataFrame, categorical: Sequence[str] = (), missing_policy: str = "error"
) -> CovEncoding:
    """Encode a samples × covariates frame of strings (NaN or an empty cell = missing).

    Non-numeric columns are categorical automatically; ``categorical`` forces the
    named columns to be categorical too (an unknown name is an error). A column of
    mostly numbers with a few non-numeric cells is an error naming the cell (add the
    string to ``missing_str``), as is a non-numeric column with very many levels that is
    not listed in ``categorical``, an infinite value, or two levels whose encoded column
    names coincide.
    """
    if missing_policy not in MISSING_POLICIES:
        raise InputError(
            f"cov_missing must be one of {', '.join(MISSING_POLICIES)}, got {missing_policy!r}"
        )
    names = [str(c) for c in raw.columns]
    unknown = [c for c in categorical if c not in names]
    if unknown:
        raise InputError(
            f"cov_categorical names unknown covariate(s) {examples(unknown)}; "
            f"available: {compact_names(names)}"
        )
    enc = CovEncoding(values=pd.DataFrame(index=raw.index), missing_policy=missing_policy)
    cols: dict[str, pd.Series] = {}
    raw_cols: dict[str, pd.Series] = {}
    forced = set(categorical)
    for name in names:
        # An empty (or blank) cell is missing, never a category level or a number.
        s = raw[name].map(
            lambda v: None if v is None or (isinstance(v, str) and not v.strip()) else v
        )
        s = s.where(s.notna(), np.nan)
        present = s.dropna()
        if present.empty:
            enc.dropped_all_missing.append(name)
            enc.groups[name] = []
            continue
        num = pd.to_numeric(present, errors="coerce")
        n_bad = int(num.isna().sum())
        if name not in forced and n_bad > 0 and n_bad * 2 < len(present):
            # Mostly numbers: a continuous covariate with a few odd cells (e.g. 'nan', '.')
            # must not silently become a categorical with one indicator per value.
            bad = present[num.isna()]
            raise InputError(
                f"non-numeric value {str(bad.iloc[0])!r} for covariate {name!r} in sample "
                f"{str(bad.index[0])!r} ({n_bad} such cell(s) in a column of numbers); add the "
                "string to missing_str if it means 'missing', or list the covariate in "
                "cov_categorical if it is categorical"
            )
        if name in forced or n_bad:
            text = pd.Series(
                [None if pd.isna(v) else str(v).strip() for v in s], index=s.index, dtype=object
            )
            levels = _sorted_levels(list(dict.fromkeys(text.dropna().tolist())))
            if (
                name not in forced
                and len(levels) > MAX_AUTO_LEVELS
                and len(levels) * 10 > len(present)
            ):
                raise InputError(
                    f"covariate {name!r} is not numeric and has {len(levels)} distinct values "
                    f"(e.g. {examples(levels)}); as a categorical it would add {len(levels) - 1} "
                    "indicator columns. Fix the values, or list it in cov_categorical to confirm"
                )
            enc.categorical[name] = {lv: i + 1 for i, lv in enumerate(levels)}
            raw_cols[name] = text
            if len(levels) == 1:
                enc.dropped_constant.append(name)
                enc.groups[name] = []
            elif len(levels) == 2:
                cols[name] = text.map(enc.categorical[name]).astype(np.float64)
                enc.groups[name] = [name]
            else:
                enc.groups[name] = []
                seen: dict[str, str] = {}
                for lv in levels[1:]:
                    cname = f"{name}_{_safe(lv)}"
                    if cname in seen:
                        raise InputError(
                            f"categorical covariate {name!r}: levels {seen[cname]!r} and {lv!r} both "
                            f"encode as column {cname!r}; rename one of them"
                        )
                    seen[cname] = lv
                    isna = text.isna().to_numpy(dtype=bool)
                    ind = np.where(isna, np.nan, np.where(text.to_numpy() == lv, 2.0, 1.0))
                    cols[cname] = pd.Series(ind, index=s.index)
                    enc.groups[name].append(cname)
        else:
            vals_num = pd.to_numeric(s, errors="coerce").astype(np.float64)
            inf = np.isinf(vals_num.to_numpy())
            if inf.any():
                i = int(np.flatnonzero(inf)[0])
                raise InputError(
                    f"infinite value {str(s.iloc[i])!r} for covariate {name!r} in sample "
                    f"{str(s.index[i])!r} ({int(inf.sum())} such cell(s)); covariates must be finite"
                )
            cols[name] = vals_num
            enc.groups[name] = [name]
    clash = [c for g, cs in enc.groups.items() for c in cs if c != g and c in names]
    if clash:
        raise InputError(
            f"encoded covariate names collide with existing covariates: {examples(clash)}"
        )
    vals = pd.DataFrame(cols, index=raw.index)
    # constant columns (ignoring missing values) are dropped before handling missing values
    vals = _drop_constant(vals, enc)
    na = vals.isna()
    if na.to_numpy().any():
        if missing_policy == "error":
            counts = na.sum()
            bad = counts[counts > 0]
            first_col = str(bad.index[0])
            first_sample = str(vals.index[np.flatnonzero(na[first_col].to_numpy())[0]])
            raise InputError(
                f"missing covariate values in {len(bad)} column(s): "
                + ", ".join(f"{c} ({int(n)})" for c, n in list(bad.items())[:5])
                + f"; e.g. {first_col!r} in sample {first_sample!r}. Set cov_missing to "
                "drop_samples or impute_mean, or fix the file"
            )
        if missing_policy == "drop_samples":
            keep = ~na.any(axis=1)
            enc.dropped_samples = [str(s) for s in vals.index[~keep]]
            enc.n_missing_handled = len(enc.dropped_samples)
            vals = vals[keep]
            if vals.empty:
                raise InputError("cov_missing=drop_samples removed every sample")
            vals = _drop_constant(vals, enc)
        else:
            enc.n_missing_handled = int(na.to_numpy().sum())
            vals = vals.fillna(vals.mean())
    enc.values = vals
    enc.raw = _raw_frame(raw_cols, vals, enc)
    return enc


def _raw_frame(
    cat_text: dict[str, pd.Series], vals: pd.DataFrame, enc: CovEncoding
) -> pd.DataFrame:
    """Covariates with categorical columns as text: numeric columns as encoded (after
    constant-column dropping and missing-value handling), categorical columns as their
    levels over the same samples. A categorical whose missing values were imputed takes
    the most frequent level; one with a single level is kept (qpgen drops it with a
    warning), so the names passed as ``--cat-cov-names`` always exist."""
    out: dict[str, pd.Series] = {}
    for name, g in enc.groups.items():
        if name in cat_text:
            t = cat_text[name].reindex(vals.index)
            if enc.missing_policy == "impute_mean" and t.isna().any():
                t = t.fillna(t.mode().iloc[0])
            out[name] = t
        elif g and name in vals.columns:
            out[name] = vals[name]
    return pd.DataFrame(out, index=vals.index)


def _drop_constant(vals: pd.DataFrame, enc: CovEncoding) -> pd.DataFrame:
    const = [c for c in vals.columns if vals[c].nunique(dropna=True) <= 1]
    if not const:
        return vals
    for c in const:
        enc.dropped_constant.append(str(c))
        for g in enc.groups.values():
            if c in g:
                g.remove(c)
    return vals.drop(columns=const)


def read_and_encode(
    path: str | os.PathLike[str],
    fmt: str,
    categorical: Sequence[str] = (),
    missing_policy: str = "error",
    missing: Sequence[str] = ("NA",),
    *,
    ignore_unknown_categorical: bool = False,
) -> CovEncoding:
    """Read a covariate file in any layout and encode it (:func:`encode_covariates`)."""
    m = read_matrix(path, fmt, missing, numeric=False, what="covariate")
    cats = list(categorical)
    if ignore_unknown_categorical:
        cats = [c for c in cats if c in m.values.columns]
    try:
        return encode_covariates(m.values, cats, missing_policy)
    except InputError as e:
        raise InputError(f"{path}: {e}") from e


def write_cov(path: str | os.PathLike[str], enc: CovEncoding) -> None:
    """Write encoded (all-numeric) covariates in regenie layout."""
    write_regenie(path, enc.values, trim_int=True)


def raw_as_text(raw: pd.DataFrame) -> pd.DataFrame:
    """String frame for writing: numbers in shortest round-trip form, ``NA`` for missing."""
    from anqchor.io.tables import fmt_float

    out = {}
    for c in raw.columns:
        s = raw[c]
        if s.dtype.kind == "f":
            out[c] = s.map(lambda v: "NA" if v != v else fmt_float(v))
        else:
            out[c] = s.map(
                lambda v: "NA" if v is None or (isinstance(v, float) and v != v) else str(v)
            )
    return pd.DataFrame(out, index=raw.index, dtype=object)


def write_cov_raw(path: str | os.PathLike[str], raw: pd.DataFrame) -> None:
    """Write covariates with categorical columns as text, in regenie layout."""
    write_regenie(path, raw_as_text(raw))


def read_cov_raw(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a file written by :func:`write_cov_raw` (all columns as text)."""
    m = read_matrix(path, "regenie", ["NA"], numeric=False, what="covariate")
    return m.values


def original_names(
    include: Sequence[str], groups: dict[str, list[str]], columns: Sequence[str]
) -> list[str]:
    """Original covariate names selected by ``include`` (original or encoded names, ranges
    expanded), in order without repeats; an encoded column maps to its covariate.
    Unknown names are an error (as in :func:`select_columns`)."""
    available = list(dict.fromkeys([*groups, *columns]))
    owner = {c: g for g, cs in groups.items() for c in cs}
    out: list[str] = []
    unknown: list[str] = []
    for token in include:
        for name in expand_name(token.strip(), available):
            if name in groups:
                out.append(name)
            elif name in owner:
                out.append(owner[name])
            else:
                unknown.append(name)
    if unknown:
        raise InputError(
            f"unknown covariate name(s) {examples(unknown)}; available: {compact_names(list(groups))}"
        )
    return list(dict.fromkeys(out))


def detect_categorical(
    path: str | os.PathLike[str], fmt: str, missing: Sequence[str], declared: Sequence[str] = ()
) -> list[str]:
    """Categorical covariates of a file: the ``declared`` ones plus every column with a
    non-numeric value (ignoring missing strings), in file order. Used when the steps are
    built, to pass ``--cat-cov-names`` to qpgen; validation happens in the prepare step."""
    return list(categorical_levels(path, fmt, missing, declared))


def categorical_levels(
    path: str | os.PathLike[str], fmt: str, missing: Sequence[str], declared: Sequence[str] = ()
) -> dict[str, int]:
    """Number of distinct levels of each categorical covariate of a file (see
    :func:`detect_categorical`), in file order. Counted over all samples of the file, an
    upper bound of what qpgen counts among the analyzed samples."""
    m = read_matrix(path, fmt, missing, numeric=False, what="covariate")
    out: dict[str, int] = {}
    for name in m.values.columns:
        s = m.values[name].dropna()
        s = s[s.astype(str).str.strip() != ""]
        if name in declared or (len(s) and pd.to_numeric(s, errors="coerce").isna().any()):
            out[str(name)] = int(s.astype(str).str.strip().nunique())
    return out


# -------------------------------------------------------------------- selection
def expand_name(token: str, available: Sequence[str]) -> list[str]:
    """Expand ``PEER1-15`` (or ``PEER1-PEER15``) to ``PEER1 … PEER15``.

    A token that is itself an available name is returned unchanged, so names that
    contain a dash are never split.
    """
    if token in available:
        return [token]
    m = _RANGE_RE.match(token)
    if not m:
        return [token]
    a, b = int(m.group("a")), int(m.group("b"))
    if b < a:
        raise InputError(f"covariate range {token!r} runs backwards")
    return [f"{m.group('prefix')}{i}" for i in range(a, b + 1)]


def select_columns(
    include: Sequence[str], groups: dict[str, list[str]], columns: Sequence[str]
) -> tuple[list[str], list[str]]:
    """Encoded columns selected by ``include`` (original or encoded names, ranges
    expanded), in the given order without repeats.

    Returns ``(columns, empty)`` where ``empty`` lists requested covariates that exist
    but have no encoded column (dropped as constant). Unknown names are an error.
    """
    available = list(dict.fromkeys([*groups, *columns]))
    out: list[str] = []
    empty: list[str] = []
    unknown: list[str] = []
    for token in include:
        for name in expand_name(token.strip(), available):
            if name in groups:
                if not groups[name]:
                    empty.append(name)
                out.extend(groups[name])
            elif name in columns:
                out.append(name)
            else:
                unknown.append(name)
    if unknown:
        pos = {c: i for i, c in enumerate(columns)}
        names = sorted(
            groups, key=lambda g: min((pos[c] for c in groups[g] if c in pos), default=len(pos))
        )
        raise InputError(
            f"unknown covariate name(s) {examples(unknown)}; available: {compact_names(names or list(columns))}"
        )
    return list(dict.fromkeys(out)), empty


def compact_names(names: Sequence[str]) -> str:
    """``sex, PC1-5, PEER1-60`` style summary of covariate names for messages."""
    parts: list[str] = []
    run: tuple[str, int, int] | None = None

    def flush() -> None:
        nonlocal run
        if run is not None:
            p, a, b = run
            parts.append(f"{p}{a}" if a == b else f"{p}{a}-{b}")
            run = None

    for n in names:
        m = _NUM_SUFFIX_RE.match(n)
        if m:
            p, k = m.group("prefix"), int(m.group("n"))
            if run is not None and run[0] == p and run[2] + 1 == k:
                run = (p, run[1], k)
                continue
            flush()
            run = (p, k, k)
        else:
            flush()
            parts.append(n)
    flush()
    return ", ".join(parts)
