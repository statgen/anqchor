"""Phenotype and covariate matrices in the four layouts qpgen understands.

* ``regenie``: rows = samples, ``FID IID trait1 …``
* ``tensorqtl``: rows = traits, ``#chr start end phenotype_id s1 s2 …`` (header names
  matched case-insensitively, optional ``#``; ``chrom`` accepted for ``chr``)
* ``tsv-sample-row``: rows = samples, ``ID trait1 …``
* ``tsv-sample-col``: rows = traits, ``ID s1 s2 …``

:func:`read_matrix` returns a samples × features frame; :func:`write_regenie` writes
the canonical regenie layout; :func:`convert_to_regenie` does both, copying the input
bytes unchanged when it already is a valid regenie file.
"""

from __future__ import annotations

import csv
import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from anqchor.io.genotypes import ChromStyle, normalize_chrom
from anqchor.io.maps import RULES, trait_key
from anqchor.io.tables import (
    InputError,
    atomic_copy_gz,
    atomic_text_writer,
    detect_sep,
    duplicates,
    examples,
    first_line,
    fmt_float,
    is_int,
    read_user_csv,
    split_fields,
)

FORMATS = ("regenie", "tensorqtl", "tsv-sample-row", "tsv-sample-col")
SAMPLE_ROWS = ("regenie", "tsv-sample-row")


@dataclass
class Matrix:
    """A matrix oriented samples × features.

    Attributes:
        values: index = sample IDs, columns = feature (trait/covariate) IDs; float64
            when read numerically, otherwise strings with NaN for missing.
        fids: regenie FIDs (same order as ``values.index``), else ``None``.
        loci: tensorqtl only: ``trait chrom start end`` as in the file.
        sep: the separator of the input (``"\\t"`` or ``r"\\s+"``).
        header: the header fields of the input.
    """

    values: pd.DataFrame
    fids: list[str] | None
    loci: pd.DataFrame | None
    sep: str
    header: list[str]

    @property
    def samples(self) -> list[str]:
        return [str(s) for s in self.values.index]

    @property
    def features(self) -> list[str]:
        return [str(c) for c in self.values.columns]

    def n_missing(self) -> int:
        return int(self.values.isna().to_numpy().sum())


def _norm(name: str) -> str:
    return name.strip().lstrip("#").lower()


def _check_format(fmt: str) -> None:
    if fmt not in FORMATS:
        raise InputError(f"unknown matrix format {fmt!r}; use one of {', '.join(FORMATS)}")


def _id_layout(fmt: str, header: list[str], path: Path, what: str) -> tuple[int, list[str]]:
    """Number of leading index columns and validation of their names."""
    if fmt == "regenie":
        if len(header) < 2 or _norm(header[0]) != "fid" or _norm(header[1]) != "iid":
            raise InputError(
                f"{path}: a regenie {what} file starts with 'FID IID', found {header[:3]}; "
                "is it tsv-sample-row (one ID column) or another layout?"
            )
        return 2, header[2:]
    if fmt == "tensorqtl":
        names = [_norm(h) for h in header[:3]]
        if len(header) < 4 or names[0] not in ("chr", "chrom") or names[1:] != ["start", "end"]:
            raise InputError(
                f"{path}: a tensorqtl {what} file starts with '#chr start end ID', found {header[:4]}"
            )
        return 4, header[4:]
    return 1, header[1:]


def read_matrix(
    path: str | os.PathLike[str],
    fmt: str,
    missing: Sequence[str] = ("NA",),
    *,
    numeric: bool = True,
    what: str = "phenotype",
    exact: bool = True,
) -> Matrix:
    """Read a matrix in any supported layout, oriented samples × features.

    Args:
        missing: strings meaning "missing" (become NaN).
        numeric: parse values as floats and fail on the first non-numeric value
            (naming the feature and the sample); otherwise keep strings.
        what: ``phenotype`` or ``covariate``, used in messages.
        exact: parse decimals with correct rounding (slower); with ``False`` pandas'
            fast parser may differ from the decimal by one unit in the last place.

    Raises:
        InputError: unknown layout, bad header, duplicate sample or feature IDs,
            non-numeric values, empty matrix.
    """
    _check_format(fmt)
    path = Path(path)
    line = first_line(path)
    if not line:
        raise InputError(f"{path}: empty {what} file")
    sep = detect_sep(line)
    header = split_fields(line, sep)
    n_idx, rest = _id_layout(fmt, header, path, what)
    if not rest:
        raise InputError(
            f"{path}: no {'traits' if fmt in SAMPLE_ROWS else 'samples'} in the header"
        )
    dup = duplicates(rest)
    if dup:
        kind = "trait/covariate" if fmt in SAMPLE_ROWS else "sample"
        raise InputError(f"{path}: duplicate {kind} IDs in the header: {examples(dup)}")
    if len(set(header)) != len(header):
        raise InputError(f"{path}: an index column name repeats a data column name")
    _check_names(rest, path, "trait/covariate" if fmt in SAMPLE_ROWS else "sample", "header")
    idx_cols = header[:n_idx]
    dtype: Any = str if not numeric else {c: str for c in idx_cols}
    if fmt == "tensorqtl" and numeric:
        dtype = {header[0]: str, header[3]: str}
    # QUOTE_NONE: read exactly the bytes qpgen will see (it does not strip quotes).
    df = read_user_csv(
        path,
        f"{fmt} {what} matrix",
        sep=sep,
        header=0,
        names=header,
        dtype=dtype,
        na_values=list(missing),
        keep_default_na=False,
        low_memory=False,
        float_precision="round_trip" if exact else None,
        engine="c",
        quoting=csv.QUOTE_NONE,
    )
    if df.empty:
        raise InputError(f"{path}: no data rows")
    id_col = (
        header[1] if fmt == "regenie" else header[n_idx - 1] if fmt == "tensorqtl" else header[0]
    )
    ids = df[id_col]
    if ids.isna().any():
        row = int(np.flatnonzero(ids.isna().to_numpy())[0]) + 2
        raise InputError(f"{path}: missing ID in column {id_col!r} on line {row}")
    ids = ids.astype(str)
    _check_names(
        ids.tolist(),
        path,
        "sample" if fmt in SAMPLE_ROWS else "trait/covariate",
        f"column {id_col!r}",
    )
    if fmt == "regenie":
        _check_names(df[header[0]].astype(str).tolist(), path, "family", f"column {header[0]!r}")
    dup = duplicates(ids.tolist())
    if dup:
        kind = "sample" if fmt in SAMPLE_ROWS else "trait/covariate"
        raise InputError(f"{path}: duplicate {kind} IDs: {examples(dup)}")
    vals = df[rest]
    vals.index = pd.Index(ids.tolist())
    loci = None
    if fmt == "tensorqtl":
        loci = pd.DataFrame(
            {
                "trait": ids.tolist(),
                "chrom": df[header[0]].astype(str).tolist(),
                "start": df[header[1]].tolist(),
                "end": df[header[2]].tolist(),
            }
        )
    if numeric:
        vals = _to_numeric(vals, path, fmt, what)
    fids = df[header[0]].astype(str).tolist() if fmt == "regenie" else None
    if fmt not in SAMPLE_ROWS:
        vals = vals.T
    vals.columns = pd.Index([str(c) for c in vals.columns])
    vals.index = pd.Index([str(i) for i in vals.index])
    return Matrix(values=vals, fids=fids, loci=loci, sep=sep, header=header)


def _check_names(names: Sequence[str], path: Path, kind: str, where: str) -> None:
    """Reject IDs qpgen would read differently: it splits fields on any whitespace and
    does not strip quotes, so ``A B``, ``" A"`` or ``"A"`` would not match anqchor's IDs."""
    for n in names:
        if any(ch.isspace() for ch in n) or '"' in n:
            raise InputError(
                f"{path}: {kind} ID {n!r} in {where} contains whitespace or a quote; qpgen splits "
                "fields on any whitespace and does not strip quotes. Rename it (e.g. replace "
                "spaces with '_') and unquote the file"
            )


def _to_numeric(vals: pd.DataFrame, path: Path, fmt: str, what: str) -> pd.DataFrame:
    """Convert to float64, failing on the first non-numeric cell (in file order)."""
    obj_cols = [c for c in vals.columns if not pd.api.types.is_numeric_dtype(vals[c].dtype)]
    first: tuple[int, int, Any] | None = None
    converted: dict[str, pd.Series] = {}
    for j, c in enumerate(vals.columns):
        if c not in obj_cols:
            continue
        col = vals[c]
        num = pd.to_numeric(col, errors="coerce")
        bad = np.flatnonzero((num.isna() & col.notna()).to_numpy(dtype=bool))
        if len(bad):
            cand = (int(bad[0]), j, col.iloc[int(bad[0])])
            if first is None or cand[:2] < first[:2]:
                first = cand
        converted[c] = num
    if first is not None:
        i, j, v = first
        row_id, col_id = str(vals.index[i]), str(vals.columns[j])
        trait, sample = (col_id, row_id) if fmt in SAMPLE_ROWS else (row_id, col_id)
        label = "trait" if what == "phenotype" else what
        col = vals.iloc[:, j]
        if (
            fmt not in SAMPLE_ROWS
            and len(col) >= 2
            and col.notna().all()
            and converted[vals.columns[j]].isna().all()
        ):
            pos = j + (4 if fmt == "tensorqtl" else 1) + 1
            raise InputError(
                f"{path}: header column {pos} {col_id!r} has no numeric value; it looks like an "
                f"annotation column, not a sample. A {fmt} file has "
                + (
                    "exactly 4 leading columns (#chr start end phenotype_id) "
                    if fmt == "tensorqtl"
                    else "one leading ID column "
                )
                + "followed by samples: remove the extra column(s)"
            )
        raise InputError(
            f"{path}: non-numeric value {v!r} for {label} {trait!r} in sample {sample!r}"
            + (
                "; declare it in cov_categorical or let anqchor detect categorical covariates"
                if what == "covariate"
                else "; add the string to missing_str if it means 'missing'"
            )
        )
    if converted:
        vals = vals.assign(**converted)
    out = vals.astype(np.float64)
    inf = np.isinf(out.to_numpy())
    if inf.any():
        i, j = (int(x[0]) for x in np.nonzero(inf))
        row_id, col_id = str(out.index[i]), str(out.columns[j])
        trait, sample = (col_id, row_id) if fmt in SAMPLE_ROWS else (row_id, col_id)
        label = "trait" if what == "phenotype" else what
        raise InputError(
            f"{path}: infinite value ({out.iat[i, j]}) for {label} {trait!r} in sample {sample!r} "
            f"({int(inf.sum())} such cell(s)); {what} values must be finite"
        )
    return out


# ---------------------------------------------------------------------- writing
def write_regenie(
    path: str | os.PathLike[str],
    values: pd.DataFrame,
    *,
    trim_int: bool = False,
    block: int = 64,
) -> None:
    """Write a samples × features frame in regenie layout (``FID IID …``, FID = IID).

    Floats are written with the shortest text that round-trips to the same double
    (``NA`` for NaN); string frames are written as is. ``trim_int`` writes ``2``
    instead of ``2.0``.
    """
    cols = [str(c) for c in values.columns]
    samples = [str(s) for s in values.index]
    is_float = all(dt.kind == "f" for dt in values.dtypes)
    arr = values.to_numpy(dtype=np.float64 if is_float else object)
    with atomic_text_writer(path) as fh:
        fh.write("\t".join(["FID", "IID", *cols]) + "\n")
        buf: list[str] = []
        for i, s in enumerate(samples):
            row = arr[i].tolist()
            if is_float:
                if trim_int or any(v != v for v in row):
                    txt = (
                        "\t".join(map(fmt_float, row))
                        if trim_int
                        else "\t".join("NA" if v != v else repr(v) for v in row)
                    )
                else:
                    txt = "\t".join(map(repr, row))
            else:
                txt = "\t".join("NA" if pd.isna(v) else str(v) for v in row)
            buf.append(f"{s}\t{s}\t{txt}\n")
            if len(buf) >= block:
                fh.write("".join(buf))
                buf = []
        if buf:
            fh.write("".join(buf))


@dataclass
class ConvertResult:
    """What :func:`convert_to_regenie` did."""

    matrix: Matrix
    fast_path: bool
    n_missing: int


def convert_to_regenie(
    src: str | os.PathLike[str],
    fmt: str,
    dest: str | os.PathLike[str],
    missing: Sequence[str] = ("NA",),
) -> ConvertResult:
    """Validate a phenotype matrix and write it in regenie layout to ``dest`` (.gz).

    Fast path: a tab-separated regenie file with FID = IID, numeric values and no
    missing-value spelling other than ``NA`` is copied byte for byte (compressed when it
    was plain), so values are preserved exactly. Otherwise the values are parsed with
    correct rounding and written with round-trip precision (the same doubles).
    """
    if fmt == "regenie":
        m = read_matrix(src, fmt, missing, exact=False)
        n_missing = m.n_missing()
        ok = (
            m.sep == "\t"
            and m.fids == m.samples
            and (n_missing == 0 or list(dict.fromkeys(missing)) == ["NA"])
            and "\r" not in first_line(src)
        )
        if ok:
            atomic_copy_gz(src, dest)
            return ConvertResult(matrix=m, fast_path=True, n_missing=n_missing)
    m = read_matrix(src, fmt, missing, exact=True)
    write_regenie(dest, m.values)
    return ConvertResult(matrix=m, fast_path=False, n_missing=m.n_missing())


# ------------------------------------------------------------- trait annotation
_TRAIT_NAMES = ("trait", "id", "phenotype_id", "gene_id", "gene", "pid", "name")
_CHROM_NAMES = ("chr", "chrom", "chromosome", "seqnames", "seqname")


_CHROM_TOKENS = {*(str(i) for i in range(1, 23)), "X", "Y", "M", "MT"}


def _looks_like_chrom(x: str) -> bool:
    """``chr``-prefixed names, 1-22, X, Y, M, MT."""
    t = x.strip()
    return t.lower().startswith("chr") or t.upper() in _CHROM_TOKENS


def read_trait_annotation(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read trait loci as ``trait chrom start end`` (strings / int64).

    Accepts a BED-like file (``chrom start end id`` + anything, such as a tensorQTL
    matrix; header optional) or a ``trait chrom start end`` table with a header. Rows
    with a missing or non-integer position are dropped; a repeated trait keeps its first
    row.
    """
    path = Path(path)
    line = first_line(path)
    if not line:
        raise InputError(f"{path}: empty trait annotation")
    sep = detect_sep(line)
    f = split_fields(line, sep)
    if len(f) < 4:
        raise InputError(f"{path}: trait annotation needs at least 4 columns (chrom start end id)")
    names = [_norm(x) for x in f]
    bed_ok = is_int(f[1]) and is_int(f[2])  # chrom start end id
    tcse_ok = is_int(f[2]) and is_int(f[3])  # trait chrom start end
    has_header = line.startswith("#") or not (bed_ok or tcse_ok)
    if has_header and names[0] not in _CHROM_NAMES and "start" in names and "end" in names:
        tcol = next((names.index(n) for n in _TRAIT_NAMES if n in names), None)
        ccol = next((names.index(n) for n in _CHROM_NAMES if n in names), None)
        if tcol is None or ccol is None:
            raise InputError(f"{path}: cannot find trait and chrom columns in header {f}")
        cols = [ccol, names.index("start"), names.index("end"), tcol]
    elif has_header or (bed_ok and not tcse_ok):
        cols = [0, 1, 2, 3]
    elif tcse_ok and not bed_ok:
        cols = [1, 2, 3, 0]
    else:
        # Both layouts fit (numeric trait IDs): decide by which column looks like a
        # chromosome; refuse to guess when both or neither do.
        c0, c1 = _looks_like_chrom(f[0]), _looks_like_chrom(f[1])
        if c0 and not c1:
            cols = [0, 1, 2, 3]
        elif c1 and not c0:
            cols = [1, 2, 3, 0]
        else:
            raise InputError(
                f"{path}: cannot tell whether the headerless trait annotation is BED "
                f"(chrom start end id) or trait chrom start end from its first line {f[:4]}; "
                "add a header line (e.g. '#chr start end phenotype_id')"
            )
    df = read_user_csv(
        path,
        "trait annotation",
        sep=sep,
        header=0 if has_header else None,
        usecols=cols,
        dtype=str,
        keep_default_na=False,
        na_values=["NA", ""],
        low_memory=False,
        engine="c",
    )
    # usecols keeps file order: map each wanted position to its column label
    label = {c: df.columns[k] for k, c in enumerate(sorted(cols))}
    ordered = [label[c] for c in cols]
    df = df[df[ordered[3]].notna()]
    out = pd.DataFrame(
        {
            "trait": df[ordered[3]].astype(str),
            "chrom": df[ordered[0]].astype(str),
            "start": pd.to_numeric(df[ordered[1]], errors="coerce"),
            "end": pd.to_numeric(df[ordered[2]], errors="coerce"),
        }
    )
    out = out.dropna(subset=["start", "end"])
    out = out.drop_duplicates("trait", keep="first")
    out["start"] = out["start"].astype(np.int64)
    out["end"] = out["end"].astype(np.int64)
    return out.reset_index(drop=True)


@dataclass
class LociResult:
    """Trait loci for every phenotype trait (NaN/None when unknown)."""

    table: pd.DataFrame
    """``trait chrom start end`` in phenotype trait order."""
    rule: str | None
    """The trait-ID rule that matched annotation IDs to phenotype IDs."""
    n_with_loci: int


def annotate_traits(
    traits: Sequence[str],
    annot: pd.DataFrame | None,
    style: ChromStyle,
    preferred_rule: str = "as_is",
) -> LociResult:
    """Attach loci to phenotype traits.

    IDs are compared directly first; when that leaves traits unmatched, the dataset's
    ``preferred_rule`` and then the other trait-ID rules are tried and the one matching
    the most traits wins (annotation keys that are ambiguous under a rule are ignored).
    Chromosome names are rewritten to the genotype ``style``.
    """
    table = pd.DataFrame({"trait": list(traits), "chrom": None, "start": np.nan, "end": np.nan})
    if annot is None or annot.empty:
        return LociResult(table=table, rule=None, n_with_loci=0)
    order = ["as_is", preferred_rule, *[r for r in RULES if r not in ("as_is", preferred_rule)]]
    best: tuple[int, str, dict[str, int]] = (-1, "as_is", {})
    ann_ids = annot["trait"].tolist()
    for rule in dict.fromkeys(order):
        keys: dict[str, int] = {}
        amb: set[str] = set()
        for i, t in enumerate(ann_ids):
            k = trait_key(t, rule)
            if k in keys:
                amb.add(k)
            keys[k] = i
        for k in amb:
            keys.pop(k, None)
        hits = {t: keys[trait_key(t, rule)] for t in traits if trait_key(t, rule) in keys}
        if len(hits) > best[0]:
            best = (len(hits), rule, hits)
        if len(hits) == len(traits):
            break
    n, rule, hits = best
    if n <= 0:
        return LociResult(table=table, rule=None, n_with_loci=0)
    rows = [hits.get(t, -1) for t in traits]
    idx = np.array(rows)
    has = idx >= 0
    sub = annot.iloc[idx[has]]
    chroms = np.array([None] * len(traits), dtype=object)
    chroms[has] = [normalize_chrom(c, style) for c in sub["chrom"].tolist()]
    starts = np.full(len(traits), np.nan)
    ends = np.full(len(traits), np.nan)
    starts[has] = sub["start"].to_numpy(dtype=float)
    ends[has] = sub["end"].to_numpy(dtype=float)
    table = pd.DataFrame({"trait": list(traits), "chrom": chroms, "start": starts, "end": ends})
    return LociResult(table=table, rule=rule, n_with_loci=int(has.sum()))


def write_traits(path: str | os.PathLike[str], table: pd.DataFrame) -> None:
    """Write ``trait chrom start end`` with ``NA`` for unknown loci."""
    with atomic_text_writer(path) as fh:
        fh.write("trait\tchrom\tstart\tend\n")
        for t, c, s, e in zip(
            table["trait"], table["chrom"], table["start"], table["end"], strict=True
        ):
            if c is None or pd.isna(c) or pd.isna(s):
                fh.write(f"{t}\tNA\tNA\tNA\n")
            else:
                fh.write(f"{t}\t{c}\t{int(s)}\t{int(e)}\n")
