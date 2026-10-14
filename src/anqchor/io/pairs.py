"""Reference pairs, CPRA variant IDs and presence checks against indexed pvar files.

A variant is addressed as ``CHROM:POS:REF:ALT`` (CPRA). qpgen looks variants up by
exact CHROM, POS, REF and ALT string equality; :func:`check_presence` reproduces that
lookup with ``tabix -R`` region queries or, without tabix, by streaming the pvar.
"""

from __future__ import annotations

import gzip
import math
import os
import subprocess
import tempfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import pandas as pd

from anqchor.io.genotypes import ChromStyle, PgenEntry, find_tool, normalize_chrom
from anqchor.io.tables import (
    InputError,
    atomic_text_writer,
    detect_sep,
    first_line,
    read_user_csv,
    split_fields,
)
from anqchor.util.files import open_text


class Variant(NamedTuple):
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def cpra(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt)


_NUCLEOTIDES = set("ACGTNacgtn")


def parse_cpra(text: str) -> Variant | None:
    """Parse ``CHROM:POS:REF:ALT``; ``None`` when malformed (POS must be a positive
    integer, REF and ALT non-empty without whitespace). Nucleotide alleles are
    upper-cased (``a`` = ``A``), as in the pvar files qpgen compares them with."""
    f = text.strip().split(":")
    if len(f) != 4 or not f[0] or not f[1].isdigit() or not f[2] or not f[3]:
        return None
    pos = int(f[1])
    if pos <= 0 or any(c.isspace() for c in f[2] + f[3]):
        return None
    ref, alt = (a.upper() if set(a) <= _NUCLEOTIDES else a for a in (f[2], f[3]))
    return Variant(f[0], pos, ref, alt)


@dataclass
class RefPairs:
    """Normalised reference pairs and the counts of what was dropped."""

    table: pd.DataFrame
    """``trait variant beta se pval`` (strings; pval ``NA`` when absent)."""
    variants: dict[str, Variant] = field(default_factory=dict)
    """CPRA (normalised) → parsed variant, unique, in first-seen order."""
    n_input: int = 0
    n_invalid: int = 0
    n_duplicates: int = 0
    n_renamed: int = 0
    """Unique variants whose chromosome name was rewritten to the genotype style."""
    invalid_examples: list[str] = field(default_factory=list)


def _finite(s: str) -> bool:
    try:
        return math.isfinite(float(s))
    except (TypeError, ValueError):
        return False


def read_reference_pairs(
    path: str | os.PathLike[str], columns: dict[str, str | None], style: ChromStyle
) -> RefPairs:
    """Read a pair file with a header, mapping user column names (``columns`` keys:
    trait, variant, beta, se, optional pval) and normalising chromosome names. Without a
    ``pval`` mapping, a column named ``pval`` is read when present.

    Rows with a malformed variant, a missing trait, or non-finite beta/se are dropped
    (``n_invalid``); repeated ``(trait, variant)`` keep the first row (``n_duplicates``).
    """
    path = Path(path)
    line = first_line(path)
    if not line:
        raise InputError(f"{path}: empty pair file")
    sep = detect_sep(line)
    header = split_fields(line, sep)
    need = {k: columns.get(k) or k for k in ("trait", "variant", "beta", "se")}
    # The p-value column is optional: an explicit name must exist; otherwise a column
    # called 'pval' is used when the file has one (the '[pval]' of '--pairs').
    pcol = columns.get("pval") or ("pval" if "pval" in header else None)
    missing = [f"{k} (column {v!r})" for k, v in need.items() if v not in header]
    if pcol and pcol not in header:
        missing.append(f"pval (column {pcol!r})")
    if missing:
        raise InputError(
            f"{path}: missing column(s) {', '.join(missing)}; header has {header[:12]}. "
            "Set references.<name>.columns to the file's column names"
        )
    use = [need["trait"], need["variant"], need["beta"], need["se"], *([pcol] if pcol else [])]
    df = read_user_csv(
        path,
        "reference pair table",
        sep=sep,
        usecols=use,
        dtype=str,
        keep_default_na=False,
        engine="c",
    )
    res = RefPairs(table=pd.DataFrame(), n_input=len(df))
    traits = df[need["trait"]].str.strip().tolist()
    var_text = df[need["variant"]].tolist()
    betas = df[need["beta"]].str.strip().tolist()
    ses = df[need["se"]].str.strip().tolist()
    pvals = df[pcol].str.strip().tolist() if pcol else ["NA"] * len(df)
    out_rows: list[tuple[str, str, str, str, str]] = []
    seen: set[tuple[str, str]] = set()
    renamed: set[str] = set()
    for t, v, b, s, p in zip(traits, var_text, betas, ses, pvals, strict=True):
        var = parse_cpra(v)
        if var is None or not t or not _finite(b) or not _finite(s):
            res.n_invalid += 1
            if len(res.invalid_examples) < 5:
                res.invalid_examples.append(f"{t}\t{v}\t{b}\t{s}")
            continue
        nchrom = normalize_chrom(var.chrom, style)
        if nchrom != var.chrom:
            var = var._replace(chrom=nchrom)
            renamed.add(var.cpra)
        cpra = var.cpra
        key = (t, cpra)
        if key in seen:
            res.n_duplicates += 1
            continue
        seen.add(key)
        res.variants.setdefault(cpra, var)
        out_rows.append((t, cpra, b, s, p if p and _is_number_or_na(p) else "NA"))
    res.n_renamed = len(renamed)
    res.table = pd.DataFrame(out_rows, columns=["trait", "variant", "beta", "se", "pval"])
    return res


def _is_number_or_na(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_variant_list(
    path: str | os.PathLike[str], style: ChromStyle
) -> tuple[dict[str, Variant], int, int]:
    """Read a CPRA list (first column of each line; ``#`` lines and a non-CPRA first
    line treated as a header are skipped). Returns (variants, n_invalid, n_renamed)."""
    out: dict[str, Variant] = {}
    n_invalid = n_renamed = 0
    with open_text(path) as fh:
        first = True
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()[0]
            var = parse_cpra(tok)
            if var is None:
                if first and ":" not in tok:
                    first = False
                    continue
                n_invalid += 1
                continue
            first = False
            nchrom = normalize_chrom(var.chrom, style)
            if nchrom != var.chrom:
                var = var._replace(chrom=nchrom)
                n_renamed += 1
            out.setdefault(var.cpra, var)
    return out, n_invalid, n_renamed


# --------------------------------------------------------------------- presence
def _targets(variants: Iterable[Variant], entry: PgenEntry) -> list[Variant]:
    return [v for v in variants if v.chrom == entry.chrom and entry.beg <= v.pos <= entry.end]


def _present_tabix(tabix: str, entry: PgenEntry, want: list[Variant]) -> set[Variant]:
    wanted = set(want)
    found: set[Variant] = set()
    with tempfile.NamedTemporaryFile("w", suffix=".regions.tsv", delete=False) as fh:
        for pos in sorted({v.pos for v in want}):
            fh.write(f"{entry.chrom}\t{pos}\n")
        regions = fh.name
    try:
        proc = subprocess.Popen(
            [tabix, "-R", regions, str(entry.pivar)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            f = line.split("\t", 5)
            if len(f) < 5:
                continue
            try:
                v = Variant(f[0], int(f[1]), f[3], f[4].rstrip("\n"))
            except ValueError:
                continue
            if v in wanted:
                found.add(v)
        err = proc.stderr.read() if proc.stderr else ""
        if proc.wait() != 0:
            raise InputError(f"tabix failed on {entry.pivar}: {err.strip()}")
    finally:
        os.unlink(regions)
    return found


def _present_stream(entry: PgenEntry, want: list[Variant]) -> set[Variant]:
    chrom_b = entry.chrom.encode()
    by_pos: dict[bytes, set[tuple[bytes, bytes]]] = {}
    for v in want:
        by_pos.setdefault(str(v.pos).encode(), set()).add((v.ref.encode(), v.alt.encode()))
    found: set[Variant] = set()
    with gzip.open(entry.pivar, "rb") as fh:
        for line in fh:
            if line[:1] == b"#":
                continue
            c, p, rest = line.split(b"\t", 2)
            if c != chrom_b or p not in by_pos:
                continue
            f = rest.split(b"\t", 3)
            if len(f) < 3:
                continue
            ref, alt = f[1], f[2].rstrip(b"\r\n")
            if (ref, alt) in by_pos[p]:
                found.add(Variant(entry.chrom, int(p), ref.decode(), alt.decode()))
    return found


def check_presence(
    variants: Sequence[Variant], entries: Sequence[PgenEntry], method: str = "auto"
) -> tuple[set[Variant], str]:
    """Variants found in the genotype pvars with identical CHROM, POS, REF and ALT.

    Each variant is looked up in the pgen-list entries whose CHROM matches and whose
    BEG..END covers its position. ``method``: ``auto`` (tabix when available),
    ``tabix`` or ``python`` (stream the bgzipped pvar). Returns (present, method used).
    """
    tabix = find_tool("tabix") if method in ("auto", "tabix") else None
    if method == "tabix" and tabix is None:
        raise InputError("tabix requested but not found")
    used: set[str] = set()
    present: set[Variant] = set()
    for e in entries:
        want = _targets(variants, e)
        if not want:
            continue
        if tabix and Path(str(e.pivar) + ".tbi").exists():
            present |= _present_tabix(tabix, e, want)
            used.add("tabix")
        else:
            present |= _present_stream(e, want)
            used.add("python")
    return present, "+".join(sorted(used)) or ("tabix" if tabix else "python")


def chrom_order(entries: Sequence[PgenEntry]) -> dict[str, int]:
    """Rank of each chromosome in pgen-list order."""
    order: dict[str, int] = {}
    for e in entries:
        order.setdefault(e.chrom, len(order))
    return order


def sort_variants(variants: Iterable[Variant], order: dict[str, int]) -> list[Variant]:
    """Sort by pgen-list chromosome order (unknown chromosomes last), position, alleles."""
    n = len(order)
    return sorted(
        variants,
        key=lambda v: (
            order.get(v.chrom, n),
            v.chrom if v.chrom not in order else "",
            v.pos,
            v.ref,
            v.alt,
        ),
    )


def write_pairs(path: str | os.PathLike[str], table: pd.DataFrame, columns: Sequence[str]) -> None:
    """Write a string table with a header (tab separated)."""
    with atomic_text_writer(path) as fh:
        fh.write("\t".join(columns) + "\n")
        for row in table[list(columns)].itertuples(index=False, name=None):
            fh.write("\t".join(str(x) for x in row) + "\n")


def read_pairs_table(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a pairs/variants table written by anqchor (all columns as strings)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, engine="c")
