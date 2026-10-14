"""Identifier maps: trait-ID rules, reference-to-dataset trait maps, sample maps.

The four trait-ID rules of :data:`anqchor.config.model.TraitIdRule` are implemented
once here (:func:`trait_key`) and used for both reference-to-dataset mapping and for
matching trait annotations to phenotype IDs.
"""

from __future__ import annotations

import os
import re
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field

from anqchor.io.tables import InputError, examples
from anqchor.util.files import open_text

RULES = ("as_is", "before_colon", "strip_version", "before_colon_strip_version")
_VERSION_RE = re.compile(r"\.\d+$")


def _before_colon(s: str) -> str:
    return s.split(":", 1)[0]


def _strip_version(s: str) -> str:
    return _VERSION_RE.sub("", s)


_RULE_FUNCS: dict[str, Callable[[str], str]] = {
    "as_is": lambda s: s,
    "before_colon": _before_colon,
    "strip_version": _strip_version,
    "before_colon_strip_version": lambda s: _strip_version(_before_colon(s)),
}


def trait_key(trait_id: str, rule: str) -> str:
    """Apply a trait-ID rule: ``as_is``, ``before_colon`` (``ENSG…5:WASH7P`` →
    ``ENSG…5``), ``strip_version`` (drop a trailing ``.<digits>``) or
    ``before_colon_strip_version`` (both)."""
    try:
        return _RULE_FUNCS[rule](trait_id)
    except KeyError:
        raise InputError(f"unknown trait_id_rule {rule!r}; use one of {', '.join(RULES)}") from None


def key_index(ids: Iterable[str], rule: str) -> dict[str, list[str]]:
    """Map each rule key to the IDs that produce it (in input order)."""
    idx: dict[str, list[str]] = {}
    for i in ids:
        idx.setdefault(trait_key(i, rule), []).append(i)
    return idx


@dataclass
class TraitMapping:
    """Result of mapping reference trait IDs onto dataset trait IDs."""

    mapping: dict[str, str] = field(default_factory=dict)
    """reference trait → dataset trait (one to one)."""
    unmapped: list[str] = field(default_factory=list)
    """reference traits without a dataset trait."""
    rule: str = "as_is"


def count_matches(ref_ids: Sequence[str], ds_ids: Sequence[str], rule: str) -> int:
    """Number of reference IDs whose key matches some dataset key (ambiguity ignored)."""
    keys = {trait_key(d, rule) for d in ds_ids}
    return sum(1 for r in ref_ids if trait_key(r, rule) in keys)


def best_rule(ref_ids: Sequence[str], ds_ids: Sequence[str]) -> tuple[str, int]:
    """The rule matching the most reference IDs (earliest rule on ties) and its count."""
    best = ("as_is", -1)
    for rule in RULES:
        n = count_matches(ref_ids, ds_ids, rule)
        if n > best[1]:
            best = (rule, n)
    return best


def _check_one_to_one(pairs: dict[str, str], what: str) -> None:
    rev: dict[str, list[str]] = {}
    for r, d in pairs.items():
        rev.setdefault(d, []).append(r)
    many = {d: rs for d, rs in rev.items() if len(rs) > 1}
    if many:
        ex = "; ".join(f"{d!r} <- {examples(rs, 3)}" for d, rs in list(many.items())[:3])
        raise InputError(
            f"{what}: {len(many)} dataset trait(s) receive several reference traits, e.g. {ex}. "
            "Use a stricter trait_id_rule or an explicit trait_map file."
        )


def map_traits_by_rule(ref_ids: Sequence[str], ds_ids: Sequence[str], rule: str) -> TraitMapping:
    """Map reference traits onto dataset traits by comparing rule keys.

    Raises :class:`InputError` when a key is ambiguous: several dataset traits share
    the key of a reference trait, or several reference traits land on one dataset trait.
    """
    idx = key_index(ds_ids, rule)
    res = TraitMapping(rule=rule)
    ambiguous: list[tuple[str, list[str]]] = []
    for r in dict.fromkeys(ref_ids):
        cands = idx.get(trait_key(r, rule))
        if not cands:
            res.unmapped.append(r)
        elif len(cands) > 1:
            ambiguous.append((r, cands))
        else:
            res.mapping[r] = cands[0]
    if ambiguous:
        ex = "; ".join(f"{r!r} -> {examples(c, 3)}" for r, c in ambiguous[:3])
        raise InputError(
            f"trait_id_rule {rule}: {len(ambiguous)} reference trait(s) match several dataset "
            f"traits, e.g. {ex}. Use a stricter rule (e.g. as_is) or an explicit trait_map file."
        )
    _check_one_to_one(res.mapping, f"trait_id_rule {rule}")
    return res


def read_trait_map(path: str | os.PathLike[str]) -> list[tuple[str, str]]:
    """Read a two-column ``ref_trait pheno_trait`` file (no header; a header line whose
    fields are ``ref_trait``/``pheno_trait`` is tolerated)."""
    rows: list[tuple[str, str]] = []
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            f = s.split("\t") if "\t" in s else s.split()
            if len(f) < 2:
                raise InputError(f"{path}:{lineno}: expected two columns 'ref_trait pheno_trait'")
            if (
                lineno == 1
                and f[0].lower() in ("ref_trait", "ref_trait_id")
                and f[1].lower().startswith("pheno")
            ):
                continue
            rows.append((f[0], f[1]))
    return rows


def map_traits_by_file(
    ref_ids: Sequence[str],
    ds_ids: Sequence[str],
    rows: Sequence[tuple[str, str]],
    path: str = "trait_map",
) -> TraitMapping:
    """Map reference traits with explicit ``(ref_trait, pheno_trait)`` rows, keeping rows
    whose IDs exist on both sides. A reference trait listed with several dataset traits
    (or vice versa) is an error."""
    ref_set, ds_set = set(ref_ids), set(ds_ids)
    res = TraitMapping(rule="file")
    multi: dict[str, set[str]] = {}
    for r, d in rows:
        if r in ref_set and d in ds_set:
            multi.setdefault(r, set()).add(d)
    bad = {r: sorted(ds) for r, ds in multi.items() if len(ds) > 1}
    if bad:
        ex = "; ".join(f"{r!r} -> {examples(d, 3)}" for r, d in list(bad.items())[:3])
        raise InputError(
            f"{path}: {len(bad)} reference trait(s) map to several dataset traits, e.g. {ex}"
        )
    for r in dict.fromkeys(ref_ids):
        if r in multi:
            res.mapping[r] = next(iter(multi[r]))
        else:
            res.unmapped.append(r)
    _check_one_to_one(res.mapping, str(path))
    return res


# ------------------------------------------------------------------ sample maps
_GENO_NAMES = ("geno_id", "nwd_id", "genotype_id", "geno")
_PHENO_NAMES = ("pheno_id", "sample_id", "phenotype_id", "pheno")


def read_sample_map(path: str | os.PathLike[str]) -> list[tuple[str, str]]:
    """Read ``(geno_id, pheno_id)`` rows from a sample map with a header.

    Recognised headers: ``geno_id pheno_id`` and the older pipeline's ``NWD_ID
    SAMPLE_ID`` (genotype ID, omics sample ID), in either column order. Exact duplicate
    rows are dropped.
    """
    rows: list[tuple[str, str]] = []
    gi = pi = -1
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            s = line.rstrip("\n").rstrip("\r")
            if not s.strip():
                continue
            f = s.split("\t") if "\t" in s else s.split()
            if gi < 0:
                names = [x.lstrip("#").lower() for x in f]
                gi = next((names.index(n) for n in _GENO_NAMES if n in names), -1)
                pi = next((names.index(n) for n in _PHENO_NAMES if n in names), -1)
                if gi < 0 or pi < 0:
                    raise InputError(
                        f"{path}: expected a header with columns 'geno_id pheno_id' "
                        f"(or 'NWD_ID SAMPLE_ID'), got {s!r}"
                    )
                continue
            if len(f) <= max(gi, pi):
                raise InputError(f"{path}:{lineno}: too few fields")
            rows.append((f[gi], f[pi]))
    if gi < 0:
        raise InputError(f"{path}: empty sample map")
    return list(dict.fromkeys(rows))


@dataclass
class SampleMapResult:
    """Sample map restricted to IDs present in both the genotypes and the phenotypes."""

    pairs: list[tuple[str, str]]
    n_rows: int
    n_rows_geno_missing: int
    n_rows_pheno_missing: int
    n_geno_with_several_pheno: int


def filter_sample_map(
    rows: Sequence[tuple[str, str]], geno_ids: Iterable[str], pheno_ids: Sequence[str]
) -> SampleMapResult:
    """Keep rows whose genotype ID is genotyped and whose phenotype ID is phenotyped.

    A phenotype sample mapped to more than one genotype is an error. Pairs are ordered
    like ``pheno_ids``.
    """
    gset = set(geno_ids)
    pset = set(pheno_ids)
    kept: dict[str, list[str]] = {}
    n_g = n_p = 0
    for g, p in rows:
        if g not in gset:
            n_g += 1
            continue
        if p not in pset:
            n_p += 1
            continue
        kept.setdefault(p, []).append(g)
    many = {p: gs for p, gs in kept.items() if len(gs) > 1}
    if many:
        ex = "; ".join(f"{p!r} -> {examples(gs, 3)}" for p, gs in list(many.items())[:3])
        raise InputError(
            f"sample map: {len(many)} phenotype sample(s) map to several genotyped samples, e.g. {ex}"
        )
    pairs = [(kept[p][0], p) for p in pheno_ids if p in kept]
    per_geno: dict[str, int] = {}
    for g, _ in pairs:
        per_geno[g] = per_geno.get(g, 0) + 1
    return SampleMapResult(
        pairs=pairs,
        n_rows=len(rows),
        n_rows_geno_missing=n_g,
        n_rows_pheno_missing=n_p,
        n_geno_with_several_pheno=sum(1 for v in per_geno.values() if v > 1),
    )
