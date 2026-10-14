"""``match.*`` step implementations: run ``match-prs-pheno`` and summarise its calls.

Two step kinds live here:

* ``match_qpgen`` (step ``match.<ds>.qpgen``) restricts the trait and sample maps to what
  the PRS matrix and the prepared phenotypes contain, decides λ, runs
  ``qpgentools match-prs-pheno`` and records the resolved parameters.
* ``summarize_match`` (step ``match.<ds>.summarize``) turns qpgen's per-sample calls
  into ``samples.tsv``, ``swaps.tsv``, ``candidates.tsv``, ``duplicates.tsv``, ``related.tsv``,
  ``weights.tsv`` and ``summary.json``.

Status rules
------------
:func:`classify` is a line-by-line port of qpgen's ``MatchStatus`` in
``match.all.tsv.gz`` (``cmd_match_prs_pheno.cpp``; the methods page,
section 3). With ``Z(1) >= ... >= Z(5)`` the top five Z-scores, ``d`` = ``z_diff`` and
``t`` = ``z_threshold``, and checked in this order:

1. ``NO_OBS_TRAITS``: no observed trait with non-zero weight (``N.Traits == 0``).
2. ``SELF_BEST``: mapped and the labelled genotype ranks first (``Rank.self == 1``).
3. ``SINGLE_NEW_BEST``: ``Z(1) > Z(2) + d``.
4. ``MULTI_NEW_BEST``: ``Z(k) > Z(k+1) + d`` for some k in 2..4.
5. ``SELF_LENIENT``: mapped and ``Z.self > t`` (strictly).
6. ``UNCLEAR`` when mapped, ``NO_SELF`` otherwise.

Missing top-k Z-scores (fewer than five genotyped samples) count as -9999, as in qpgen.

Groups (:func:`group_of`): ``PASS`` = ``SELF_BEST``; ``WARN`` = ``SELF_LENIENT``,
``UNCLEAR``; ``FAIL`` = ``SINGLE_NEW_BEST``, ``MULTI_NEW_BEST`` of a *mapped* sample;
``NO_TEST`` = ``NO_SELF``, ``NO_OBS_TRAITS`` and every unmapped sample. An unmapped
sample with ``SINGLE_NEW_BEST`` is a *rescue candidate*: it stays ``NO_TEST`` (there is
no label to confirm or reject) and is listed in ``candidates.tsv`` with its best
genotype as ``proposed_geno_id``. A mapped ``FAIL`` sample gets ``proposed_geno_id`` =
its best genotype when ``Z(1) - Z(2) >= z_diff``. Neither gets one when that genotype is
the confirmed (``PASS``) match of another phenotyped sample, since two samples would then
map to one genotype; the note says so (:func:`proposal_conflict`), and such an unmapped
sample is not a rescue candidate.

``duplicates.tsv`` lists samples that *match two genotype IDs about equally*
(:func:`is_duplicate`): the top two candidates stand clear of the third by ``z_diff``
(qpgen's ``MULTI_NEW_BEST`` test at k = 2), the second scores above ``z_threshold``,
and the second scores at least :data:`DUPLICATE_RATIO` of the first. The two genotypes
then most likely come from one person (a duplicated genotype sample or identical twins).
Relatives do not qualify: a relative's Z is roughly half the self Z. A duplicate below
the clear-lead power cannot be told from noise and is not flagged. A sample listed here
is not listed in ``related.tsv``.

``related.tsv`` lists *possible relatives in the genotype set* (1000 Genomes holds trios):
samples whose top five start with a leading block of k = 2..4 candidates (the best
included) that stands out from the rest. The block ends at the largest drop
``Z(k) - Z(k+1)`` for k in 2..4 (a candidate k+1 must exist); it counts when that drop is at
least :data:`RELATED_GAP` and ``Z(k) >= z_threshold``. A relative's Z is roughly half the
self Z, so "within z_diff of the best" does not describe it, and with thousands of
genotypes the background candidates reach Z ~ 3 (MAGE against the full 1000 Genomes
genotype set: ``Z.3rd`` has median 3.2; against the 729 genotypes of the example data, 2.7),
so a Z threshold alone does not either; the drop to the background does.

The report (:mod:`anqchor.report.match`) recomputes these rules in the browser with a
JavaScript port; :data:`CLASSIFY_TEST_VECTORS` pins both implementations.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from collections.abc import Collection, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from anqchor.io.tables import examples, read_lines, write_lines
from anqchor.pipelines.layout import Layout
from anqchor.qpgen import commands as qc
from anqchor.qpgen.outputs import (
    MATCH_ORDINALS,
    QpgenOutputError,
    parse_auto_lambda,
    read_match_all,
    read_match_weights,
    read_prs_traits,
)
from anqchor.steps import StepContext, StepError, register
from anqchor.util.files import atomic_write, read_header, read_json, write_json

SCHEMA_VERSION = 1
MISSING_Z = -9999.0
"""qpgen's placeholder for an absent top-k candidate."""
FALLBACK_LAMBDA = 0.5
"""λ used with ``mahalanobis`` and ``lambda: auto`` when too few samples are mapped."""
RELATED_GAP = 1.0
"""Minimum Z drop between the leading block of the top five and the next candidate for
:func:`close_candidates` (possible relatives)."""
DUPLICATE_RATIO = 0.8
""":func:`is_duplicate`: the second candidate's Z must be at least this fraction of the
best's. Observed on MAGE: relatives 0.47-0.71 (one outlier at 0.98 fails the clear-lead
test), duplicated genotypes 0.96-1.00, every sample with Z(1) > 4 at most 0.85."""
QPGEN_MIN_BEST_DEFAULT = 100
"""qpgen's default ``--auto-lambda-min-best`` (and ``--auto-lambda-min-self``)."""

STATUS_ORDER = (
    "SELF_BEST",
    "SELF_LENIENT",
    "UNCLEAR",
    "SINGLE_NEW_BEST",
    "MULTI_NEW_BEST",
    "NO_SELF",
    "NO_OBS_TRAITS",
)
GROUP_ORDER = ("PASS", "WARN", "FAIL", "NO_TEST")
STATUS_GROUP = {
    "SELF_BEST": "PASS",
    "SELF_LENIENT": "WARN",
    "UNCLEAR": "WARN",
    "SINGLE_NEW_BEST": "FAIL",
    "MULTI_NEW_BEST": "FAIL",
    "NO_SELF": "NO_TEST",
    "NO_OBS_TRAITS": "NO_TEST",
}
NOT_SCORED = "NOT_SCORED"
"""anqchor status of a phenotyped sample missing from qpgen's output (no covariates)."""

SAMPLE_COLUMNS = [
    "pheno_id",
    "geno_id_assigned",
    "qpgen_status",
    "group",
    "z_self",
    "cor_self",
    "rank_self",
    "n_traits",
    "best_id",
    "z_best",
    "cor_best",
    "second_id",
    "z_second",
    "z_gap",
    "proposed_geno_id",
    "note",
    *[f"id_{k}" for k in range(1, 6)],
    *[f"z_{k}" for k in range(1, 6)],
]
SWAP_COLUMNS = [
    "kind",
    "pheno_id",
    "geno_id_assigned",
    "group",
    "best_id",
    "z_best",
    "z_self",
    "rank_self",
    "partner_pheno_id",
    "partner_geno_id_assigned",
    "partner_group",
    "partner_best_id",
]
CANDIDATE_COLUMNS = [
    "pheno_id",
    "candidate_geno_id",
    "z_best",
    "second_id",
    "z_second",
    "z_gap",
    "n_traits",
    "candidate_assigned_to",
]
DUPLICATE_COLUMNS = [
    "pheno_id",
    "geno_id_assigned",
    "group",
    "geno_id_1",
    "z_1",
    "geno_id_2",
    "z_2",
    "z_3",
    "assigned_in_pair",
]
RELATED_COLUMNS = [
    "pheno_id",
    "geno_id_assigned",
    "group",
    "n_close",
    "close_ids",
    "close_z",
    "assigned_in_close",
]
WEIGHT_COLUMNS = ["trait", "ref_trait", "weight", "n_obs", "weight_le0"]
SUMMARY_FILES = (
    "samples.tsv",
    "swaps.tsv",
    "candidates.tsv",
    "duplicates.tsv",
    "related.tsv",
    "weights.tsv",
    "summary.json",
)


# =========================================================================== rules
def qpgen_threads(p: dict[str, Any]) -> int:
    """Threads for ``match-prs-pheno``: ``$ANQCHOR_THREADS``, else ``params['threads']``
    (older signatures), else 1."""
    for v in (os.environ.get(qc.THREADS_ENV), p.get("threads")):
        try:
            n = int(v)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            continue
        if n >= 1:
            return n
    return 1


def _z(v: object) -> float:
    """A Z-score with missing values (None/NaN/NA) mapped to qpgen's -9999."""
    if v is None or v is pd.NA:
        return MISSING_Z
    try:
        f = float(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return MISSING_Z
    return MISSING_Z if math.isnan(f) else f


def classify(
    *,
    mapped: bool,
    n_traits: float | None,
    rank_self: float | None,
    z_self: float | None,
    z_top: Sequence[float | None],
    z_threshold: float,
    z_diff: float,
) -> str:
    """qpgen's ``match.all`` status for one phenotyped sample (see the module docstring).

    Args:
        mapped: the sample has a labelled genotype (``ID.self`` is not NA).
        n_traits: observed traits with non-zero weight (``N.Traits``).
        rank_self: ``Rank.self`` (ignored when not mapped).
        z_self: ``Z.self`` (ignored when not mapped).
        z_top: ``Z.1st`` ... ``Z.5th``; missing entries count as -9999.
        z_threshold: ``--z-threshold``.
        z_diff: ``--z-diff``.
    """
    if n_traits is None or (isinstance(n_traits, float) and math.isnan(n_traits)) or n_traits == 0:
        return "NO_OBS_TRAITS"
    z = [_z(v) for v in list(z_top)[:5]]
    z += [MISSING_Z] * (5 - len(z))
    if (
        mapped
        and rank_self is not None
        and not (isinstance(rank_self, float) and math.isnan(rank_self))
        and int(rank_self) == 1
    ):
        return "SELF_BEST"
    if z[0] > z[1] + z_diff:
        return "SINGLE_NEW_BEST"
    for k in range(2, 5):
        if z[k - 1] > z[k] + z_diff:
            return "MULTI_NEW_BEST"
    if mapped and _z(z_self) > z_threshold:
        return "SELF_LENIENT"
    return "UNCLEAR" if mapped else "NO_SELF"


def group_of(status: str, mapped: bool) -> str:
    """Coarse verdict group of a status; unmapped samples are always ``NO_TEST``."""
    if not mapped or status == NOT_SCORED:
        return "NO_TEST"
    return STATUS_GROUP.get(status, "NO_TEST")


def is_rescue(status: str, mapped: bool) -> bool:
    """An unmapped sample with one clearly leading genotype."""
    return (not mapped) and status == "SINGLE_NEW_BEST"


def close_candidates(
    ids: Sequence[str | None],
    zs: Sequence[float | None],
    z_threshold: float,
    gap: float = RELATED_GAP,
) -> list[tuple[str, float]]:
    """The leading block of the top five (best included) that stands out from the rest.

    The block is the first k candidates, k in 2..4, where ``Z(k) - Z(k+1)`` is largest
    (first k on ties; candidate k+1 must exist). Returns it when that drop is ``>= gap``
    and ``Z(k) >= z_threshold``, else an empty list.
    """
    pairs = [
        (i, _z(z)) for i, z in zip(ids, zs, strict=False) if i is not None and _z(z) > MISSING_Z
    ]
    best_k, best_drop = 0, -math.inf
    for k in range(2, min(len(pairs), 5)):
        drop = pairs[k - 1][1] - pairs[k][1]
        if drop > best_drop:
            best_k, best_drop = k, drop
    if best_k == 0 or best_drop < gap or pairs[best_k - 1][1] < z_threshold:
        return []
    return pairs[:best_k]


def is_related(n_close: int) -> bool:
    """Possible relatives: a leading block of at least two candidates."""
    return n_close >= 2


def two_genotype_ids(
    assigned: str,
    ids: Sequence[str | None],
    z_top: Sequence[float | None],
    z_threshold: float,
    z_diff: float,
) -> bool:
    """The sample matches two genotype IDs about equally (:func:`is_duplicate`) and has no
    assigned ID or its assigned ID is one of the two: the two genotypes are probably one
    person, so the sample is not a mismatch and gets no proposed ID."""
    if not is_duplicate(z_top, z_threshold, z_diff):
        return False
    return not assigned or assigned in (ids[0], ids[1] if len(ids) > 1 else None)


def is_duplicate(
    z_top: Sequence[float | None], z_threshold: float, z_diff: float, ratio: float = DUPLICATE_RATIO
) -> bool:
    """The sample matches two genotype IDs about equally (see the module docstring):
    ``Z(2) > Z(3) + z_diff``, ``Z(2) > z_threshold`` and ``Z(2) >= ratio * Z(1)``.
    Missing candidates count as -9999, as in :func:`classify`."""
    z = [_z(v) for v in list(z_top)[:3]]
    z += [MISSING_Z] * (3 - len(z))
    z1, z2, z3 = z
    if z1 <= 0 or z2 <= MISSING_Z:
        return False
    return z2 > z3 + z_diff and z2 > z_threshold and z2 >= ratio * z1


CLASSIFY_TEST_VECTORS: list[dict[str, Any]] = [
    # each: kwargs of classify() and the expected status; also used by the JS port's test
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 1,
        "z_self": 5.0,
        "z_top": [5.0, 1.0, 0.5, 0.4, 0.3],
        "expected": "SELF_BEST",
    },
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 1,
        "z_self": 1.0,
        "z_top": [1.0, 0.9, 0.8, 0.7, 0.6],
        "expected": "SELF_BEST",
    },
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 76,
        "z_self": 2.0,
        "z_top": [5.98, 3.5, 3.4, 3.3, 3.2],
        "expected": "SINGLE_NEW_BEST",
    },
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 2,
        "z_self": 4.0,
        "z_top": [4.1, 4.0, 1.0, 0.9, 0.8],
        "expected": "MULTI_NEW_BEST",
    },
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 2,
        "z_self": 3.0,
        "z_top": [3.2, 3.0, 2.9, 2.8, 2.7],
        "expected": "SELF_LENIENT",
    },
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 3,
        "z_self": 1.96,
        "z_top": [2.5, 2.2, 1.96, 1.9, 1.8],
        "expected": "UNCLEAR",
    },
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 2,
        "z_self": 1.0,
        "z_top": [1.5, 1.0, 0.9, 0.8, 0.7],
        "expected": "UNCLEAR",
    },
    # strict inequalities: a gap of exactly z_diff is not "clear"
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 2,
        "z_self": 3.0,
        "z_top": [5.0, 3.0, 2.9, 2.8, 2.7],
        "expected": "SELF_LENIENT",
    },
    {
        "mapped": False,
        "n_traits": 10,
        "rank_self": None,
        "z_self": None,
        "z_top": [5.0, 3.0, 2.9, 2.8, 2.7],
        "expected": "NO_SELF",
    },
    {
        "mapped": False,
        "n_traits": 10,
        "rank_self": None,
        "z_self": None,
        "z_top": [5.01, 3.0, 2.9, 2.8, 2.7],
        "expected": "SINGLE_NEW_BEST",
    },
    {
        "mapped": False,
        "n_traits": 10,
        "rank_self": None,
        "z_self": None,
        "z_top": [5.0, 4.9, 4.8, 2.7, 2.6],
        "expected": "MULTI_NEW_BEST",
    },
    {
        "mapped": False,
        "n_traits": 10,
        "rank_self": None,
        "z_self": None,
        "z_top": [5.0, 4.9, 4.8, 4.7, 2.6],
        "expected": "MULTI_NEW_BEST",
    },
    {
        "mapped": False,
        "n_traits": 10,
        "rank_self": None,
        "z_self": None,
        "z_top": [3.0, 2.9, 2.8, 2.7, 2.6],
        "expected": "NO_SELF",
    },
    # MULTI_NEW_BEST takes precedence over SELF_LENIENT (qpgen checks it last and overrides)
    {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 3,
        "z_self": 4.0,
        "z_top": [4.5, 4.2, 4.0, 1.0, 0.5],
        "expected": "MULTI_NEW_BEST",
    },
    {
        "mapped": True,
        "n_traits": 0,
        "rank_self": None,
        "z_self": None,
        "z_top": [None] * 5,
        "expected": "NO_OBS_TRAITS",
    },
    {
        "mapped": False,
        "n_traits": 0,
        "rank_self": None,
        "z_self": None,
        "z_top": [None] * 5,
        "expected": "NO_OBS_TRAITS",
    },
    # fewer than five genotyped samples: missing Z count as -9999
    {
        "mapped": False,
        "n_traits": 3,
        "rank_self": None,
        "z_self": None,
        "z_top": [0.7, None, None, None, None],
        "expected": "SINGLE_NEW_BEST",
    },
    {
        "mapped": True,
        "n_traits": 3,
        "rank_self": 2,
        "z_self": 0.5,
        "z_top": [0.7, 0.5, None, None, None],
        "expected": "MULTI_NEW_BEST",
    },
    # threshold on Z.self is strict
    {
        "mapped": True,
        "n_traits": 5,
        "rank_self": 2,
        "z_self": 1.96,
        "z_top": [2.0, 1.96, 1.5, 1.4, 1.3],
        "expected": "UNCLEAR",
    },
    {
        "mapped": True,
        "n_traits": 5,
        "rank_self": 2,
        "z_self": 1.97,
        "z_top": [2.0, 1.97, 1.5, 1.4, 1.3],
        "expected": "SELF_LENIENT",
    },
    # negative Z values
    {
        "mapped": True,
        "n_traits": 5,
        "rank_self": 900,
        "z_self": -1.0,
        "z_top": [1.0, 0.9, 0.8, 0.7, 0.6],
        "expected": "UNCLEAR",
    },
]
"""Edge cases at ``z_threshold=1.96, z_diff=2.0`` shared by the Python and JS tests."""


# ===================================================================== qpgen step
def _restrict_trait_map(
    pairs: list[tuple[str, str]], prs_traits: Sequence[str], pheno_traits: Sequence[str]
) -> tuple[list[tuple[str, str]], dict[str, int]]:
    """Keep pairs whose reference trait has a PRS and whose dataset trait is in the
    phenotype matrix; each side used once (first occurrence wins)."""
    prs_set, ph_set = set(prs_traits), set(pheno_traits)
    seen_ref: set[str] = set()
    seen_ph: set[str] = set()
    kept: list[tuple[str, str]] = []
    n = Counter[str]()
    for r, d in pairs:
        if r not in prs_set:
            n["no_prs"] += 1
        elif d not in ph_set:
            n["not_in_pheno"] += 1
        elif r in seen_ref or d in seen_ph:
            n["duplicate"] += 1
        else:
            kept.append((r, d))
            seen_ref.add(r)
            seen_ph.add(d)
    return kept, {
        "n_trait_map_rows": len(pairs),
        "n_dropped_no_prs": n["no_prs"],
        "n_dropped_not_in_pheno": n["not_in_pheno"],
        "n_dropped_duplicate": n["duplicate"],
    }


def _read_ids(path: Path, column: str = "IID") -> list[str]:
    """One ID column of a regenie-layout matrix, without parsing the values."""
    return [str(x) for x in pd.read_csv(path, sep="\t", usecols=[column], dtype=str)[column]]


def _read_two_columns(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for line in read_lines(path):
        f = line.split("\t") if "\t" in line else line.split()
        if len(f) >= 2:
            out.append((f[0], f[1]))
    return out


def _log_excerpt(lines: Sequence[str], limit: int = 60) -> list[str]:
    """Informative qpgen log lines (no option banner, no progress noise)."""
    keep = []
    for line in lines:
        s = line.rstrip()
        if not s or "Expanding phenotype matrix" in s or "Standardizing columns of matrix" in s:
            continue
        if s.startswith(("NOTICE", "WARNING", "ERROR")) or "error" in s.lower():
            keep.append(s)
    return keep[-limit:]


def _run_streaming(argv: Sequence[str], cwd: Path) -> tuple[int, list[str]]:
    """Run ``argv``; echo every output line to stdout (the step log) and keep a copy."""
    print("+ " + qc.to_shell(argv), flush=True)
    try:
        proc = subprocess.Popen(
            [str(a) for a in argv],
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            bufsize=1,
        )
    except OSError as e:
        raise StepError(
            f"cannot run qpgentools ({argv[0]}): {e}; check 'anqchor qpgen check'"
        ) from e
    lines: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(line)
        lines.append(line)
    sys.stdout.flush()
    return proc.wait(), lines


def _failure_message(rc: int, lines: Sequence[str], p: dict[str, Any]) -> str:
    text = "".join(lines)
    raw = [s.rstrip() for s in lines if s.strip() and "auto-lambda evaluation" not in s]
    tail = "\n".join(raw[-6:])
    if "rank their own PRS first at the tuned lambda" in text:
        return (
            "match-prs-pheno --auto-lambda stopped: too few mapped samples rank their own PRS first, "
            "which almost always means the sample mapping is wrong (IDs that do not correspond). "
            "Check datasets.<ds>.sample_map and the sample IDs; to proceed anyway set "
            "pipelines.match.lambda to a number (e.g. 0.5) or use the default independent scoring "
            f"(mahalanobis: false).\n{tail}"
        )
    if "fewer than --auto-lambda-min-self" in text:
        return (
            "match-prs-pheno --auto-lambda: too few mapped samples with observed traits to tune "
            "lambda; set a numeric pipelines.match.lambda (e.g. 0.5) or check the sample mapping."
            f"\n{tail}"
        )
    return f"qpgentools match-prs-pheno failed (exit {rc}) for dataset {p.get('dataset')}:\n{tail}"


@register("match_qpgen")
def match_qpgen(ctx: StepContext) -> None:
    """Run ``match-prs-pheno`` for one dataset (see the module docstring)."""
    p = ctx.params
    ds, ref = p["dataset"], p["ref"]
    L = Layout
    t0 = time.time()
    warnings: list[str] = []
    prs = L.prs_matrix(ref)
    pheno = L.ds_pheno(ds)
    cov = L.ds_cov(ds) if L.ds_cov(ds) in ctx.inputs else None
    ds_info = read_json(ctx.path(L.ds_json(ds)))

    # ---- traits: reference traits with a PRS × dataset traits in the matrix
    prs_traits = read_prs_traits(ctx.path(prs))
    pheno_traits = read_header(ctx.path(pheno))[2:]
    tmap = _read_two_columns(ctx.path(L.dsref_trait_map(ds, ref)))
    used, tcounts = _restrict_trait_map(tmap, prs_traits, pheno_traits)
    if not used:
        raise StepError(
            f"dataset {ds}: none of the {len(tmap)} mapped traits of reference {ref} has a PRS "
            f"and a phenotype column (PRS traits look like {examples(prs_traits, 3)}, dataset traits "
            f"like {examples(pheno_traits, 3)}). Were all reference variants absent from the genotypes?"
        )
    write_lines(ctx.path(L.match_trait_map(ds)), (f"{r}\t{d}" for r, d in used))
    identity = all(r == d for r, d in used) and {r for r, _ in used} == set(prs_traits) & set(
        pheno_traits
    )
    trait_tsv = None if identity else L.match_trait_map(ds)

    # ---- samples: the map qpgen sees must only name IDs present in both matrices
    prs_ids = set(_read_ids(ctx.path(prs)))
    pheno_ids = read_lines(ctx.path(L.ds_samples(ds)))
    eligible = set(pheno_ids)
    if cov is not None:
        eligible &= set(_read_ids(ctx.path(cov)))
    sample_tsv: Path | None = None
    n_map_dropped = 0
    if p.get("has_sample_map"):
        rows = _read_two_columns(ctx.path(L.ds_sample_map(ds)))
        keep = [(g, ph) for g, ph in rows if g in prs_ids and ph in eligible]
        n_map_dropped = len(rows) - len(keep)
        if n_map_dropped:
            sample_tsv = L.match_prefix(ds).with_name("qpgen.sample_map.used.tsv")
            write_lines(ctx.path(sample_tsv), (f"{g}\t{ph}" for g, ph in keep))
            warnings.append(
                f"{n_map_dropped} sample-map rows dropped (genotype without PRS or phenotype "
                "sample without covariates)"
            )
        else:
            sample_tsv = L.ds_sample_map(ds)
        n_mapped = len(keep)
    else:
        n_mapped = sum(1 for s in pheno_ids if s in eligible and s in prs_ids)
    if n_mapped == 0:
        warnings.append(
            "no phenotyped sample maps to a genotyped sample; only rescue calls are possible"
        )

    # ---- lambda
    mahalanobis = bool(p.get("mahalanobis"))
    lam_cfg = p.get("lambda", "auto")
    lam: float | None = None
    auto = False
    min_self = int(p.get("min_self_for_auto", QPGEN_MIN_BEST_DEFAULT))
    if not mahalanobis:
        lambda_mode = "independent"
    elif lam_cfg == "auto":
        if n_mapped >= min_self:
            auto, lambda_mode = True, "auto"
        else:
            lam, lambda_mode = FALLBACK_LAMBDA, "fallback"
            msg = (
                f"lambda auto: only {n_mapped} mapped samples (< min_self_for_auto {min_self}); "
                f"using --lambda {FALLBACK_LAMBDA}"
            )
            warnings.append(msg)
            if ctx.log is not None:
                ctx.log.warning(msg)
    else:
        lam, lambda_mode = float(lam_cfg), "fixed"
    auto_min_self = auto_min_best = None
    if auto and min_self != QPGEN_MIN_BEST_DEFAULT:
        auto_min_self = min_self
        auto_min_best = min(QPGEN_MIN_BEST_DEFAULT, min_self)

    # ---- run qpgen into a temporary prefix, then move the outputs into place
    mdir = ctx.path(L.match_dir(ds))
    mdir.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(tempfile.mkdtemp(prefix=".qpgen.", dir=mdir))
    try:
        tmp_prefix = tmpdir / "qpgen"
        missing = list(dict.fromkeys(["NA", *(p.get("missing_str") or [])]))

        def build(out_prefix: Path) -> list[str]:
            return qc.match_prs_pheno(
                p["qpgen"],
                prs=prs,
                pheno=pheno,
                out_prefix=out_prefix,
                cov=cov,
                cat_cov_names=list(p.get("cat_cov_names") or []) if cov is not None else [],
                max_cat_levels=p.get("max_cat_levels") if cov is not None else None,
                sample_tsv=sample_tsv,
                trait_tsv=trait_tsv,
                missing_str=missing,
                rint=bool(p.get("rint", True)),
                mahalanobis=mahalanobis,
                lam=lam,
                auto_lambda=auto,
                min_weight=p.get("min_weight") or None,
                z_threshold=float(p["z_threshold"]),
                z_diff=float(p["z_diff"]),
                threads=qpgen_threads(p),
                auto_lambda_min_self=auto_min_self,
                auto_lambda_min_best=auto_min_best,
            )

        argv = build(tmp_prefix)
        rc, lines = _run_streaming(argv, ctx.run_dir)
        if rc != 0:
            raise StepError(_failure_message(rc, lines, p))
        finals = L.match_qpgen_outputs(ds)
        for final in finals:
            src = tmpdir / final.name
            if not src.exists():
                raise StepError(f"qpgentools match-prs-pheno did not write {final.name}")
            os.replace(src, ctx.path(final))
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    lambda_used: float | None
    if lambda_mode == "auto":
        lambda_used = parse_auto_lambda("".join(lines))
        if lambda_used is None:
            warnings.append("could not parse the tuned lambda from the qpgen log")
    elif lambda_mode == "independent":
        lambda_used = None
    else:
        lambda_used = lam
    write_json(
        ctx.path(L.match_params(ds)),
        {
            "schema_version": SCHEMA_VERSION,
            "dataset": ds,
            "reference": ref,
            "argv": [str(a) for a in build(L.match_prefix(ds))],
            "qpgen": str(p["qpgen"]),
            "mahalanobis": mahalanobis,
            "lambda_config": lam_cfg,
            "lambda_mode": lambda_mode,
            "lambda_used": lambda_used,
            "rint": bool(p.get("rint", True)),
            "use_cov": cov is not None,
            "z_threshold": float(p["z_threshold"]),
            "z_diff": float(p["z_diff"]),
            "min_weight": p.get("min_weight"),
            "n_traits_used": len(used),
            **tcounts,
            "n_prs_samples": len(prs_ids),
            "n_samples_pheno": len(pheno_ids),
            "n_samples_mapped": n_mapped,
            "n_samples_mapped_dataset": ds_info.get("n_samples_mapped"),
            "n_sample_map_rows_dropped": n_map_dropped,
            "trait_tsv_used": trait_tsv is not None,
            "sample_tsv_used": sample_tsv is not None,
            "warnings": warnings,
            "log_excerpt": _log_excerpt(lines),
            "seconds": round(time.time() - t0, 2),
        },
    )
    if ctx.log is not None:
        ctx.log.info(
            "match %s: %d traits, %d mapped samples, lambda %s (%s)",
            ds,
            len(used),
            n_mapped,
            lambda_used,
            lambda_mode,
        )


# ================================================================== summarize step
def _num(v: Any) -> float | None:
    if v is None or v is pd.NA:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def _str(v: Any) -> str:
    if v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v)):
        return ""
    return str(v)


def build_samples(
    match_all: pd.DataFrame,
    pheno_ids: Sequence[str],
    z_threshold: float,
    z_diff: float,
    assigned_ids: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, int]:
    """``samples.tsv`` rows from qpgen's ``match.all`` table.

    ``assigned_ids`` (phenotype ID to genotype ID, from the prepared sample map) gives
    the assigned genotype of the phenotyped samples qpgen did not score (``NOT_SCORED``,
    e.g. no covariates), so that they are not counted as samples without an ID.

    Returns the table and the number of rows where :func:`classify` disagrees with
    qpgen's status (expected 0; non-zero means the qpgen rules changed).
    """
    recs: list[dict[str, Any]] = []
    n_disagree = 0
    for r in match_all.to_dict("records"):
        assigned = _str(r.get("ID.self"))
        mapped = bool(assigned)
        status = str(r["MatchStatus"])
        ids = [_str(r.get(f"ID.{o}")) or None for o in MATCH_ORDINALS]
        zs = [_num(r.get(f"Z.{o}")) for o in MATCH_ORDINALS]
        mine = classify(
            mapped=mapped,
            n_traits=_num(r.get("N.Traits")),
            rank_self=_num(r.get("Rank.self")),
            z_self=_num(r.get("Z.self")),
            z_top=zs,
            z_threshold=z_threshold,
            z_diff=z_diff,
        )
        if mine != status:
            n_disagree += 1
        z1, z2 = zs[0], zs[1]
        gap = z1 - z2 if z1 is not None and z2 is not None else None
        group = group_of(status, mapped)
        proposed = ""
        if (
            (group == "FAIL" or is_rescue(status, mapped))
            and gap is not None
            and gap >= z_diff
            and ids[0]
        ):
            proposed = ids[0]
        if proposed and two_genotype_ids(assigned, ids, zs, z_threshold, z_diff):
            # one person under two genotype IDs: neither ID is a relabel
            proposed = ""
        rank = _num(r.get("Rank.self"))
        recs.append(
            {
                "pheno_id": str(r["ID.Pheno"]),
                "geno_id_assigned": assigned,
                "qpgen_status": status,
                "group": group,
                "z_self": _num(r.get("Z.self")),
                "cor_self": _num(r.get("COR.self")),
                "rank_self": int(rank) if rank is not None else None,
                "n_traits": int(_num(r.get("N.Traits")) or 0),
                "best_id": ids[0] or "",
                "z_best": z1,
                "cor_best": _num(r.get("COR.1st")),
                "second_id": ids[1] or "",
                "z_second": z2,
                "z_gap": gap,
                "proposed_geno_id": proposed,
                "note": "",
                **{f"id_{k + 1}": ids[k] or "" for k in range(5)},
                **{f"z_{k + 1}": zs[k] for k in range(5)},
            }
        )
    seen = {r["pheno_id"] for r in recs}
    for s in pheno_ids:
        if s not in seen:
            recs.append(
                {c: "" for c in SAMPLE_COLUMNS}
                | {
                    "pheno_id": s,
                    "geno_id_assigned": (assigned_ids or {}).get(s, ""),
                    "qpgen_status": NOT_SCORED,
                    "group": "NO_TEST",
                    "n_traits": 0,
                    "z_self": None,
                    "cor_self": None,
                    "rank_self": None,
                    "z_best": None,
                    "cor_best": None,
                    "z_second": None,
                    "z_gap": None,
                    **{f"z_{k}": None for k in range(1, 6)},
                    "note": "not in qpgen output (no covariates for this sample?)",
                }
            )
    df = pd.DataFrame.from_records(recs, columns=SAMPLE_COLUMNS)
    df["rank_self"] = df["rank_self"].astype("Int64")
    df["n_traits"] = df["n_traits"].astype("Int64")
    _annotate_notes(df, z_threshold, z_diff)
    _drop_conflicting_proposals(df)
    return df, n_disagree


def _assigned_index(df: pd.DataFrame) -> dict[str, list[int]]:
    idx: dict[str, list[int]] = {}
    for i, g in enumerate(df["geno_id_assigned"].tolist()):
        if g:
            idx.setdefault(g, []).append(i)
    return idx


def proposal_conflict(
    best_id: str, i: int, by_geno: dict[str, list[int]], groups: Sequence[str]
) -> bool:
    """Whether ``best_id`` is the confirmed (``PASS``) genotype of another phenotyped
    sample. Relabelling sample ``i`` to it would map two samples to one genotype: the
    pattern points to a duplicated or contaminated library, not a mislabel."""
    return any(groups[j] == "PASS" for j in by_geno.get(best_id, []) if j != i)


CONFLICT_NOTE = "no relabel proposed: best genotype is the confirmed match of another sample"


def _drop_conflicting_proposals(df: pd.DataFrame) -> None:
    """Clear ``proposed_geno_id`` where :func:`proposal_conflict` holds and say why."""
    by_geno = _assigned_index(df)
    groups = df["group"].tolist()
    prop = df["proposed_geno_id"].tolist()
    notes = df["note"].tolist()
    for i, best in enumerate(prop):
        if best and proposal_conflict(best, i, by_geno, groups):
            prop[i] = ""
            notes[i] = "; ".join([n for n in (notes[i], CONFLICT_NOTE) if n])
    df["proposed_geno_id"] = prop
    df["note"] = notes


def _annotate_notes(df: pd.DataFrame, z_threshold: float, z_diff: float) -> None:
    """Fill the ``note`` column: who else holds the best genotype, rescue, two genotype IDs,
    relatives."""
    by_geno = _assigned_index(df)
    notes = df["note"].tolist()
    pheno = df["pheno_id"].tolist()
    groups = df["group"].tolist()
    for i, row in enumerate(df.itertuples(index=False)):
        parts = [notes[i]] if notes[i] else []
        status, mapped = row.qpgen_status, bool(row.geno_id_assigned)
        if (
            is_rescue(status, mapped)
            and row.proposed_geno_id
            and not proposal_conflict(row.best_id, i, by_geno, groups)
        ):
            parts.append(f"rescue candidate: {row.best_id}")
        if row.group == "FAIL" and row.best_id:
            holders = [j for j in by_geno.get(row.best_id, []) if j != i]
            if holders:
                parts.append(
                    "best genotype is assigned to "
                    + ", ".join(f"{pheno[j]} ({groups[j]})" for j in holders)
                )
            else:
                parts.append("best genotype is not assigned to any phenotyped sample")
        zs = [getattr(row, f"z_{k}") for k in range(1, 6)]
        dup = status != NOT_SCORED and is_duplicate(zs, z_threshold, z_diff)
        if dup:
            parts.append(
                f"matches two genotype IDs about equally: {row.id_1}, {row.id_2} (same person?)"
            )
        elif status == "MULTI_NEW_BEST":
            parts.append("several genotypes lead about equally (mixed sample?)")
        close = close_candidates(
            [getattr(row, f"id_{k}") or None for k in range(1, 6)], zs, z_threshold
        )
        if not dup and is_related(len(close)):
            parts.append("possible relatives: " + ", ".join(i_ for i_, _ in close))
        notes[i] = "; ".join(parts)
    df["note"] = notes


def find_swaps(df: pd.DataFrame, two_ids: Collection[str] = ()) -> pd.DataFrame:
    """Reciprocal swaps (A's best is B's genotype and B's best is A's) and the other
    ``FAIL`` samples as ``one-way``.

    A reciprocal pair is a ``reciprocal`` swap only when at least one side is a ``FAIL``
    whose best genotype leads by at least ``z_diff`` (``SINGLE_NEW_BEST``). When neither side is (e.g. two
    ``SELF_LENIENT`` samples whose own genotype is second by a hair, typical of weak
    signal or relatives), the pair is listed as ``reciprocal-weak``: not significant
    evidence of a swap, and not counted in ``n_reciprocal_swaps``.

    ``two_ids``: phenotype IDs of samples that match two genotype IDs about equally, one of
    them their own (:func:`two_id_samples`); they are not listed as ``one-way``."""
    by_geno = _assigned_index(df)
    rows = df.to_dict("records")
    out: list[dict[str, Any]] = []
    in_pair: set[int] = set()
    for i, a in enumerate(rows):
        if not a["geno_id_assigned"] or not a["best_id"] or a["best_id"] == a["geno_id_assigned"]:
            continue
        for j in by_geno.get(a["best_id"], []):
            b = rows[j]
            if j > i and b["best_id"] == a["geno_id_assigned"]:
                in_pair.update((i, j))
                kind = (
                    "reciprocal" if _clear_mismatch(a) or _clear_mismatch(b) else "reciprocal-weak"
                )
                for x, y in ((a, b), (b, a)):
                    out.append(_swap_row(kind, x, y))
    for i, a in enumerate(rows):
        if a["group"] == "FAIL" and i not in in_pair and a["pheno_id"] not in two_ids:
            holders = [rows[j] for j in by_geno.get(a["best_id"], []) if j != i]
            out.append(_swap_row("one-way", a, holders[0] if holders else None))
    return pd.DataFrame.from_records(out, columns=SWAP_COLUMNS)


def _clear_mismatch(r: dict[str, Any]) -> bool:
    """A mapped sample whose best genotype is not its own and leads the second by at least
    ``z_diff`` (``SINGLE_NEW_BEST``, a ``FAIL``)."""
    return bool(r["group"] == "FAIL" and r["qpgen_status"] == "SINGLE_NEW_BEST")


def _swap_row(kind: str, a: dict[str, Any], b: dict[str, Any] | None) -> dict[str, Any]:
    return {
        "kind": kind,
        "pheno_id": a["pheno_id"],
        "geno_id_assigned": a["geno_id_assigned"],
        "group": a["group"],
        "best_id": a["best_id"],
        "z_best": a["z_best"],
        "z_self": a["z_self"],
        "rank_self": a["rank_self"],
        "partner_pheno_id": b["pheno_id"] if b else "",
        "partner_geno_id_assigned": b["geno_id_assigned"] if b else "",
        "partner_group": b["group"] if b else "",
        "partner_best_id": b["best_id"] if b else "",
    }


def find_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """Rescue candidates: unmapped samples with one clearly leading genotype that is not
    the confirmed match of another sample (:func:`proposal_conflict`; such a sample has no
    ``proposed_geno_id``)."""
    by_geno = _assigned_index(df)
    pheno = df["pheno_id"].tolist()
    out = []
    for r in df.to_dict("records"):
        if is_rescue(r["qpgen_status"], bool(r["geno_id_assigned"])) and r["proposed_geno_id"]:
            holders = by_geno.get(r["best_id"], [])
            out.append(
                {
                    "pheno_id": r["pheno_id"],
                    "candidate_geno_id": r["best_id"],
                    "z_best": r["z_best"],
                    "second_id": r["second_id"],
                    "z_second": r["z_second"],
                    "z_gap": r["z_gap"],
                    "n_traits": r["n_traits"],
                    "candidate_assigned_to": ",".join(pheno[j] for j in holders),
                }
            )
    return pd.DataFrame.from_records(out, columns=CANDIDATE_COLUMNS)


def find_duplicates(df: pd.DataFrame, z_threshold: float, z_diff: float) -> pd.DataFrame:
    """Samples that match two genotype IDs about equally (:func:`is_duplicate`)."""
    out = []
    for r in df.to_dict("records"):
        if r["qpgen_status"] == NOT_SCORED or not is_duplicate(
            [r["z_1"], r["z_2"], r["z_3"]], z_threshold, z_diff
        ):
            continue
        out.append(
            {
                "pheno_id": r["pheno_id"],
                "geno_id_assigned": r["geno_id_assigned"],
                "group": r["group"],
                "geno_id_1": r["id_1"],
                "z_1": r["z_1"],
                "geno_id_2": r["id_2"],
                "z_2": r["z_2"],
                "z_3": r["z_3"],
                "assigned_in_pair": bool(r["geno_id_assigned"])
                and r["geno_id_assigned"] in (r["id_1"], r["id_2"]),
            }
        )
    return pd.DataFrame.from_records(out, columns=DUPLICATE_COLUMNS)


def two_id_samples(duplicates: pd.DataFrame) -> set[str]:
    """Phenotype IDs of the :func:`find_duplicates` rows whose assigned ID is one of the
    two genotypes (``assigned_in_pair``): not mismatches."""
    if not len(duplicates):
        return set()
    return {
        str(p)
        for p, a in zip(duplicates["pheno_id"], duplicates["assigned_in_pair"], strict=True)
        if bool(a)
    }


def find_related(df: pd.DataFrame, z_threshold: float, z_diff: float | None = None) -> pd.DataFrame:
    """Samples whose top five start with a block of possible relatives
    (see :func:`close_candidates`); with ``z_diff``, samples that match two genotype IDs
    about equally (:func:`find_duplicates`) are left out."""
    out = []
    for r in df.to_dict("records"):
        ids = [r[f"id_{k}"] or None for k in range(1, 6)]
        zs = [r[f"z_{k}"] for k in range(1, 6)]
        if (
            z_diff is not None
            and r["qpgen_status"] != NOT_SCORED
            and is_duplicate(zs, z_threshold, z_diff)
        ):
            continue
        close = close_candidates(ids, zs, z_threshold)
        if is_related(len(close)):
            out.append(
                {
                    "pheno_id": r["pheno_id"],
                    "geno_id_assigned": r["geno_id_assigned"],
                    "group": r["group"],
                    "n_close": len(close),
                    "close_ids": ",".join(i for i, _ in close),
                    "close_z": ",".join(f"{z:.3f}" for _, z in close),
                    "assigned_in_close": bool(r["geno_id_assigned"])
                    and r["geno_id_assigned"] in {i for i, _ in close},
                }
            )
    return pd.DataFrame.from_records(out, columns=RELATED_COLUMNS)


def build_weights(weights: pd.DataFrame, trait_map: Sequence[tuple[str, str]]) -> pd.DataFrame:
    """qpgen weights (dataset trait IDs) with the reference trait they came from."""
    back = {d: r for r, d in trait_map}
    w = pd.to_numeric(weights["Weight"], errors="coerce")
    return pd.DataFrame(
        {
            "trait": weights["Trait"].astype(str),
            "ref_trait": [back.get(t, t) for t in weights["Trait"].astype(str)],
            "weight": w,
            "n_obs": pd.to_numeric(weights["N.Obs"], errors="coerce").astype("Int64"),
            "weight_le0": (w <= 0).astype(int),
        },
        columns=WEIGHT_COLUMNS,
    )


def fail_phrase(pheno_id: str, assigned: str, best: str) -> str:
    """``"A best matches B"``, naming the labelled genotype when it differs from the ID."""
    label = f" (labelled {assigned})" if assigned and assigned != pheno_id else ""
    return f"{pheno_id}{label} best matches {best}"


def make_verdict(
    df: pd.DataFrame, max_listed: int = 3, duplicates: pd.DataFrame | None = None
) -> tuple[str, str]:
    """One-sentence verdict and its status (good, warning, critical).

    ``duplicates`` (:func:`find_duplicates`) adds a part naming the samples that match two
    genotype IDs about equally; those whose assigned ID is one of the two
    (:func:`two_id_samples`) are not counted as mismatches. Only duplicates give a
    ``warning``."""
    mapped = df[df["geno_id_assigned"] != ""]
    n_mapped = len(mapped)
    st = Counter(mapped["qpgen_status"])
    confirmed = st["SELF_BEST"] + st["SELF_LENIENT"]
    fail = mismatches(df, duplicates)
    parts: list[str] = []
    if n_mapped:
        lead = f"{confirmed} of {n_mapped} mapped samples confirmed"
        if st["SELF_LENIENT"]:
            lead += f" ({st['SELF_LENIENT']} only leniently)"
        parts.append(lead)
    else:
        parts.append("no phenotyped sample is mapped to a genotype")
    if len(fail):
        listed = [
            fail_phrase(r.pheno_id, r.geno_id_assigned, r.best_id)
            for r in fail.head(max_listed).itertuples()
        ]
        more = f"; and {len(fail) - max_listed} more" if len(fail) > max_listed else ""
        parts.append(
            f"{len(fail)} mismatch{'es' if len(fail) != 1 else ''} ({'; '.join(listed)}{more})"
        )
    if st["UNCLEAR"]:
        parts.append(f"{st['UNCLEAR']} unclear")
    n_noobs = int((df["qpgen_status"] == "NO_OBS_TRAITS").sum())
    if n_noobs:
        parts.append(f"{n_noobs} without observed traits")
    n_not_scored = int((df["qpgen_status"] == NOT_SCORED).sum())
    if n_not_scored:
        parts.append(f"{n_not_scored} not scored (no covariates)")
    if duplicates is not None and len(duplicates):
        listed = [
            f"{r.pheno_id}: {r.geno_id_1}, {r.geno_id_2}"
            for r in duplicates.head(max_listed).itertuples()
        ]
        more = f"; and {len(duplicates) - max_listed} more" if len(duplicates) > max_listed else ""
        n_dup = len(duplicates)
        parts.append(
            f"{n_dup} sample{'s' if n_dup != 1 else ''} matching two genotype IDs about equally "
            f"({'; '.join(listed)}{more})"
        )
    n_unmapped = len(df) - n_mapped
    if n_unmapped:
        rescue = (df["geno_id_assigned"] == "") & (df["qpgen_status"] == "SINGLE_NEW_BEST")
        if "proposed_geno_id" in df.columns:
            rescue &= df["proposed_geno_id"] != ""
        n_rescue = int(rescue.sum())
        parts.append(
            f"{n_unmapped} unmapped sample{'s' if n_unmapped != 1 else ''}, "
            f"{n_rescue} rescue candidate{'s' if n_rescue != 1 else ''}"
        )
    n_dups = len(duplicates) if duplicates is not None else 0
    status = (
        "critical"
        if len(fail)
        else ("warning" if st["UNCLEAR"] or not n_mapped or n_dups else "good")
    )
    return "; ".join(parts), status


def mismatches(df: pd.DataFrame, duplicates: pd.DataFrame | None = None) -> pd.DataFrame:
    """Samples with an assigned ID in group ``FAIL``, without those that match two genotype
    IDs about equally, one of them their own (:func:`two_id_samples`)."""
    fail = df[(df["geno_id_assigned"] != "") & (df["group"] == "FAIL")]
    two = two_id_samples(duplicates) if duplicates is not None else set()
    return fail[~fail["pheno_id"].isin(two)] if two else fail


def _quantiles(x: pd.Series) -> dict[str, float | None]:
    v = pd.to_numeric(x, errors="coerce").dropna().to_numpy(dtype=float)
    if not len(v):
        return {k: None for k in ("min", "q05", "q25", "median", "mean", "q75", "q95", "max")}
    q = np.quantile(v, [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
    return {
        "min": float(q[0]),
        "q05": float(q[1]),
        "q25": float(q[2]),
        "median": float(q[3]),
        "mean": float(v.mean()),
        "q75": float(q[4]),
        "q95": float(q[5]),
        "max": float(q[6]),
    }


def _write_tsv(path: Path, df: pd.DataFrame) -> None:
    with atomic_write(path) as fh:
        df.to_csv(fh, sep="\t", index=False, na_rep="NA", float_format="%.6g")


@register("summarize_match")
def summarize_match(ctx: StepContext) -> None:
    """Write the per-dataset match summaries (see :data:`SUMMARY_FILES`)."""
    p = ctx.params
    ds = p["dataset"]
    L = Layout
    _, _, all_path = L.match_qpgen_outputs(ds)
    w_path = L.match_qpgen_outputs(ds)[0]
    params = read_json(ctx.path(L.match_params(ds)))
    ds_info = read_json(ctx.path(L.ds_json(ds)))
    zt, zd = float(params["z_threshold"]), float(params["z_diff"])
    try:
        match_all = read_match_all(ctx.path(all_path))
        weights_raw = read_match_weights(ctx.path(w_path))
    except QpgenOutputError as e:
        raise StepError(str(e)) from e
    pheno_ids = read_lines(ctx.path(L.ds_samples(ds))) if L.ds_samples(ds) in ctx.inputs else []
    assigned_ids: dict[str, str] = {}
    if L.ds_sample_map(ds) in ctx.inputs and ctx.path(L.ds_sample_map(ds)).exists():
        assigned_ids = {ph: g for g, ph in _read_two_columns(ctx.path(L.ds_sample_map(ds)))}
    samples, n_disagree = build_samples(match_all, pheno_ids, zt, zd, assigned_ids)
    if n_disagree and ctx.log is not None:
        ctx.log.warning("%d qpgen statuses differ from anqchor's port of the rules", n_disagree)
    duplicates = find_duplicates(samples, zt, zd)
    swaps = find_swaps(samples, two_id_samples(duplicates))
    cands = find_candidates(samples)
    related = find_related(samples, zt, zd)
    tmap = (
        _read_two_columns(ctx.path(L.match_trait_map(ds)))
        if L.match_trait_map(ds) in ctx.inputs
        else []
    )
    weights = build_weights(weights_raw, tmap)

    _write_tsv(ctx.path(L.match_summary(ds, "samples.tsv")), samples)
    _write_tsv(ctx.path(L.match_summary(ds, "swaps.tsv")), swaps)
    _write_tsv(ctx.path(L.match_summary(ds, "candidates.tsv")), cands)
    _write_tsv(ctx.path(L.match_summary(ds, "duplicates.tsv")), duplicates)
    _write_tsv(ctx.path(L.match_summary(ds, "related.tsv")), related)
    _write_tsv(ctx.path(L.match_summary(ds, "weights.tsv")), weights)

    mapped = samples["geno_id_assigned"] != ""
    status_counts = {
        s: int((samples["qpgen_status"] == s).sum()) for s in (*STATUS_ORDER, NOT_SCORED)
    }
    status_counts = {k: v for k, v in status_counts.items() if v or k != NOT_SCORED}
    group_counts = {g: int((samples["group"] == g).sum()) for g in GROUP_ORDER}
    fail = mismatches(samples, duplicates)
    verdict, verdict_status = make_verdict(samples, duplicates=duplicates)
    n_not_scored = int((samples["qpgen_status"] == NOT_SCORED).sum())
    warnings = list(params.get("warnings") or [])
    if n_not_scored:
        warnings.append(
            f"{n_not_scored} sample{'s' if n_not_scored != 1 else ''} not scored: no covariates"
        )
    failing = [
        {
            "pheno_id": r.pheno_id,
            "geno_id_assigned": r.geno_id_assigned,
            "qpgen_status": r.qpgen_status,
            "best_id": r.best_id,
            "z_best": _num(r.z_best),
            "z_self": _num(r.z_self),
            "rank_self": None if pd.isna(r.rank_self) else int(r.rank_self),
            "z_gap": _num(r.z_gap),
            "proposed_geno_id": r.proposed_geno_id,
            "note": r.note,
        }
        for r in fail.itertuples()
    ]
    wv = weights["weight"]
    write_json(
        ctx.path(L.match_summary(ds, "summary.json")),
        {
            "schema_version": SCHEMA_VERSION,
            "pipeline": "match",
            "dataset": ds,
            "omics": ds_info.get("omics"),
            "reference": params.get("reference"),
            "n_pheno": len(samples),
            "n_scored": len(match_all),
            "n_not_scored": n_not_scored,
            "n_mapped": int(mapped.sum()),
            "n_unmapped": int((~mapped).sum()),
            "n_prs_samples": params.get("n_prs_samples"),
            "status_counts": status_counts,
            "group_counts": group_counts,
            "n_fail": len(fail),
            "failing": failing,
            "n_rescue_candidates": len(cands),
            "n_reciprocal_swaps": int((swaps["kind"] == "reciprocal").sum() // 2),
            "n_duplicates": len(duplicates),
            "n_fail_two_ids": int(
                (
                    (samples["group"] == "FAIL")
                    & samples["pheno_id"].isin(two_id_samples(duplicates))
                ).sum()
            ),
            "n_related": len(related),
            "verdict": verdict,
            "verdict_status": verdict_status,
            "lambda_mode": params.get("lambda_mode"),
            "lambda_used": params.get("lambda_used"),
            "mahalanobis": params.get("mahalanobis"),
            "rint": params.get("rint"),
            "use_cov": params.get("use_cov"),
            "z_threshold": zt,
            "z_diff": zd,
            "n_traits_used": params.get("n_traits_used"),
            "n_traits_weight_le0": int((wv <= 0).sum()),
            "weights": _quantiles(wv),
            "n_status_disagreements": n_disagree,
            "warnings": warnings,
            "qpgen_argv": params.get("argv"),
            "qpgen_log_excerpt": params.get("log_excerpt") or [],
            "trait_tsv_used": params.get("trait_tsv_used"),
            "sample_tsv_used": params.get("sample_tsv_used"),
            "dataset_counts": {
                k: ds_info.get(k)
                for k in (
                    "n_samples_pheno",
                    "n_traits",
                    "n_samples_geno",
                    "n_samples_mapped",
                    "n_pheno_unmapped",
                    "n_cov_columns",
                    "pheno_format",
                    "sample_map",
                )
            },
        },
    )
    if ctx.log is not None:
        ctx.log.info("match %s: %s", ds, verdict)
