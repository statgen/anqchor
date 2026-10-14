"""Builders of small synthetic match runs shared by the ``test_match_*`` modules.

The scenario (thresholds 1.96 / 2.0), one row per phenotyped sample:

* S1-S6: SELF_BEST (PASS); S6 has a close second candidate (possible relative).
* S7/S8: reciprocal swap (each best matches the other's genotype; FAIL).
* S9: one-way mismatch, best genotype G1 which S1 holds (FAIL).
* S10: SELF_LENIENT, S11: UNCLEAR (WARN); S12: MULTI_NEW_BEST (FAIL, mapped) whose two
  leading genotypes (Y4, Y5, not its own) score about equally: it matches two genotype IDs.
* U1: unmapped rescue candidate (SINGLE_NEW_BEST, best G20); U2: unmapped NO_SELF;
  U3: unmapped NO_OBS_TRAITS.
* S13: phenotyped but absent from qpgen's output (NOT_SCORED); the sample map assigns it G13.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

from anqchor.pipelines.layout import Layout as L
from anqchor.qpgen.outputs import MATCH_ORDINALS
from anqchor.steps import StepContext
from anqchor.util.files import write_json
from anqchor.util.log import get_logger

HEADER = ["ID.Pheno", "N.Traits", "MatchStatus", "ID.self", "Z.self", "COR.self", "Rank.self"]
for _o in MATCH_ORDINALS:
    HEADER += [f"ID.{_o}", f"Z.{_o}", f"COR.{_o}"]


def row(
    pheno: str,
    status: str,
    self_id: str | None,
    z_self: float | None,
    rank: int | None,
    top: list[tuple[str, float]],
    n_traits: int = 50,
) -> list[str]:
    """One ``match.all`` line (``NA`` for missing values)."""

    def f(v: Any) -> str:
        return "NA" if v is None else (f"{v:.6f}" if isinstance(v, float) else str(v))

    out = [
        pheno,
        str(n_traits),
        status,
        f(self_id),
        f(z_self),
        f(None if z_self is None else z_self / 20),
        f(rank),
    ]
    for k in range(5):
        if k < len(top):
            out += [top[k][0], f(top[k][1]), f(top[k][1] / 20)]
        else:
            out += ["NA", "NA", "NA"]
    return out


def _tail(start: float) -> list[tuple[str, float]]:
    return [(f"X{i}", start - 0.1 * i) for i in range(4)]


def scenario_rows() -> list[list[str]]:
    rows = [
        row(f"S{i}", "SELF_BEST", f"G{i}", 6.0 + i / 10, 1, [(f"G{i}", 6.0 + i / 10), *_tail(3.0)])
        for i in range(1, 6)
    ]
    rows.append(row("S6", "SELF_BEST", "G6", 6.0, 1, [("G6", 6.0), ("G6P", 5.0), *_tail(3.0)[:3]]))
    rows.append(row("S7", "SINGLE_NEW_BEST", "G7", 0.5, 900, [("G8", 6.5), *_tail(3.0)]))
    rows.append(row("S8", "SINGLE_NEW_BEST", "G8", 0.2, 1000, [("G7", 6.0), *_tail(3.0)]))
    rows.append(row("S9", "SINGLE_NEW_BEST", "G9", 1.0, 400, [("G1", 6.2), *_tail(3.0)]))
    rows.append(
        row("S10", "SELF_LENIENT", "G10", 3.0, 2, [("Y1", 3.2), ("G10", 3.0), *_tail(2.9)[:3]])
    )
    rows.append(row("S11", "UNCLEAR", "G11", 1.0, 30, [("Y2", 2.5), ("Y3", 2.4), *_tail(2.3)[:3]]))
    rows.append(
        row(
            "S12",
            "MULTI_NEW_BEST",
            "G12",
            0.1,
            800,
            [("Y4", 6.0), ("Y5", 5.8), ("Y6", 1.0), ("Y7", 0.9), ("Y8", 0.8)],
        )
    )
    rows.append(row("U1", "SINGLE_NEW_BEST", None, None, None, [("G20", 7.0), *_tail(3.0)]))
    rows.append(row("U2", "NO_SELF", None, None, None, [("Y9", 3.0), *_tail(2.9)]))
    rows.append(row("U3", "NO_OBS_TRAITS", None, None, None, [], n_traits=0))
    return rows


PHENO_IDS = [f"S{i}" for i in range(1, 14)] + ["U1", "U2", "U3"]


def _write_gz(path: Path, lines: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for r in lines:
            fh.write("\t".join(r) + "\n")


def write_qpgen_outputs(run: Path, ds: str, rows: list[list[str]] | None = None) -> None:
    """The three qpgen match outputs, the params json and trait map of a finished qpgen step."""
    rows = scenario_rows() if rows is None else rows
    weights, assigned, match_all = L.match_qpgen_outputs(ds)
    _write_gz(run / match_all, [HEADER, *rows])
    _write_gz(run / assigned, [HEADER[:7], *[r[:7] for r in rows if r[3] != "NA"]])
    _write_gz(
        run / weights,
        [
            ["Trait", "Weight", "N.Obs"],
            ["T1", "0.300000", "12"],
            ["T2", "0.100000", "12"],
            ["T3", "-0.050000", "12"],
            ["T4", "0.000000", "11"],
        ],
    )
    (run / L.match_trait_map(ds)).write_text("R1:T1\tT1\nR2:T2\tT2\nR3:T3\tT3\nR4:T4\tT4\n")
    write_json(
        run / L.match_params(ds),
        {
            "schema_version": 1,
            "dataset": ds,
            "reference": "ref",
            "argv": ["qpgentools", "match-prs-pheno"],
            "mahalanobis": False,
            "lambda_mode": "independent",
            "lambda_used": None,
            "rint": True,
            "use_cov": True,
            "z_threshold": 1.96,
            "z_diff": 2.0,
            "n_traits_used": 4,
            "n_prs_samples": 40,
            "trait_tsv_used": True,
            "sample_tsv_used": False,
            "warnings": [],
            "log_excerpt": ["NOTICE - Analysis finished"],
        },
    )
    (run / L.ds_dir(ds)).mkdir(parents=True, exist_ok=True)
    write_json(
        run / L.ds_json(ds),
        {
            "schema_version": 1,
            "dataset": ds,
            "omics": "rna",
            "n_samples_pheno": 16,
            "n_samples_mapped": 13,
            "n_traits": 4,
        },
    )
    (run / L.ds_samples(ds)).write_text("\n".join(PHENO_IDS) + "\n")
    (run / L.ds_sample_map(ds)).write_text("".join(f"G{i}\tS{i}\n" for i in range(1, 14)))


def summarize_ctx(run: Path, ds: str) -> StepContext:
    """A StepContext for ``summarize_match`` matching the pipeline's declared inputs."""
    return StepContext(
        name=f"match.{ds}.summarize",
        run_dir=run,
        params={"dataset": ds},
        inputs=[
            *L.match_qpgen_outputs(ds),
            L.match_params(ds),
            L.match_trait_map(ds),
            L.ds_json(ds),
            L.ds_samples(ds),
            L.ds_sample_map(ds),
        ],
        outputs=[],
        log=get_logger("test.summarize"),
    )


def summarized_run(run: Path, ds: str = "rna") -> Path:
    """A run directory with the summaries of the synthetic scenario."""
    from anqchor.steps.match import summarize_match

    write_qpgen_outputs(run, ds)
    summarize_match(summarize_ctx(run, ds))
    return run
