"""match_qpgen and summarize_match steps, and the match pipeline's step graph."""

from __future__ import annotations

import gzip
import json
import os
import stat
import sys
import textwrap
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from test_match_helpers import PHENO_IDS, summarize_ctx, summarized_run, write_qpgen_outputs

from anqchor.config.model import Config
from anqchor.pipelines import build_graph
from anqchor.pipelines import match as match_pipeline
from anqchor.pipelines.context import BuildContext
from anqchor.pipelines.layout import Layout as L
from anqchor.steps import StepContext, StepError
from anqchor.steps.match import (
    SAMPLE_COLUMNS,
    build_samples,
    make_verdict,
    match_qpgen,
    summarize_match,
)
from anqchor.util.log import get_logger

DS = "rna"


def tsv(run: Path, ds: str, name: str) -> pd.DataFrame:
    return pd.read_csv(
        run / L.match_summary(ds, name), sep="\t", keep_default_na=False, na_values=["NA"]
    )


# ====================================================================== summarize
def test_summarize_outputs(tmp_path: Path) -> None:
    run = summarized_run(tmp_path)
    s = tsv(run, DS, "samples.tsv")
    assert list(s.columns) == SAMPLE_COLUMNS
    assert s["pheno_id"].tolist() == [*PHENO_IDS[:12], "U1", "U2", "U3", "S13"]
    by = s.set_index("pheno_id")
    assert (
        by.loc["S1", "group"] == "PASS"
        and by.loc["S10", "group"] == "WARN"
        and by.loc["S11", "group"] == "WARN"
    )
    assert by.loc["S7", "group"] == "FAIL" and by.loc["S12", "group"] == "FAIL"
    assert by.loc["U1", "group"] == "NO_TEST" and by.loc["U1", "proposed_geno_id"] == "G20"
    assert by.loc["S7", "proposed_geno_id"] == "G8" and by.loc["S12", "proposed_geno_id"] == ""
    assert by.loc["S13", "qpgen_status"] == "NOT_SCORED" and by.loc["S13", "group"] == "NO_TEST"
    # not scored (no covariates) but genotyped: keeps its assigned ID from the sample map
    assert by.loc["S13", "geno_id_assigned"] == "G13"
    assert by.loc["S9", "note"].startswith("best genotype is assigned to S1 (PASS)")
    # G1 is S1's confirmed genotype: relabelling S9 to it would map two samples to one genotype
    assert (
        by.loc["S9", "proposed_geno_id"] == ""
        and "confirmed match of another sample" in by.loc["S9", "note"]
    )
    assert "duplicate" not in by.loc["S9", "note"]
    assert "rescue candidate: G20" in by.loc["U1", "note"]
    assert "possible relatives: G6, G6P" in by.loc["S6", "note"]
    # S12's two leaders score about equally: two genotype IDs, not relatives
    assert "matches two genotype IDs about equally: Y4, Y5" in by.loc["S12", "note"]
    assert "relatives" not in by.loc["S12", "note"]
    assert by.loc["S7", "z_gap"] == pytest.approx(3.5)
    assert by.loc["S1", "id_1"] == "G1" and by.loc["U3", "id_1"] == ""

    sw = tsv(run, DS, "swaps.tsv")
    recip = sw[sw["kind"] == "reciprocal"]
    assert sorted(recip["pheno_id"]) == ["S7", "S8"]
    assert set(recip["partner_pheno_id"]) == {"S7", "S8"}
    one = sw[sw["kind"] == "one-way"].set_index("pheno_id")
    assert sorted(one.index) == ["S12", "S9"]
    assert one.loc["S9", "partner_pheno_id"] == "S1" and one.loc["S9", "partner_group"] == "PASS"

    c = tsv(run, DS, "candidates.tsv")
    assert c["pheno_id"].tolist() == ["U1"] and c["candidate_geno_id"].tolist() == ["G20"]
    assert c["candidate_assigned_to"].tolist() == [""]

    d = tsv(run, DS, "duplicates.tsv")
    assert d["pheno_id"].tolist() == ["S12"]
    assert d[["geno_id_1", "geno_id_2"]].values.tolist() == [["Y4", "Y5"]]
    assert d["assigned_in_pair"].tolist() == [False]

    r = tsv(run, DS, "related.tsv")
    assert sorted(r["pheno_id"]) == [
        "S6"
    ]  # U2's flat top five is not relatives; S12 is a duplicate

    w = tsv(run, DS, "weights.tsv")
    assert w["ref_trait"].tolist() == ["R1:T1", "R2:T2", "R3:T3", "R4:T4"]
    assert w["weight_le0"].tolist() == [0, 0, 1, 1]

    sm = json.loads((run / L.match_summary(DS, "summary.json")).read_text())
    assert sm["schema_version"] == 1 and sm["dataset"] == DS
    assert (
        sm["n_pheno"],
        sm["n_scored"],
        sm["n_not_scored"],
        sm["n_mapped"],
        sm["n_unmapped"],
    ) == (16, 15, 1, 13, 3)
    assert sm["warnings"] == ["1 sample not scored: no covariates"]
    assert sm["group_counts"] == {"PASS": 6, "WARN": 2, "FAIL": 4, "NO_TEST": 4}
    assert sm["status_counts"]["SINGLE_NEW_BEST"] == 4 and sm["status_counts"]["NOT_SCORED"] == 1
    assert sm["n_fail"] == 4 and [f["pheno_id"] for f in sm["failing"]] == ["S7", "S8", "S9", "S12"]
    assert sm["n_rescue_candidates"] == 1 and sm["n_reciprocal_swaps"] == 1 and sm["n_related"] == 1
    assert sm["n_duplicates"] == 1
    assert sm["n_traits_weight_le0"] == 2 and sm["weights"]["max"] == pytest.approx(0.3)
    assert sm["n_status_disagreements"] == 0
    assert sm["verdict_status"] == "critical"
    assert sm["verdict"].startswith(
        "7 of 13 mapped samples confirmed (1 only leniently); 4 mismatches (S7 (labelled G7) best matches G8; "
    )
    assert (
        "1 unclear" in sm["verdict"] and "3 unmapped samples, 1 rescue candidate" in sm["verdict"]
    )
    assert "1 not scored (no covariates)" in sm["verdict"]
    assert "1 sample matching two genotype IDs about equally (S12: Y4, Y5)" in sm["verdict"]


def test_weak_reciprocal_is_not_a_swap() -> None:
    """Two SELF_LENIENT samples pointing at each other's genotype by a 0.1 Z gap: listed as
    'reciprocal-weak', not counted as a swap, and the verdict stays good."""
    from test_match_helpers import HEADER, row

    from anqchor.steps.match import find_swaps

    rows = [
        row(
            "A",
            "SELF_LENIENT",
            "gA",
            3.0,
            2,
            [("gB", 3.1), ("gA", 3.0), *[(f"X{i}", 1.0) for i in range(3)]],
        ),
        row(
            "B",
            "SELF_LENIENT",
            "gB",
            3.0,
            2,
            [("gA", 3.1), ("gB", 3.0), *[(f"X{i}", 1.0) for i in range(3)]],
        ),
    ]
    df = pd.DataFrame(rows, columns=HEADER)
    for c in df.columns:
        if c.startswith(("Z.", "COR.", "Rank.", "N.")):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    s, _ = build_samples(df, ["A", "B"], 1.96, 2.0)
    sw = find_swaps(s)
    assert sw["kind"].tolist() == ["reciprocal-weak", "reciprocal-weak"]
    assert int((sw["kind"] == "reciprocal").sum()) == 0
    assert make_verdict(s)[1] == "good"


def test_summarize_counts_disagreement(tmp_path: Path) -> None:
    from test_match_helpers import row, scenario_rows

    rows = scenario_rows()
    rows[0] = row(
        "S1", "UNCLEAR", "G1", 6.1, 1, [("G1", 6.1), ("X", 3.0)]
    )  # qpgen says UNCLEAR, rules say SELF_BEST
    write_qpgen_outputs(tmp_path, DS, rows)
    summarize_match(summarize_ctx(tmp_path, DS))
    sm = json.loads((tmp_path / L.match_summary(DS, "summary.json")).read_text())
    assert sm["n_status_disagreements"] == 1


def test_verdict_variants() -> None:
    def frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
        return pd.DataFrame(
            rows, columns=["pheno_id", "geno_id_assigned", "qpgen_status", "group", "best_id"]
        )

    good = frame(
        [
            {
                "pheno_id": "A",
                "geno_id_assigned": "A",
                "qpgen_status": "SELF_BEST",
                "group": "PASS",
                "best_id": "A",
            }
        ]
    )
    assert make_verdict(good) == ("1 of 1 mapped samples confirmed", "good")
    relabelled = frame(
        [
            {
                "pheno_id": "A",
                "geno_id_assigned": "B",
                "qpgen_status": "SINGLE_NEW_BEST",
                "group": "FAIL",
                "best_id": "A",
            }
        ]
    )
    text, st = make_verdict(relabelled)
    assert st == "critical" and "A (labelled B) best matches A" in text and "1 mismatch (" in text
    many = frame(
        [
            {
                "pheno_id": f"P{i}",
                "geno_id_assigned": f"P{i}",
                "qpgen_status": "SINGLE_NEW_BEST",
                "group": "FAIL",
                "best_id": "Z",
            }
            for i in range(5)
        ]
    )
    assert "; and 2 more)" in make_verdict(many)[0]
    none = frame(
        [
            {
                "pheno_id": "A",
                "geno_id_assigned": "",
                "qpgen_status": "NO_SELF",
                "group": "NO_TEST",
                "best_id": "Z",
            }
        ]
    )
    text, st = make_verdict(none)
    assert (
        st == "warning"
        and text.startswith("no phenotyped sample is mapped")
        and "0 rescue candidates" in text
    )


def test_build_samples_handles_empty_output() -> None:
    df = pd.DataFrame(
        columns=[
            "ID.Pheno",
            "N.Traits",
            "MatchStatus",
            "ID.self",
            "Z.self",
            "COR.self",
            "Rank.self",
        ]
        + [f"{p}.{o}" for o in ("1st", "2nd", "3rd", "4th", "5th") for p in ("ID", "Z", "COR")]
    )
    s, n = build_samples(df, ["A", "B"], 1.96, 2.0)
    assert n == 0 and s["qpgen_status"].tolist() == ["NOT_SCORED", "NOT_SCORED"]


# ======================================================================= qpgen step
FAKE_QPGEN = textwrap.dedent("""\
    import gzip, json, os, sys
    argv = sys.argv[1:]
    with open(os.environ["FAKE_ARGV"], "w") as fh:
        json.dump(argv, fh)
    print("NOTICE [x] - Loading things")
    print("NOTICE [x] - Expanding phenotype matrix to 200 samples")
    if "--auto-lambda" in argv:
        print("NOTICE [x] - Using automatically tuned lambda = 0.6008 after 13 evaluations")
    rc = int(os.environ.get("FAKE_RC", "0"))
    if rc:
        print("FATAL ERROR -")
        print(os.environ.get("FAKE_MSG", "something broke"))
        sys.exit(rc)
    out = argv[argv.index("--out") + 1]
    src = os.environ["FAKE_SRC"]
    for suffix in (".weights.tsv.gz", ".match.assigned.tsv.gz", ".match.all.tsv.gz"):
        with open(src + suffix, "rb") as a, open(out + suffix, "wb") as b:
            b.write(a.read())
""")


def _gz(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


@pytest.fixture
def qpgen_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """A run directory with the inputs of ``match.rna.qpgen`` and a fake qpgentools."""
    run = tmp_path / "run"
    ref = "ref"
    # template outputs for the fake binary
    tmpl = tmp_path / "tmpl"
    write_qpgen_outputs(tmpl, "t")
    src = tmp_path / "src" / "qpgen"
    src.parent.mkdir()
    for p in L.match_qpgen_outputs("t"):
        os.replace(tmpl / p, Path(str(src) + p.name[len("qpgen") :]))
    fake = tmp_path / "qpgentools"
    fake.write_text(f"#!{sys.executable}\n" + FAKE_QPGEN)
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("FAKE_ARGV", str(tmp_path / "argv.json"))
    monkeypatch.setenv("FAKE_SRC", str(src))
    monkeypatch.delenv("FAKE_RC", raising=False)
    # inputs
    geno = [f"G{i}" for i in range(1, 9)]
    pheno = ["G1", "G2", "G3", "G4", "P5"]
    _gz(
        run / L.prs_matrix(ref),
        ["FID\tIID\tA:x\tB:y\tC:z"] + [f"{g}\t{g}\t0.1\t0.2\t0.3" for g in geno],
    )
    _gz(run / L.ds_pheno(DS), ["FID\tIID\tA\tB\tD"] + [f"{s}\t{s}\t1\t2\t3" for s in pheno])
    _gz(run / L.ds_cov(DS), ["FID\tIID\tc1"] + [f"{s}\t{s}\t1" for s in pheno if s != "G4"])
    (run / L.ds_samples(DS)).write_text("\n".join(pheno) + "\n")
    (run / L.ds_sample_map(DS)).write_text("G1\tG1\nG2\tG2\nG3\tG3\nG4\tG4\nG9\tP5\n")
    (run / L.dsref_dir(DS, ref)).mkdir(parents=True, exist_ok=True)
    (run / L.dsref_trait_map(DS, ref)).write_text("A:x\tA\nB:y\tB\nC:z\tC\nE:w\tD\n")
    (run / L.ds_json(DS)).write_text(json.dumps({"n_samples_mapped": 4}))
    (run / L.dsref_json(DS, ref)).write_text("{}")
    return {"run": run, "fake": fake, "argv": tmp_path / "argv.json", "ref": ref}


def _qctx(env: dict[str, Any], *, cov: bool = True, **params: Any) -> StepContext:
    run, ref = env["run"], env["ref"]
    p = {
        "qpgen": str(env["fake"]),
        "ref": ref,
        "dataset": DS,
        "has_sample_map": False,
        "use_cov": cov,
        "rint": True,
        "mahalanobis": False,
        "lambda": "auto",
        "min_self_for_auto": 100,
        "z_threshold": 1.96,
        "z_diff": 2.0,
        "min_weight": 0.0,
        "threads": 2,
        "missing_str": ["NA", "."],
    }
    p.update(params)
    inputs = (
        [L.prs_matrix(ref), L.ds_pheno(DS)]
        + ([L.ds_cov(DS)] if cov else [])
        + [
            L.ds_samples(DS),
            L.ds_sample_map(DS),
            L.dsref_trait_map(DS, ref),
            L.ds_json(DS),
            L.dsref_json(DS, ref),
        ]
    )
    return StepContext(
        name=f"match.{DS}.qpgen",
        run_dir=run,
        params=p,
        inputs=inputs,
        outputs=[],
        log=get_logger("test.qpgen"),
    )


def _argv(env: dict[str, Any]) -> list[str]:
    return json.loads(env["argv"].read_text())


def _params(env: dict[str, Any]) -> dict[str, Any]:
    return json.loads((env["run"] / L.match_params(DS)).read_text())


def test_qpgen_step_default(qpgen_run: dict[str, Any]) -> None:
    match_qpgen(_qctx(qpgen_run))
    argv = _argv(qpgen_run)
    assert argv[0] == "match-prs-pheno"
    assert "--trait-tsv" in argv  # A:x -> A is not an identity map
    assert "--sample-tsv" not in argv and "--mahalanobis" not in argv and "--lambda" not in argv
    assert argv[argv.index("--cov") + 1] == str(L.ds_cov(DS))
    assert argv[argv.index("--missing-str") + 1] == "NA,."
    assert argv[argv.index("--threads") + 1] == "2"
    run = qpgen_run["run"]
    for p in L.match_qpgen_outputs(DS):
        assert (run / p).exists()
    assert not list((run / L.match_dir(DS)).glob(".qpgen.*"))  # temporary directory removed
    # E:w has no PRS; C:z -> C is not a phenotype column
    assert (run / L.match_trait_map(DS)).read_text() == "A:x\tA\nB:y\tB\n"
    p = _params(qpgen_run)
    assert p["lambda_mode"] == "independent" and p["lambda_used"] is None
    assert (
        p["n_traits_used"] == 2 and p["n_dropped_no_prs"] == 1 and p["n_dropped_not_in_pheno"] == 1
    )
    assert p["n_samples_mapped"] == 3  # G4 has no covariates; P5 has no genotype
    assert p["trait_tsv_used"] and not p["sample_tsv_used"] and p["use_cov"]
    assert p["argv"][p["argv"].index("--out") + 1] == str(L.match_prefix(DS))
    assert all("Expanding" not in x for x in p["log_excerpt"])


def test_qpgen_step_identity_traits_without_cov(qpgen_run: dict[str, Any]) -> None:
    run, ref = qpgen_run["run"], qpgen_run["ref"]
    _gz(
        run / L.prs_matrix(ref), ["FID\tIID\tA\tB"] + [f"G{i}\tG{i}\t0.1\t0.2" for i in range(1, 9)]
    )
    (run / L.dsref_trait_map(DS, ref)).write_text("A\tA\nB\tB\n")
    match_qpgen(_qctx(qpgen_run, cov=False))
    argv = _argv(qpgen_run)
    assert "--trait-tsv" not in argv and "--cov" not in argv
    assert _params(qpgen_run)["n_samples_mapped"] == 4


def test_qpgen_step_sample_map_filtered(qpgen_run: dict[str, Any]) -> None:
    match_qpgen(_qctx(qpgen_run, has_sample_map=True))
    argv = _argv(qpgen_run)
    used = argv[argv.index("--sample-tsv") + 1]
    assert used != str(L.ds_sample_map(DS))
    # G4 lacks covariates, G9 has no PRS
    assert (qpgen_run["run"] / used).read_text() == "G1\tG1\nG2\tG2\nG3\tG3\n"
    p = _params(qpgen_run)
    assert p["n_sample_map_rows_dropped"] == 2 and p["sample_tsv_used"] and p["warnings"]


def test_qpgen_step_sample_map_unfiltered(qpgen_run: dict[str, Any]) -> None:
    run = qpgen_run["run"]
    (run / L.ds_sample_map(DS)).write_text("G1\tG1\nG2\tG2\n")
    match_qpgen(_qctx(qpgen_run, has_sample_map=True))
    argv = _argv(qpgen_run)
    assert argv[argv.index("--sample-tsv") + 1] == str(L.ds_sample_map(DS))


@pytest.mark.parametrize(
    ("lam", "min_self", "expect_flags", "mode", "used"),
    [
        (
            "auto",
            3,
            ["--auto-lambda", "--auto-lambda-min-self", "3", "--auto-lambda-min-best", "3"],
            "auto",
            0.6008,
        ),
        ("auto", 100, ["--lambda", "0.5"], "fallback", 0.5),
        (0.3, 100, ["--lambda", "0.3"], "fixed", 0.3),
    ],
)
def test_qpgen_step_lambda(
    qpgen_run: dict[str, Any],
    lam: Any,
    min_self: int,
    expect_flags: list[str],
    mode: str,
    used: float,
) -> None:
    match_qpgen(
        _qctx(qpgen_run, mahalanobis=True, **{"lambda": lam, "min_self_for_auto": min_self})
    )
    argv = _argv(qpgen_run)
    assert "--mahalanobis" in argv
    joined = " ".join(argv)
    assert " ".join(expect_flags) in joined
    p = _params(qpgen_run)
    assert p["lambda_mode"] == mode and p["lambda_used"] == pytest.approx(used)
    assert bool(p["warnings"]) == (mode == "fallback")


def test_qpgen_step_min_best_failure(
    qpgen_run: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_RC", "1")
    monkeypatch.setenv(
        "FAKE_MSG",
        "--auto-lambda: only 0 of 729 mapped samples rank their own PRS first at the "
        "tuned lambda, fewer than --auto-lambda-min-best 100.",
    )
    with pytest.raises(StepError, match="sample mapping is wrong") as ei:
        match_qpgen(_qctx(qpgen_run, mahalanobis=True, min_self_for_auto=2))
    assert "fewer than --auto-lambda-min-best" in str(ei.value)
    assert not (qpgen_run["run"] / L.match_params(DS)).exists()
    assert not list((qpgen_run["run"] / L.match_dir(DS)).glob(".qpgen.*"))


def test_qpgen_step_min_self_and_generic_failures(
    qpgen_run: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_RC", "1")
    monkeypatch.setenv(
        "FAKE_MSG",
        "--auto-lambda: only 5 phenotyped samples ... fewer than --auto-lambda-min-self 100",
    )
    with pytest.raises(StepError, match="too few mapped samples with observed traits"):
        match_qpgen(_qctx(qpgen_run, mahalanobis=True, min_self_for_auto=2))
    monkeypatch.setenv("FAKE_MSG", "boom")
    with pytest.raises(StepError, match=r"failed \(exit 1\).*\n(.*\n)*boom"):
        match_qpgen(_qctx(qpgen_run))


def test_qpgen_step_missing_binary(qpgen_run: dict[str, Any]) -> None:
    with pytest.raises(StepError, match="cannot run qpgentools"):
        match_qpgen(_qctx(qpgen_run, qpgen=str(qpgen_run["run"] / "nope")))


def test_qpgen_step_no_traits(qpgen_run: dict[str, Any]) -> None:
    (qpgen_run["run"] / L.dsref_trait_map(DS, qpgen_run["ref"])).write_text("Q:q\tQ\n")
    with pytest.raises(StepError, match="none of the 1 mapped traits"):
        match_qpgen(_qctx(qpgen_run))


# ======================================================================== pipeline
def _cfg(tmp_path: Path, **match: Any) -> Config:
    return Config.model_validate(
        {
            "project": "p",
            "genotypes": {"pgen_list": str(tmp_path / "g.tsv")},
            "references": {"r": {"pairs": str(tmp_path / "pairs.tsv")}},
            "datasets": {
                "a": {"pheno": str(tmp_path / "a.tsv"), "cov": str(tmp_path / "c.tsv")},
                "b": {
                    "pheno": str(tmp_path / "b.tsv"),
                    "sample_map": str(tmp_path / "m.tsv"),
                    "omics": "prot",
                },
            },
            "pipelines": {"match": match, "cis": {"enabled": False}, "trans": {"enabled": False}},
        }
    )


def test_pipeline_steps(tmp_path: Path) -> None:
    ctx = BuildContext(cfg=_cfg(tmp_path), qpgen=Path("/opt/qpgentools"))
    steps = {s.name: s for s in match_pipeline.build(ctx)}
    for name in (
        "prepare.genotypes",
        "prepare.reference.r",
        "prepare.dataset.a",
        "prepare.mapping.b.r",
        "prs.r",
        "match.a.qpgen",
        "match.a.summarize",
        "match.b.qpgen",
        "match.b.summarize",
        "report.match",
    ):
        assert name in steps, name
    qa, qb = steps["match.a.qpgen"], steps["match.b.qpgen"]
    assert (
        qa.kind == "match_qpgen" and L.ds_cov("a") in qa.inputs and L.ds_cov("b") not in qb.inputs
    )
    assert L.prs_matrix("r") in qa.inputs and L.dsref_trait_map("a", "r") in qa.inputs
    assert qa.outputs == (*L.match_qpgen_outputs("a"), L.match_trait_map("a"), L.match_params("a"))
    assert qa.params["use_cov"] and not qa.params["has_sample_map"] and qb.params["has_sample_map"]
    assert qa.params["qpgen"] == "/opt/qpgentools" and qa.params["lambda"] == "auto"
    assert set(qa.params) >= {
        "ref",
        "dataset",
        "rint",
        "mahalanobis",
        "min_self_for_auto",
        "z_threshold",
        "z_diff",
        "min_weight",
        "missing_str",
    }
    assert "threads" not in qa.params  # passed via $ANQCHOR_THREADS: changing it re-runs nothing
    sa = steps["match.a.summarize"]
    assert sa.kind == "summarize_match" and L.match_params("a") in sa.inputs
    assert len(sa.outputs) == 7 and L.match_summary("a", "summary.json") in sa.outputs
    assert (
        L.match_summary("a", "duplicates.tsv") in sa.outputs and L.ds_sample_map("a") in sa.inputs
    )
    rep = steps["report.match"]
    assert rep.outputs == (L.report("match"),) and rep.params["datasets"] == ["a", "b"]
    assert rep.params["omics"] == {"a": "omics", "b": "prot"}
    assert L.match_summary("b", "samples.tsv") in rep.inputs


def test_pipeline_use_cov_off_and_dataset_subset(tmp_path: Path) -> None:
    ctx = BuildContext(cfg=_cfg(tmp_path, use_cov=False, datasets=["a"]), qpgen=Path("q"))
    names = [s.name for s in match_pipeline.build(ctx)]
    assert "match.b.qpgen" not in names and "prepare.dataset.b" not in names
    q = next(s for s in match_pipeline.build(ctx) if s.name == "match.a.qpgen")
    assert L.ds_cov("a") not in q.inputs and q.params["use_cov"] is False


def test_threshold_change_touches_only_match_steps(tmp_path: Path) -> None:
    a = {
        s.name: s.signature()
        for s in build_graph(BuildContext(cfg=_cfg(tmp_path), qpgen=Path("q")), ["match"])[0]
    }
    b = {
        s.name: s.signature()
        for s in build_graph(
            BuildContext(cfg=_cfg(tmp_path, z_threshold=2.5), qpgen=Path("q")), ["match"]
        )[0]
    }
    changed = sorted(n for n in a if a[n] != b[n])
    assert changed == ["match.a.qpgen", "match.b.qpgen", "report.match"]


def test_threads_not_in_signature_and_read_from_env(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from anqchor.steps.match import qpgen_threads

    a = {
        s.name: s.signature()
        for s in build_graph(BuildContext(cfg=_cfg(tmp_path), qpgen=Path("q")), ["match"])[0]
    }
    cfg = _cfg(tmp_path)
    cfg.qpgen.threads = 8
    b = {
        s.name: s.signature()
        for s in build_graph(BuildContext(cfg=cfg, qpgen=Path("q")), ["match"])[0]
    }
    assert a == b
    monkeypatch.delenv("ANQCHOR_THREADS", raising=False)
    assert qpgen_threads({}) == 1 and qpgen_threads({"threads": 3}) == 3
    monkeypatch.setenv("ANQCHOR_THREADS", "6")
    assert qpgen_threads({"threads": 3}) == 6
    monkeypatch.setenv("ANQCHOR_THREADS", "x")
    assert qpgen_threads({}) == 1


def test_rescue_conflict_and_duplicates() -> None:
    """An unassigned sample whose best genotype is another sample's confirmed match is not a
    rescue candidate; samples matching two genotype IDs about equally are listed in
    duplicates.tsv (assigned or not, own ID first or not) and not as relatives."""
    from test_match_helpers import HEADER, row

    from anqchor.steps.match import find_candidates, find_duplicates, find_related

    bg = [(f"X{i}", 2.6 - 0.1 * i) for i in range(3)]
    rows = [
        row("A", "SELF_BEST", "gA", 6.0, 1, [("gA", 6.0), *bg, ("X9", 2.0)]),
        row("_A", "SINGLE_NEW_BEST", None, None, None, [("gA", 5.0), *bg, ("X9", 2.0)]),
        row("_B", "SINGLE_NEW_BEST", None, None, None, [("gB", 5.0), *bg, ("X9", 2.0)]),
        row("C", "SELF_BEST", "gC", 6.99, 1, [("gC", 6.99), ("gC_dup", 6.99), *bg]),
        row("D", "MULTI_NEW_BEST", "gD", 6.83, 2, [("gD_twin", 7.0), ("gD", 6.83), *bg]),
        row("_E", "MULTI_NEW_BEST", None, None, None, [("gE", 6.9), ("gE_dup", 6.9), *bg]),
        row("F", "SELF_BEST", "gF", 9.2, 1, [("gF", 9.2), ("gP", 5.1), *bg]),
    ]
    df = pd.DataFrame(rows, columns=HEADER).replace({"NA": None})
    for c in df.columns:
        if c.startswith(("Z.", "COR.", "Rank.", "N.")):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    s, n = build_samples(df, [r[0] for r in rows], 1.96, 2.0)
    assert n == 0
    by = s.set_index("pheno_id")
    assert by.loc["_A", "proposed_geno_id"] == "" and "rescue candidate" not in by.loc["_A", "note"]
    assert "confirmed match of another sample" in by.loc["_A", "note"]
    assert find_candidates(s)["pheno_id"].tolist() == ["_B"]
    d = find_duplicates(s, 1.96, 2.0).set_index("pheno_id")
    assert d.index.tolist() == ["C", "D", "_E"]
    assert d["assigned_in_pair"].tolist() == [True, True, False]
    assert d.loc["D", ["geno_id_1", "geno_id_2"]].tolist() == ["gD_twin", "gD"]
    assert (
        "several genotypes" not in by.loc["D", "note"]
        and "about equally: gD_twin, gD" in by.loc["D", "note"]
    )
    assert find_related(s, 1.96, 2.0)["pheno_id"].tolist() == ["F"]
    text, _ = make_verdict(s, duplicates=find_duplicates(s, 1.96, 2.0))
    assert (
        "3 samples matching two genotype IDs about equally (C: gC, gC_dup; D: gD_twin, gD; _E: gE, gE_dup)"
        in text
    )
    assert "1 rescue candidate" in text


def _frame(rows: list[Any]) -> pd.DataFrame:
    from test_match_helpers import HEADER

    df = pd.DataFrame(rows, columns=HEADER).replace({"NA": None})
    for c in df.columns:
        if c.startswith(("Z.", "COR.", "Rank.", "N.")):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def test_two_genotype_ids_are_not_mismatches_or_relabels() -> None:
    """A sample that matches its own ID and a second genotype ID about equally is one person
    under two IDs: not a mismatch (verdict, swaps, n_fail) and no proposed ID, even when
    the second ID leads clearly (strong signal: 12.4 vs 10)."""
    from test_match_helpers import row

    from anqchor.steps.match import (
        find_candidates,
        find_duplicates,
        find_swaps,
        mismatches,
        two_id_samples,
    )

    bg = [(f"X{i}", 2.9 - 0.1 * i) for i in range(3)]
    rows = [
        # the dupnoisy shape: own ID second, its noisy copy first (qpgen: MULTI_NEW_BEST)
        row(
            "D",
            "MULTI_NEW_BEST",
            "gD",
            6.83,
            2,
            [("gD_twin", 7.0), ("gD", 6.83), ("X9", 2.56), *bg[:2]],
        ),
        # strong signal: the copy leads clearly (SINGLE_NEW_BEST) and both are about equal
        row(
            "A",
            "SINGLE_NEW_BEST",
            "G2",
            10.0,
            2,
            [("G1", 12.4), ("G2", 10.0), ("X8", 3.0), *bg[:2]],
        ),
        row(
            "_U",
            "SINGLE_NEW_BEST",
            None,
            None,
            None,
            [("H1", 12.4), ("H2", 10.0), ("X7", 3.0), *bg[:2]],
        ),
        # a real mismatch stays one
        row("M", "SINGLE_NEW_BEST", "gM", 2.0, 40, [("gQ", 6.0), ("X6", 3.0), *bg]),
    ]
    s, n = build_samples(_frame(rows), [r[0] for r in rows], 1.96, 2.0)
    assert n == 0
    by = s.set_index("pheno_id")
    assert by.loc["A", "proposed_geno_id"] == "" and by.loc["_U", "proposed_geno_id"] == ""
    assert by.loc["M", "proposed_geno_id"] == "gQ"
    assert "rescue candidate" not in by.loc["_U", "note"]
    assert find_candidates(s).empty
    dups = find_duplicates(s, 1.96, 2.0)
    assert two_id_samples(dups) == {"D", "A"}
    assert mismatches(s, dups)["pheno_id"].tolist() == ["M"]
    sw = find_swaps(s, two_id_samples(dups))
    assert sw["pheno_id"].tolist() == ["M"]
    text, status = make_verdict(s, duplicates=dups)
    assert "1 mismatch (M (labelled gM) best matches gQ)" in text and status == "critical"
    assert "0 rescue candidates" in text
    # duplicates alone give a warning, not a critical verdict
    only = s[s["pheno_id"] != "M"]
    assert make_verdict(only, duplicates=find_duplicates(only, 1.96, 2.0))[1] == "warning"
