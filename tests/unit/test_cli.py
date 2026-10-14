"""The ``anqchor`` command line: help, version, exit codes, dry runs and partial runs.

Most tests replace the pipeline modules with a small fake pipeline (shell steps plus the
real shared prepare steps), so they need neither qpgen nor the pipeline implementations.
Tests marked ``integration`` use the real pipelines on the MAGE example data.
"""

from __future__ import annotations

import json
import os
import shutil
import stat
import time
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from typer.testing import CliRunner

import anqchor.pipelines as pipelines_pkg
from anqchor import __version__
from anqchor.cli.execute import RunOptions, build_config, run_lock
from anqchor.cli.main import app
from anqchor.config.load import ConfigError
from anqchor.engine.minimake import stamp_path
from anqchor.engine.step import Step
from anqchor.pipelines.common import prs, shared_prepare_steps
from anqchor.pipelines.context import BuildContext
from anqchor.pipelines.layout import Layout

runner = CliRunner()
needs_make = pytest.mark.skipif(shutil.which("make") is None, reason="make not installed")

COMMANDS = [
    [],
    ["match"],
    ["match", "run"],
    ["cis"],
    ["cis", "run"],
    ["trans"],
    ["trans", "run"],
    ["trans", "all"],
    ["trans", "compare"],
    ["trans", "report"],
    ["run"],
    ["steps"],
    ["report"],
    ["qpgen"],
    ["qpgen", "build"],
    ["qpgen", "check"],
    ["prepare"],
    ["_step"],
]


@pytest.mark.parametrize("cmd", COMMANDS, ids=lambda c: " ".join(c) or "root")
def test_help(cmd: list[str]) -> None:
    r = runner.invoke(app, [*cmd, "--help"])
    assert r.exit_code == 0, r.output
    assert "Usage" in r.output


def test_pairs_help_keeps_optional_pval_column() -> None:
    """``[pval]`` must not be swallowed as rich markup."""
    r = runner.invoke(app, ["cis", "run", "--help"], env={"COLUMNS": "200"})
    assert r.exit_code == 0, r.output
    assert "beta se [pval]" in r.output


def test_threads_help_names_identity_matching() -> None:
    """--threads only affects match-prs-pheno; cis/trans help must say so."""
    for cmd in (["cis", "run"], ["trans", "run"], ["match", "run"]):
        r = runner.invoke(app, [*cmd, "--help"], env={"COLUMNS": "200"})
        text = " ".join(r.output.replace("│", " ").split())
        assert "only identity matching runs" in text, cmd


def test_version() -> None:
    r = runner.invoke(app, ["--version"])
    assert r.exit_code == 0
    assert r.output.strip() == f"anqchor {__version__}"


# ----------------------------------------------------------------- fake inputs
@pytest.fixture
def inputs(tmp_path: Path) -> SimpleNamespace:
    d = tmp_path / "data"
    d.mkdir()
    geno = d / "list.tsv"
    geno.write_text("22\t1\t1000000000\tgeno/chr22\n")
    pheno = d / "pheno.tsv"
    pheno.write_text("FID\tIID\tt1\nA\tA\t1\n")
    pairs = d / "pairs.tsv"
    pairs.write_text("trait\tvariant\tbeta\tse\nt1\t22:1:A:G\t0.1\t0.01\n")
    qpgen = tmp_path / "bin" / "qpgentools"
    qpgen.parent.mkdir()
    qpgen.write_text("#!/bin/sh\nexit 0\n")
    qpgen.chmod(qpgen.stat().st_mode | stat.S_IXUSR)
    return SimpleNamespace(
        dir=d, geno=geno, pheno=pheno, pairs=pairs, qpgen=qpgen, out=tmp_path / "run"
    )


def flag_args(i: SimpleNamespace, *extra: str) -> list[str]:
    return [
        "--qpgen",
        str(i.qpgen),
        "match",
        "run",
        "--geno-list",
        str(i.geno),
        "--pheno",
        str(i.pheno),
        "--pairs",
        str(i.pairs),
        "--out",
        str(i.out),
        *extra,
    ]


def fake_match_build(fail: bool = False, report: bool = False):  # type: ignore[no-untyped-def]
    """A match pipeline made of the real shared steps plus shell steps."""

    def build(ctx: BuildContext) -> list[Step]:
        L = ctx.layout
        ref = ctx.cfg.reference_for("match")
        out: list[Step] = list(shared_prepare_steps(ctx, ctx.cfg.datasets_for("match"), [ref]))
        out.append(prs(ctx, ref))
        for ds in ctx.cfg.datasets_for("match"):
            s1 = Path(f"match/{ds}/s1.txt")
            s2 = Path(f"match/{ds}/s2.txt")
            out += [
                Step(
                    name=f"match.{ds}.s1",
                    outputs=(s1,),
                    commands=(f"echo one > {s1}",),
                    pipeline="match",
                ),
                Step(
                    name=f"match.{ds}.s2",
                    inputs=(s1,),
                    outputs=(s2,),
                    pipeline="match",
                    commands=(
                        f"cat {s1} > {s2} && echo s2 >> ledger.txt" + ("; exit 7" if fail else ""),
                    ),
                ),
                Step(
                    name=f"match.{ds}.summary",
                    inputs=(s2,),
                    pipeline="match",
                    outputs=(L.match_summary(ds, "summary.json"), L.report("match")),
                    commands=(
                        f"echo '{{}}' > {L.match_summary(ds, 'summary.json')} && "
                        f"mkdir -p report && echo '<html></html>' > {L.report('match')}",
                    ),
                ),
            ]
            if report:
                out.append(
                    Step(
                        name="report.match",
                        inputs=(L.match_summary(ds, "summary.json"),),
                        outputs=(Path("report/extra.html"),),
                        pipeline="match",
                        commands=("mkdir -p report && echo '<html></html>' > report/extra.html",),
                    )
                )
        return out

    return build


@pytest.fixture
def fake_pipelines(monkeypatch: pytest.MonkeyPatch) -> Iterator[dict[str, object]]:
    """Replace the pipeline modules; ``state['fail']`` makes step s2 fail."""
    state: dict[str, object] = {"fail": False, "report": False}

    def module(name: str) -> SimpleNamespace:
        if name != "match":
            raise AssertionError(f"unexpected pipeline {name}")
        return SimpleNamespace(build=fake_match_build(bool(state["fail"]), bool(state["report"])))

    monkeypatch.setattr(pipelines_pkg, "pipeline_module", module)
    yield state


# ------------------------------------------------------------------ exit codes
def test_missing_flag_input_is_exit_2(inputs: SimpleNamespace) -> None:
    args = flag_args(inputs)
    args[args.index("--pheno") + 1] = str(inputs.dir / "absent.tsv")
    r = runner.invoke(app, args)
    assert r.exit_code == 2


def test_missing_config_input_is_exit_2(inputs: SimpleNamespace, tmp_path: Path) -> None:
    cfg = {
        "genotypes": {"pgen_list": str(inputs.geno)},
        "references": {"ref": {"pairs": "missing_pairs.tsv"}},
        "datasets": {"rna": {"pheno": str(inputs.pheno)}},
        "pipelines": {"cis": {"enabled": False}, "trans": {"enabled": False}},
    }
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(cfg))
    r = runner.invoke(app, ["--qpgen", str(inputs.qpgen), "run", "--config", str(p), "--dry-run"])
    assert r.exit_code == 2
    assert "input files not found" in r.output and "missing_pairs.tsv" in r.output


def test_invalid_config_is_exit_2(tmp_path: Path) -> None:
    p = tmp_path / "cfg.yaml"
    p.write_text("genotypes: {pgen_list: x}\nreferences: {}\ndatasets: {}\n")
    r = runner.invoke(app, ["run", "--config", str(p)])
    assert r.exit_code == 2 and "invalid configuration" in r.output
    r = runner.invoke(app, ["run"])
    assert r.exit_code == 2 and "--config is required" in r.output


def test_missing_qpgen_is_exit_4(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object]
) -> None:
    args = flag_args(inputs, "--dry-run")
    args[args.index("--qpgen") + 1] = str(inputs.dir / "no" / "qpgentools")
    r = runner.invoke(app, args)
    assert r.exit_code == 4, r.output
    assert "qpgentools not found" in r.output


def test_bad_lambda_is_exit_2(inputs: SimpleNamespace) -> None:
    r = runner.invoke(app, flag_args(inputs, "--mahalanobis", "--lambda", "fast", "--dry-run"))
    assert r.exit_code == 2
    r = runner.invoke(app, flag_args(inputs, "--mahalanobis", "--lambda", "3", "--dry-run"))
    assert r.exit_code == 2 and "lambda" in r.output


# -------------------------------------------------------------------- dry runs
@needs_make
def test_flag_only_dry_run(inputs: SimpleNamespace, fake_pipelines: dict[str, object]) -> None:
    r = runner.invoke(app, flag_args(inputs, "--dry-run", "--z-threshold", "2.5"))
    assert r.exit_code == 0, r.output
    run = inputs.out
    assert (run / "Makefile").exists()
    resolved = yaml.safe_load((run / Layout.resolved_config).read_text())
    assert resolved["genotypes"]["pgen_list"] == str(inputs.geno)
    assert resolved["datasets"]["dataset"]["pheno"] == str(inputs.pheno)
    assert resolved["references"]["reference"]["pairs"] == str(inputs.pairs)
    assert resolved["pipelines"]["match"]["enabled"] is True
    assert resolved["pipelines"]["cis"]["enabled"] is False
    assert resolved["pipelines"]["match"]["z_threshold"] == 2.5
    assert resolved["qpgen"]["binary"] == str(inputs.qpgen)
    assert "would run: prepare.genotypes" in r.output
    assert "match.dataset.s1" in r.output
    assert "report.index" not in r.output  # one pipeline: no index page
    # nothing ran
    assert not (run / "match").exists() and not (run / "run.json").exists()
    assert not list((run / "steps").glob("*.done"))
    mk = (run / "Makefile").read_text()
    assert "pair-prs" in mk and str(inputs.qpgen) in mk


def test_flags_override_single_config_dataset(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object], tmp_path: Path
) -> None:
    """--pheno with a config file overrides its only dataset instead of adding one."""
    cfg = {
        "genotypes": {"pgen_list": str(inputs.geno)},
        "references": {"mage": {"pairs": str(inputs.pairs)}},
        "datasets": {"rna": {"pheno": "nonexistent.tsv"}},
        "pipelines": {"cis": {"enabled": False}, "trans": {"enabled": False}},
        "out_dir": str(inputs.out),
    }
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(cfg))
    r = runner.invoke(
        app,
        [
            "--qpgen",
            str(inputs.qpgen),
            "match",
            "run",
            "--config",
            str(p),
            "--pheno",
            str(inputs.pheno),
            "--dry-run",
        ],
    )
    assert r.exit_code == 0, r.output
    resolved = yaml.safe_load((inputs.out / Layout.resolved_config).read_text())
    assert list(resolved["datasets"]) == ["rna"]
    assert resolved["datasets"]["rna"]["pheno"] == str(inputs.pheno)


# -------------------------------------------------------------- real make runs
@needs_make
def test_until_only_and_rerun(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object], monkeypatch: pytest.MonkeyPatch
) -> None:
    run = inputs.out
    r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.s2"))
    assert r.exit_code == 0, r.output
    assert "ran: match.dataset.s1, match.dataset.s2" in r.output
    assert (run / "ledger.txt").read_text().split() == ["s2"]
    assert not (run / stamp_path("match.dataset.summary")).exists()
    run_json = json.loads((run / "run.json").read_text())
    assert run_json["invocations"][-1]["returncode"] == 0

    # --only re-runs exactly that step; a dry run of it changes nothing
    stamp = run / stamp_path("match.dataset.s2")
    before = stamp.stat().st_mtime_ns
    r = runner.invoke(app, flag_args(inputs, "--only", "match.dataset.s2", "--dry-run"))
    assert r.exit_code == 0, r.output
    assert "would run: match.dataset.s2" in r.output
    assert stamp.exists() and stamp.stat().st_mtime_ns == before
    r = runner.invoke(app, flag_args(inputs, "--only", "match.dataset.s2"))
    assert r.exit_code == 0, r.output
    assert "ran: match.dataset.s2" in r.output
    assert (run / "ledger.txt").read_text().split() == ["s2", "s2"]

    # --only on a step whose own upstream (s1, s2) is done runs just that step, even
    # though unrelated prepare steps have never run
    r = runner.invoke(app, flag_args(inputs, "--only", "match.dataset.summary"))
    assert r.exit_code == 0, r.output
    assert "ran: match.dataset.summary" in r.output
    assert (run / Layout.report("match")).exists()
    # unknown steps are config errors
    r = runner.invoke(app, flag_args(inputs, "--until", "nope"))
    assert r.exit_code == 2 and "unknown step" in r.output
    r = runner.invoke(
        app, flag_args(inputs, "--until", "match.dataset.s1", "--only", "match.dataset.s2")
    )
    assert r.exit_code == 2


@needs_make
def test_dry_run_and_steps_keep_the_runs_settings(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object], tmp_path: Path
) -> None:
    """anqchor.resolved.yaml describes the settings the outputs were made with: a dry run
    or 'anqchor steps' with other settings must not overwrite it (it is used to reproduce
    the run); the next real run does."""
    run = inputs.out
    r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.s2"))
    assert r.exit_code == 0, r.output
    resolved = run / Layout.resolved_config

    def z() -> float:
        return float(yaml.safe_load(resolved.read_text())["pipelines"]["match"]["z_threshold"])

    assert z() == 1.96
    r = runner.invoke(app, flag_args(inputs, "--z-threshold", "2.5", "--dry-run"))
    assert r.exit_code == 0, r.output
    assert z() == 1.96
    cfg = tmp_path / "other.yaml"
    other = yaml.safe_load(resolved.read_text())
    other["pipelines"]["match"]["z_threshold"] = 3.0
    cfg.write_text(yaml.safe_dump(other))
    r = runner.invoke(
        app, ["--qpgen", str(inputs.qpgen), "steps", "--config", str(cfg), "--out", str(run)]
    )
    assert r.exit_code == 0, r.output
    assert z() == 1.96
    r = runner.invoke(app, flag_args(inputs, "--z-threshold", "2.5", "--until", "match.dataset.s2"))
    assert r.exit_code == 0, r.output
    assert z() == 2.5


@needs_make
def test_run_flags_an_out_of_date_report(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object]
) -> None:
    """'match run' does not render the report; when one from an earlier render is older
    than the results, the run says so instead of pointing at it as current."""
    run = inputs.out
    r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.summary"))
    assert r.exit_code == 0, r.output
    rep = run / Layout.report("match")
    assert rep.exists()  # the fake summary step writes it
    os.utime(rep, (0, 0))
    r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.summary"))
    assert r.exit_code == 0, r.output
    assert "is older than the results" in r.output and "anqchor match report --out" in r.output
    future = time.time() + 1000
    os.utime(rep, (future, future))
    r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.summary"))
    assert r.exit_code == 0, r.output
    assert f"match report: {rep}\n" in r.output


@needs_make
def test_only_refuses_when_upstream_missing(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object]
) -> None:
    r = runner.invoke(app, flag_args(inputs, "--only", "match.dataset.s2"))
    assert r.exit_code == 2
    assert "upstream steps not run" in r.output and "match.dataset.s1" in r.output


@needs_make
def test_failing_step_is_exit_3(inputs: SimpleNamespace, fake_pipelines: dict[str, object]) -> None:
    fake_pipelines["fail"] = True
    r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.s2"))
    assert r.exit_code == 3, r.output
    assert "run failed in step match.dataset.s2" in r.output
    assert str(inputs.out / "logs" / "match.dataset.s2.log") in r.output


@needs_make
def test_steps_command(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object], tmp_path: Path
) -> None:
    cfg = {
        "genotypes": {"pgen_list": str(inputs.geno)},
        "references": {"ref": {"pairs": str(inputs.pairs)}},
        "datasets": {"rna": {"pheno": str(inputs.pheno)}},
        "pipelines": {"cis": {"enabled": False}, "trans": {"enabled": False}},
        "out_dir": str(inputs.out),
    }
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(cfg))
    base = ["--qpgen", str(inputs.qpgen)]
    r = runner.invoke(app, [*base, "run", "--config", str(p), "--until", "match.rna.s1"])
    assert r.exit_code == 0, r.output
    r = runner.invoke(app, [*base, "steps", "--config", str(p)])
    assert r.exit_code == 0, r.output
    lines = {
        ln.split()[0]: ln.split()[1]
        for ln in r.output.splitlines()
        if ln.startswith(("match.", "prepare."))
    }
    assert lines["match.rna.s1"] == "done"
    assert lines["match.rna.s2"] == "missing"
    assert lines["prepare.genotypes"] == "missing"


@needs_make
def test_report_command_single_pipeline(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object], tmp_path: Path
) -> None:
    """Top-level 'anqchor report' on a one-pipeline run renders that report; there is no
    report.index step to ask for (it exists only when several pipelines run)."""
    fake_pipelines["report"] = True
    cfg = {
        "genotypes": {"pgen_list": str(inputs.geno)},
        "references": {"ref": {"pairs": str(inputs.pairs)}},
        "datasets": {"rna": {"pheno": str(inputs.pheno)}},
        "pipelines": {"cis": {"enabled": False}, "trans": {"enabled": False}},
        "out_dir": str(inputs.out),
    }
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(cfg))
    base = ["--qpgen", str(inputs.qpgen)]
    r = runner.invoke(app, [*base, "run", "--config", str(p), "--until", "match.rna.summary"])
    assert r.exit_code == 0, r.output
    for extra in ([], ["--force"]):
        r = runner.invoke(app, [*base, "report", "--config", str(p), "--pipeline", "match", *extra])
        assert r.exit_code == 0, r.output
        assert "report.index" not in r.output
        assert (inputs.out / "report" / "extra.html").exists()


def _multi_config(inputs: SimpleNamespace, tmp_path: Path) -> Path:
    cfg = {
        "genotypes": {"pgen_list": str(inputs.geno)},
        "references": {"ref": {"kind": "cis", "pairs": str(inputs.pairs)}},
        "datasets": {"rnaseq": {"pheno": str(inputs.pheno)}, "other": {"pheno": str(inputs.pheno)}},
        "pipelines": {
            "match": {"datasets": ["rnaseq"]},
            "cis": {
                "datasets": ["other"],
                "strategies": [{"name": "peer60", "cov": None}, {"name": "peer15", "cov": None}],
            },
            "trans": {"strategies": [{"name": "peer60", "cov": None}]},
        },
        "out_dir": str(inputs.out),
    }
    p = tmp_path / "multi.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return p


def test_dataset_filter_skips_pipelines_without_it(inputs: SimpleNamespace, tmp_path: Path) -> None:
    """'anqchor run --dataset X' must not force X into a pipeline configured without it."""
    p = _multi_config(inputs, tmp_path)
    cfg, active = build_config(RunOptions(config=p, dataset=["rnaseq"]), ["match", "cis", "trans"])
    assert active == ["match", "trans"]
    assert cfg.pipelines.match.datasets == ["rnaseq"] and cfg.pipelines.cis.datasets == ["other"]
    # a single, explicitly named pipeline still accepts the dataset
    cfg, active = build_config(RunOptions(config=p, dataset=["rnaseq"]), ["cis"])
    assert active == ["cis"] and cfg.pipelines.cis.datasets == ["rnaseq"]
    with pytest.raises(ConfigError, match="unknown dataset"):
        build_config(RunOptions(config=p, dataset=["nope"]), ["match", "cis"])


def test_strategy_filter_skips_pipelines_without_it(
    inputs: SimpleNamespace, tmp_path: Path
) -> None:
    p = _multi_config(inputs, tmp_path)
    cfg, active = build_config(RunOptions(config=p, strategy=["peer15"]), ["match", "cis", "trans"])
    assert active == ["match", "cis"]
    assert [s.name for s in cfg.pipelines.cis.strategies] == ["peer15"]
    cfg, active = build_config(RunOptions(config=p, strategy=["peer60"]), ["cis", "trans"])
    assert active == ["cis", "trans"]
    with pytest.raises(ConfigError, match="no trans strategy"):
        build_config(RunOptions(config=p, strategy=["peer15"]), ["trans"])
    with pytest.raises(ConfigError, match="no cis/trans strategy"):
        build_config(RunOptions(config=p, strategy=["nope"]), ["match", "cis", "trans"])


@needs_make
def test_concurrent_run_is_refused(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object]
) -> None:
    inputs.out.mkdir(parents=True)
    with run_lock(inputs.out):
        r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.s1"))
        assert r.exit_code == 2, r.output
        assert "another anqchor command is running" in r.output
    r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.s1"))
    assert r.exit_code == 0, r.output


@needs_make
def test_only_refuses_when_upstream_stale(
    inputs: SimpleNamespace, fake_pipelines: dict[str, object]
) -> None:
    run = inputs.out
    r = runner.invoke(app, flag_args(inputs, "--until", "match.dataset.s2"))
    assert r.exit_code == 0, r.output
    # an upstream stamp older than one of its prerequisites (here its signature): stale, not done
    sig = run / "steps" / "match.dataset.s1.json"
    old = sig.stat().st_mtime
    (run / stamp_path("match.dataset.s1")).touch()
    os.utime(run / stamp_path("match.dataset.s1"), (old - 100, old - 100))
    r = runner.invoke(app, flag_args(inputs, "--only", "match.dataset.s2"))
    assert r.exit_code == 2, r.output
    assert "out of date" in r.output and "match.dataset.s1" in r.output


# ------------------------------------------------------- real pipelines (MAGE)
def mage_config(out: Path) -> dict[str, object]:
    """A two-pipeline configuration over the in-repo example data (paths relative to the
    YAML, which sits next to a link to ``examples/mage/data``)."""
    return {
        "project": "mage_chr21_22",
        "out_dir": str(out),
        "genotypes": {"pgen_list": "data/geno/mage.pgen_list.tsv"},
        "references": {
            "mage_cis": {"kind": "cis", "pairs": "data/reference/gtex_lcl_eqtl.chr21_22.tsv.gz"}
        },
        "datasets": {
            "rnaseq": {
                "omics": "rnaseq",
                "pheno": "data/pheno/mage_expr.chr21_22.regenie.tsv.gz",
                "trait_annot": "data/pheno/mage_genes.chr21_22.bed",
                "cov": "data/pheno/mage_covariates.regenie.tsv.gz",
            }
        },
        "pipelines": {
            "cis": {
                "strategies": [{"name": "peer60", "cov": "default"}, {"name": "nocov", "cov": None}]
            },
            "trans": {"strategies": [{"name": "peer60", "cov": "default"}]},
        },
    }


@pytest.mark.integration
@needs_make
def test_run_dry_run_mage(mage_data: Path, qpgen_bin: Path, tmp_path: Path) -> None:
    cfg_dir = tmp_path / "cfg"
    cfg_dir.mkdir()
    # relative paths in the YAML resolve against the YAML's directory: link the data there
    os.symlink(mage_data, cfg_dir / "data")
    p = cfg_dir / "anqchor.yaml"
    out = tmp_path / "run"
    p.write_text(yaml.safe_dump(mage_config(out)))
    r = runner.invoke(app, ["--qpgen", str(qpgen_bin), "run", "--config", str(p), "--dry-run"])
    assert r.exit_code == 0, r.output
    assert (out / "Makefile").exists() and (out / Layout.resolved_config).exists()
    resolved = yaml.safe_load((out / Layout.resolved_config).read_text())
    assert Path(resolved["genotypes"]["pgen_list"]).is_absolute()
    assert "prepare.genotypes" in r.output
    assert not list((out / "steps").glob("*.done"))


@pytest.mark.integration
@needs_make
def test_match_flag_only_dry_run_mage(mage_data: Path, qpgen_bin: Path, tmp_path: Path) -> None:
    out = tmp_path / "run"
    r = runner.invoke(
        app,
        [
            "--qpgen",
            str(qpgen_bin),
            "match",
            "run",
            "--dry-run",
            "--geno-list",
            str(mage_data / "geno/mage.pgen_list.tsv"),
            "--pheno",
            str(mage_data / "pheno/mage_expr.chr21_22.regenie.tsv.gz"),
            "--pairs",
            str(mage_data / "reference/gtex_lcl_eqtl.chr21_22.tsv.gz"),
            "--out",
            str(out),
        ],
    )
    assert r.exit_code == 0, r.output
    assert (out / "Makefile").exists() and (out / Layout.resolved_config).exists()
    assert "pair-prs" in (out / "Makefile").read_text()


def test_match_report_requires_a_run(tmp_path: Path) -> None:
    """'match report' renders from an existing run only; a plain directory is exit 2."""
    r = runner.invoke(app, ["match", "report", "--out", str(tmp_path)])
    assert r.exit_code == 2
    assert "not an anqchor run directory" in r.output


def test_match_has_run_all_report_commands() -> None:
    r = runner.invoke(app, ["match", "--help"])
    assert r.exit_code == 0
    for cmd in ("run", "all", "report"):
        assert cmd in r.output
