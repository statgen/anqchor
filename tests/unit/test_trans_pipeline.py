"""pipelines.trans: step graph, qpgen command lines and step parameters."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from anqchor.config.model import Config
from anqchor.engine.step import StepGraph
from anqchor.pipelines import BuildContext, trans
from anqchor.pipelines.layout import Layout as L


def make_cfg(tmp_path: Path, **trans_cfg: Any) -> Config:
    t = {"strategies": [{"name": "all", "cov": "default"}], **trans_cfg}
    return Config.model_validate(
        {
            "out_dir": str(tmp_path / "run"),
            "genotypes": {"pgen_list": str(tmp_path / "geno.tsv")},
            "references": {"ref": {"pairs": str(tmp_path / "pairs.tsv")}},
            "datasets": {"ds": {"pheno": str(tmp_path / "p.tsv"), "cov": str(tmp_path / "c.tsv")}},
            "pipelines": {"match": {"enabled": False}, "cis": {"enabled": False}, "trans": t},
        }
    )


def build(tmp_path: Path, **trans_cfg: Any) -> StepGraph:
    ctx = BuildContext(cfg=make_cfg(tmp_path, **trans_cfg), qpgen=Path("/opt/qpgentools"))
    return StepGraph(trans.build(ctx))


def test_steps_and_rect_commands(tmp_path: Path) -> None:
    g = build(tmp_path)
    names = g.names()
    for n in (
        "prepare.genotypes",
        "prepare.reference.ref",
        "prepare.dataset.ds",
        "prepare.mapping.ds.ref",
        "trans.ds.lists",
        "trans.ds.all.inputs",
        "trans.ds.all.rect.raw",
        "trans.ds.summarize",
        "report.trans",
    ):
        assert n in names
    assert "trans.ds.all.rect.rint" not in names  # one transform per run
    raw = g.get("trans.ds.all.rect.raw").commands[0]
    assert "rect-assoc" in raw and "--rint-after-adj" not in raw
    assert f"--var-list {L.trans_variants('ds')}" in raw
    assert f"--pheno-list {L.trans_traits('ds')}" in raw  # traits: reference
    assert f"--cov {L.strategy_cov('trans', 'ds', 'all')}" in raw
    assert f"--out {L.rect_assoc('trans', 'ds', 'all', 'raw')}" in raw
    assert "--min-maf 0.01" in raw
    assert set(g.upstream("trans.ds.all.rect.raw")) >= {
        "trans.ds.lists",
        "trans.ds.all.inputs",
        "prepare.genotypes",
    }


def test_rint_is_one_run_with_rint_after_adj(tmp_path: Path) -> None:
    g = build(tmp_path, rint=True)
    assert "trans.ds.all.rect.raw" not in g.names()
    step = g.get("trans.ds.all.rect.rint")
    assert "--rint-after-adj" in step.commands[0]
    assert step.outputs == (L.rect_assoc("trans", "ds", "all", "rint"),)
    s = g.get("trans.ds.summarize")
    assert s.params["rint"] is True and s.params["rect"] == str(
        L.rect_assoc("trans", "ds", "all", "rint")
    )
    assert L.rect_assoc("trans", "ds", "all", "rint") in s.inputs


def test_lists_step_declares_restriction_inputs(tmp_path: Path) -> None:
    g = build(tmp_path)
    s = g.get("trans.ds.lists")
    assert s.kind == "trans_lists"
    assert L.ds_traits("ds") in s.inputs and L.geno_json in s.inputs
    assert L.ref_variants_present("ref") in s.inputs and L.dsref_traits("ds", "ref") in s.inputs
    assert s.outputs == (L.trans_variants("ds"), L.trans_traits("ds"), L.trans_lists_json("ds"))
    assert s.params["trait_chroms"] == "genotyped" and s.params["traits"] == "reference"


def test_trait_chroms_controls_pheno_list(tmp_path: Path) -> None:
    # traits: all + genotyped (default): restricted, so qpgen gets --pheno-list
    g = build(tmp_path, traits="all")
    cmd = g.get("trans.ds.all.rect.raw").commands[0]
    assert "--pheno-list" in cmd
    assert L.trans_traits("ds") in g.get("trans.ds.all.rect.raw").inputs
    assert L.ds_pheno("ds") in g.get("trans.ds.lists").inputs
    # traits: all + trait_chroms: all: the whole matrix, no --pheno-list
    g = build(tmp_path, traits="all", trait_chroms="all")
    step = g.get("trans.ds.all.rect.raw")
    assert "--pheno-list" not in step.commands[0]
    assert L.trans_traits("ds") not in step.inputs
    # an explicit list (YAML integers) arrives as strings
    g = build(tmp_path, traits="all", trait_chroms=[21, "chr22"])
    assert g.get("trans.ds.lists").params["trait_chroms"] == ["21", "chr22"]
    assert g.get("trans.ds.summarize").params["trait_chroms"] == ["21", "chr22"]
    assert "--pheno-list" in g.get("trans.ds.all.rect.raw").commands[0]


def test_summarize_and_report_steps(tmp_path: Path) -> None:
    g = build(tmp_path, hotspot_min_traits=3, hotspot_bin_bp=500_000, levels={"nominal": 1e-6})
    s = g.get("trans.ds.summarize")
    assert s.kind == "summarize_trans"
    assert s.outputs == tuple(L.trans_summary("ds", n) for n in trans.SUMMARY_FILES)
    assert L.trans_lists_json("ds") in s.inputs and L.ds_traits("ds") in s.inputs
    assert L.dsref_pairs("ds", "ref") in s.inputs
    p = s.params
    assert (
        p["run"] == "run" and p["strategy"] == "all" and p["rint"] is False
    )  # run = directory name
    assert p["hotspot_min_traits"] == 3 and p["hotspot_bin_bp"] == 500_000
    assert p["levels"] == {"nominal": 1e-6, "strict": None, "bonferroni": True}
    assert p["known_pairs"] is True and p["schema_version"] == 2
    r = g.get("report.trans")
    assert r.kind == "report_trans" and r.outputs == (L.report("trans"),)
    assert r.params["reference_file"] == "pairs.tsv" and r.params["covariates"] == "all"
    assert set(g.upstream("report.trans")) >= {"trans.ds.summarize"}


def test_known_pairs_off_drops_pairs_input(tmp_path: Path) -> None:
    g = build(tmp_path, known_pairs=None)
    s = g.get("trans.ds.summarize")
    assert s.params["known_pairs"] is False
    assert L.dsref_pairs("ds", "ref") not in s.inputs


def test_trans_cli_flags_and_report_command(tmp_path: Path) -> None:
    """--rint and --hotspot-bin-bp reach the configuration; 'trans report' needs a run directory."""
    from typer.testing import CliRunner

    from anqchor.cli.main import app

    runner = CliRunner()
    r = runner.invoke(app, ["trans", "run", "--help"], env={"COLUMNS": "200"})
    assert (
        "--rint" in r.output
        and "--hotspot-bin-bp" in r.output
        and "--hotspot-min-traits" in r.output
    )
    r = runner.invoke(app, ["trans", "report", "--out", str(tmp_path)])
    assert r.exit_code == 2 and "not an anqchor run directory" in r.output
