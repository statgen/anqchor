"""pipelines.cis: step graph for the pair and rect engines."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from anqchor.config.model import Config
from anqchor.engine.step import StepGraph
from anqchor.pipelines import BuildContext, cis
from anqchor.pipelines.layout import Layout as L


def make_cfg(
    tmp_path: Path,
    engine: str = "pair",
    strategies: list[dict[str, Any]] | None = None,
    rint: bool = False,
) -> Config:
    return Config.model_validate(
        {
            "out_dir": str(tmp_path / "run"),
            "genotypes": {"pgen_list": str(tmp_path / "geno.tsv")},
            "references": {"ref": {"pairs": str(tmp_path / "pairs.tsv")}},
            "datasets": {"ds": {"pheno": str(tmp_path / "p.tsv"), "cov": str(tmp_path / "c.tsv")}},
            "pipelines": {
                "match": {"enabled": False},
                "trans": {"enabled": False},
                "cis": {
                    "assoc_engine": engine,
                    "rint": rint,
                    "strategies": strategies
                    or [
                        {"name": "all", "cov": "default"},
                        {"name": "few", "cov": {"include": ["PC1"]}},
                        {"name": "nocov", "cov": None},
                    ],
                },
            },
        }
    )


def build(tmp_path: Path, engine: str = "pair", **kw: Any) -> StepGraph:
    ctx = BuildContext(cfg=make_cfg(tmp_path, engine, **kw), qpgen=Path("/opt/qpgentools"))
    return StepGraph(cis.build(ctx))


def test_pair_engine_steps(tmp_path: Path) -> None:
    g = build(tmp_path)
    names = g.names()
    for n in (
        "prepare.genotypes",
        "prepare.reference.ref",
        "prepare.dataset.ds",
        "prepare.mapping.ds.ref",
        "cis.ds.all.inputs",
        "cis.ds.all.assoc",
        "cis.ds.few.assoc",
        "cis.ds.nocov.assoc",
        "cis.ds.summarize",
        "report.cis",
    ):
        assert n in names
    a = g.get("cis.ds.all.assoc")
    cmd = a.commands[0]
    # phenotypes as given by default: pair-assoc skips its RINT columns
    assert "pair-assoc" in cmd and "--skip-rint" in cmd
    assert f"--cov {L.strategy_cov('cis', 'ds', 'all')}" in cmd
    assert f"--pairs {L.dsref_pairs('ds', 'ref')}" in cmd
    assert f"--out {L.cis_pair_assoc('ds', 'all')}" in cmd
    assert "--pgen-list prepare/genotypes/pgen_list.tsv" in cmd
    nocov = g.get("cis.ds.nocov.assoc")
    assert "--cov" not in nocov.commands[0]
    # the association depends on the strategy covariates and the prepared pairs
    assert set(g.upstream("cis.ds.all.assoc")) >= {
        "cis.ds.all.inputs",
        "prepare.mapping.ds.ref",
        "prepare.genotypes",
    }
    assert "cis.ds.nocov.inputs" not in g.upstream("cis.ds.nocov.assoc")


def test_rint_run(tmp_path: Path) -> None:
    """--rint: pair-assoc writes its RINT columns and the summary reads them."""
    g = build(tmp_path, rint=True)
    for s in ("all", "few", "nocov"):
        assert "--skip-rint" not in g.get(f"cis.ds.{s}.assoc").commands[0]
    assert g.get("cis.ds.summarize").params["rint"] is True
    # the flag changes the association and summary signatures, not the inputs steps
    g0 = build(tmp_path)
    changed = {n for n in g0.names() if g0.get(n).signature() != g.get(n).signature()}
    assert changed == {
        "cis.ds.all.assoc",
        "cis.ds.few.assoc",
        "cis.ds.nocov.assoc",
        "cis.ds.summarize",
    }


def test_summarize_and_report_steps(tmp_path: Path) -> None:
    g = build(tmp_path)
    s = g.get("cis.ds.summarize")
    assert s.kind == "summarize_cis"
    assert s.outputs == tuple(L.cis_summary("ds", n) for n in cis.SUMMARY_FILES)
    assert (
        L.cis_pair_assoc("ds", "nocov") in s.inputs
        and L.strategy_json("cis", "ds", "nocov") in s.inputs
    )
    assert L.ds_traits("ds") not in s.inputs  # only the rect engine needs trait loci
    p = s.params
    assert (
        p["engine"] == "pair"
        and p["levels"]["strict"] == 1e-5
        and p["max_unresolved_fraction"] == 0.2
    )
    assert p["rint"] is False
    assert [x["name"] for x in p["strategies"]] == ["all", "few", "nocov"]
    # several strategies (multi-pipeline configuration): labelled by strategy name
    assert [x["label"] for x in p["strategies"]] == ["all", "few", "nocov"]
    assert p["strategies"][2]["file"] == str(L.cis_pair_assoc("ds", "nocov"))
    assert "rint" in p["transform_semantics"] and "as_given" in p["transform_semantics"]
    # shown in summary.json and the report: refers to the qpgen release, not dev commits
    assert "dev" not in json.dumps(p["transform_semantics"])
    r = g.get("report.cis")
    assert r.kind == "report_cis" and r.outputs == (L.report("cis"),)
    assert set(g.upstream("report.cis")) == {"cis.ds.summarize"}


def test_single_strategy_label_is_run_directory(tmp_path: Path) -> None:
    g = build(tmp_path, strategies=[{"name": "default"}])
    assert g.get("cis.ds.summarize").params["strategies"][0]["label"] == "run"  # out_dir name


@pytest.mark.parametrize("rint", [False, True])
def test_rect_engine_steps(tmp_path: Path, rint: bool) -> None:
    g = build(tmp_path, "rect", rint=rint)
    t, other = ("rint", "raw") if rint else ("raw", "rint")
    names = g.names()
    assert f"cis.ds.all.rect.{t}" in names and f"cis.ds.all.rect.{other}" not in names
    cmd = g.get(f"cis.ds.all.rect.{t}").commands[0]
    assert "rect-assoc" in cmd and ("--rint-after-adj" in cmd) == rint
    assert f"--var-list {L.dsref_variants('ds', 'ref')}" in cmd
    assert f"--pheno-list {L.dsref_traits('ds', 'ref')}" in cmd
    assert f"--out {L.rect_assoc('cis', 'ds', 'all', t)}" in cmd
    s = g.get("cis.ds.summarize")
    assert L.ds_traits("ds") in s.inputs
    assert s.params["strategies"][0]["file"] == str(L.rect_assoc("cis", "ds", "all", t))


def test_change_of_levels_changes_only_summary_signature(tmp_path: Path) -> None:
    g1 = build(tmp_path)
    cfg = make_cfg(tmp_path)
    cfg.pipelines.cis.levels.strict = 1e-6
    g2 = StepGraph(cis.build(BuildContext(cfg=cfg, qpgen=Path("/opt/qpgentools"))))
    changed = [n for n in g1.names() if g1.get(n).signature() != g2.get(n).signature()]
    assert changed == ["cis.ds.summarize"]


@pytest.mark.parametrize("engine", ["pair", "rect"])
def test_steps_dedupe_with_shared_prepare(tmp_path: Path, engine: str) -> None:
    ctx = BuildContext(cfg=make_cfg(tmp_path, engine), qpgen=Path("/opt/qpgentools"))
    g = StepGraph(cis.build(ctx))
    g.extend(cis.build(ctx))  # identical definitions are a no-op
    assert len(g) == len(cis.build(ctx))
