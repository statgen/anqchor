"""Run-directory layout and the shared step definitions built from it."""

from __future__ import annotations

from pathlib import Path

import pytest

from anqchor.config import load_config
from anqchor.config.model import StrategyConfig
from anqchor.engine.step import StepGraph
from anqchor.pipelines.common import (
    prepare_dataset,
    prepare_genotypes,
    prs,
    read_pgen_list_files,
    shared_prepare_steps,
    strategy_files,
    strategy_inputs,
)
from anqchor.pipelines.context import BuildContext
from anqchor.pipelines.layout import Layout as L
from anqchor.pipelines.report import index_step


def test_paths_are_relative_and_distinct() -> None:
    paths = [
        L.resolved_config,
        L.run_json,
        L.makefile,
        L.pgen_list,
        L.geno_samples,
        L.geno_json,
        L.ref_pairs("r"),
        L.ref_variants("r"),
        L.ref_variants_present("r"),
        L.ref_json("r"),
        L.ds_pheno("d"),
        L.ds_cov("d"),
        L.ds_cov_json("d"),
        L.ds_traits("d"),
        L.ds_samples("d"),
        L.ds_sample_map("d"),
        L.ds_json("d"),
        L.dsref_trait_map("d", "r"),
        L.dsref_pairs("d", "r"),
        L.dsref_variants("d", "r"),
        L.dsref_traits("d", "r"),
        L.dsref_json("d", "r"),
        L.prs_matrix("r"),
        L.prs_se("r"),
        *L.match_qpgen_outputs("d"),
        L.match_trait_map("d"),
        L.match_params("d"),
        L.strategy_cov("cis", "d", "s"),
        L.strategy_pheno("cis", "d", "s"),
        L.strategy_json("cis", "d", "s"),
        L.cis_pair_assoc("d", "s"),
        L.rect_assoc("trans", "d", "s", "raw"),
        L.trans_variants("d"),
        L.trans_traits("d"),
        L.report("index"),
    ]
    assert all(not p.is_absolute() for p in paths)
    assert len(set(paths)) == len(paths)


def test_qpgen_prefix_naming() -> None:
    """qpgen appends suffixes to --out prefixes; the layout must agree with it."""
    assert L.prs_matrix("r") == Path(str(L.prs_prefix("r")) + ".prs.tsv.gz")
    assert L.prs_se("r") == Path(str(L.prs_prefix("r")) + ".se.tsv.gz")
    w, a, al = L.match_qpgen_outputs("d")
    pre = str(L.match_prefix("d"))
    assert (w, a, al) == (
        Path(pre + ".weights.tsv.gz"),
        Path(pre + ".match.assigned.tsv.gz"),
        Path(pre + ".match.all.tsv.gz"),
    )
    assert L.strategy_dir("trans", "d", "s") == Path("trans/d/strategies/s")
    assert L.report("cis") == Path("report/cis.html")


@pytest.fixture
def ctx(tmp_path: Path) -> BuildContext:
    data = tmp_path / "data"
    (data / "geno").mkdir(parents=True)
    for suffix in (".pgen", ".psam", ".pvar.idx.gz", ".pvar.idx.gz.tbi"):
        (data / "geno" / f"chr22{suffix}").write_text("")
    (data / "list.tsv").write_text("22\t1\t1000000000\tgeno/chr22\n")
    cfg = load_config(
        None,
        {
            "genotypes": {"pgen_list": str(data / "list.tsv")},
            "references": {"ref": {"pairs": str(data / "pairs.tsv")}},
            "datasets": {
                "rna": {
                    "pheno": str(data / "p.tsv"),
                    "cov": str(data / "c.tsv"),
                    "cov_categorical": ["sex"],
                },
                "prot": {"pheno": str(data / "q.tsv")},
            },
            "pipelines": {"match": {"prs_samples": str(data / "ids.txt")}},
            "out_dir": str(tmp_path / "run"),
        },
    )
    return BuildContext(cfg=cfg, qpgen=Path("/opt/qpgentools"))


def test_read_pgen_list_files(ctx: BuildContext) -> None:
    files = read_pgen_list_files(ctx.cfg.genotypes.pgen_list)  # type: ignore[arg-type]
    names = sorted(p.name for p in files)
    assert names == [
        "chr22.pgen",
        "chr22.psam",
        "chr22.pvar.idx.gz",
    ]  # .pgen too: genotype changes re-run
    assert all(p.is_absolute() for p in files)
    assert read_pgen_list_files(Path("/nonexistent/list.tsv")) == []


def test_shared_steps(ctx: BuildContext) -> None:
    steps = shared_prepare_steps(ctx, ["rna", "prot"], ["ref"])
    names = [s.name for s in steps]
    assert names == [
        "prepare.genotypes",
        "prepare.reference.ref",
        "prepare.dataset.rna",
        "prepare.dataset.prot",
        "prepare.mapping.rna.ref",
        "prepare.mapping.prot.ref",
    ]
    g = StepGraph(steps)
    g.extend(shared_prepare_steps(ctx, ["rna"], ["ref"]))  # shared steps de-duplicate
    assert len(g) == 6
    assert g.upstream("prepare.mapping.rna.ref") == ["prepare.dataset.rna", "prepare.reference.ref"]
    geno = prepare_genotypes(ctx)
    assert geno.kind == "prepare_genotypes" and geno.outputs == (
        L.pgen_list,
        L.geno_samples,
        L.geno_json,
    )
    assert ctx.cfg.genotypes.pgen_list in geno.inputs
    rna = prepare_dataset(ctx, "rna")
    assert L.ds_cov("rna") in rna.outputs and rna.params["cov_categorical"] == ["sex"]
    assert L.ds_cov("prot") not in prepare_dataset(ctx, "prot").outputs


def test_prs_step(ctx: BuildContext) -> None:
    s = prs(ctx, "ref")
    assert s.kind is None and len(s.commands) == 1
    cmd = s.commands[0]
    assert cmd.startswith("/opt/qpgentools pair-prs --pgen-list prepare/genotypes/pgen_list.tsv")
    assert "--out prs/ref/prs" in cmd and "--sample " in cmd
    assert s.outputs == (L.prs_matrix("ref"), L.prs_se("ref"))
    assert ctx.cfg.pipelines.match.prs_samples in s.inputs


def test_strategy_inputs_and_files(ctx: BuildContext) -> None:
    default = StrategyConfig(name="all")
    s = strategy_inputs(ctx, "cis", "rna", default)
    assert s.name == "cis.rna.all.inputs"
    assert strategy_files(ctx, "cis", "rna", default) == (
        L.ds_pheno("rna"),
        L.strategy_cov("cis", "rna", "all"),
    )
    nocov = StrategyConfig(name="nocov", cov=None)
    assert strategy_files(ctx, "cis", "rna", nocov) == (L.ds_pheno("rna"), None)
    # 'default' on a dataset without covariates means none
    assert strategy_files(ctx, "trans", "prot", default) == (L.ds_pheno("prot"), None)
    sub = StrategyConfig.model_validate({"name": "pcs", "cov": {"include": ["PC1-2"]}})
    assert strategy_files(ctx, "cis", "rna", sub)[1] == L.strategy_cov("cis", "rna", "pcs")
    with pytest.raises(ValueError, match=r"cov\.include needs covariates"):
        strategy_inputs(ctx, "cis", "prot", sub)
    alt = StrategyConfig.model_validate({"name": "alt", "pheno": "/x/alt.tsv", "cov": None})
    assert strategy_files(ctx, "cis", "rna", alt) == (L.strategy_pheno("cis", "rna", "alt"), None)


def test_index_step(ctx: BuildContext) -> None:
    s = index_step(ctx, ["match"])
    assert s.outputs == (L.report("index"),)
    assert L.match_summary("rna", "summary.json") in s.inputs and L.report("match") in s.inputs
    assert s.params["datasets"] == {"match": ["rna", "prot"]}


def test_qpgen_identity_in_signatures(ctx: BuildContext, tmp_path: Path) -> None:
    """A rebuilt qpgentools at the same path must make every qpgen step stale."""
    import dataclasses

    from anqchor import __version__
    from anqchor.engine.minimake import write_pipeline_makefile
    from anqchor.pipelines import QPGEN_ID_PARAM, build_graph

    def graph(sha: str) -> StepGraph:
        return build_graph(dataclasses.replace(ctx, qpgen_sha256=sha), ["match", "cis", "trans"])[0]

    g1, g2 = graph("a" * 64), graph("b" * 64)
    tagged = {s.name for s in g1 if QPGEN_ID_PARAM in s.params}
    assert tagged == {s.name for s in g1 if s.kind is None or "qpgen" in s.params}
    assert {"prs.ref", "match.rna.qpgen"} <= tagged
    assert any(n.startswith("cis.") and n.endswith(".assoc") for n in tagged)
    assert any(n.startswith("trans.") and ".rect." in n for n in tagged)
    assert not any(
        n.startswith(("prepare.", "report.")) or n.endswith(".summarize") for n in tagged
    )
    changed = {s.name for s in g1 if s.signature() != g2.get(s.name).signature()}
    assert changed == tagged
    # Python steps record the anqchor version (an upgrade re-runs them)
    assert g1.get("report.index").signature()["anqchor_version"] == __version__
    assert "anqchor_version" not in g1.get("prs.ref").signature()
    # end to end: rewriting the Makefile with the new hash drops the qpgen and downstream stamps
    run = tmp_path / "mk"
    write_pipeline_makefile(g1, run, {})
    for s in g1:
        (run / "steps" / f"{s.name}.done").write_text("")
    write_pipeline_makefile(g2, run, {})
    gone = {s.name for s in g1 if not (run / "steps" / f"{s.name}.done").exists()}
    assert "prs.ref" in gone and "match.rna.summarize" in gone and "report.match" in gone
    assert "prepare.genotypes" not in gone
