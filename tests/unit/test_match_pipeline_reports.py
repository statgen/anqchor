"""``build_graph(reports=False)`` stops at the summaries (``anqchor match run``)."""

from __future__ import annotations

from pathlib import Path

from anqchor.config.model import Config
from anqchor.pipelines import BuildContext, build_graph, is_report_step


def _cfg(tmp_path: Path) -> Config:
    for name in ("list.tsv", "pairs.tsv", "pheno.tsv"):
        (tmp_path / name).write_text("x\n")
    return Config.model_validate(
        {
            "out_dir": str(tmp_path / "run"),
            "genotypes": {"pgen_list": str(tmp_path / "list.tsv")},
            "references": {"ref": {"pairs": str(tmp_path / "pairs.tsv")}},
            "datasets": {"ds": {"pheno": str(tmp_path / "pheno.tsv")}},
            "pipelines": {"cis": {"enabled": False}, "trans": {"enabled": False}},
        }
    )


def test_run_stage_has_no_report_steps(tmp_path: Path) -> None:
    ctx = BuildContext(cfg=_cfg(tmp_path), qpgen=Path("/bin/true"))
    graph, groups = build_graph(ctx, ["match"], reports=False)
    assert "match.ds.summarize" in graph
    assert not [s.name for s in graph if is_report_step(s)]
    assert all(not n.startswith("report.") for n in groups["match"])


def test_all_stage_has_the_match_report_but_no_index_for_one_pipeline(tmp_path: Path) -> None:
    ctx = BuildContext(cfg=_cfg(tmp_path), qpgen=Path("/bin/true"))
    graph, _ = build_graph(ctx, ["match"], reports=True, index=False)
    assert "report.match" in graph and "report.index" not in graph
    assert "has_index" not in graph.get("report.match").params


def test_reports_link_to_the_index_only_when_the_run_has_one(tmp_path: Path) -> None:
    ctx = BuildContext(cfg=_cfg(tmp_path), qpgen=Path("/bin/true"))
    graph, _ = build_graph(ctx, ["match"], reports=True, index=True)
    assert "report.index" in graph
    assert graph.get("report.match").params["has_index"] is True
    assert "has_index" not in graph.get("report.index").params
