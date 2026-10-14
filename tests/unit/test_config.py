"""Configuration loading, path resolution, overrides and validation."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import pytest
import yaml

from anqchor.config import ConfigError, CovSpec, load_config
from anqchor.config.load import (
    config_to_dict,
    deep_merge,
    load_resolved_config,
    write_resolved_config,
)
from anqchor.config.model import Config


def _write(path: Path, data: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def _base(**extra: Any) -> dict[str, Any]:
    data: dict[str, Any] = {
        "project": "t",
        "genotypes": {"pgen_list": "geno/list.tsv"},
        "references": {"ref": {"pairs": "pair/ref.tsv.gz"}},
        "datasets": {"rna": {"pheno": "exprs/p.tsv.gz", "cov": "exprs/c.tsv.gz"}},
    }
    data.update(extra)
    return data


@pytest.fixture
def chdir(tmp_path: Path):
    """Run the test from a separate working directory and restore the old one."""
    old = Path.cwd()
    work = tmp_path / "work"
    work.mkdir()
    os.chdir(work)
    yield work
    os.chdir(old)


# ----------------------------------------------------------------- path resolution
def test_relative_paths_resolve_against_yaml_dir(tmp_path: Path, chdir: Path) -> None:
    cfg_path = _write(tmp_path / "proj" / "cfg.yaml", _base(out_dir="results"))
    cfg = load_config(cfg_path)
    proj = (tmp_path / "proj").resolve()
    assert cfg.genotypes.pgen_list == proj / "geno/list.tsv"
    assert cfg.references["ref"].pairs == proj / "pair/ref.tsv.gz"
    assert cfg.datasets["rna"].pheno == proj / "exprs/p.tsv.gz"
    assert cfg.datasets["rna"].cov == proj / "exprs/c.tsv.gz"
    assert cfg.out_dir == proj / "results"


def test_default_out_dir_is_next_to_config(tmp_path: Path, chdir: Path) -> None:
    cfg = load_config(_write(tmp_path / "proj" / "cfg.yaml", _base()))
    assert cfg.out_dir == (tmp_path / "proj" / "anqchor_out").resolve()


def test_absolute_and_home_paths_kept(tmp_path: Path, chdir: Path) -> None:
    data = _base()
    data["genotypes"] = {"pfile": "/abs/geno/chr"}
    data["datasets"]["rna"]["trait_annot"] = "~/annot.bed"
    cfg = load_config(_write(tmp_path / "cfg.yaml", data))
    assert cfg.genotypes.pfile == Path("/abs/geno/chr")
    assert cfg.datasets["rna"].trait_annot == Path("~/annot.bed").expanduser()


def test_flag_overrides_resolve_against_cwd_and_win(tmp_path: Path, chdir: Path) -> None:
    cfg_path = _write(tmp_path / "proj" / "cfg.yaml", _base(out_dir="results"))
    overrides = {
        "out_dir": "here_out",
        "datasets": {"rna": {"pheno": "flag_pheno.tsv"}},
        "pipelines": {"match": {"z_threshold": 3.0}},
    }
    cfg = load_config(cfg_path, overrides)
    assert cfg.out_dir == chdir.resolve() / "here_out"
    assert cfg.datasets["rna"].pheno == chdir.resolve() / "flag_pheno.tsv"
    # untouched file values survive the merge, still relative to the YAML directory
    assert cfg.datasets["rna"].cov == (tmp_path / "proj" / "exprs/c.tsv.gz").resolve()
    assert cfg.pipelines.match.z_threshold == 3.0


def test_none_overrides_do_not_clobber(tmp_path: Path, chdir: Path) -> None:
    cfg = load_config(
        _write(tmp_path / "cfg.yaml", _base(out_dir="r")),
        {"out_dir": None, "qpgen": {"threads": None}},
    )
    assert cfg.out_dir == (tmp_path / "r").resolve()
    assert cfg.qpgen.threads == 1


def test_genotype_flag_replaces_file_section(tmp_path: Path, chdir: Path) -> None:
    """pgen_list and pfile are alternatives: a --pfile flag must not merge with a
    pgen_list from the file (that would fail 'exactly one')."""
    cfg = load_config(_write(tmp_path / "cfg.yaml", _base()), {"genotypes": {"pfile": "g/chr"}})
    assert cfg.genotypes.pgen_list is None
    assert cfg.genotypes.pfile == chdir.resolve() / "g/chr"


def test_flag_only_resolves_against_cwd(chdir: Path) -> None:
    cfg = load_config(None, _base())
    assert cfg.genotypes.pgen_list == chdir.resolve() / "geno/list.tsv"
    assert cfg.out_dir == chdir.resolve() / "anqchor_out"


def test_deep_merge() -> None:
    base = {"a": {"b": 1, "c": [1]}, "d": 1}
    out = deep_merge(base, {"a": {"b": 2}, "d": None, "e": {"f": 1}})
    assert out == {"a": {"b": 2, "c": [1]}, "d": 1, "e": {"f": 1}}
    assert base == {"a": {"b": 1, "c": [1]}, "d": 1}  # not mutated


# ------------------------------------------------------------------- validation
def test_unknown_keys_fail(tmp_path: Path) -> None:
    data = _base()
    data["datasets"]["rna"]["phenoo"] = "x"
    with pytest.raises(ConfigError, match=r"datasets\.rna\.phenoo"):
        load_config(_write(tmp_path / "cfg.yaml", data))
    with pytest.raises(ConfigError, match=r"pipelines\.match\.zthreshold"):
        load_config(_write(tmp_path / "c2.yaml", _base(pipelines={"match": {"zthreshold": 1}})))


def test_missing_file_and_bad_yaml(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="not found"):
        load_config(tmp_path / "nope.yaml")
    bad = tmp_path / "bad.yaml"
    bad.write_text("a: [1, 2\n")
    with pytest.raises(ConfigError, match="cannot parse YAML"):
        load_config(bad)
    lst = tmp_path / "list.yaml"
    lst.write_text("- 1\n- 2\n")
    with pytest.raises(ConfigError, match="mapping"):
        load_config(lst)
    with pytest.raises(ConfigError, match="cannot read"):
        load_config(tmp_path)  # a directory


def test_genotypes_exactly_one(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="exactly one"):
        load_config(_write(tmp_path / "c.yaml", _base(genotypes={"pgen_list": "a", "pfile": "b"})))
    with pytest.raises(ConfigError, match="exactly one"):
        load_config(_write(tmp_path / "c2.yaml", _base(genotypes={})))


def test_invalid_names(tmp_path: Path) -> None:
    data = _base()
    data["datasets"] = {"bad name": {"pheno": "p"}}
    with pytest.raises(ConfigError, match="invalid dataset name"):
        load_config(_write(tmp_path / "c.yaml", data))
    with pytest.raises(ConfigError, match="at least one entry"):
        load_config(_write(tmp_path / "c2.yaml", _base(references={})))


def test_cross_reference_errors(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="unknown dataset 'nope'"):
        load_config(_write(tmp_path / "a.yaml", _base(pipelines={"cis": {"datasets": ["nope"]}})))
    with pytest.raises(ConfigError, match="unknown reference 'nope'"):
        load_config(_write(tmp_path / "b.yaml", _base(pipelines={"match": {"reference": "nope"}})))
    two = _base()
    two["references"]["ref2"] = {"pairs": "r2.tsv"}
    with pytest.raises(ConfigError, match="reference must be set"):
        load_config(_write(tmp_path / "c.yaml", two))
    # a cis/match pipeline needs a reference with pairs
    trans_only = _base(references={"tr": {"kind": "trans", "variants": "v.txt"}})
    with pytest.raises(ConfigError, match=r"pipelines\.cis needs a reference with 'pairs'"):
        load_config(_write(tmp_path / "d.yaml", trans_only))


def test_reference_choice_only_required_for_enabled_pipelines(tmp_path: Path) -> None:
    data = _base()
    data["references"]["ref2"] = {"kind": "trans", "variants": "v.txt"}
    data["pipelines"] = {
        "match": {"reference": "ref"},
        "cis": {"reference": "ref"},
        "trans": {"enabled": False},
    }
    cfg = load_config(_write(tmp_path / "c.yaml", data))
    assert cfg.reference_for("cis") == "ref"
    assert cfg.enabled_pipelines() == ["match", "cis"]


def test_reference_inputs_required(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="a cis reference needs 'pairs'"):
        load_config(_write(tmp_path / "a.yaml", _base(references={"r": {"kind": "cis"}})))
    with pytest.raises(ConfigError, match="needs 'pairs' or 'variants'"):
        load_config(_write(tmp_path / "b.yaml", _base(references={"r": {"kind": "trans"}})))


@pytest.mark.parametrize(
    ("value", "expected"),
    [("auto", "auto"), (0.5, 0.5), (0, 0.0), (1, 1.0), ("0.25", 0.25)],
)
def test_lambda_auto_or_float(tmp_path: Path, value: Any, expected: Any) -> None:
    cfg = load_config(_write(tmp_path / "c.yaml", _base(pipelines={"match": {"lambda": value}})))
    assert cfg.pipelines.match.lambda_ == expected


@pytest.mark.parametrize("value", [1.5, -0.1, "fast", True, float("nan")])
def test_lambda_invalid(tmp_path: Path, value: Any) -> None:
    with pytest.raises(ConfigError, match="lambda"):
        load_config(None, _base(pipelines={"match": {"lambda": value}}))


def test_lambda_alias_round_trip(tmp_path: Path) -> None:
    cfg = load_config(None, _base(pipelines={"match": {"lambda": 0.3, "mahalanobis": True}}))
    d = config_to_dict(cfg)
    assert d["pipelines"]["match"]["lambda"] == 0.3
    assert "lambda_" not in d["pipelines"]["match"]


def test_prs_samples_forms(tmp_path: Path, chdir: Path) -> None:
    cfg = load_config(None, _base(pipelines={"match": {"prs_samples": "all"}}))
    assert cfg.pipelines.match.prs_samples is None
    cfg = load_config(None, _base(pipelines={"match": {"prs_samples": "file:ids.txt"}}))
    assert cfg.pipelines.match.prs_samples == chdir.resolve() / "ids.txt"
    cfg = load_config(None, _base(pipelines={"match": {"prs_samples": "ids.txt"}}))
    assert cfg.pipelines.match.prs_samples == chdir.resolve() / "ids.txt"


# ------------------------------------------------------------------- strategies
def test_strategy_cov_forms(tmp_path: Path) -> None:
    strategies = [
        {"name": "peer60", "cov": "default"},
        {"name": "peer15", "cov": {"include": ["sex", "PC1-5", "PEER1-15"]}},
        {"name": "nocov", "cov": None},
        {"name": "other", "cov": {"file": "covs/other.tsv", "format": "tsv-sample-col"}},
        {"name": "short", "cov": "covs/short.tsv"},
        {"name": "prefixed", "cov": "file:covs/pre.tsv"},
        {"name": "implicit"},
    ]
    cfg = load_config(
        _write(tmp_path / "c.yaml", _base(pipelines={"cis": {"strategies": strategies}}))
    )
    s = {x.name: x for x in cfg.pipelines.cis.strategies}
    assert s["peer60"].cov == "default"
    assert isinstance(s["peer15"].cov, CovSpec)
    assert s["peer15"].cov.include == ["sex", "PC1-5", "PEER1-15"]
    assert s["peer15"].cov.file is None
    assert s["nocov"].cov is None
    assert isinstance(s["other"].cov, CovSpec)
    assert s["other"].cov.file == (tmp_path / "covs/other.tsv").resolve()
    assert s["other"].cov.format == "tsv-sample-col"
    assert isinstance(s["short"].cov, CovSpec)
    assert s["short"].cov.file == (tmp_path / "covs/short.tsv").resolve()
    assert isinstance(s["prefixed"].cov, CovSpec)
    assert s["prefixed"].cov.file == (tmp_path / "covs/pre.tsv").resolve()
    assert s["implicit"].cov == "default"
    assert not hasattr(s["implicit"], "transforms")  # a cis run tests one transform
    assert cfg.pipelines.cis.rint is False


def test_strategy_errors(tmp_path: Path) -> None:
    def load(strategies: Any) -> Config:
        return load_config(None, _base(pipelines={"cis": {"strategies": strategies}}))

    with pytest.raises(ConfigError, match="duplicate strategy names: a"):
        load([{"name": "a"}, {"name": "a"}])
    with pytest.raises(ConfigError, match="at least one strategy"):
        load([])
    with pytest.raises(ConfigError, match="exactly one of 'include' or 'file'"):
        load([{"name": "a", "cov": {"include": ["x"], "file": "f"}}])
    with pytest.raises(ConfigError, match="exactly one of 'include' or 'file'"):
        load([{"name": "a", "cov": {}}])
    with pytest.raises(ConfigError, match="invalid strategy name"):
        load([{"name": "a/b"}])


@pytest.mark.parametrize("transforms", [["raw", "rint"], ["rint"], ["raw"]])
def test_cis_transforms_list_is_rejected_with_guidance(transforms: list[str]) -> None:
    """Old cis configurations listed the transforms; a cis run now tests one (pipelines.cis.rint)."""
    with pytest.raises(
        ConfigError, match=r"no longer accepted.*pipelines\.cis\.rint: true.*cis compare"
    ):
        load_config(
            None,
            _base(pipelines={"cis": {"strategies": [{"name": "a", "transforms": transforms}]}}),
        )


def test_cis_and_trans_rint_flags() -> None:
    cfg = load_config(None, _base(pipelines={"cis": {"rint": True}}))
    assert cfg.pipelines.cis.rint is True
    with pytest.raises(ConfigError, match="rint"):
        load_config(None, _base(pipelines={"cis": {"rint": "maybe"}}))
    # trans: one transform per run too (pipelines.trans.rint, default as given)
    cfg = load_config(None, _base(pipelines={"trans": {"strategies": [{"name": "a"}]}}))
    assert cfg.pipelines.trans.rint is False
    assert not hasattr(cfg.pipelines.trans.strategies[0], "transforms")
    assert (
        load_config(None, _base(pipelines={"trans": {"rint": True}})).pipelines.trans.rint is True
    )


@pytest.mark.parametrize("transforms", [["raw", "rint"], ["rint"]])
def test_trans_transforms_list_is_rejected_with_guidance(transforms: list[str]) -> None:
    with pytest.raises(
        ConfigError, match=r"no longer accepted.*pipelines\.trans\.rint: true.*two runs"
    ):
        load_config(
            None,
            _base(pipelines={"trans": {"strategies": [{"name": "a", "transforms": transforms}]}}),
        )


def test_trans_one_strategy_and_new_keys() -> None:
    with pytest.raises(ConfigError, match="a trans run evaluates one set of inputs, got 2"):
        load_config(
            None, _base(pipelines={"trans": {"strategies": [{"name": "a"}, {"name": "b"}]}})
        )
    t = load_config(None, _base()).pipelines.trans
    assert (t.hotspot_bin_bp, t.hotspot_min_traits, t.cis_window_bp) == (1_000_000, 5, 1_000_000)
    assert (t.levels.nominal, t.levels.strict, t.levels.bonferroni) == (1e-5, None, True)
    t = load_config(None, _base(pipelines={"trans": {"hotspot_bin_bp": 250_000}})).pipelines.trans
    assert t.hotspot_bin_bp == 250_000
    with pytest.raises(ConfigError, match="hotspot_bin_bp"):
        load_config(None, _base(pipelines={"trans": {"hotspot_bin_bp": 0}}))
    with pytest.raises(ConfigError, match=r"levels\.fdr is no longer used"):
        load_config(None, _base(pipelines={"trans": {"levels": {"fdr": 0.05}}}))


# ------------------------------------------------------------------------ trans
def test_trans_variants_traits_literal_vs_path(tmp_path: Path, chdir: Path) -> None:
    cfg = load_config(None, _base())
    t = cfg.pipelines.trans
    assert t.variants == "reference" and t.traits == "reference" and t.known_pairs is True

    cfg = load_config(None, _base(pipelines={"trans": {"variants": "vars.txt", "traits": "all"}}))
    assert cfg.pipelines.trans.variants == chdir.resolve() / "vars.txt"
    assert cfg.pipelines.trans.traits == "all"

    cfg = load_config(None, _base(pipelines={"trans": {"traits": "file:traits.txt"}}))
    assert cfg.pipelines.trans.traits == chdir.resolve() / "traits.txt"

    yaml_path = _write(
        tmp_path / "p" / "c.yaml", _base(pipelines={"trans": {"traits": "t/traits.txt"}})
    )
    cfg = load_config(yaml_path)
    assert cfg.pipelines.trans.traits == (tmp_path / "p" / "t/traits.txt").resolve()


def test_trans_reference_prefix_forms(tmp_path: Path) -> None:
    data = _base()
    data["references"]["other"] = {"pairs": "o.tsv"}
    data["pipelines"] = {
        "match": {"reference": "ref"},
        "cis": {"reference": "ref"},
        "trans": {
            "variants": "reference:other",
            "traits": "reference:other",
            "known_pairs": "reference:other",
        },
    }
    cfg = load_config(_write(tmp_path / "c.yaml", data))
    t = cfg.pipelines.trans
    assert t.reference == "other"
    assert t.variants == "reference" and t.traits == "reference" and t.known_pairs is True
    assert cfg.reference_for("trans") == "other"

    data["pipelines"]["trans"] = {"reference": "ref", "variants": "reference:other"}
    with pytest.raises(ConfigError, match="different references"):
        load_config(_write(tmp_path / "d.yaml", data))

    data["pipelines"]["trans"] = {"reference": "other", "known_pairs": None}
    assert load_config(_write(tmp_path / "e.yaml", data)).pipelines.trans.known_pairs is False


def test_trans_traits_reference_needs_pairs(tmp_path: Path) -> None:
    data = _base(references={"tr": {"kind": "trans", "variants": "v.txt"}})
    data["pipelines"] = {"match": {"enabled": False}, "cis": {"enabled": False}}
    with pytest.raises(ConfigError, match="traits: 'reference' needs pairs"):
        load_config(None, data)
    data["pipelines"]["trans"] = {"traits": "all"}
    cfg = load_config(None, data)
    assert cfg.references["tr"].variants is not None


def test_trait_id_rule_file_shorthand(tmp_path: Path, chdir: Path) -> None:
    data = _base()
    data["datasets"]["rna"]["trait_id_rule"] = "file:maps/traits.tsv"
    cfg = load_config(_write(tmp_path / "c.yaml", data))
    d = cfg.datasets["rna"]
    assert d.trait_id_rule == "as_is"
    assert d.trait_map == (tmp_path / "maps/traits.tsv").resolve()
    with pytest.raises(ConfigError, match="trait_id_rule"):
        load_config(None, _base(datasets={"rna": {"pheno": "p", "trait_id_rule": "lowercase"}}))


# ----------------------------------------------------------- resolved round trip
def test_resolved_config_round_trip(tmp_path: Path) -> None:
    strategies = [{"name": "a", "cov": {"include": ["PC1-2"]}}, {"name": "b", "cov": None}]
    cfg = load_config(
        _write(
            tmp_path / "c.yaml",
            _base(pipelines={"cis": {"strategies": strategies}, "match": {"lambda": 0.4}}),
        )
    )
    run = tmp_path / "run"
    run.mkdir()
    write_resolved_config(cfg, run / "anqchor.resolved.yaml")
    mtime = (run / "anqchor.resolved.yaml").stat().st_mtime_ns
    write_resolved_config(cfg, run / "anqchor.resolved.yaml")  # unchanged → not rewritten
    assert (run / "anqchor.resolved.yaml").stat().st_mtime_ns == mtime
    back = load_resolved_config(run)
    assert back == cfg
    text = (run / "anqchor.resolved.yaml").read_text()
    assert str((tmp_path / "geno/list.tsv").resolve()) in text


def test_helpers(tmp_path: Path) -> None:
    data = _base()
    data["datasets"]["prot"] = {"pheno": "q.tsv"}
    data["pipelines"] = {"cis": {"datasets": ["prot"]}}
    cfg = load_config(None, data)
    assert cfg.datasets_for("cis") == ["prot"]
    assert cfg.datasets_for("match") == ["rna", "prot"]
    assert cfg.reference_for("match") == "ref"


def test_negative_z_diff_rejected(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="z_diff"):
        load_config(_write(tmp_path / "c.yaml", _base(pipelines={"match": {"z_diff": -1}})))
