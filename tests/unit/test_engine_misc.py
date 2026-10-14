"""Python step dispatch (``anqchor _step``), provenance and file helpers."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

import anqchor.steps as steps_mod
from anqchor.engine.minimake import write_pipeline_makefile
from anqchor.engine.provenance import write_run_json
from anqchor.engine.step import Step, StepGraph
from anqchor.steps import REGISTRY, StepContext, StepError, register, run_step
from anqchor.util.files import atomic_write, open_text, read_header, resolve_near, write_json


@pytest.fixture
def no_step_modules(monkeypatch: pytest.MonkeyPatch) -> None:
    """Do not import the real step modules (other packages may be mid-edit)."""
    monkeypatch.setattr(steps_mod, "STEP_MODULES", ())


def test_run_step_dispatch(tmp_path: Path, no_step_modules: None) -> None:
    seen: dict[str, object] = {}

    @register("unit_test_kind")
    def impl(ctx: StepContext) -> None:
        seen["cwd"] = Path.cwd()
        seen["params"] = ctx.params
        seen["inputs"] = ctx.inputs
        with atomic_write(ctx.outputs[0]) as fh:
            fh.write("done\n")
        seen["abs"] = ctx.path("x/y.txt")

    try:
        run = tmp_path / "run"
        step = Step(
            name="t.py",
            kind="unit_test_kind",
            inputs=(Path("in.txt"),),
            outputs=(Path("out/o.txt"),),
            params={"alpha": 0.05, "names": ["a"]},
        )
        write_pipeline_makefile(StepGraph([step]), run, {})
        cwd = Path.cwd()
        run_step("t.py", run)
        assert Path.cwd() == cwd  # restored
        assert seen["cwd"] == run.resolve()
        assert seen["params"] == {"alpha": 0.05, "names": ["a"]}
        assert seen["inputs"] == [Path("in.txt")]
        assert seen["abs"] == run.resolve() / "x/y.txt"
        assert (run / "out/o.txt").read_text() == "done\n"
    finally:
        REGISTRY.pop("unit_test_kind", None)


def test_run_step_errors(tmp_path: Path, no_step_modules: None) -> None:
    run = tmp_path / "run"
    g = StepGraph(
        [
            Step(name="shell", commands=("true",), outputs=(Path("a"),)),
            Step(name="nokind", kind="no_such_kind", outputs=(Path("b"),)),
        ]
    )
    write_pipeline_makefile(g, run, {})
    with pytest.raises(StepError, match="regenerate the Makefile"):
        run_step("missing", run)
    with pytest.raises(StepError, match="not a Python step"):
        run_step("shell", run)
    with pytest.raises(StepError, match="no implementation registered"):
        run_step("nokind", run)


def test_register_twice() -> None:
    def a(ctx: StepContext) -> None: ...

    def b(ctx: StepContext) -> None: ...

    try:
        register("dup_kind")(a)
        register("dup_kind")(a)  # same function: fine
        with pytest.raises(ValueError, match="registered twice"):
            register("dup_kind")(b)
    finally:
        REGISTRY.pop("dup_kind", None)


def test_write_run_json(tmp_path: Path) -> None:
    binary = tmp_path / "qpgentools"
    binary.write_text("#!/bin/sh\n")
    g = StepGraph([Step(name="a", commands=("true",), outputs=(Path("a"),))])
    (tmp_path / "steps").mkdir()
    (tmp_path / "steps" / "a.seconds").write_text("7\n")
    for i in range(3):
        write_run_json(
            tmp_path,
            config_text="x: 1\n",
            qpgen_binary=binary,
            argv=["anqchor", str(i)],
            started="s",
            finished="f",
            returncode=0,
            graph=g,
        )
    data = json.loads((tmp_path / "run.json").read_text())
    assert data["schema_version"] == 1
    assert data["qpgen"]["binary"] == str(binary) and len(data["qpgen"]["sha256"]) == 64
    assert [inv["argv"][1] for inv in data["invocations"]] == ["0", "1", "2"]
    assert data["step_seconds"] == {"a": 7.0}
    assert len(data["config_sha256"]) == 64


def test_file_helpers(tmp_path: Path) -> None:
    gz = tmp_path / "t.tsv.gz"
    with atomic_write(gz) as fh:
        fh.write("a\tb\n1\t2\n")
    assert read_header(gz) == ["a", "b"]
    with open_text(gz) as fh:
        assert fh.read().splitlines()[1] == "1\t2"
    assert not [p for p in tmp_path.iterdir() if p.name.endswith(".tmp")]
    with pytest.raises(RuntimeError), atomic_write(tmp_path / "x.txt") as fh:
        fh.write("partial")
        raise RuntimeError("boom")
    assert not (tmp_path / "x.txt").exists()
    assert not [p for p in tmp_path.iterdir() if p.name.endswith(".tmp")]
    write_json(tmp_path / "j.json", {"nan": float("nan"), "p": Path("/a"), "s": {1}})
    assert json.loads((tmp_path / "j.json").read_text()) == {"nan": None, "p": "/a", "s": [1]}
    c = tmp_path / "hdr.tsv"
    c.write_text("## comment\nx\ty\r\n")
    assert read_header(c, comment="##") == ["x", "y"]


def test_resolve_near(tmp_path: Path) -> None:
    proj = tmp_path / "proj"
    (proj / "geno" / "sub").mkdir(parents=True)
    (proj / "geno" / "sub" / "chr1.pgen").write_text("")
    anchor = proj / "geno" / "list.tsv"
    anchor.write_text("")
    # relative to an ancestor of the list (as in the MAGE example)
    assert resolve_near("geno/sub/chr1", anchor) == (proj / "geno/sub/chr1").resolve()
    assert resolve_near("sub/chr1", anchor) == (proj / "geno/sub/chr1").resolve()
    assert resolve_near("/abs/x", anchor) == Path("/abs/x")
    assert resolve_near("nothing", anchor) == (proj / "geno/nothing").resolve()
    assert os.path.isabs(resolve_near("nothing", anchor))
