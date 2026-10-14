"""Step and StepGraph: validation, de-duplication, producers, closure, cycles."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from anqchor.engine.step import Step, StepGraph, python_step_command


def sh(
    name: str, inputs: tuple[str, ...], outputs: tuple[str, ...], cmd: str = "true", **kw: object
) -> Step:
    return Step(
        name=name,
        inputs=tuple(Path(i) for i in inputs),
        outputs=tuple(Path(o) for o in outputs),
        commands=(cmd,),
        **kw,
    )  # type: ignore[arg-type]


def test_step_validation() -> None:
    with pytest.raises(ValueError, match="invalid step name"):
        sh("bad name", (), ("a",))
    with pytest.raises(ValueError, match="invalid step name"):
        sh("", (), ("a",))
    with pytest.raises(ValueError, match="no outputs"):
        sh("a", (), ())
    with pytest.raises(ValueError, match="neither commands nor a Python kind"):
        Step(name="a", outputs=(Path("x"),))
    with pytest.raises(ValueError, match="single-line"):
        sh("a", (), ("x",), cmd="echo 1\necho 2")


def test_step_normalises_fields() -> None:
    s = Step(name="a", inputs=["in.txt"], outputs="out.txt", commands="echo hi")  # type: ignore[arg-type]
    assert s.inputs == (Path("in.txt"),)
    assert s.outputs == (Path("out.txt"),)
    assert s.commands == ("echo hi",)
    hash(s)  # frozen and hashable


def test_python_step_command_and_signature() -> None:
    s = Step(
        name="prepare.genotypes", kind="prepare_genotypes", outputs=(Path("o"),), params={"x": 1}
    )
    cmd = python_step_command("prepare.genotypes")
    assert cmd.endswith("-m anqchor _step prepare.genotypes --run-dir .")
    assert sys.executable in cmd
    assert s.all_commands() == (cmd,)
    sig = s.signature()
    assert sig["params"] == {"x": 1} and sig["kind"] == "prepare_genotypes"
    assert sig["outputs"] == ["o"]


def test_graph_dedup_and_conflicts() -> None:
    g = StepGraph()
    a = sh("a", (), ("a.txt",))
    g.add(a)
    g.add(sh("a", (), ("a.txt",)))  # identical → no-op
    assert len(g) == 1
    with pytest.raises(ValueError, match="conflicting definitions for step a"):
        g.add(sh("a", (), ("a.txt",), cmd="false"))
    with pytest.raises(ValueError, match="conflicting definitions"):
        g.add(Step(name="a", outputs=(Path("a.txt"),), commands=("true",), params={"k": 1}))
    with pytest.raises(ValueError, match="produced by both a and b"):
        g.add(sh("b", (), ("a.txt",)))
    assert "a" in g and "b" not in g
    with pytest.raises(KeyError, match="unknown step 'zz'"):
        g.get("zz")


def chain() -> StepGraph:
    #   a → b → c
    #   a ------→ d (d also reads an external file)
    return StepGraph(
        [
            sh("a", ("ext/in.txt",), ("a.txt",)),
            sh("b", ("a.txt",), ("b.txt", "b2.txt")),
            sh("c", ("b2.txt",), ("c.txt",)),
            sh("d", ("a.txt", "/abs/ext.txt"), ("d.txt",), pipeline="other"),
        ]
    )


def test_producer_lookup_and_upstream() -> None:
    g = chain()
    assert g.producer_of(Path("a.txt")) == "a"
    assert g.producer_of(Path("b2.txt")) == "b"
    assert g.producer_of(Path("ext/in.txt")) is None
    assert g.upstream("c") == ["b"]
    assert g.upstream("d") == ["a"]
    assert g.upstream("a") == []


def test_closure_and_downstream() -> None:
    g = chain()
    assert g.closure(["c"]) == ["a", "b", "c"]
    assert g.closure(["d", "c"]) == ["a", "d", "b", "c"]
    assert sorted(g.downstream("a")) == ["b", "c", "d"]
    assert g.downstream("c") == []
    assert [s.name for s in g.by_pipeline("other")] == ["d"]
    assert g.names() == ["a", "b", "c", "d"]


def test_cycle_detection() -> None:
    g = StepGraph([sh("x", ("y.txt",), ("x.txt",)), sh("y", ("x.txt",), ("y.txt",))])
    with pytest.raises(ValueError, match="dependency cycle"):
        g.closure(["x"])


def test_self_input_is_not_a_dependency() -> None:
    g = StepGraph([sh("x", ("x.txt",), ("x.txt",))])
    assert g.upstream("x") == []
    assert g.closure(["x"]) == ["x"]
