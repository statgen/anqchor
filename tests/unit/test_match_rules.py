"""Match status rules: classify(), its JavaScript port, and agreement with qpgen."""

from __future__ import annotations

import json
import random
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from anqchor.qpgen import commands as qc
from anqchor.qpgen.outputs import MATCH_ORDINALS, read_match_all
from anqchor.report.match import MATCH_RULES_JS
from anqchor.steps.match import (
    CLASSIFY_TEST_VECTORS,
    classify,
    close_candidates,
    group_of,
    is_duplicate,
    is_related,
    is_rescue,
)

ZT, ZD = 1.96, 2.0
needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")


def _kw(v: dict[str, Any]) -> dict[str, Any]:
    return {k: v[k] for k in ("mapped", "n_traits", "rank_self", "z_self", "z_top")}


@pytest.mark.parametrize("vec", CLASSIFY_TEST_VECTORS, ids=lambda v: v["expected"])
def test_classify_vectors(vec: dict[str, Any]) -> None:
    assert classify(**_kw(vec), z_threshold=ZT, z_diff=ZD) == vec["expected"]


def test_vectors_cover_every_status() -> None:
    assert {v["expected"] for v in CLASSIFY_TEST_VECTORS} == {
        "SELF_BEST",
        "SELF_LENIENT",
        "UNCLEAR",
        "SINGLE_NEW_BEST",
        "MULTI_NEW_BEST",
        "NO_SELF",
        "NO_OBS_TRAITS",
    }


def test_classify_thresholds_move_calls() -> None:
    kw = {
        "mapped": True,
        "n_traits": 10,
        "rank_self": 2,
        "z_self": 3.0,
        "z_top": [3.5, 3.0, 2.0, 1.9, 1.8],
    }
    assert classify(**kw, z_threshold=1.96, z_diff=2.0) == "SELF_LENIENT"
    assert classify(**kw, z_threshold=3.0, z_diff=2.0) == "UNCLEAR"
    assert classify(**kw, z_threshold=1.96, z_diff=0.4) == "SINGLE_NEW_BEST"
    assert classify(**kw, z_threshold=1.96, z_diff=0.6) == "MULTI_NEW_BEST"  # 3.0 > 2.0 + 0.6


def test_classify_nan_handling() -> None:
    nan = float("nan")
    assert (
        classify(
            mapped=True,
            n_traits=nan,
            rank_self=1,
            z_self=1.0,
            z_top=[1.0] * 5,
            z_threshold=ZT,
            z_diff=ZD,
        )
        == "NO_OBS_TRAITS"
    )
    assert (
        classify(
            mapped=False,
            n_traits=5,
            rank_self=nan,
            z_self=nan,
            z_top=[3.0, nan],
            z_threshold=ZT,
            z_diff=ZD,
        )
        == "SINGLE_NEW_BEST"
    )
    # rank NaN for a mapped sample: not SELF_BEST
    assert (
        classify(
            mapped=True,
            n_traits=5,
            rank_self=nan,
            z_self=0.0,
            z_top=[1.0, 0.9, 0.8, 0.7, 0.6],
            z_threshold=ZT,
            z_diff=ZD,
        )
        == "UNCLEAR"
    )


def test_groups_and_rescue() -> None:
    assert group_of("SELF_BEST", True) == "PASS"
    assert group_of("SELF_LENIENT", True) == group_of("UNCLEAR", True) == "WARN"
    assert group_of("SINGLE_NEW_BEST", True) == group_of("MULTI_NEW_BEST", True) == "FAIL"
    assert group_of("NO_OBS_TRAITS", True) == "NO_TEST"
    for s in ("SINGLE_NEW_BEST", "MULTI_NEW_BEST", "NO_SELF", "NO_OBS_TRAITS"):
        assert group_of(s, False) == "NO_TEST"
    assert group_of("NOT_SCORED", True) == "NO_TEST"
    assert is_rescue("SINGLE_NEW_BEST", False) and not is_rescue("SINGLE_NEW_BEST", True)
    assert not is_rescue("MULTI_NEW_BEST", False)


def test_close_candidates_and_related() -> None:
    ids = ["A", "B", "C", "D", "E"]
    # one relative: a block of two that drops to the background
    close = close_candidates(ids, [6.0, 5.0, 3.0, 2.9, 2.8], ZT)
    assert close == [("A", 6.0), ("B", 5.0)] and is_related(len(close))
    # a relative far below the best (Z about half the self Z) still counts
    assert close_candidates(ids, [9.2, 5.1, 3.4, 3.3, 3.2], ZT) == [("A", 9.2), ("B", 5.1)]
    # trio child: both parents; the block ends at the largest drop
    assert [i for i, _ in close_candidates(ids, [9.0, 6.0, 4.8, 3.4, 3.3], ZT)] == ["A", "B", "C"]
    assert [i for i, _ in close_candidates(ids, [9.0, 5.1, 5.0, 3.4, 3.3], ZT)] == ["A", "B", "C"]
    # a smaller later drop does not extend the block
    assert [i for i, _ in close_candidates(ids, [9.0, 5.0, 3.4, 3.3, 2.2], ZT)] == ["A", "B"]
    # flat background: no block
    assert close_candidates(ids, [3.0, 2.9, 2.8, 2.7, 2.6], ZT) == []
    assert close_candidates(ids, [7.0, 3.3, 3.1, 3.0, 2.9], ZT) == []
    # drop exactly at the gap counts; just below does not
    assert len(close_candidates(ids, [6.0, 4.0, 3.0, 2.9, 2.8], ZT, gap=1.0)) == 2
    assert close_candidates(ids, [6.0, 3.99, 3.0, 2.95, 2.9], ZT, gap=1.0) == []
    # block members must reach z_threshold
    assert close_candidates(ids, [3.0, 1.9, 0.5, 0.4, 0.3], ZT) == []
    assert len(close_candidates(ids, [3.0, 1.9, 0.5, 0.4, 0.3], 1.5)) == 2
    # candidate k+1 must exist: two genotyped samples cannot show a block
    assert close_candidates(["A", "B", None, None, None], [6.0, 5.0, None, None, None], ZT) == []
    assert close_candidates(["A", "B", "C", None, None], [6.0, 5.0, 1.0, None, None], ZT) == [
        ("A", 6.0),
        ("B", 5.0),
    ]
    assert close_candidates([None] * 5, [None] * 5, ZT) == []
    assert not is_related(1) and not is_related(0)


# ------------------------------------------------------------- against qpgen output
def _check_against_qpgen(df: pd.DataFrame, zt: float, zd: float) -> int:
    n = 0
    for r in df.to_dict("records"):
        mapped = not pd.isna(r["ID.self"])
        got = classify(
            mapped=mapped,
            n_traits=r["N.Traits"],
            rank_self=r["Rank.self"],
            z_self=r["Z.self"],
            z_top=[r[f"Z.{o}"] for o in MATCH_ORDINALS],
            z_threshold=zt,
            z_diff=zd,
        )
        assert got == r["MatchStatus"], (r["ID.Pheno"], got, r["MatchStatus"])
        n += 1
    return n


@pytest.mark.parametrize("name", ["match.match.all.tsv.gz", "matchauto.match.all.tsv.gz"])
def test_classify_reproduces_qpgen_fixtures(qpgen_fixtures: Path, name: str) -> None:
    df = read_match_all(qpgen_fixtures / name)
    assert _check_against_qpgen(df, ZT, ZD) == len(df) > 0


@pytest.mark.integration
def test_classify_reproduces_qpgen_on_mage(
    qpgen_bin: Path, mage_data: Path, mage_pgen_list: Path, tmp_path: Path
) -> None:
    """Every row of real MAGE match-prs-pheno output, at the default and other thresholds."""
    prs = tmp_path / "prs"
    subprocess.run(
        qc.pair_prs(
            qpgen_bin,
            genotypes=qc.Genotypes(pgen_list=mage_pgen_list),
            pairs=mage_data / "reference" / "gtex_lcl_eqtl.chr21_22.tsv.gz",
            out_prefix=prs,
        ),
        check=True,
        capture_output=True,
    )
    pheno = mage_data / "pheno" / "mage_expr.chr21_22.regenie.tsv.gz"
    statuses: dict[tuple[float, float], dict[str, int]] = {}
    for zt, zd in ((1.96, 2.0), (3.0, 1.0), (1.0, 0.5), (5.0, 3.0)):
        out = tmp_path / f"m_{zt}_{zd}"
        subprocess.run(
            qc.match_prs_pheno(
                qpgen_bin,
                prs=str(prs) + ".prs.tsv.gz",
                pheno=pheno,
                out_prefix=out,
                rint=True,
                z_threshold=zt,
                z_diff=zd,
                threads=4,
            ),
            check=True,
            capture_output=True,
        )
        df = read_match_all(Path(str(out) + ".match.all.tsv.gz"))
        assert _check_against_qpgen(df, zt, zd) == 729
        statuses[(zt, zd)] = df["MatchStatus"].value_counts().to_dict()
    assert statuses[(1.96, 2.0)]["SELF_BEST"] > 650  # 678 without covariates (GTEx reference)
    assert (
        len({json.dumps(v, sort_keys=True) for v in statuses.values()}) > 1
    )  # thresholds mattered


# ------------------------------------------------------------------ JavaScript port
def _run_js(tmp_path: Path, payload: dict[str, Any], body: str) -> Any:
    (tmp_path / "rules.js").write_text(MATCH_RULES_JS, encoding="utf-8")
    (tmp_path / "in.json").write_text(json.dumps(payload), encoding="utf-8")
    script = (
        "const M = require('./rules.js'); const fs = require('fs');\n"
        "const P = JSON.parse(fs.readFileSync('./in.json', 'utf8'));\n"
        + body
        + "\nprocess.stdout.write(JSON.stringify(out));\n"
    )
    (tmp_path / "run.js").write_text(script, encoding="utf-8")
    proc = subprocess.run(
        ["node", "run.js"], cwd=tmp_path, capture_output=True, text=True, check=False
    )
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


def _random_vectors(n: int, seed: int = 7) -> list[dict[str, Any]]:
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        mapped = rng.random() < 0.8
        z = sorted((round(rng.gauss(3, 1.5), 2) for _ in range(5)), reverse=True)
        if rng.random() < 0.1:  # fewer than five genotyped samples
            k = rng.randrange(2, 5)
            z[k:] = [None] * (5 - k)
        if rng.random() < 0.05:
            z = [z[0], None, None, None, None]
        rank = rng.choice([1, 1, 1, 2, 3, 50]) if mapped else None
        zself = (z[0] if rank == 1 else round(rng.gauss(2, 1.5), 2)) if mapped else None
        out.append(
            {
                "mapped": mapped,
                "n_traits": rng.choice([0, 5, 100, 100, 100]),
                "rank_self": rank,
                "z_self": zself,
                "z_top": z,
                "zt": rng.choice([1.0, 1.96, 2.5, 3.0]),
                "zd": rng.choice([0.5, 1.0, 2.0, 3.0]),
            }
        )
    return out


@needs_node
def test_js_classify_matches_python(tmp_path: Path) -> None:
    vecs = [dict(v, zt=ZT, zd=ZD) for v in CLASSIFY_TEST_VECTORS] + _random_vectors(400)
    out = _run_js(
        tmp_path,
        {"v": vecs},
        (
            "const out = P.v.map(v => [M.classify(v.mapped, v.n_traits, v.rank_self, v.z_self, v.z_top, v.zt, v.zd),"
            " M.groupOf(M.classify(v.mapped, v.n_traits, v.rank_self, v.z_self, v.z_top, v.zt, v.zd), v.mapped)]);"
        ),
    )
    for v, (js_status, js_group) in zip(vecs, out, strict=True):
        py = classify(**_kw(v), z_threshold=v["zt"], z_diff=v["zd"])
        assert js_status == py, v
        assert js_group == group_of(py, v["mapped"]), v


@needs_node
def test_js_close_candidates_and_verdict(tmp_path: Path) -> None:
    from anqchor.steps.match import fail_phrase

    vecs = _random_vectors(200, seed=3)
    rng = random.Random(5)
    for v in vecs:
        v["gap"] = rng.choice([0.3, 0.5, 1.0])
    ids = ["A", "B", "C", "D", "E"]
    out = _run_js(
        tmp_path,
        {"v": vecs, "ids": ids},
        (
            "const out = P.v.map(v => { const c = M.closeCandidates(P.ids, v.z_top, v.zt, v.gap);"
            " return [c, M.isRelated(c.length)]; });"
        ),
    )
    for v, (js_close, js_rel) in zip(vecs, out, strict=True):
        py = close_candidates(ids, v["z_top"], v["zt"], v["gap"])
        assert [tuple(x) for x in js_close] == [(i, pytest.approx(z)) for i, z in py]
        assert js_rel == is_related(len(py))
    assert sum(1 for _, rel in out if rel) > 10  # the vectors exercise both outcomes
    rows = [{"pheno_id": "S1", "geno_id_assigned": "G2", "id_1": "G1"}]
    evals = [{"status": "SINGLE_NEW_BEST", "group": "FAIL", "mapped": True}]
    v = _run_js(tmp_path, {"rows": rows, "evals": evals}, "const out = M.verdict(P.rows, P.evals);")
    assert v["status"] == "critical"
    assert fail_phrase("S1", "G2", "G1") in v["text"]


# ------------------------------------------------------------- two genotype IDs
def test_is_duplicate() -> None:
    # MAGE: a genotype duplicated (exact copy, noisy copy): equal scores far above the rest
    assert is_duplicate([6.99, 6.99, 2.56, 2.5, 2.4], ZT, ZD)
    assert is_duplicate([7.00, 6.83, 2.56, 2.5, 2.4], ZT, ZD)
    # parent/child: the relative scores about half the self score
    assert not is_duplicate([9.2, 5.1, 3.4, 3.3, 3.2], ZT, ZD)
    assert not is_duplicate([15.0, 8.0, 3.0, 2.9, 2.8], ZT, ZD)
    # a weak duplicate does not stand clear of the background: not flagged (honest limit)
    assert not is_duplicate([3.84, 3.84, 3.17, 3.0, 2.9], ZT, ZD)
    # each condition at its boundary: clear lead (strict), threshold (strict), ratio (>=)
    assert not is_duplicate([6.0, 5.0, 3.0, 2.0, 1.0], ZT, ZD)
    assert is_duplicate([6.0, 5.01, 3.0, 2.0, 1.0], ZT, ZD)
    assert not is_duplicate([2.0, 1.96, -0.1, -0.2, -0.3], ZT, ZD)
    assert is_duplicate([2.0, 1.97, -0.1, -0.2, -0.3], ZT, ZD)
    assert is_duplicate([10.0, 8.0, 3.0, 2.0, 1.0], ZT, ZD) and not is_duplicate(
        [10.0, 7.99, 3.0, 2.0, 1.0], ZT, ZD
    )
    # fewer than three genotypes: the missing third counts as -9999
    assert is_duplicate([5.0, 4.9, None, None, None], ZT, ZD)
    assert not is_duplicate([5.0, None, None, None, None], ZT, ZD)
    assert not is_duplicate([None] * 5, ZT, ZD)
    # a duplicate is always MULTI_NEW_BEST or SINGLE_NEW_BEST unless the own ID ranks first
    for z in ([6.99, 6.99, 2.56, 2.5, 2.4], [12.0, 9.8, 3.0, 2.9, 2.8]):
        st = classify(
            mapped=False,
            n_traits=9,
            rank_self=None,
            z_self=None,
            z_top=z,
            z_threshold=ZT,
            z_diff=ZD,
        )
        assert st in ("MULTI_NEW_BEST", "SINGLE_NEW_BEST")


@needs_node
def test_js_is_duplicate_matches_python(tmp_path: Path) -> None:
    vecs = _random_vectors(400, seed=11)
    rng = random.Random(13)
    for v in vecs:  # add near-ties at strong scores so that both outcomes occur
        if rng.random() < 0.3 and v["z_top"][0] is not None:
            top = round(rng.uniform(4, 10), 2)
            v["z_top"] = [
                top,
                round(top * rng.uniform(0.7, 1.0), 2),
                *[round(rng.gauss(2.5, 0.5), 2)] * 3,
            ]
    out = _run_js(
        tmp_path, {"v": vecs}, "const out = P.v.map(v => M.isDuplicate(v.z_top, v.zt, v.zd));"
    )
    py = [is_duplicate(v["z_top"], v["zt"], v["zd"]) for v in vecs]
    assert out == py
    assert 10 < sum(py) < len(py) - 10
