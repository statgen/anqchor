"""Exact argv of the typed qpgentools command builders, and binary discovery."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from anqchor.qpgen import commands as qc
from anqchor.qpgen import locate as loc

Q = "/opt/qpgen/bin/qpgentools"
GL = qc.Genotypes(pgen_list=Path("prepare/genotypes/pgen_list.tsv"))


def test_genotypes_args() -> None:
    assert GL.args() == ["--pgen-list", "prepare/genotypes/pgen_list.tsv"]
    pf = qc.Genotypes(pfile=Path("/g/chr22"))
    assert pf.args() == [
        "--pgen",
        "/g/chr22.pgen",
        "--psam",
        "/g/chr22.psam",
        "--pivar",
        "/g/chr22.pvar.idx.gz",
    ]
    with pytest.raises(ValueError):
        qc.Genotypes().args()


def test_pair_prs() -> None:
    assert qc.pair_prs(Q, genotypes=GL, pairs="r.tsv.gz", out_prefix="prs/ref/prs") == [
        Q,
        "pair-prs",
        "--pgen-list",
        "prepare/genotypes/pgen_list.tsv",
        "--pairs",
        "r.tsv.gz",
        "--out",
        "prs/ref/prs",
    ]
    argv = qc.pair_prs(
        Q,
        genotypes=GL,
        pairs="r",
        out_prefix="o",
        sample=Path("ids.txt"),
        colname_trait="gene",
        colname_variant="snp",
        colname_beta="b",
        colname_se="s",
        max_chunk_vars=500,
    )
    assert argv[8:] == [
        "--sample",
        "ids.txt",
        "--colname-pair-trait",
        "gene",
        "--colname-pair-variant",
        "snp",
        "--colname-pair-beta",
        "b",
        "--colname-pair-se",
        "s",
        "--max-chunk-vars",
        "500",
    ]


def test_match_prs_pheno_defaults() -> None:
    assert qc.match_prs_pheno(
        Q, prs="p.prs.tsv.gz", pheno="ph.tsv.gz", out_prefix="match/rna/qpgen"
    ) == [
        Q,
        "match-prs-pheno",
        "--prs",
        "p.prs.tsv.gz",
        "--pheno",
        "ph.tsv.gz",
        "--out",
        "match/rna/qpgen",
        "--rint",
        "--z-threshold",
        "1.96",
        "--z-diff",
        "2",
        "--threads",
        "1",
    ]


def test_match_prs_pheno_full() -> None:
    argv = qc.match_prs_pheno(
        Q,
        prs="p",
        pheno="ph",
        out_prefix="o",
        cov="c",
        sample_tsv="s.tsv",
        trait_tsv="t.tsv",
        pheno_format="tsv-sample-col",
        cov_format="tsv-sample-row",
        missing_str=("NA", "."),
        rint=False,
        mahalanobis=True,
        auto_lambda=True,
        auto_lambda_min_self=50,
        auto_lambda_min_best=0,
        min_weight=0.01,
        z_threshold=2.5,
        z_diff=1.5,
        threads=8,
        cov_impute_mean=True,
    )
    assert argv == [
        Q,
        "match-prs-pheno",
        "--prs",
        "p",
        "--pheno",
        "ph",
        "--out",
        "o",
        "--cov",
        "c",
        "--sample-tsv",
        "s.tsv",
        "--trait-tsv",
        "t.tsv",
        "--pheno-format",
        "tsv-sample-col",
        "--cov-format",
        "tsv-sample-row",
        "--missing-str",
        "NA,.",
        "--mahalanobis",
        "--auto-lambda",
        "--auto-lambda-min-self",
        "50",
        "--auto-lambda-min-best",
        "0",
        "--cov-impute-mean",
        "--min-weight",
        "0.01",
        "--z-threshold",
        "2.5",
        "--z-diff",
        "1.5",
        "--threads",
        "8",
    ]
    fixed = qc.match_prs_pheno(Q, prs="p", pheno="ph", out_prefix="o", mahalanobis=True, lam=0.5)
    assert fixed[fixed.index("--mahalanobis") + 1 : fixed.index("--mahalanobis") + 3] == [
        "--lambda",
        "0.5",
    ]
    # cov_format is only passed with a covariate file
    assert "--cov-format" not in qc.match_prs_pheno(
        Q, prs="p", pheno="ph", out_prefix="o", cov_format="tsv-sample-col"
    )


@pytest.mark.parametrize(
    ("kw", "msg"),
    [
        ({"auto_lambda": True}, "auto_lambda requires mahalanobis"),
        ({"lam": 0.5}, "lam requires mahalanobis"),
        ({"mahalanobis": True, "lam": 0.5, "auto_lambda": True}, "either lam or auto_lambda"),
        ({"mahalanobis": True, "lam": 1.5}, "between 0 and 1"),
        ({"mahalanobis": True, "auto_lambda_min_self": 10}, "require auto_lambda"),
    ],
)
def test_match_prs_pheno_invalid(kw: dict[str, object], msg: str) -> None:
    with pytest.raises(ValueError, match=msg):
        qc.match_prs_pheno(Q, prs="p", pheno="ph", out_prefix="o", **kw)  # type: ignore[arg-type]


def test_pair_assoc() -> None:
    assert qc.pair_assoc(Q, genotypes=GL, pheno="ph", pairs="pairs.tsv.gz", out="cis/x.tsv.gz") == [
        Q,
        "pair-assoc",
        "--pgen-list",
        "prepare/genotypes/pgen_list.tsv",
        "--pheno",
        "ph",
        "--pairs",
        "pairs.tsv.gz",
        "--out",
        "cis/x.tsv.gz",
    ]
    argv = qc.pair_assoc(
        Q,
        genotypes=GL,
        pheno="ph",
        pairs="p",
        out="o",
        cov="c",
        cov_format="tsv-sample-col",
        pheno_format="tensorqtl",
        sample="s",
        skip_rint=True,
        colname_trait="g",
        colname_variant="v",
        max_chunk_vars=10,
    )
    assert argv[10:] == [
        "--cov",
        "c",
        "--pheno-format",
        "tensorqtl",
        "--cov-format",
        "tsv-sample-col",
        "--sample",
        "s",
        "--skip-rint",
        "--colname-pair-trait",
        "g",
        "--colname-pair-variant",
        "v",
        "--max-chunk-vars",
        "10",
    ]


def test_rect_assoc() -> None:
    assert qc.rect_assoc(Q, genotypes=GL, pheno="ph", var_list="v.txt", out="r.tsv.gz") == [
        Q,
        "rect-assoc",
        "--pgen-list",
        "prepare/genotypes/pgen_list.tsv",
        "--pheno",
        "ph",
        "--var-list",
        "v.txt",
        "--out",
        "r.tsv.gz",
    ]
    argv = qc.rect_assoc(
        Q,
        genotypes=GL,
        pheno="ph",
        var_list="v",
        out="o",
        cov="c",
        pheno_list="t.txt",
        sample="s",
        rint_after_adj=True,
        min_maf=0.01,
        min_mac=5,
        max_chunk_vars=100,
    )
    assert argv[10:] == [
        "--cov",
        "c",
        "--pheno-list",
        "t.txt",
        "--sample",
        "s",
        "--rint-after-adj",
        "--min-maf",
        "0.01",
        "--min-mac",
        "5",
        "--max-chunk-vars",
        "100",
    ]
    assert "--rint-before-adj" in qc.rect_assoc(
        Q, genotypes=GL, pheno="p", var_list="v", out="o", rint_before_adj=True
    )


def test_to_shell_quotes() -> None:
    assert qc.to_shell([Q, "pair-prs", "--out", "a b/$x"]) == f"{Q} pair-prs --out 'a b/$x'"


# --------------------------------------------------------------------------- locate
def _exe(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\necho '   - pair-prs : x'\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


@pytest.fixture
def clean_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """No env var, an empty PATH, and no submodule build."""
    monkeypatch.delenv(loc.ENV_VAR, raising=False)
    empty = tmp_path / "emptybin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    monkeypatch.setattr(loc, "submodule_dir", lambda: None)
    return empty


def test_locate_order(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_env: Path) -> None:
    explicit = _exe(tmp_path / "a" / "qpgentools")
    env = _exe(tmp_path / "b" / "qpgentools")
    on_path = _exe(clean_env / "qpgentools")
    sub = tmp_path / "sub"
    _exe(sub / "bin" / "qpgentools")
    monkeypatch.setattr(loc, "submodule_dir", lambda: sub)

    monkeypatch.setenv(loc.ENV_VAR, str(env))
    assert loc.locate(explicit) == loc.Located(explicit.resolve(), "explicit")
    assert loc.locate() == loc.Located(env.resolve(), "env")
    monkeypatch.delenv(loc.ENV_VAR)
    assert loc.locate() == loc.Located(on_path.resolve(), "path")
    on_path.unlink()
    assert loc.locate() == loc.Located(sub / "bin" / "qpgentools", "submodule")
    assert loc.find_qpgentools() == sub / "bin" / "qpgentools"


def test_locate_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_env: Path) -> None:
    with pytest.raises(loc.QpgenNotFound, match="not found or not executable"):
        loc.locate(tmp_path / "missing")
    plain = tmp_path / "plain"
    plain.write_text("")
    os.chmod(plain, 0o644)
    with pytest.raises(loc.QpgenNotFound, match="not found or not executable"):
        loc.locate(plain)
    # an explicitly set env var that is wrong is an error, not a silent fallback
    _exe(clean_env / "qpgentools")
    monkeypatch.setenv(loc.ENV_VAR, str(tmp_path / "gone"))
    with pytest.raises(loc.QpgenNotFound, match=loc.ENV_VAR):
        loc.locate()
    monkeypatch.delenv(loc.ENV_VAR)
    (clean_env / "qpgentools").unlink()
    with pytest.raises(loc.QpgenNotFound, match="anqchor qpgen build"):
        loc.locate()


def test_check_parses_help(tmp_path: Path) -> None:
    fake = tmp_path / "qpgentools"
    fake.write_text(
        "#!/bin/sh\ncat <<'EOF'\nAvailable Commands:\n"
        "   - pair-assoc      : Perform pairwise association\n"
        "   - pair-prs        : Perform pairwise PRS\n"
        "   - rect-assoc      : rect\nEOF\n"
    )
    fake.chmod(0o755)
    assert loc.available_commands(fake) == ["pair-assoc", "pair-prs", "rect-assoc"]
    assert loc.check(fake) == ["match-prs-pheno"]


def test_check_real_binary(qpgen_bin: Path) -> None:
    assert loc.check(qpgen_bin) == []
    assert loc.missing_flags(qpgen_bin) == []


def test_missing_flags(tmp_path: Path) -> None:
    # an older qpgen 0.2.0 without categorical covariates: only pair-assoc knows the flags
    fake = tmp_path / "qpgentools"
    fake.write_text(
        '#!/bin/sh\nif [ "$1" = pair-assoc ]; then echo "--cat-cov-names --max-cat-levels"; fi\n'
    )
    fake.chmod(0o755)
    assert loc.missing_flags(fake) == [
        "match-prs-pheno --cat-cov-names",
        "match-prs-pheno --max-cat-levels",
        "rect-assoc --cat-cov-names",
        "rect-assoc --max-cat-levels",
    ]


def test_cat_cov_flags() -> None:
    for build in (
        lambda **k: qc.pair_assoc(Q, genotypes=GL, pheno="ph", pairs="p", out="o", **k),
        lambda **k: qc.rect_assoc(Q, genotypes=GL, pheno="ph", var_list="v", out="o", **k),
        lambda **k: qc.match_prs_pheno(Q, prs="prs", pheno="ph", out_prefix="o", **k),
    ):
        argv = build(cov="c", cat_cov_names=["sex", "batch"], max_cat_levels=12)
        i = argv.index("--cat-cov-names")
        assert argv[i : i + 4] == ["--cat-cov-names", "sex,batch", "--max-cat-levels", "12"]
        assert "--max-cat-levels" not in build(cov="c", cat_cov_names=["sex"])
        # without covariates neither flag is passed
        assert "--cat-cov-names" not in build(cat_cov_names=["sex"], max_cat_levels=12)
        assert "--max-cat-levels" not in build(cat_cov_names=["sex"], max_cat_levels=12)
