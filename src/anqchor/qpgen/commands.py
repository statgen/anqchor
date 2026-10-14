"""Typed builders for the qpgentools command lines anqchor uses.

This is the only module that knows qpgen flag names. Each builder returns an argv
list; ``None`` options are omitted and booleans become bare flags. Use
:func:`to_shell` to turn an argv into a command string for a :class:`~anqchor.engine.Step`.
"""

from __future__ import annotations

import shlex
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

PathLike = str | Path


THREADS_ENV = "ANQCHOR_THREADS"
"""Environment variable through which ``--threads`` / ``qpgen.threads`` reaches the
``match-prs-pheno`` step. The thread count does not change results, so it is not a step
parameter (changing it re-runs nothing)."""


@dataclass(frozen=True)
class Genotypes:
    """Genotype input: a pgen list, or one prefix with ``.pgen/.psam/.pvar.idx.gz``."""

    pgen_list: Path | None = None
    pfile: Path | None = None

    def args(self) -> list[str]:
        """The genotype flags of a qpgen command line."""
        if self.pgen_list is not None:
            return ["--pgen-list", str(self.pgen_list)]
        if self.pfile is not None:
            p = str(self.pfile)
            return ["--pgen", p + ".pgen", "--psam", p + ".psam", "--pivar", p + ".pvar.idx.gz"]
        raise ValueError("Genotypes needs pgen_list or pfile")


def _opt(argv: list[str], flag: str, value: object | None) -> None:
    if value is None or value is False:
        return
    if value is True:
        argv.append(flag)
    else:
        argv += [flag, str(value)]


def to_shell(argv: Sequence[PathLike]) -> str:
    """Quote ``argv`` as one POSIX shell command line."""
    return shlex.join(str(a) for a in argv)


def pair_prs(
    qpgen: PathLike,
    *,
    genotypes: Genotypes,
    pairs: PathLike,
    out_prefix: PathLike,
    sample: PathLike | None = None,
    colname_trait: str = "trait",
    colname_variant: str = "variant",
    colname_beta: str = "beta",
    colname_se: str = "se",
    max_chunk_vars: int | None = None,
) -> list[str]:
    """``pair-prs``: writes ``<out_prefix>.prs.tsv.gz`` and ``<out_prefix>.se.tsv.gz``."""
    argv = [
        str(qpgen),
        "pair-prs",
        *genotypes.args(),
        "--pairs",
        str(pairs),
        "--out",
        str(out_prefix),
    ]
    _opt(argv, "--sample", sample)
    if colname_trait != "trait":
        _opt(argv, "--colname-pair-trait", colname_trait)
    if colname_variant != "variant":
        _opt(argv, "--colname-pair-variant", colname_variant)
    if colname_beta != "beta":
        _opt(argv, "--colname-pair-beta", colname_beta)
    if colname_se != "se":
        _opt(argv, "--colname-pair-se", colname_se)
    _opt(argv, "--max-chunk-vars", max_chunk_vars)
    return argv


def match_prs_pheno(
    qpgen: PathLike,
    *,
    prs: PathLike,
    pheno: PathLike,
    out_prefix: PathLike,
    cov: PathLike | None = None,
    sample_tsv: PathLike | None = None,
    trait_tsv: PathLike | None = None,
    pheno_format: str = "regenie",
    cov_format: str = "regenie",
    cat_cov_names: Sequence[str] = (),
    max_cat_levels: int | None = None,
    missing_str: Sequence[str] = ("NA",),
    rint: bool = True,
    mahalanobis: bool = False,
    lam: float | None = None,
    auto_lambda: bool = False,
    min_weight: float | None = None,
    z_threshold: float = 1.96,
    z_diff: float = 2.0,
    threads: int = 1,
    auto_lambda_min_self: int | None = None,
    auto_lambda_min_best: int | None = None,
    cov_impute_mean: bool = False,
) -> list[str]:
    """``match-prs-pheno``: writes ``<out>.weights.tsv.gz``, ``<out>.match.assigned.tsv.gz``
    and ``<out>.match.all.tsv.gz``.

    ``lam`` and ``auto_lambda`` need ``mahalanobis`` (without it qpgen scores traits as
    independent, which equals ``lam = 1``). ``auto_lambda_min_self`` and
    ``auto_lambda_min_best`` (qpgen default 100 each; ``0`` disables the latter) set the
    sample counts below which ``--auto-lambda`` is an error."""
    if auto_lambda and not mahalanobis:
        raise ValueError("auto_lambda requires mahalanobis")
    if lam is not None and not mahalanobis:
        raise ValueError("lam requires mahalanobis")
    if auto_lambda and lam is not None:
        raise ValueError("give either lam or auto_lambda, not both")
    if lam is not None and not 0.0 <= lam <= 1.0:
        raise ValueError(f"lam must be between 0 and 1, got {lam}")
    if (auto_lambda_min_self is not None or auto_lambda_min_best is not None) and not auto_lambda:
        raise ValueError("auto_lambda_min_self/auto_lambda_min_best require auto_lambda")
    argv = [
        str(qpgen),
        "match-prs-pheno",
        "--prs",
        str(prs),
        "--pheno",
        str(pheno),
        "--out",
        str(out_prefix),
    ]
    _opt(argv, "--cov", cov)
    _opt(argv, "--sample-tsv", sample_tsv)
    _opt(argv, "--trait-tsv", trait_tsv)
    if pheno_format != "regenie":
        _opt(argv, "--pheno-format", pheno_format)
    if cov is not None and cov_format != "regenie":
        _opt(argv, "--cov-format", cov_format)
    if cov is not None and cat_cov_names:
        _opt(argv, "--cat-cov-names", ",".join(cat_cov_names))
        _opt(argv, "--max-cat-levels", max_cat_levels)
    if list(missing_str) != ["NA"]:
        _opt(argv, "--missing-str", ",".join(missing_str))
    _opt(argv, "--rint", rint)
    _opt(argv, "--mahalanobis", mahalanobis)
    if mahalanobis:
        if auto_lambda:
            _opt(argv, "--auto-lambda", True)
            _opt(argv, "--auto-lambda-min-self", auto_lambda_min_self)
            _opt(argv, "--auto-lambda-min-best", auto_lambda_min_best)
        elif lam is not None:
            _opt(argv, "--lambda", f"{lam:g}")
    _opt(argv, "--cov-impute-mean", cov_impute_mean)
    if min_weight:
        _opt(argv, "--min-weight", f"{min_weight:g}")
    _opt(argv, "--z-threshold", f"{z_threshold:g}")
    _opt(argv, "--z-diff", f"{z_diff:g}")
    _opt(argv, "--threads", max(1, threads))
    return argv


def pair_assoc(
    qpgen: PathLike,
    *,
    genotypes: Genotypes,
    pheno: PathLike,
    pairs: PathLike,
    out: PathLike,
    cov: PathLike | None = None,
    pheno_format: str = "regenie",
    cov_format: str = "regenie",
    cat_cov_names: Sequence[str] = (),
    max_cat_levels: int | None = None,
    sample: PathLike | None = None,
    skip_rint: bool = False,
    colname_trait: str = "trait",
    colname_variant: str = "variant",
    max_chunk_vars: int | None = None,
) -> list[str]:
    """``pair-assoc``: ``out`` is a file name (bgzipped when it ends in ``.gz``)."""
    argv = [
        str(qpgen),
        "pair-assoc",
        *genotypes.args(),
        "--pheno",
        str(pheno),
        "--pairs",
        str(pairs),
        "--out",
        str(out),
    ]
    _opt(argv, "--cov", cov)
    if pheno_format != "regenie":
        _opt(argv, "--pheno-format", pheno_format)
    if cov is not None and cov_format != "regenie":
        _opt(argv, "--cov-format", cov_format)
    if cov is not None and cat_cov_names:
        _opt(argv, "--cat-cov-names", ",".join(cat_cov_names))
        _opt(argv, "--max-cat-levels", max_cat_levels)
    _opt(argv, "--sample", sample)
    _opt(argv, "--skip-rint", skip_rint)
    if colname_trait != "trait":
        _opt(argv, "--colname-pair-trait", colname_trait)
    if colname_variant != "variant":
        _opt(argv, "--colname-pair-variant", colname_variant)
    _opt(argv, "--max-chunk-vars", max_chunk_vars)
    return argv


def rect_assoc(
    qpgen: PathLike,
    *,
    genotypes: Genotypes,
    pheno: PathLike,
    var_list: PathLike,
    out: PathLike,
    cov: PathLike | None = None,
    pheno_list: PathLike | None = None,
    pheno_format: str = "regenie",
    cov_format: str = "regenie",
    cat_cov_names: Sequence[str] = (),
    max_cat_levels: int | None = None,
    sample: PathLike | None = None,
    rint_before_adj: bool = False,
    rint_after_adj: bool = False,
    min_maf: float | None = None,
    min_mac: float | None = None,
    max_chunk_vars: int | None = None,
) -> list[str]:
    """``rect-assoc``: ``out`` is a file name (bgzipped when it ends in ``.gz``)."""
    argv = [
        str(qpgen),
        "rect-assoc",
        *genotypes.args(),
        "--pheno",
        str(pheno),
        "--var-list",
        str(var_list),
        "--out",
        str(out),
    ]
    _opt(argv, "--cov", cov)
    _opt(argv, "--pheno-list", pheno_list)
    if pheno_format != "regenie":
        _opt(argv, "--pheno-format", pheno_format)
    if cov is not None and cov_format != "regenie":
        _opt(argv, "--cov-format", cov_format)
    if cov is not None and cat_cov_names:
        _opt(argv, "--cat-cov-names", ",".join(cat_cov_names))
        _opt(argv, "--max-cat-levels", max_cat_levels)
    _opt(argv, "--sample", sample)
    _opt(argv, "--rint-before-adj", rint_before_adj)
    _opt(argv, "--rint-after-adj", rint_after_adj)
    if min_maf is not None:
        _opt(argv, "--min-maf", f"{min_maf:g}")
    if min_mac is not None:
        _opt(argv, "--min-mac", f"{min_mac:g}")
    _opt(argv, "--max-chunk-vars", max_chunk_vars)
    return argv
