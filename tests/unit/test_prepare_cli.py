"""``anqchor prepare`` commands."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

from typer.testing import CliRunner

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fixtures" / "io"))
from anq_io_fixtures import make_geno, needs_htslib, tsv, write_text

from anqchor.cli.prepare import app, expand_chroms

runner = CliRunner()


def test_expand_chroms() -> None:
    assert expand_chroms("1-3,X, Y") == ["1", "2", "3", "X", "Y"]


def test_pheno_cmd(tmp_path: Path) -> None:
    src = write_text(
        tmp_path / "x.bed",
        tsv([["#Chr", "start", "end", "ID", "A", "B"], ["chr1", 1, 2, "T1", 0.5, "."]]),
    )
    out = tmp_path / "o.tsv.gz"
    r = runner.invoke(
        app,
        [
            "pheno",
            "--in",
            str(src),
            "--format",
            "tensorqtl",
            "--out",
            str(out),
            "--missing-str",
            ".",
        ],
    )
    assert r.exit_code == 0, r.output
    assert gzip.decompress(out.read_bytes()).decode() == "FID\tIID\tT1\nA\tA\t0.5\nB\tB\tNA\n"
    r = runner.invoke(app, ["pheno", "--in", str(src), "--format", "tensorqtl", "--out", str(out)])
    assert r.exit_code == 2 and "non-numeric" in r.output


def test_cov_cmd(tmp_path: Path) -> None:
    src = write_text(
        tmp_path / "c.tsv",
        tsv(
            [
                ["ID", "S1", "S2", "S3"],
                ["sex", "XX", "XY", "XY"],
                ["batch", 1, 2, 3],
                ["PC1", 1, "NA", 3],
            ]
        ),
    )
    out = tmp_path / "c.out.tsv.gz"
    args = [
        "cov",
        "--in",
        str(src),
        "--format",
        "tsv-sample-col",
        "--out",
        str(out),
        "--categorical",
        "batch",
    ]
    r = runner.invoke(app, args)
    assert r.exit_code == 2 and "missing covariate values" in r.output
    r = runner.invoke(app, [*args, "--missing", "impute_mean", "--json", str(tmp_path / "c.json")])
    assert r.exit_code == 0, r.output
    lines = gzip.decompress(out.read_bytes()).decode().splitlines()
    # default: categorical covariates keep their text levels (qpgen --cat-cov-names)
    assert lines[0] == "FID\tIID\tsex\tbatch\tPC1" and lines[2] == "S2\tS2\tXY\t2\t2"
    assert "--cat-cov-names sex,batch" in r.output
    assert (tmp_path / "c.json").exists()
    # --encode writes the numeric indicator columns instead
    r = runner.invoke(app, [*args, "--missing", "impute_mean", "--encode"])
    assert r.exit_code == 0, r.output
    lines = gzip.decompress(out.read_bytes()).decode().splitlines()
    assert lines[0] == "FID\tIID\tsex\tbatch_2\tbatch_3\tPC1" and lines[2] == "S2\tS2\t2\t2\t1\t2"


@needs_htslib
def test_pgen_list_cmd(tmp_path: Path) -> None:
    make_geno(tmp_path, style="chr")
    out = tmp_path / "l.tsv"
    r = runner.invoke(
        app,
        [
            "pgen-list",
            "--prefix-template",
            str(tmp_path / "geno" / "set.chr{chrom}"),
            "--chroms",
            "1-2",
            "--out",
            str(out),
        ],
    )
    assert r.exit_code == 0, r.output
    lines = [x.split("\t") for x in out.read_text().splitlines()]
    assert [x[0] for x in lines] == ["chr1", "chr2"] and lines[0][3] == str(
        (tmp_path / "geno" / "set.chr1").resolve()
    )
    r = runner.invoke(
        app,
        [
            "pgen-list",
            "--prefix-template",
            str(tmp_path / "geno" / "set.chr{chrom}"),
            "--chroms",
            "1-3",
            "--out",
            str(out),
        ],
    )
    assert r.exit_code == 2 and "set.chr3" in r.output
    r = runner.invoke(
        app,
        [
            "pgen-list",
            "--prefix-template",
            "nowhere/c{chrom}",
            "--chroms",
            "X",
            "--out",
            str(out),
            "--no-check",
        ],
    )
    assert r.exit_code == 0 and out.read_text().startswith("X\t1\t")


@needs_htslib
def test_index_pvar_cmd(tmp_path: Path) -> None:
    (tmp_path / "x.pvar").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t10\ta\tA\tG\t50\tPASS\tAC=1\n"
    )
    r = runner.invoke(app, ["index-pvar", "--pfile", str(tmp_path / "x")])
    assert r.exit_code == 0, r.output
    lines = gzip.decompress((tmp_path / "x.pvar.idx.gz").read_bytes()).decode().splitlines()
    assert lines[1] == "1\t10\ta\tA\tG\t50\tPASS\tAC=1\t1"
    r = runner.invoke(app, ["index-pvar", "--pfile", str(tmp_path / "x")])
    assert r.exit_code == 2 and "--force" in r.output
