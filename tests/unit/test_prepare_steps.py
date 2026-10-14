"""prepare_genotypes / prepare_reference / prepare_dataset / prepare_mapping steps."""

from __future__ import annotations

import gzip
import json
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fixtures" / "io"))
from anq_io_fixtures import SAMPLES, make_geno, make_prefix, needs_htslib, step_ctx, tsv, write_text

from anqchor.pipelines.layout import Layout as L
from anqchor.steps import StepError
from anqchor.steps.prepare import (
    prepare_dataset,
    prepare_genotypes,
    prepare_mapping,
    prepare_reference,
)

pytestmark = needs_htslib

PAIRS = [
    ["trait", "variant", "beta", "se", "pval"],
    ["G1.5:AAA", "chr1:100:A:G", 0.5, 0.1, 1e-4],
    ["G1.5:AAA", "chr1:200:C:CA", 0.2, 0.1, 0.01],
    ["G2.1:BBB", "chr2:50:T:C", -0.3, 0.1, 1e-3],
    ["G3.2:CCC", "chr2:999:G:A", 0.3, 0.1, 1e-3],
    ["G3.2:CCC", "chr2:150:AT:A", 0.3, 0.1, 1e-3],
]
TRAITS = ["G1.5:AAA", "G2.1:BBB", "G3.2:CCC", "G4.1:DDD"]


def j(run: Path, rel: Path) -> Any:
    return json.loads((run / rel).read_text())


def run_genotypes(tmp: Path, style: str = "plain") -> Path:
    run = tmp / "run"
    lst = make_geno(tmp, style=style)
    prepare_genotypes(step_ctx(run, "prepare.genotypes", {"pgen_list": str(lst), "pfile": None}))
    return run


def ref_params(tmp: Path, **kw: Any) -> dict[str, Any]:
    pairs = write_text(tmp / "ref.tsv.gz", tsv(PAIRS))
    p = {
        "reference": "r",
        "kind": "cis",
        "pairs": str(pairs),
        "variants": None,
        "columns": {
            "trait": "trait",
            "variant": "variant",
            "beta": "beta",
            "se": "se",
            "pval": "pval",
        },
    }
    p.update(kw)
    return p


def regenie(
    tmp: Path, samples: list[str], traits: list[str] = TRAITS, name: str = "pheno.tsv.gz"
) -> Path:
    rows = [["FID", "IID", *traits]] + [
        [s, s, *(0.1 * (i + k) for k in range(len(traits)))] for i, s in enumerate(samples)
    ]
    return write_text(tmp / name, tsv(rows))


def ds_params(tmp: Path, pheno: Path, **kw: Any) -> dict[str, Any]:
    p = {
        "dataset": "d",
        "omics": "rna",
        "pheno": str(pheno),
        "pheno_format": "regenie",
        "trait_annot": None,
        "cov": None,
        "cov_format": "regenie",
        "cov_categorical": [],
        "cov_missing": "error",
        "missing_str": ["NA"],
        "sample_map": None,
        "trait_id_rule": "as_is",
        "trait_map": None,
    }
    p.update(kw)
    return p


# ------------------------------------------------------------------- genotypes
def test_genotypes_pgen_list(tmp_path: Path) -> None:
    run = run_genotypes(tmp_path)
    lines = (run / L.pgen_list).read_text().splitlines()
    assert len(lines) == 2 and all(
        len(x.split("\t")) == 4 and x.split("\t")[3].startswith("/") for x in lines
    )
    assert lines[0].split("\t")[3] == str((tmp_path / "geno" / "set.chr1").resolve())
    assert (run / L.geno_samples).read_text().split() == SAMPLES
    g = j(run, L.geno_json)
    assert (g["n_samples"], g["chrom_style"], g["chromosomes"], g["n_entries"]) == (
        6,
        "plain",
        ["1", "2"],
        2,
    )
    assert g["schema_version"] == 1 and set(g["entries"][0]) == {
        "chrom",
        "beg",
        "end",
        "pgen",
        "psam",
        "pivar",
    }


@pytest.mark.parametrize("order", ["psam_pivar", "pivar_psam"])
def test_genotypes_six_columns_nonstandard_names(tmp_path: Path, order: str) -> None:
    make_geno(tmp_path)
    g = tmp_path / "geno"
    for c in ("1", "2"):
        for suf in (".pgen", ".psam", ".pvar.idx.gz", ".pvar.idx.gz.tbi"):
            (g / f"set.chr{c}{suf}").rename(g / f"odd{c}{suf.replace('.pvar.idx.gz', '.idx.gz')}")
    rows = []
    for c in ("1", "2"):
        files = [f"geno/odd{c}.pgen", f"geno/odd{c}.psam", f"geno/odd{c}.idx.gz"]
        if order == "pivar_psam":
            files = [files[0], files[2], files[1]]
        rows.append("\t".join([c, "1", "1000", *files]))
    lst = tmp_path / "six.tsv"
    lst.write_text("\n".join(rows) + "\n")
    run = tmp_path / "run"
    prepare_genotypes(step_ctx(run, "g", {"pgen_list": str(lst), "pfile": None}))
    line = (run / L.pgen_list).read_text().splitlines()[0].split("\t")
    assert line[:3] == ["1", "1", "1000"] and len(line) == 4
    assert Path(line[3] + ".psam").resolve() == (g / "odd1.psam").resolve()
    assert Path(line[3] + ".pvar.idx.gz.tbi").resolve() == (g / "odd1.idx.gz.tbi").resolve()
    assert j(run, L.geno_json)["entries"][0]["psam"] == str((g / "odd1.psam").resolve())


def test_genotypes_pfile(tmp_path: Path) -> None:
    prefix = make_prefix(
        tmp_path / "all",
        [("chr1", 5, "A", "C"), ("chr1", 9, "A", "C"), ("chr2", 3, "G", "T")],
        psam_header="#FID\tIID\tSEX",
    )
    run = tmp_path / "run"
    prepare_genotypes(step_ctx(run, "g", {"pgen_list": None, "pfile": str(prefix)}))
    g = j(run, L.geno_json)
    assert (
        g["chromosomes"] == ["chr1", "chr2"] and g["chrom_style"] == "chr" and g["n_samples"] == 6
    )


def test_genotypes_style_mismatch(tmp_path: Path) -> None:
    lst = make_geno(tmp_path, style="plain")
    lst.write_text(lst.read_text().replace("1\t1\t", "chr1\t1\t", 1))
    with pytest.raises(StepError, match=r"without a 'chr' prefix.*'chr1'"):
        prepare_genotypes(step_ctx(tmp_path / "run", "g", {"pgen_list": str(lst), "pfile": None}))


def test_genotypes_sample_mismatch_and_missing_files(tmp_path: Path) -> None:
    lst = make_geno(tmp_path)
    psam = tmp_path / "geno" / "set.chr2.psam"
    psam.write_text("#IID\nS1\nS2\nZ9\n")
    with pytest.raises(StepError, match=r"share one sample list.*Z9"):
        prepare_genotypes(step_ctx(tmp_path / "run", "g", {"pgen_list": str(lst), "pfile": None}))
    (tmp_path / "geno" / "set.chr2.pvar.idx.gz.tbi").unlink()
    with pytest.raises(StepError, match=r"not found.*index-pvar"):
        prepare_genotypes(step_ctx(tmp_path / "run", "g", {"pgen_list": str(lst), "pfile": None}))


# ------------------------------------------------------------------ references
@pytest.mark.parametrize("no_tabix", [False, True])
def test_reference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, no_tabix: bool) -> None:
    run = run_genotypes(tmp_path)
    if no_tabix:
        monkeypatch.setenv("ANQCHOR_NO_TABIX", "1")
    prepare_reference(step_ctx(run, "r", ref_params(tmp_path)))
    pairs = pd.read_csv(run / L.ref_pairs("r"), sep="\t", dtype=str)
    assert list(pairs.columns) == ["trait", "variant", "beta", "se", "pval"]
    assert pairs["variant"].tolist()[0] == "1:100:A:G"
    v = pd.read_csv(run / L.ref_variants("r"), sep="\t", dtype=str)
    assert list(v.columns) == ["variant", "chrom", "pos", "ref", "alt", "is_indel", "present"]
    assert dict(zip(v["variant"], v["present"], strict=True)) == {
        "1:100:A:G": "1",
        "1:200:C:CA": "1",
        "2:50:T:C": "1",
        "2:150:AT:A": "1",
        "2:999:G:A": "0",
    }
    assert (run / L.ref_variants_present("r")).read_text().split() == [
        "1:100:A:G",
        "1:200:C:CA",
        "2:50:T:C",
        "2:150:AT:A",
    ]
    r = j(run, L.ref_json("r"))
    assert (
        r["n_pairs"],
        r["n_traits"],
        r["n_variants"],
        r["n_present"],
        r["n_absent"],
        r["n_absent_snv"],
        r["n_absent_indel"],
        r["chrom_renamed"],
    ) == (5, 3, 5, 4, 1, 1, 0, 5)
    assert r["presence_method"] == ("python" if no_tabix else "tabix")
    assert r["per_chrom"]["2"] == {"n_pairs": 3, "n_variants": 3, "n_present": 2}


def test_reference_none_present(tmp_path: Path) -> None:
    run = run_genotypes(tmp_path)
    bad = write_text(tmp_path / "bad.tsv", tsv([PAIRS[0], ["T", "1:101:A:G", 0.1, 0.1, 0.1]]))
    with pytest.raises(StepError, match=r"chromosome naming.*genome"):
        prepare_reference(step_ctx(run, "r", ref_params(tmp_path, pairs=str(bad))))


def test_trans_reference_variants_only(tmp_path: Path) -> None:
    run = run_genotypes(tmp_path, style="chr")
    vl = tmp_path / "v.txt"
    vl.write_text("1:100:A:G\n2:50:T:C\n2:51:T:C\n")
    prepare_reference(
        step_ctx(run, "r", ref_params(tmp_path, kind="trans", pairs=None, variants=str(vl)))
    )
    r = j(run, L.ref_json("r"))
    assert (r["n_pairs"], r["n_variants"], r["n_present"]) == (0, 3, 2)
    assert (run / L.ref_variants_present("r")).read_text().split() == [
        "chr1:100:A:G",
        "chr2:50:T:C",
    ]
    assert (
        gzip.decompress((run / L.ref_pairs("r")).read_bytes()).decode()
        == "trait\tvariant\tbeta\tse\tpval\n"
    )


# -------------------------------------------------------------------- datasets
def test_dataset_regenie_cov_annot(tmp_path: Path) -> None:
    run = run_genotypes(tmp_path)
    pheno = regenie(tmp_path, ["S1", "S2", "S3", "NOGENO"])
    cov = write_text(
        tmp_path / "cov.tab.gz",
        tsv(
            [
                ["id", "S1", "S2", "S3", "NOGENO"],
                ["sex", "XX", "XY", "XY", "XX"],
                ["PC1", 0.1, 0.2, 0.3, 0.4],
                ["k", 1, 1, 1, 1],
            ]
        ),
    )
    annot = write_text(
        tmp_path / "a.bed",
        tsv([["#chr", "start", "end", "id"], ["chr1", 10, 11, "G1.5"], ["chr2", 20, 21, "G2.1"]]),
    )
    p = ds_params(
        tmp_path, pheno, cov=str(cov), cov_format="tsv-sample-col", trait_annot=str(annot)
    )
    prepare_dataset(step_ctx(run, "d", p))
    assert (run / L.ds_pheno("d")).read_bytes() == pheno.read_bytes()
    assert (run / L.ds_samples("d")).read_text().split() == ["S1", "S2", "S3", "NOGENO"]
    assert (run / L.ds_sample_map("d")).read_text() == "S1\tS1\nS2\tS2\nS3\tS3\n"
    traits = pd.read_csv(run / L.ds_traits("d"), sep="\t", dtype=str, keep_default_na=False)
    assert traits.values.tolist() == [
        ["G1.5:AAA", "1", "10", "11"],
        ["G2.1:BBB", "2", "20", "21"],
        ["G3.2:CCC", "NA", "NA", "NA"],
        ["G4.1:DDD", "NA", "NA", "NA"],
    ]
    c = j(run, L.ds_cov_json("d"))
    assert c["columns"] == ["sex", "PC1"] and c["categorical"] == {"sex": {"XX": 1, "XY": 2}}
    assert c["dropped_constant"] == ["k"] and c["groups"]["sex"] == ["sex"]
    # categorical covariates keep their text levels in the prepared file
    assert (
        gzip.decompress((run / L.ds_cov("d")).read_bytes()).decode().splitlines()[1]
        == "S1\tS1\tXX\t0.1"
    )
    assert c["categorical_names"] == ["sex"]
    d = j(run, L.ds_json("d"))
    assert (
        d["n_samples_pheno"],
        d["n_traits"],
        d["n_samples_mapped"],
        d["n_pheno_unmapped"],
        d["n_cov_columns"],
        d["n_traits_with_loci"],
        d["trait_annot_rule"],
        d["pheno_fast_path"],
    ) == (4, 4, 3, 1, 2, 2, "before_colon", True)
    assert d["n_samples_analysable"] == 3 and d["warnings"] == []


def test_dataset_tensorqtl_with_sample_map(tmp_path: Path) -> None:
    run = run_genotypes(tmp_path, style="chr")
    rows = [
        ["#Chr", "start", "end", "ID", "p1", "p2", "p3"],
        ["1", 5, 6, "G1.5", 1, 2, 3],
        ["2", 7, 8, "G2.1", 4, "NA", 6],
    ]
    pheno = write_text(tmp_path / "x.bed.gz", tsv(rows))
    smap = tmp_path / "map.tsv"
    smap.write_text("NWD_ID\tSAMPLE_ID\nS1\tp1\nS2\tp2\nGX\tp3\nS3\tpX\n")
    prepare_dataset(
        step_ctx(
            run, "d", ds_params(tmp_path, pheno, pheno_format="tensorqtl", sample_map=str(smap))
        )
    )
    out = gzip.decompress((run / L.ds_pheno("d")).read_bytes()).decode().splitlines()
    assert out == [
        "FID\tIID\tG1.5\tG2.1",
        "p1\tp1\t1.0\t4.0",
        "p2\tp2\t2.0\tNA",
        "p3\tp3\t3.0\t6.0",
    ]
    assert (run / L.ds_sample_map("d")).read_text() == "S1\tp1\nS2\tp2\n"
    assert (run / L.ds_traits("d")).read_text().splitlines()[1] == "G1.5\tchr1\t5\t6"
    d = j(run, L.ds_json("d"))
    assert (
        d["n_samples_mapped"],
        d["n_pheno_unmapped"],
        d["n_missing_values"],
        d["n_sample_map_rows_geno_missing"],
    ) == (2, 1, 1, 1)
    # cis/trans ignore the sample map: no phenotype ID is a genotype ID here
    assert d["n_samples_analysable"] == 0 and any(
        "sample_map is used by match only" in w for w in d["warnings"]
    )


def test_dataset_reports_samples_lost_to_covariates(tmp_path: Path) -> None:
    run = run_genotypes(tmp_path)
    pheno = regenie(tmp_path, ["S1", "S2", "S3", "NOGENO"])
    cov = write_text(
        tmp_path / "cov.tsv",
        tsv([["FID", "IID", "PC1"], ["S1", "S1", 0.1], ["S2", "S2", 0.2], ["S3", "S3", "NA"]]),
    )
    p = ds_params(tmp_path, pheno, cov=str(cov), cov_missing="drop_samples")
    prepare_dataset(step_ctx(run, "d", p))
    d = j(run, L.ds_json("d"))
    assert d["n_samples_mapped"] == 3 and d["n_samples_analysable"] == 2
    assert any("no covariates" in w for w in d["warnings"])
    assert any("dropped (cov_missing)" in w for w in d["warnings"])


def test_dataset_errors(tmp_path: Path) -> None:
    run = run_genotypes(tmp_path)
    bad = write_text(tmp_path / "bad.tsv", tsv([["FID", "IID", "T1"], ["S1", "S1", "x"]]))
    with pytest.raises(StepError, match="non-numeric value 'x' for trait 'T1' in sample 'S1'"):
        prepare_dataset(step_ctx(run, "d", ds_params(tmp_path, bad)))
    dup = write_text(
        tmp_path / "dup.tsv", tsv([["FID", "IID", "T1"], ["S1", "S1", 1], ["S1", "S1", 2]])
    )
    with pytest.raises(StepError, match="duplicate sample"):
        prepare_dataset(step_ctx(run, "d", ds_params(tmp_path, dup)))
    pheno = regenie(tmp_path, ["S1", "S2", "S3"])
    cov = write_text(
        tmp_path / "c.tsv",
        tsv([["FID", "IID", "x"], ["S1", "S1", 1], ["S2", "S2", "NA"], ["S3", "S3", 2]]),
    )
    with pytest.raises(StepError, match="cov_missing"):
        prepare_dataset(step_ctx(run, "d", ds_params(tmp_path, pheno, cov=str(cov))))
    smap = tmp_path / "m.tsv"
    smap.write_text("geno_id\tpheno_id\nS1\tS1\nS2\tS1\n")
    with pytest.raises(StepError, match="several genotyped samples"):
        prepare_dataset(step_ctx(run, "d", ds_params(tmp_path, pheno, sample_map=str(smap))))


# -------------------------------------------------------------------- mappings
def _mapping_run(tmp: Path, traits: list[str]) -> Path:
    run = run_genotypes(tmp)
    prepare_reference(step_ctx(run, "r", ref_params(tmp)))
    prepare_dataset(step_ctx(run, "d", ds_params(tmp, regenie(tmp, ["S1", "S2"], traits))))
    return run


def map_params(**kw: Any) -> dict[str, Any]:
    return {"dataset": "d", "reference": "r", "trait_id_rule": "as_is", "trait_map": None, **kw}


def test_mapping_as_is(tmp_path: Path) -> None:
    run = _mapping_run(tmp_path, ["G2.1:BBB", *TRAITS[:1], "G3.2:CCC"])
    prepare_mapping(step_ctx(run, "m", map_params()))
    pairs = pd.read_csv(run / L.dsref_pairs("d", "r"), sep="\t", dtype=str)
    assert list(pairs.columns) == ["trait", "variant", "beta", "se", "pval", "ref_trait", "present"]
    assert pairs["present"].tolist() == ["1", "1", "1", "0", "1"]
    assert (run / L.dsref_traits("d", "r")).read_text().split() == [
        "G2.1:BBB",
        "G1.5:AAA",
        "G3.2:CCC",
    ]
    assert (run / L.dsref_variants("d", "r")).read_text().split() == [
        "1:100:A:G",
        "1:200:C:CA",
        "2:50:T:C",
        "2:150:AT:A",
    ]
    m = j(run, L.dsref_json("d", "r"))
    assert (
        m["n_ref_traits"],
        m["n_mapped_traits"],
        m["n_pairs_mapped"],
        m["n_pairs_present"],
        m["rule"],
    ) == (3, 3, 5, 4, "as_is")


def test_mapping_rule_and_zero_mapped(tmp_path: Path) -> None:
    run = _mapping_run(tmp_path, ["G1.5", "G2.1", "G9.9"])
    with pytest.raises(StepError, match="trait_id_rule 'before_colon' would map 2"):
        prepare_mapping(step_ctx(run, "m", map_params()))
    prepare_mapping(step_ctx(run, "m", map_params(trait_id_rule="before_colon")))
    assert (run / L.dsref_trait_map("d", "r")).read_text() == "G1.5:AAA\tG1.5\nG2.1:BBB\tG2.1\n"
    m = j(run, L.dsref_json("d", "r"))
    assert (m["n_mapped_traits"], m["n_unmapped_traits"], m["unmapped_examples"]) == (
        2,
        1,
        ["G3.2:CCC"],
    )


def test_mapping_trait_map_file(tmp_path: Path) -> None:
    run = _mapping_run(tmp_path, ["A", "B"])
    tm = tmp_path / "tm.tsv"
    tm.write_text("G1.5:AAA\tB\nG3.2:CCC\tA\nG2.1:BBB\tZ\n")
    prepare_mapping(step_ctx(run, "m", map_params(trait_map=str(tm))))
    assert (run / L.dsref_traits("d", "r")).read_text().split() == ["A", "B"]
    assert j(run, L.dsref_json("d", "r"))["rule"] == "file"
    shutil.rmtree(run)
