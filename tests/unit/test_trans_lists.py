"""steps.trans: ``trans_lists`` and the ``trait_chroms`` restriction (no qpgen needed)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from anqchor.pipelines.layout import Layout as L
from anqchor.steps import StepContext, StepError
from anqchor.steps.trans import (
    no_trait_kept_message,
    select_trait_chroms,
    trait_chrom_map,
    trans_lists,
)

TRAIT_CHROM = {"t1": "1", "t2": "2", "t21": "21", "t22": "22", "tX": "X"}
TRAITS = ["t1", "t21", "tnone", "t22", "tX", "t2"]


def test_genotyped_keeps_genotype_chromosomes() -> None:
    sel = select_trait_chroms(TRAITS, TRAIT_CHROM, "genotyped", ["21", "22"], "plain")
    assert sel.kept == ["t21", "t22"]
    assert sel.allowed == ["21", "22"]
    assert (sel.n_candidates, sel.n_off_chrom, sel.n_no_locus) == (6, 3, 1)
    assert sel.examples_no_locus == ["tnone"]
    j = sel.as_json()
    assert j["n_kept"] == 2 and j["n_dropped_off_chrom"] == 3 and j["n_dropped_no_locus"] == 1


def test_all_keeps_everything_including_traits_without_locus() -> None:
    sel = select_trait_chroms(TRAITS, TRAIT_CHROM, "all", ["21", "22"], "plain")
    assert sel.kept == TRAITS and sel.allowed is None
    assert sel.n_off_chrom == 0 and sel.n_no_locus == 0


def test_explicit_list_mixed_chromosome_styles() -> None:
    # genotype style 'chr': the trait loci and the list are normalised to it
    tc = {"a": "chr1", "b": "2", "c": "chrX", "d": "22"}
    sel = select_trait_chroms(["a", "b", "c", "d", "e"], tc, ["1", "chrX"], ["chr21"], "chr")
    assert sel.allowed == ["chr1", "chrX"]
    assert sel.kept == ["a", "c"]  # raw loci are compared as given; trait_chrom_map normalises them
    tcn = trait_chrom_map(pd.DataFrame({"trait": list(tc), "chrom": list(tc.values())}), "chr")
    assert tcn == {"a": "chr1", "b": "chr2", "c": "chrX", "d": "chr22"}
    sel = select_trait_chroms(["a", "b", "c", "d", "e"], tcn, ["1", "chrX"], ["chr21"], "chr")
    assert sel.kept == ["a", "c"] and sel.n_off_chrom == 2 and sel.n_no_locus == 1
    # plain genotype style with 'chr' names in the list
    tcp = trait_chrom_map(pd.DataFrame({"trait": list(tc), "chrom": list(tc.values())}), "plain")
    sel = select_trait_chroms(["a", "b", "c", "d"], tcp, ["chr22", "chr1"], ["21"], "plain")
    assert sel.kept == ["a", "d"] and sel.allowed == ["22", "1"]


def test_trait_chrom_map_skips_missing_loci() -> None:
    loci = pd.DataFrame({"trait": ["a", "b", "c"], "chrom": ["chr21", "NA", ""]})
    assert trait_chrom_map(loci, "plain") == {"a": "21"}


# ------------------------------------------------------------------ the step
def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)


def make_run(tmp: Path, style: str = "plain") -> Path:
    """A tiny prepared run directory: two genotyped chromosomes (21, 22)."""
    c = "chr" if style == "chr" else ""
    run = tmp / "run"
    geno = {
        "chrom_style": style,
        "chromosomes": [f"{c}21", f"{c}22"],
        "entries": [{"chrom": f"{c}21"}, {"chrom": f"{c}22"}],
        "n_samples": 3,
    }
    _write(run / L.geno_json, json.dumps(geno))
    _write(run / L.ref_variants_present("ref"), f"{c}22:300:A:G\n{c}21:100:C:T\n{c}21:100:C:T\n")
    _write(run / L.dsref_traits("ds", "ref"), "g21\ng1\n")
    _write(run / L.ds_pheno("ds"), "FID\tIID\tg1\tg21\tg22\tgnone\ns1\ts1\t1\t2\t3\t4\n")
    # trait loci in the other chromosome style than the genotypes, to exercise normalisation
    o = "" if style == "chr" else "chr"
    _write(
        run / L.ds_traits("ds"),
        f"trait\tchrom\tstart\tend\ng1\t{o}1\t10\t11\ng21\t{o}21\t50\t51\ng22\t{o}22\t60\t61\ngnone\tNA\tNA\tNA\n",
    )
    return run


def run_lists(run: Path, **params: Any) -> dict[str, Any]:
    p = {
        "dataset": "ds",
        "reference": "ref",
        "variants": "reference",
        "traits": "reference",
        "trait_chroms": "genotyped",
        **params,
    }
    ctx = StepContext(name="trans.ds.lists", run_dir=run, params=p, inputs=[], outputs=[])
    trans_lists(ctx)
    return {
        "variants": (run / L.trans_variants("ds")).read_text().split(),
        "traits": (run / L.trans_traits("ds")).read_text().split(),
        "json": json.loads((run / L.trans_lists_json("ds")).read_text()),
    }


@pytest.mark.parametrize("style", ["plain", "chr"])
def test_lists_reference_genotyped(tmp_path: Path, style: str) -> None:
    run = make_run(tmp_path, style)
    out = run_lists(run)
    c = "chr" if style == "chr" else ""
    assert out["variants"] == [f"{c}21:100:C:T", f"{c}22:300:A:G"]  # sorted, deduplicated
    assert out["traits"] == ["g21"]  # g1 is on chromosome 1, not genotyped
    t = out["json"]["traits"]
    assert t["n_candidates"] == 2 and t["n_kept"] == 1 and t["n_dropped_off_chrom"] == 1
    assert t["allowed_chroms"] == [f"{c}21", f"{c}22"]
    assert out["json"]["variants"]["n_kept"] == 2


def test_lists_all_traits_restricted_and_unrestricted(tmp_path: Path) -> None:
    run = make_run(tmp_path)
    out = run_lists(run, traits="all")
    assert out["traits"] == ["g21", "g22"]
    assert out["json"]["traits"]["n_dropped_no_locus"] == 1
    out = run_lists(run, traits="all", trait_chroms="all")
    assert out["traits"] == ["g1", "g21", "g22", "gnone"]
    assert out["json"]["traits"]["allowed_chroms"] is None
    out = run_lists(run, traits="all", trait_chroms=["chr1", "22"])
    assert out["traits"] == ["g1", "g22"]


def test_lists_trait_file_intersected_with_dataset(tmp_path: Path) -> None:
    run = make_run(tmp_path)
    f = tmp_path / "traits.txt"
    f.write_text("g22\nmissing\ng1\n")
    out = run_lists(run, traits=str(f), trait_chroms="all")
    assert out["traits"] == ["g1", "g22"]  # dataset order
    assert out["json"]["traits"]["n_missing_from_dataset"] == 1


def test_lists_restriction_selecting_nothing_is_an_actionable_error(tmp_path: Path) -> None:
    run = make_run(tmp_path)
    with pytest.raises(StepError, match="trait_chroms"):
        run_lists(run, traits="all", trait_chroms=["5"])
    # no trait has a locus at all: point at trait_annot / trait_chroms: all
    _write(run / L.ds_traits("ds"), "trait\tchrom\tstart\tend\ng21\tNA\tNA\tNA\n")
    with pytest.raises(StepError, match="trait_annot"):
        run_lists(run)


def test_no_trait_kept_message_has_no_empty_example() -> None:
    """Without off-chromosome traits the message must not print an empty '(e.g. )'."""
    sel = select_trait_chroms(["T1", "T2"], {}, "genotyped", ["21", "22"], "plain")
    assert not sel.kept and sel.n_no_locus == 2
    msg = no_trait_kept_message("rnaseq", "genotyped", "reference", 2, sel)
    assert "(e.g. )" not in msg and "0 traits on other chromosomes," in msg
    assert "no trait has a locus" in msg
    sel = select_trait_chroms(["T1"], {"T1": "1"}, "genotyped", ["21"], "plain")
    msg = no_trait_kept_message("rnaseq", "genotyped", "reference", 1, sel)
    assert "1 traits on other chromosomes (e.g. 'T1')" in msg
