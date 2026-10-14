# anqchor

**Genetically anchored quality control for large-scale multi-omics xQTL analysis.**

anqchor uses known xQTLs from reference studies to check a multi-omics dataset before an
expensive xQTL mapping run:

- **match**: is every omics sample paired with the right genotype? (PRS of known xQTLs vs
  measured traits; sample swaps, mislabels, rescue candidates)
- **cis**: which covariate set / normalisation replicates the most known cis-xQTLs?
- **trans**: do known cis and trans pairs replicate, are there hotspots, and is a trans analysis setup calibrated (λ_GC, QQ)?

The statistics run in the [qpgen](https://github.com/hyunminkang/qpgen) C++ toolkit (a git
submodule). anqchor prepares inputs, runs the steps through a generated `Makefile` (incremental,
parallel), and writes TSV/JSON summaries and interactive HTML reports.

## Install

```bash
git clone --recursive https://github.com/hyunminkang/anqchor.git
cd anqchor
uv sync                       # or: pip install -e .
uv run anqchor qpgen build    # builds extern/qpgen/bin/qpgentools (macOS: brew install libomp llvm)
uv run anqchor qpgen check
```

Requires Python ≥ 3.10, GNU make, a C++14 compiler with OpenMP and cmake (for qpgen). The
pinned qpgen submodule (`qpgentools --version`: `qpgentools 0.2.0 (git 5b811e2)`) is needed for
categorical covariates.

## Quickstart

The repository ships a small example in `examples/mage/data/`: 729 individuals with MAGE RNA-seq
data, their genotypes on chromosomes 21 and 22 (256,645 variants with MAF ≥ 5%), covariates,
and 367 known cis-eQTLs from an independent study (GTEx, same cell type). Each analysis runs
on its own:

```bash
D=examples/mage/data
IN=(--geno-list $D/geno/mage.pgen_list.tsv --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz
    --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz)
COV=(--cov $D/pheno/mage_covariates.regenie.tsv.gz)

# sample identity: analysis, then the HTML report
anqchor match run "${IN[@]}" "${COV[@]}" --out results/match
anqchor match report --out results/match

# cis evaluation: one run per setup, compare, re-render the page when needed
anqchor cis run "${IN[@]}"                --out results/cis/nocov
anqchor cis run "${IN[@]}" "${COV[@]}"    --out results/cis/all_cov
anqchor cis compare results/cis/nocov results/cis/all_cov --out results/cis/compare
anqchor cis report --out results/cis/compare     # re-render the page only

# trans evaluation: one run per setup, compare, re-render the page when needed
GENES=(--trait-annot $D/pheno/mage_genes.chr21_22.bed)
anqchor trans run "${IN[@]}" "${GENES[@]}"            --out results/trans/nocov
anqchor trans run "${IN[@]}" "${GENES[@]}" "${COV[@]}" --cov-include sex,PC1-5,PEER1-10 \
  --out results/trans/peer10
anqchor trans compare results/trans/nocov results/trans/peer10 --out results/trans/compare
anqchor trans report --out results/trans/compare   # re-render the page only
```

`cis compare` and `trans compare` write `report/cis.html` and `report/trans.html` in the
comparison directory; `cis report --out DIR` and `trans report --out DIR` re-render them.

(`IN`, `COV` and `GENES` are shell arrays, so the commands work in both bash and zsh.)

Identity matching confirms 727 of 729 samples; `NA11919` clearly does not match its own
genotypes, but no other individual is a clear match. In the cis evaluation, 178 of the 367
known pairs replicate at p < 1e-5 without covariates and 251 with all 66 covariates. In the trans evaluation, the cross-chromosome λ_GC
is 1.74 without covariates (37 hotspots: 1-Mb bins with trans signals for ≥ 5 genes at
p < 1e-5) and 1.00 with sex, 5 PCs and 10 PEER factors (no hotspot).

Covariates with text values, such as `sex` (`XX`/`XY`) in the example, are categorical:
anqchor keeps the text and passes `--cat-cov-names sex` to qpgen, which encodes them as
indicator columns. For categories coded as numbers (a batch number), add
`--cat-cov-names batch`.

## Documentation

<https://hyunminkang.github.io/anqchor/> (sources in `docs/`; `uv sync --extra docs && uv run mkdocs serve`).

## License

Apache-2.0. See [LICENSE](LICENSE).
