# Tools in qpgentools

`qpgentools` is a collection of C++ tools that are designed to help process modified version of PLINK2 files to rapidly perform quality control, normalization, and evaluation of analysis pipelines.
Currently, the following tools are available

* xQTL Analysis Tools
  * [`qpgentools pair-assoc`](tools/pair_assoc.md): Perform rapid association tests between selected variant-trait pairs.
  * [`qpgentools pair-prs`](tools/pair_prs.md): Perform rapid calculation of polygenic risk scores (PRS) based on selected variant-trait pairs
  * [`qpgentools rect-assoc`](tools/rect_assoc.md): Perform rapid association tests between all pairs between selected subset of variants and selected subset of traits 
  * [`qpgentools region-assoc`](tools/region_assoc.md): Test every variant in a genomic region against one or more traits, with optional SuSiE fine-mapping (credible sets, PIPs, log Bayes factors)
  * [`qpgentools match-prs-pheno`](tools/match_prs_pheno.md): Perform sample identity matching based on PRS and phenotype data.
* Sequence Data Tools
  * [`qpgentools match-plp-geno`](tools/match_plp_geno.md): Identify the best matching samples based on sequence pileups
  * [`qpgentools join-plp-geno`](tools/join_plp_geno.md): Join the sequence pileups and genotypes for visual inspection
* Genotype Extraction Tools
  * [`qpgentools geno2tsv`](tools/geno2tsv.md): Extract genotypes from PLINK2 or PLINK 1.9 files to TSV format.
  * [`qpgentools pgen2tsv`](tools/pgen2tsv.md): Extract PGEN genotypes (0/1/2 ALT coding) to TSV, with region streaming, AF/AC filters, and sparse output.

## Choosing between the association tools

| Tool | Variants tested | Traits tested | Typical use |
|---|---|---|---|
| [`pair-assoc`](tools/pair_assoc.md) | An explicit list of (trait, variant) pairs | One per pair | Re-testing or validating specific hits |
| [`rect-assoc`](tools/rect_assoc.md) | An explicit variant list, streamed in chunks | All, or a selected subset | A "rectangular" slice of summary statistics across many traits |
| [`region-assoc`](tools/region_assoc.md) | Every variant in one region, held in memory at once | A selected subset | Region-based analysis and SuSiE fine-mapping |

!!! note
    This documentation contains only a subset of primary tools implemented in `qpgentools`.
    Other tools not documented here (e.g. `pair-assoc-v0`, `test-qpgen`) are considered non-primary
    tools that have limited support, and they are NOT supported officially.
