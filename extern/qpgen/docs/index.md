# Welcome to qpgentools documentation

This is the documentation for the `qpgentools` toolkit. 

!!! warning
    This documentation is for an actively developed software package.
    
    Note that features and usage may change frequently.

    This documentation provides information for the tools available in the `main` branch
    of the repository. However, the `dev` branch may contain newer tools or features that are not yet documented here.

## What is qpgentools?

`qpgentools` is a collection of C++ tools that are designed to help process modified version of PLINK2 files to
rapidly perform quality control, normalization, and evaluation of analysis pipelines.
These tools are under active development, so they may change frequently. 

The toolkit is built around an "indexed pvar" file (`[prefix].pvar.idx.gz`), a bgzipped and
tabix-indexed variant file that carries the variant's row index in the `.pgen` file. This lets
`qpgentools` jump directly to an arbitrary variant or genomic region without scanning the whole
genotype file. See [Genotype file formats](formats/genotypes.md) for the format definition and
for a recipe to build these files from existing PLINK2 datasets.

## Documentation Overview

This documentation provides the following information:

* [Quickstart](quickstart.md): A quick guide to get started with `qpgentools`.
* [Install](install.md): How to install `qpgentools`.
* [Formats](formats/genotypes.md): Input file formats accepted by `qpgentools`.
* [Tools](tools.md): Documentation of individual tools implemented in `qpgentools`.

## Recent additions

* [`qpgentools region-assoc`](tools/region_assoc.md) : region-wide association analysis for one or
  more traits, with optional [SuSiE fine-mapping](tools/region_assoc.md#susie-fine-mapping)
  (credible sets, PIPs, and log Bayes factors), including the `SuSiE-inf` and `SuSiE-ash`
  unmappable-effects models.
* [`qpgentools pgen2tsv`](tools/pgen2tsv.md) : PGEN-only genotype extraction with region streaming,
  allele frequency/count filters, and a sparse output mode.
* [`qpgentools rect-assoc`](tools/rect_assoc.md) now residualizes genotypes against the covariates
  (Frisch-Waugh-Lovell), so effect sizes are proper partial effects when `--cov` is used.
