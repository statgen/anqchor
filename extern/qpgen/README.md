# qpgen
Libraries and Command Line Tools for PLINK files

## Overview

`qpgen` is a collection of C++ libraries and tools to help analysis involving PLINK-formatted files. This repository is under development and is not fully documented yet.

Full documentation is available at [https://hyunminkang.github.io/qpgen/](https://hyunminkang.github.io/qpgen/).

## Available commands

* xQTL analysis
  * `pair-assoc` : association tests for an explicit list of (trait, variant) pairs
  * `pair-prs` : polygenic risk scores from per-trait summary statistics
  * `rect-assoc` : association tests between a variant list and a set of traits
  * `region-assoc` : association tests for every variant in a region, with optional **SuSiE fine-mapping** (`--susie`), including the SuSiE-inf and SuSiE-ash unmappable-effects models
  * `match-prs-pheno` : sample identity matching between PRS and phenotype matrices
* Sequence data
  * `match-plp-geno` : concordance between sequence pileups and genotypes
  * `join-plp-geno` : join pileups and individual genotypes for inspection
* Genotype extraction
  * `geno2tsv` : extract genotypes from PLINK2/PLINK 1.9 files into TSV
  * `pgen2tsv` : extract PGEN genotypes (0/1/2 ALT coding) into TSV, with region streaming, AF/AC filters, and sparse output

## Installing qpgen

You can install `qpgen` by following the instructions below.

```bash
## clone the repository
git clone --recursive https://github.com/hyunminkang/qpgen.git
cd qpgen

## build the submodules
cd submodules
sh -x build.sh
cd ..

## build qpgen
mkdir build
cd build
cmake ..
make

## list available package
../bin/qpgentools --help
```

To see the usage of individual commands, type:

```bash
../bin/qpgentools [command] --help
```

To see the compiled library, type:

```bash
ls -l ../lib/libqpgen.a
```

## Preparing input files

`qpgen` reads standard PLINK2 `.pgen` / `.psam` files, plus an "indexed pvar" file
(`[prefix].pvar.idx.gz`): a bgzipped, tabix-indexed variant file with the variant's row index
appended as a 9th column. Create one with:

```bash
PREFIX=/path/to/plink/prefix

( echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tINDEX';
  grep -v ^# ${PREFIX}.pvar | awk -v OFS='\t' '{ print $0, FNR }'; ) \
  | bgzip -c > ${PREFIX}.pvar.idx.gz
tabix -pvcf ${PREFIX}.pvar.idx.gz
```

See [the genotype format documentation](https://hyunminkang.github.io/qpgen/formats/genotypes/)
for `.pvar.zst` inputs and for indexing many chromosomes in parallel.

## Advanced Options for Installation

In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

<pre>
For libqgen:
  - $ cmake -DQGEN_INCLUDE_DIRS=/qgenlib_absolute_path/include/ -DQGEN_LIBRARIES=/qgenlib_absolute_path/lib/libqgen.a ..

For libhts:
  - $ cmake -DHTS_INCLUDE_DIRS=/htslib_absolute_path/include/  -DHTS_LIBRARIES=/htslib_absolute_path/lib/libhts.a ..

For bzip2:
  - $ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - $ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
</pre>

A system LAPACK/BLAS is optional but recommended for `region-assoc --susie --unmappable-effects inf/ash`: it is detected automatically (Accelerate on macOS). Choose a specific library with `-DBLA_VENDOR=OpenBLAS` (or `Intel10_64lp`, etc.), or disable it with `-DQPGEN_USE_LAPACK=OFF`.

Finally, to build the binary, run

<pre>
$ make
</pre>
