# Quickstart for qpgen

## Installing qpgen

Please follow the instruction below to install `qpgentools`

```sh
git clone --recursive https://github.com/hyunminkang/qpgen.git
cd qpgen
cd submodules
sh -x build.sh
cd ..
mkdir build
cd build
cmake ..
make
```

## List available tools

To list the available tools, run the following command:

```sh
../bin/qpgentools --help
```

To see the usage of an individual command, type:

```sh
../bin/qpgentools [command] --help
```

If you encounter any difficulties, see [Install](install.md) for more details.

## Prepare the indexed pvar files

`qpgentools` reads standard PLINK2 `.pgen` and `.psam` files, but it needs the variant information
in an "indexed pvar" file (`[prefix].pvar.idx.gz`) that carries each variant's row index and is
bgzipped and tabix-indexed. Build one from an existing PLINK2 dataset with:

```bash
PREFIX=/path/to/plink/prefix

( echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tINDEX';
  grep -v ^# ${PREFIX}.pvar | awk -v OFS='\t' '{ print $0, FNR }'; ) \
  | bgzip -c > ${PREFIX}.pvar.idx.gz
tabix -pvcf ${PREFIX}.pvar.idx.gz
```

See [Genotype file formats](formats/genotypes.md) for Zstandard-compressed `.pvar` files, for
indexing many chromosomes in parallel, and for how to describe a chromosome-split dataset with
`--pgen-list`.

## Run an association analysis

Test every variant near a gene against one trait, and fine-map the region with SuSiE:

```bash
../bin/qpgentools region-assoc \
    --pgen ${PREFIX}.pgen --psam ${PREFIX}.psam --pivar ${PREFIX}.pvar.idx.gz \
    --pheno expression.bed.gz --pheno-format tensorqtl \
    --cov covariates.tsv \
    --region chr12:6534517-6538371 \
    --traits ENSG00000111640 \
    --susie --output-lbf \
    --out out/ENSG00000111640
```

This writes `out/ENSG00000111640.assoc.tsv.gz` (marginal statistics),
`out/ENSG00000111640.susie.cs.tsv.gz` (credible sets and PIPs), and
`out/ENSG00000111640.susie.lbf.tsv.gz` (per-variant log Bayes factors).
See [`region-assoc`](tools/region_assoc.md) for the full column definitions.

## Where to go next

* [Genotype file formats](formats/genotypes.md) and [Phenotype file formats](formats/phenotypes.md)
* [Tools](tools.md) for the list of available commands
