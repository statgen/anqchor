# Genotype File Format for `qpgentools`

## Overview

`qpgentools` uses PLINK2 format genotype files as the primary input format, but it requires a modified version of `pvar` file, which is called as "indexed pvar" format, to allow fast partial access.

## Indexed PVAR format

The indexed PVAR format is a VCF-like format to store variant information, except that it has additional column to indicate the variant index as integer. The expected header is as follows (tab-delimited):

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    INDEX
```

The first 8 columns are the standard VCF columns, and the last column is the variant index as 1-based integer. The variant index is used to quickly access the variant in the genotype file.

These files are bgzipped and tabix-indexed, so that it can be accessed with standard tabix commands or libraries. In `qpgentools`, this file is called as `pivar` file, and the bgzipped files are typically named as `[prefix].pvar.idx.gz`, and the tabix index file is named as `[prefix].pvar.idx.gz.tbi`.

!!! note
    The `INDEX` column must list the variants in exactly the same order as they appear in the
    `.pgen` file, because `qpgentools` uses it as the row index for a direct seek. Do not sort,
    filter, or deduplicate the `.pvar` records before assigning the index.

    The index column is 1-based by default. Its position is configurable through the
    `--icol-pivar-idx` option (default: `9`), so extra columns can be inserted before it
    as long as the option is adjusted accordingly.

## Creating an indexed pvar file

An indexed pvar file is created from the `.pvar` (or `.pvar.zst`) file that PLINK2 produced,
by (1) writing the 9-column header, (2) appending the 1-based record number to each variant
line, and (3) bgzipping and tabix-indexing the result.

### Single PLINK2 dataset

For one dataset with a plain (uncompressed) `.pvar` file:

```bash
PREFIX=/path/to/plink/prefix   ## expects ${PREFIX}.pgen / ${PREFIX}.psam / ${PREFIX}.pvar

( echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tINDEX';
  grep -v ^# ${PREFIX}.pvar | awk -v OFS='\t' '{ print $0, FNR }'; ) \
  | bgzip -c > ${PREFIX}.pvar.idx.gz
tabix -pvcf ${PREFIX}.pvar.idx.gz
```

If the `.pvar` file is Zstandard-compressed (`.pvar.zst`, PLINK2's default when
`--out ... vzs` is used), decompress it on the fly with `zstd -dc`:

```bash
( echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tINDEX';
  zstd -dc ${PREFIX}.pvar.zst | grep -v ^# | awk -v OFS='\t' '{ print $0, FNR }'; ) \
  | bgzip -c > ${PREFIX}.pvar.idx.gz
tabix -pvcf ${PREFIX}.pvar.idx.gz
```

!!! warning
    `grep -v ^#` removes both the `##` meta lines and the original `#CHROM` header line, so that
    `FNR` counts only variant records and the index matches the `.pgen` row order.

    Some `.pvar` files contain fewer than 8 columns (e.g. only `#CHROM POS ID REF ALT`). In that
    case the pasted `INDEX` column will not land in column 9, and you must either pad the missing
    `QUAL`/`FILTER`/`INFO` columns or pass the actual column index through `--icol-pivar-idx`.
    For example, to pad a 5-column `.pvar`:

    ```bash
    awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, ".", ".", ".", FNR }'
    ```

### Many datasets in parallel (one per chromosome)

Genotype datasets are usually split by chromosome. The example below indexes chromosomes 1-22
in parallel, running 10 jobs at a time:

```bash
OUT=/path/to/genotypes/chr{}.GTEx_v8_WGS_IMPUTED.dose

seq 1 22 | xargs -I {} -P 10 bash -c "(echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tINDEX'; zstd -dc ${OUT}.pvar.zst | grep -v ^# | awk -v OFS='\t' '{ print \$0, FNR }';) | bgzip -c > ${OUT}.pvar.idx.gz; tabix -pvcf ${OUT}.pvar.idx.gz"
```

Here `{}` is substituted by `xargs` with the chromosome number, so `${OUT}` expands to a
per-chromosome path such as `.../chr7.GTEx_v8_WGS_IMPUTED.dose`. The `\$0` escape is required
because the whole `bash -c` argument is double-quoted, and `$0` must reach `awk` unexpanded.

Adjust `-P 10` to the number of parallel jobs your filesystem can sustain, and use
`seq 1 22; echo X` (or a similar list) if the sex chromosomes need to be included.

!!! tip
    On macOS and other BSD systems, `xargs -I` caps the assembled command at 255 bytes and fails
    with `xargs: command line cannot be assembled, too long`. Add `-S 4096` (right after `-I {}`)
    to raise the limit. GNU `xargs`, as found on Linux, has no such restriction.

!!! warning
    `tabix` requires the records to be sorted by position within each chromosome, and it will fail
    with `Unsorted positions on sequence #N` otherwise. PLINK2 writes sorted `.pvar` files, so this
    normally only bites when the `.pvar` was manually concatenated or edited. Note that the sort
    must be applied to the `.pvar` **before** the `INDEX` column is assigned, and the resulting
    order must still match the `.pgen` row order --- if it does not, the `.pgen` itself needs to be
    rebuilt with `plink2 --sort-vars`.

### Verifying the result

The number of variant lines in the indexed pvar file must equal the variant count in the `.pgen`
file, and the last `INDEX` value must equal that count:

```bash
bgzip -dc ${PREFIX}.pvar.idx.gz | grep -v ^# | tail -1 | cut -f 9   ## last index
bgzip -dc ${PREFIX}.pvar.idx.gz | grep -vc ^#                       ## number of variants
tabix ${PREFIX}.pvar.idx.gz chr1:1000000-1010000 | head        ## random access works
```

## Specifying genotype files

In `qpgentools`, many commands that accepted this modified PLINK2 format, including 
`qpgentools pair-assoc`, `qpgentools pair-prs`, `qpgentools rect-assoc`, and
`qpgentools region-assoc`, can accept the genotype files in two ways:

1. Using `--pgen`, `--psam`, and `--pivar` options, which takes the genotype, sample, and pvar files, respectively. 
    * `--pgen` : PLINK2 binary genotype file, typically named as `[prefix].pgen`
    * `--psam` : PLINK2 binary sample file, typically named as `[prefix].psam`
    * `--pivar` : Indexed pvar file, typically named as `[prefix].pvar.idx.gz`

2. If the genotype files are separated by chromosomes or chunks, using `--pgen-list` option, which takes a TSV file containing the list of genotype files separated by chunk. Each line is expected to have the following format

```
[CHROM] [BEG] [END] [PLINK2_PREFIX]
```
or 

```
[CHROM] [BEG] [END] [PGEN_PATH] [PSAM_PATH] [PIVAR_PATH]
```

If `[PLINK2_PREFIX]` is provided, `[PLINK2_PREFIX].pgen`, `[PLINK2_PREFIX].psam`, and `[PLINK2_PREFIX].pvar.idx.gz` are used as the genotype, sample, and pvar files, respectively.

For example, if the genotype files are separated by chromosomes, the `--pgen-list` option can be written as follows:

```
1   1   1000000000  /path/to/plink/prefix_chr1
2   1   1000000000  /path/to/plink/prefix_chr2
...
22  1   1000000000  /path/to/plink/prefix_chr22
23  1   1000000000  /path/to/plink/prefix_chrX
```

The chromosome names in the first column must match the names used in the indexed pvar files
(e.g. use `chr1` rather than `1` when the pvar file uses `chr1`).

## Dosages and missing genotypes

* When the `.pgen` file stores dosages, `qpgentools` reads the dosage values. In
  `rect-assoc` and `region-assoc`, a missing dosage is replaced by the variant's mean dosage.
* When the `.pgen` file stores hard calls, `rect-assoc` and `region-assoc` mean-impute missing
  genotypes (i.e. set them to the allele-frequency-based mean) for association and fine-mapping.
* `pair-assoc` instead excludes samples with a missing genotype or dosage from that variant's
  test (complete-case analysis), unless `--mean-impute-geno` is given. See
  [pair-assoc](../tools/pair_assoc.md#samples-in-each-test).
* Genotype columns are mean-centered before analysis, and monomorphic variants are dropped.
