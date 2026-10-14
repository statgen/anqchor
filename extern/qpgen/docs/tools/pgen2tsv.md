# qpgentools pgen2tsv

## Summary 

`qpgentools pgen2tsv` extracts genotypes from a PLINK2 `.pgen` file into TSV format.

It is a modernized counterpart to [`geno2tsv`](geno2tsv.md). Compared to `geno2tsv`, it:

* is PGEN-only (no PLINK 1.9 `.bed` support);
* emits genotypes in **ALT-dosage coding**: `0` = REF/REF, `1` = heterozygous, `2` = ALT/ALT, `NA` = missing (`geno2tsv` uses the internal REF-count coding);
* can select variants either by explicit index (`--var`) or by genomic region (`--region`);
* supports allele frequency and allele count filters;
* opportunistically uses the `.pgen` file's sparse (difflist) representation, so that reading and, in `--sparse` output mode, writing only touch the non-reference carriers.

A typical running example command is given below:

```bash
qpgentools pgen2tsv --pfile [plink2_prefix] --region chr1:1000000-2000000 --out [out].tsv.gz
```

or, selecting variants by index:

```bash
qpgentools pgen2tsv --pfile [plink2_prefix] --var [variants.tsv] --idx-begin 1 --out [out].tsv.gz
```

## Required options

* Input genotypes : either `--pfile`, or the individual `--pgen` / `--psam` (and `--pivar`) paths.
    * `--pfile` : PLINK 2.0 file prefix. Expands to `[prefix].pgen`, `[prefix].psam`, and `[prefix].pvar.idx.gz`.
    * `--pgen`, `--psam`, `--pivar` : Explicit paths; each overrides the value derived from `--pfile`.
* Variant selection : **exactly one** of `--var` or `--region` must be given.
    * `--var` : TSV file whose first column is the variant index in the `.pgen` file; remaining columns are echoed into the output.
    * `--region` : Genomic region `CHROM:BEG-END` to stream. Requires an indexed pvar file (see [Genotype file formats](../formats/genotypes.md)).
* `--out` : Output file. Use a `.gz` suffix for bgzipped output; any other suffix writes plain text.

## Additional Options

* `--space-delim` : Read input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) as delimited by single spaces instead of tabs (default: tab-delimited). See [Phenotype file formats](../formats/phenotypes.md#delimiters).
* `--sample-list` : Comma-separated list of sample IDs to include.
* `--sample-file` : File containing the list of sample IDs to include (one per line). At most one of `--sample-list` and `--sample-file` may be given; if neither is given, all samples are included.
* `--idx-begin` : Starting index used by `--var` (default: 0). Use `1` when the variant file holds 1-based indices, such as the `INDEX` column of an indexed pvar file.
* `--icol-pivar-idx` : 1-based column index for the variant index in the indexed pvar file (default: 9).
* `--min-af` / `--max-af` : Minimum / maximum ALT allele frequency (defaults: 0.0 and 1.0).
* `--min-ac` / `--max-ac` : Minimum / maximum ALT allele count (defaults: 0 and 1e18).
* `--sparse` : Emit `[SAMPLE_INDEX]:[GENOTYPE]` for non-REF/REF samples only, instead of one column per sample.

Variants with `AN = 0` (no called genotype among the selected samples) are always dropped.

## Expected Output

### With `--region`

The header is `#CHROM`, `POS`, `ID`, `REF`, `ALT`, `AC`, `AN`, followed by one column per selected
sample, named by its sample ID.

* `AC` : ALT allele count among the selected samples
* `AN` : Number of called alleles among the selected samples (2 x number of non-missing samples)

### With `--var`

The leading columns of the `--var` file (all columns *except* the first index column) are echoed
first, followed by `AC`, `AN`, and one column per selected sample.

* If the `--var` file has a header line starting with a single `#`, its column names are reused. Lines starting with `##` are skipped.
* If it has no header, the echoed columns are named `V1`, `V2`, ... (the numbering follows the original 1-based column positions, so it starts at `V1` for the second column of the file).

### Genotype encoding

* Default (dense) output : `0` = REF/REF, `1` = heterozygous, `2` = ALT/ALT, `NA` = missing.
* `--sparse` output : only non-REF/REF samples are written, as `[i]:[g]` where `i` is the 0-based index of the sample **within the selected sample subset** and `g` is `1`, `2`, or `NA`. Samples not listed are REF/REF.

## Full Usage 

The full usage of `qpgentools pgen2tsv` can be viewed with the `--help` option:

```
$ ./qpgentools pgen2tsv --help
[bin/qpgentools pgen2tsv] -- Extract PGEN genotypes (0/1/2 ALT coding) into TSV with AF/AC filters and sparse support

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input genotypes ==
   --pfile          [STR: ]             : Input PLINK 2.0 file prefix (expects .pgen/.psam/.pvar.idx.gz)
   --pgen           [STR: ]             : Input PLINK 2.0 genotype file (.pgen); overrides --pfile
   --psam           [STR: ]             : Input PLINK 2.0 sample file (.psam); overrides --pfile
   --pivar          [STR: ]             : Input tabixed/indexed pvar file (required for --region)
   --sample-list    [STR: ]             : Comma-separated list of samples to include (optional)
   --sample-file    [STR: ]             : File containing the list of samples to include (optional)

== Variant selection (use one of --var or --region) ==
   --var            [STR: ]             : Input variant file (TSV) with [variant index] [extra columns]
   --idx-begin      [INT: 0]            : Starting index for --var (default: 0). Use 1 for 1-based indices
   --region         [STR: ]             : Genomic region CHROM:BEG-END to stream (requires --pivar)
   --icol-pivar-idx [INT: 9]            : 1-based column index for the variant ID in the pvar file (default: 9)

== Variant filters ==
   --min-af         [FLT: 0.00]         : Minimum ALT allele frequency (default: 0.0)
   --max-af         [FLT: 1.00]         : Maximum ALT allele frequency (default: 1.0)
   --min-ac         [FLT: 0.00]         : Minimum ALT allele count (default: 0)
   --max-ac         [FLT: 1000000000000000000.00] : Maximum ALT allele count (default: 1e18)

== Output options ==
   --out            [STR: ]             : Output file (use .gz suffix for bgzipped output)
   --sparse         [FLG: OFF]          : Emit sparse [IDX]:[GENO] for non-REF/REF samples instead of full genotypes
   --space-delim    [FLG: OFF]          : Input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) are delimited by single spaces instead of tabs (default: tab-delimited)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
