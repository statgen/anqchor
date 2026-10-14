# qpgentools pair-prs

## Summary 

`qpgentools pair-prs` rapidly computes polygenic risk scores (PRS) for every sample, from
per-trait summary statistics given as `[trait_id] [variant_id] [beta] [se]` tuples. One score is
produced per (sample, trait) pair, together with its standard error.

A typical running example command is given below:

```bash
qpgentools pair-prs --pgen-list [list] --pairs [pairs] --out [out_prefix]
```

## Required options

* Input genotype files : Either `--pgen-list` or `--pgen`, `--psam`, `--pivar` options are required. See [Genotype file formats](../formats/genotypes.md) for more details.
* `--pairs` : Input TSV file containing `[trait_id] [variant_id] [beta] [se]` at each line to construct PRS.
* `--out` : Output prefix to store the PRS output.

## Additional Options

* `--space-delim` : Read input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) as delimited by single spaces instead of tabs (default: tab-delimited). See [Phenotype file formats](../formats/phenotypes.md#delimiters).
* `--sample` : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
* `--colname-pair-trait` : Column name of the trait ID in the pair file (default: 'trait')
* `--colname-pair-variant` : Column name of the variant ID in the pair file (default: 'variant')
* `--colname-pair-beta` : Column name of the effect sizes in the pair file (default: 'beta')
* `--colname-pair-se` : Column name of the standard errors in the pair file (default: 'se')
* `--icol-pivar-idx` : 1-based column index for the variant ID in the indexed pvar file (default: 9)
* `--jump-thres-bp` : Jump threshold in base pairs for the variant index (default: 1000000)
* `--max-chunk-vars` : Maximum number of variants to store at once in memory (default: 1000)

## Expected Output

The following two files are expected to be generated:

* `[out_prefix].prs.tsv.gz` : PRS values for each sample. Regenie-style layout: `FID`, `IID`, followed by one column per trait.
* `[out_prefix].se.tsv.gz` : Standard errors of the PRS values, in the same layout. Written only when standard errors are available from the `--pairs` file.

## Full Usage 

The full usage of `qpgentools pair-prs` can be viewed with the `--help` option:

```
$ ./qpgentools pair-prs --help
[bin/qpgentools pair-prs] -- Perform pairwise PRS generation based on summary statistics

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input Options ==
   --pgen                 [STR: ]             : Input PLINK2 genotype file
   --psam                 [STR: ]             : Input PLINK2 sample file
   --pivar                [STR: ]             : Input PLINK2 index pvar file (bgzipped and tabix)
   --pgen-list            [STR: ]             : Input file containing CHROM BEG END PGEN PSAM PIVAR
   --pairs                [STR: ]             : Input file containing TRAIT VARIANT BETA SE summary statistics
   --sample               [STR: ]             : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files

== Output options ==
   --out                  [STR: ]             : Output prefix

== Auxiliary options ==
   --jump-thres-bp        [INT: 1000000]      : Jump threshold in base pairs for the variant index (default: 1000000)
   --max-chunk-vars       [INT: 1000]         : Maximum number of variants to store at once in memory (default: 1000)
   --icol-pivar-idx       [INT: 9]            : 1-based column index for the variant ID in the pvar file (default: 9)
   --colname-pair-trait   [STR: trait]        : Column name of the trait ID in the pair file
   --colname-pair-variant [STR: variant]      : Column name of the variant ID in the pair file
   --colname-pair-beta    [STR: beta]         : Column name of the effect sizes in the pair file
   --colname-pair-se      [STR: se]           : Column name of the standard errors in the pair file
   --space-delim          [FLG: OFF]          : Input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) are delimited by single spaces instead of tabs (default: tab-delimited)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```