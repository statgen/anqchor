# qpgentools rect-assoc

## Summary 

`qpgentools rect-assoc` performs rectangular association analysis, testing a specific list of variants against a specific list of phenotypes (or all phenotypes). This is efficient for looking up associations for a subset of variants across multiple traits, producing a "rectangular" slice of the full association summary statistics.

To test *every* variant in a genomic region instead of an explicit variant list, and to fine-map
that region with SuSiE, use [`region-assoc`](region_assoc.md).

A typical running example command is given below:

```bash
qpgentools rect-assoc --pgen [pgen] --pivar [pivar] --psam [psam] --pheno [pheno] --var-list [variants] --out [out_prefix]
```

## Required options

* Input genotype files : Either `--pgen-list` or `--pgen`, `--psam`, `--pivar` options are required. See [Genotype file formats](../formats/genotypes.md) for more details.
* `--pheno` : Input phenotype matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details.
* `--var-list` : Input file containing the list of variants to be tested. The file can contain a single column of variant IDs (e.g. `chr:pos:ref:alt`) or 4 columns (CHROM, POS, REF, ALT).
* `--out` : Output **file** to store the association analysis output. Despite the help text calling it a prefix, the value is used as the file name as-is; use a `.gz` suffix for bgzipped output.

## Additional Options

* `--cov-names` : Comma-separated covariates (column names in `--cov`) to use, as REGENIE's `--covarColList`. Ranges such as `PC{1:10}` are expanded. Covariates in `--cat-cov-names` are always used. Other columns are ignored and not parsed, so they may hold text. Default: all columns. See [Categorical and selected covariates](../formats/phenotypes.md#categorical-and-selected-covariates).
* `--cat-cov-names` : Comma-separated covariates to treat as categorical, as REGENIE's `--catCovarList`. Their values may be any strings. A covariate with L levels among the analyzed samples is replaced by L-1 indicator columns.
* `--max-cat-levels` : Maximum number of levels of a categorical covariate, as REGENIE's `--maxCatLevels` (default: 10). More levels stop the run.
* `--missing-str` : Comma-separated strings marking missing phenotype or covariate values (default: `NA`). rect-assoc tests all traits on the same samples, so a missing value in a tested phenotype or a used covariate stops the run. Missing values in traits excluded by `--pheno-list`, or in covariates not selected by `--cov-names`, are ignored.
* `--loco` : REGENIE step-1 LOCO predictions of a single trait (`.loco` file). Only valid when one trait is tested. See [LOCO adjustment](#loco-adjustment).
* `--pred` : REGENIE step-1 `*_pred.list` file for one or more traits, with one `[pheno_id] [path_to_loco]` line per trait. Every tested trait must be listed. Cannot be combined with `--loco`.
* `--space-delim` : Read input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) as delimited by single spaces instead of tabs (default: tab-delimited). See [Phenotype file formats](../formats/phenotypes.md#delimiters).
* `--pheno-list` : Input file containing the list of phenotypes to be tested (single column with phenotype IDs). If not provided, all phenotypes in the phenotype matrix will be tested. Phenotype IDs that are absent from the phenotype matrix are silently dropped.
* `--cov` : Input covariate matrix in Regenie or TSV format. See [Covariate handling](#covariate-handling) below.
* `--sample` : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files.
* `--pheno-format` : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'.
* `--cov-format` : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'.
* `--rint-before-adj` : Perform rank-based inverse normal transformation before covariate adjustment (default: false).
* `--rint-after-adj` : Perform rank-based inverse normal transformation after covariate adjustment (default: false). RINT uses Blom offsets, \(\Phi^{-1}\big((r - 3/8)/(n + 1/4)\big)\) with average ranks \(r\) for ties, as REGENIE's `--apply-rint` and `--apply-rerint` do.
* `--min-maf` : Minimum minor allele frequency to include a variant in the analysis (default: 1e-10, i.e. only monomorphic variants are excluded).
* `--min-mac` : Minimum minor allele count to include a variant in the analysis (default: 1.0).
* `--max-chunk-vars` : Maximum number of variants held in memory at once (default: 100). Variants from `--var-list` are processed in chunks of this size.
* `--jump-thres-bp` : Jump threshold in base pairs for the variant index (default: 1000000).
* `--icol-pivar-idx` : 1-based column index for the variant index in the indexed pvar file (default: 9).

## Covariate handling

When `--cov` is supplied, the phenotype matrix is adjusted for the covariates by linear regression
and the genotype matrix is residualized against the *same* covariates (Frisch-Waugh-Lovell), so
that `BETA` is the partial effect of the variant. Adjusting only the phenotype -- as earlier
versions of this command did -- leaves genotype variance that is collinear with the covariates
(e.g. genotype PCs) in the design, which shrinks `BETA` and the test statistic.

!!! note
    Effect sizes and p-values from `--cov` runs will therefore differ from those produced by
    versions of `qpgentools` released before this fix. [`region-assoc`](region_assoc.md) applies
    the same treatment, so the two commands now agree on shared variants.

    `SE`, `TSTAT` and `LOG10P` use the residual degrees of freedom of the full model
    `y ~ 1 + covariates + genotype`, `df = N - 2 - k`, where `N` is the number of analyzed samples
    and `k` is the number of linearly independent covariates (the numerical rank of the intercept
    and covariates, minus 1; the number of covariates when they are full rank). The statistics are
    therefore those of an ordinary least-squares fit with the covariates in the model; a redundant
    covariate (constant in the analyzed samples, or a linear combination of others) is adjusted
    for once and costs no degree of freedom. The numerical rank is computed on the intercept and covariates scaled to unit norm, counting the pivots of a column-pivoted QR above `1e-7` times the largest one, as R's `lm()` does. A covariate that is *nearly* a linear combination of others is therefore treated as redundant too. A variant whose genotype is (numerically) a linear combination of the covariates (residual variance below `1e-8` of its total, e.g. when the variant itself is a covariate) cannot be tested: its statistics are reported as NaN. With `--rint-after-adj`, the transformed phenotype
    is no longer exactly orthogonal to the covariates, and its residual sum of squares on the
    intercept and covariates replaces `y'y` in the residual variance. Without `--cov`,
    `df = N - 2`.

    Without `--cov`, the phenotypes are centered: the model is `y ~ 1 + genotype`. Earlier versions
    did not center them, which inflated `SE` when a phenotype's mean was far from zero (14.5-fold for a simulated
    phenotype with mean 20).

    Versions before this change used `df = N - 2` with `--cov` as well, which made `SE` too small
    by a factor of `sqrt((N - 2 - k) / (N - 2))` (about 5% for 66 covariates and 729 samples);
    `--cov` results therefore differ from those versions in `SE`, `TSTAT` and `LOG10P`, but not in
    `BETA`.

### LOCO adjustment

With `--loco` (one trait) or `--pred` (a `*_pred.list` for one or more traits), the phenotypes are
adjusted for the polygenic background estimated by REGENIE step 1. The steps are the same as in
[`region-assoc`](region_assoc.md#loco-adjustment): RINT before the adjustment (`--rint-before-adj`),
covariate residualization, scaling to unit SD, subtraction of the LOCO prediction, and RINT
(`--rint-after-adj`) or rescaling to the original SD.

The prediction depends on the chromosome, so the phenotype matrix is LOCO-adjusted again for each
chromosome of the variant list. A chunk of variants never spans two chromosomes, and the output
keeps the order of the variant list. A list that returns to a chromosome after another one
works, but that chromosome's predictions are read again, and a notice says so. Group the variants
by chromosome to avoid this.

Samples absent from any tested trait's `.loco` file are excluded, with a warning. A missing value
(`NA`) in a chromosome's row stops the run with an error. `--pheno-list` is applied before the
LOCO files are read, so `--pred` needs to list only the tested traits.

A `.loco` file has a header `FID_IID <FID>_<IID> ...` and one row per chromosome (`1`-`22`, and `23` for X); the format, the chromosome naming and the sample-ID matching are the same as in [`region-assoc`](region_assoc.md#loco-adjustment). For the same chromosome, `pair-assoc`, `rect-assoc` and `region-assoc` give identical LOCO-adjusted statistics.

### Missing genotypes

`rect-assoc` tests every trait against every variant with one matrix product, so all tests share
the same samples. A missing genotype is therefore **mean-imputed**: it is set to the mean of the
called genotypes of that variant (2 × `A1FREQ`), and the sample stays in the test. For dosage
files, the missing dosage is replaced by the variant's mean dosage when it is read. As a result:

* `N` in the output counts the samples with a called genotype (`N_RR + N_RA + N_AA`), while the
  regression and its degrees of freedom use all analyzed samples. For a variant with missing
  genotypes, `df = (number of analyzed samples) - 2 - k` is larger than `N - 2 - k`.
* An imputed genotype contributes nothing to `x'y` or `x'x`, so `BETA` is approximately unbiased
  when genotypes are missing at random. The imputed samples still add residual variance and
  degrees of freedom, so `SE` and `LOG10P` differ slightly from a complete-case analysis. The
  difference is negligible at low missing rates.

For complete-case results (each test uses only the samples with a called genotype, and `N` is
the regression sample size), use [`pair-assoc`](pair_assoc.md) for the pairs of interest.

!!! warning
    `rect-assoc` does not accept missing values in the phenotypes or in the used covariates. A
    missing value (`NA`, or another `--missing-str`) stops the run with an error; earlier versions
    read `NA` as `0`. Use [`pair-assoc`](pair_assoc.md) for traits with missing values.

## Expected Output

The output file named by `--out` (bgzipped when the name ends in `.gz`) contains the following columns. It is in a wide format where statistics for each phenotype are appended horizontally.

* `#CHROM` : Chromosome
* `GENPOS` : Genomic Position
* `ID` : Variant ID
* `ALLELE0` : Reference Allele
* `ALLELE1` : Alternative Allele
* `A1FREQ` : Frequency of Allele 1 (AC / AN)
* `N` : Number of samples with a called genotype (`N_RR + N_RA + N_AA`). Samples with a missing genotype are mean-imputed and still enter the regression; see [Missing genotypes](#missing-genotypes)
* `N_RR` : Count of Reference/Reference genotypes
* `N_RA` : Count of Reference/Alternative genotypes
* `N_AA` : Count of Alternative/Alternative genotypes
* `TEST` : Test type (e.g., "Linear")

For each phenotype `[PHENO]`, the following columns are added:

* `BETA.[PHENO]` : Effect size
* `SE.[PHENO]` : Standard Error
* `TSTAT.[PHENO]` : T-statistic
* `LOG10P.[PHENO]` : Log10 p-value

## Full Usage 

The full usage of `qpgentools rect-assoc` can be viewed with the `--help` option:

```
$ ./qpgentools rect-assoc --help    
[bin/qpgentools rect-assoc] -- Perform rectangular association analysis for multiple phenotypes and variants

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input Options ==
   --pgen            [STR: ]             : Input PLINK2 genotype file
   --psam            [STR: ]             : Input PLINK2 sample file
   --pivar           [STR: ]             : Input PLINK2 index pvar file (bgzipped and tabix)
   --pgen-list       [STR: ]             : Input file containing the list of region-specific BED files
   --pheno           [STR: ]             : Input phenotype matrix
   --loco            [STR: ]             : REGENIE step-1 LOCO predictions of a single trait (.loco file). For each chromosome, the covariate-adjusted phenotype is scaled to unit SD, that chromosome's prediction is subtracted, and (unless --rint-after-adj) the result is rescaled to the original SD so BETA/SE stay in phenotype units
   --pred            [STR: ]             : REGENIE step-1 *_pred.list file for one or more traits: one '[pheno_id] [path_to_loco]' line per trait (relative paths are resolved against the list's directory). Each tested trait must be listed. Same LOCO adjustment as --loco
   --sample          [STR: ]             : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
   --cov             [STR: ]             : Input covariate matrix (optional)
   --var-list        [STR: ]             : Input file containing the list of variants to be tested
   --pheno-list      [STR: ]             : Input file containing the list of phenotypes to be tested. If not provided, all phenotypes in the phenotype matrix will be tested.
   --pheno-format    [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --cov-format      [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --missing-str     [STR: NA]           : Comma-separated strings marking missing phenotype/covariate values (default: 'NA'). rect-assoc tests all traits on the same samples, so a missing value in a tested phenotype or a used covariate stops the run
   --cov-names       [STR: ]             : Comma-separated covariates (column names in --cov) to use, as REGENIE's --covarColList; ranges such as PC{1:10} are expanded. Covariates in --cat-cov-names are always used. Default: all columns
   --cat-cov-names   [STR: ]             : Comma-separated covariates (column names in --cov) to treat as categorical, as REGENIE's --catCovarList. Their values may be any strings; each covariate with L levels among the analyzed samples is replaced by L-1 indicator columns (the first level in sorted order is the reference)
   --max-cat-levels  [INT: 10]           : Maximum number of levels of a categorical covariate, as REGENIE's --maxCatLevels (default: 10)

== Output options ==
   --out             [STR: ]             : Output prefix

== Analysis options ==
   --rint-before-adj [FLG: OFF]          : Perform rank-based inverse normal transformation before covariate adjustment (default: false)
   --rint-after-adj  [FLG: OFF]          : Perform rank-based inverse normal transformation after covariate adjustment (default: false)
   --max-chunk-vars  [INT: 100]          : Maximum number of variants to store at once in memory (default: 100)
   --min-maf         [FLT: 1.0e-10]      : Minimum minor allele frequency to include a variant in the analysis (default: 0.0)
   --min-mac         [FLT: 1.00]         : Minimum minor allele count to include a variant in the analysis (default: 0.0)

== Auxiliary options ==
   --jump-thres-bp   [INT: 1000000]      : Jump threshold in base pairs for the variant index (default: 1000000)
   --icol-pivar-idx  [INT: 9]            : 1-based column index for the variant ID in the pvar file (default: 9)
   --space-delim     [FLG: OFF]          : Input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) are delimited by single spaces instead of tabs (default: tab-delimited)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
