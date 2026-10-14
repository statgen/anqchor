# qpgentools pair-assoc

## Summary 

`qpgentools pair-assoc` performs rapid association tests between selected variant-trait pairs.

A typical running example command is given below:

```bash
qpgentools pair-assoc --pgen-list [list] --pheno [pheno] --cov [cov] --pairs [pairs] --out [out_prefix]
```

## Required options

* Input genotype files : Either `--pgen-list` or `--pgen`, `--psam`, `--pivar` options are required. See [Genotype file formats](../formats/genotypes.md) for more details.
* `--pheno` : Input phenotype matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details.
* `--pairs` : Input TSV file containing `[trait_id]` and `[variant_id]` or `[trait_id]` and `[region]` at each line to perform association tests
* `--out` : Output TSV file to store the association analysis output 

## Additional Options

* `--cov-names` : Comma-separated covariates (column names in `--cov`) to use, as REGENIE's `--covarColList`. Ranges such as `PC{1:10}` are expanded. Covariates in `--cat-cov-names` are always used. Other columns are ignored and not parsed, so they may hold text. Default: all columns. See [Categorical and selected covariates](../formats/phenotypes.md#categorical-and-selected-covariates).
* `--cat-cov-names` : Comma-separated covariates to treat as categorical, as REGENIE's `--catCovarList`. Their values may be any strings. A covariate with L levels among the analyzed samples is replaced by L-1 indicator columns.
* `--max-cat-levels` : Maximum number of levels of a categorical covariate, as REGENIE's `--maxCatLevels` (default: 10). More levels stop the run.
* `--loco` : REGENIE step-1 LOCO predictions of a single trait (`.loco` file). Only valid when one trait is tested. See [LOCO adjustment](#loco-adjustment).
* `--pred` : REGENIE step-1 `*_pred.list` file for one or more traits, with one `[pheno_id] [path_to_loco]` line per trait. Every tested trait must be listed. Cannot be combined with `--loco`.
* `--space-delim` : Read input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) as delimited by single spaces instead of tabs (default: tab-delimited). See [Phenotype file formats](../formats/phenotypes.md#delimiters).
* `--cov` : Input covariate matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details. When provided, **both** the phenotypes and the genotypes are residualized against the covariates (see [Model](#model)).
* `--missing-str` : Comma-separated strings marking missing phenotype or covariate values (default: `NA`). See [Samples in each test](#samples-in-each-test).
* `--mean-impute-geno` : Replace missing genotypes by the mean of the called genotypes instead of excluding the samples from that variant's test (default: false).
* `--sample` : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
* `--pheno-format` : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
* `--cov-format` : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'. Use 'tsv-sample-col' when each row is a covariate and each column a sample
* `--skip-rint` : Skip tests based on rank-based inverse normal transformation (default: false). With `--cov`, the RINT is applied to the covariate-adjusted phenotype, as `rect-assoc --rint-after-adj` does.
* `--icol-pivar-idx` : 1-based column index for the variant index in the pvar file (default: 9)
* `--colname-pheno-sample` : When `--sample` is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')
* `--colname-geno-sample` : When `--sample` is provided, the column name for the sample IDs in the genotype file (default: 'geno')
* `--colname-pair-trait` : Column name of the trait ID in the pair file (default: 'trait')
* `--colname-pair-variant` : Column name of the variant ID in the pair file (default: 'variant')
* `--colname-pair-region` : Column name of the region string in the pair file (default: 'region')

## Model

Each (trait, variant) pair is tested by the ordinary least-squares fit
`y ~ 1 + covariates + genotype`. Without `--cov` the model is `y ~ 1 + genotype`, so the
phenotype and the genotype are both centered.

The fit uses Frisch-Waugh-Lovell: the phenotype and the genotype are both residualized against
the intercept and the covariates, and the slope of the residuals is tested. `BETA` is therefore
the partial effect of the variant, and `BETA`, `SE`, `TSTAT` and `LOG10P` equal those of the full
OLS fit. Adjusting only the phenotype, as early versions of this command did, leaves genotype
variance that is collinear with the covariates (e.g. genotype PCs) in the design, which shrinks
`BETA` and the test statistic.

### Samples in each test

* The samples are those shared by the genotype, phenotype, and covariate files (and `--sample`,
  if given). They may appear in any order in each file; they are analyzed in genotype-file order.
* Samples with a missing covariate value (any string listed in `--missing-str`, default `NA`) are
  excluded from all tests.
* Samples with a missing phenotype value are excluded from that trait's tests. The covariates
  are fitted again on the remaining samples.
* By default, samples with a missing genotype (hard call or dosage) are excluded from that
  variant's test, and the covariates are fitted again on the remaining samples (complete-case
  analysis). `N` is then the number of samples in the regression, and
  `N = N_RR + N_RA + N_AA`.
* With `--mean-impute-geno`, a missing genotype is replaced by the mean of the called genotypes
  and every sample with an observed phenotype is kept, as in [`rect-assoc`](rect_assoc.md).
  `N` is still the number of samples in the regression. It can then exceed
  `N_RR + N_RA + N_AA`.
* `A1FREQ`, `N_RR`/`N_RA`/`N_AA` and `INFO` are computed from the called genotypes of the
  trait's samples.
* A value that is neither a number nor a missing string stops the run with an error, unless its
  column is categorical (`--cat-cov-names`, e.g. sex coded as `XX`/`XY`) or not selected by
  `--cov-names`. See [Categorical and selected covariates](../formats/phenotypes.md#categorical-and-selected-covariates).

### Degrees of freedom

`SE`, `TSTAT` and `LOG10P` (and the `_RINT` columns) use the residual degrees of freedom of the
full model: `df = N - 2 - k`, where `N` is the number of samples in the regression and `k` is the
number of linearly independent covariates. `k` is the numerical rank of the intercept and
covariates, minus 1. It equals the number of covariates when they are full rank, and it is 0
without `--cov`. The value of `k` is printed in the log.

### RINT statistics

The `_RINT` statistics use the rank-based inverse normal transform of the covariate-adjusted
phenotype of the trait's samples, tested against the same residualized genotypes. They match
[`rect-assoc`](rect_assoc.md) with `--rint-after-adj`. The transformed phenotype is no longer
exactly orthogonal to the covariates, so its residual sum of squares on the intercept and
covariates replaces `y'y` in the residual variance. The `_RINT` statistics are therefore those of
the OLS fit `RINT(adjusted y) ~ 1 + covariates + genotype`. For a variant with missing genotypes,
the transformed phenotype is not recomputed; the samples with a missing genotype are dropped
from it.

### LOCO adjustment

With `--loco` (one trait) or `--pred` (a `*_pred.list` for one or more traits), each trait is
adjusted for the polygenic background estimated by REGENIE step 1, following REGENIE step 2 for
quantitative traits. The steps are run for each trait on its samples, then again for each
chromosome of its pairs:

1. Residualize the phenotype against the intercept and the covariates.
2. Scale the residual to unit SD, with `SD = ||r|| / sqrt(N - 1 - k)` (REGENIE's `scale_Y`).
3. Subtract the LOCO prediction for the chromosome of the tested variant.
4. For `BETA`, `SE`, `TSTAT` and `LOG10P`, multiply the result back by the SD from step 2, so
   effects stay in units of the covariate-adjusted phenotype. The `_RINT` statistics use the RINT
   of the result of step 3, as REGENIE's `--apply-rerint` does.

The raw statistics therefore equal [`region-assoc`](region_assoc.md) with `--loco`/`--pred`, and
the `_RINT` statistics equal `region-assoc --rint-after-adj` with `--loco`/`--pred`. Both are
OLS fits of the LOCO-adjusted phenotype on the intercept, covariates and genotype.

Samples absent from a trait's `.loco` file are excluded from that trait's tests, with a warning.
A missing value (`NA`) in a chromosome's row stops the run with an error. Pairs are tested in
chromosome order, so each chromosome's predictions are read once per trait.

A `.loco` file has a header `FID_IID <FID>_<IID> ...` and one row per chromosome (`1`-`22`, and `23` for X); the format, the chromosome naming and the sample-ID matching are the same as in [`region-assoc`](region_assoc.md#loco-adjustment). For the same chromosome, `pair-assoc`, `rect-assoc` and `region-assoc` give identical LOCO-adjusted statistics.

```bash
qpgentools pair-assoc --pgen-list [list] --pheno [pheno] --cov [cov] --pairs [pairs] \
    --pred [regenie_step1]_pred.list --out [out]
```

### Redundant covariates and untestable variants

* A variant whose genotype is (numerically) a linear combination of the covariates cannot be
  tested. This is detected when its residual variance is below `1e-8` of its total variance, e.g.
  when the variant itself is a covariate. Its `BETA`, `SE`, `TSTAT`, and `LOG10P` (and `_RINT`) are
  reported as NaN, and the number of such tests is reported at the end of the run.
* A covariate that is constant in the analyzed samples (e.g. sex when `--sample` keeps one sex) is
  redundant, and so is one that is a linear combination of other covariates. The adjustment uses
  the numerical rank of the intercept and covariates, so a redundant covariate is adjusted for
  once and costs no degree of freedom. The numerical rank is computed on the intercept and covariates scaled to unit norm, counting the pivots of a column-pivoted QR above `1e-7` times the largest one, as R's `lm()` does. A covariate that is *nearly* a linear combination of others is therefore treated as redundant too. The command prints a warning when it
  detects redundant covariates; dropping them is still recommended.

!!! note "Changes from earlier versions"
    * Without `--cov`, earlier versions neither centered the phenotype nor centered dosages. When
      the phenotype mean was far from zero, `SE` was inflated: about 100-fold for a phenotype
      with mean 50. With dosages, `BETA` was also wrong.
    * Earlier versions read `NA` as `0`. They mean-imputed missing genotypes with `--cov` and
      dropped them without it. They also failed when the genotype-file sample order was not
      lexicographic.
    * Earlier versions used `df = n - 2` with `--cov`, ignoring the covariates. Their `SE` was too
      small by a factor of `sqrt((n - 2 - k) / (n - 2))`.
    * The first header field is `TRAIT`; earlier versions wrote `#TRAIT`.

## Expected Output

The output TSV file has a header line followed by one line per tested pair. Its columns largely follow the Regenie output format:

* `TRAIT` : Trait name
* `CHROM` : Chromosome of the variant tested
* `GENPOS` : Base position (1-based) of the variant tested
* `ID` : Variant ID in `[CHROM]:[GENPOS]:[REF]:[ALT]` format
* `ALLELE0` : Reference allele
* `ALLELE1` : Alternative allele
* `A1FREQ`: Allele frequency of the alternate allele among the called genotypes
* `N`: Number of samples in the regression (see [Samples in each test](#samples-in-each-test))
* `N_RR`: Number of RR genotypes
* `N_RA`: Number of RA genotypes
* `N_AA`: Number of AA genotypes
* `INFO`: Variant quality score (typically imputation quality)
* `BETA`: Effect size for raw phenotype association (partial effect adjusted for the covariates when `--cov` is given)
* `SE`: Standard error for raw phenotype association
* `TSTAT`: Test statistic for raw phenotype association
* `LOG10P`: Log10 p-value for raw phenotype association
* `BETA_RINT`: Effect size for rank-based inverse normal transformed phenotype
* `SE_RINT`: Standard error for rank-based inverse normal transformed phenotype
* `TSTAT_RINT`: Test statistic for rank-based inverse normal transformed phenotype
* `LOG10P_RINT`: Log10 p-value for rank-based inverse normal transformed phenotype


## Full Usage 

The full usage of `qpgentools pair-assoc` can be viewed with the `--help` option:

```
$ ./qpgentools pair-assoc --help
[bin/qpgentools pair-assoc] -- Perform pairwise association analysis for specific pairs of phenotype variant pairs

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input Options ==
   --pgen                 [STR: ]             : Input PLINK2 genotype file
   --psam                 [STR: ]             : Input PLINK2 sample file
   --pivar                [STR: ]             : Input PLINK2 index pvar file (bgzipped and tabix)
   --pgen-list            [STR: ]             : Input file containing CHROM BEG END PGEN PSAM PIVAR
   --pheno                [STR: ]             : Input phenotype matrix
   --pairs                [STR: ]             : Input file containing [trait_id] [variant_id] pairs to be tested
   --sample               [STR: ]             : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
   --cov                  [STR: ]             : Input covariate matrix (optional)
   --pheno-format         [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --cov-format           [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --cov-names            [STR: ]             : Comma-separated covariates (column names in --cov) to use, as REGENIE's --covarColList; ranges such as PC{1:10} are expanded. Covariates in --cat-cov-names are always used. Default: all columns
   --cat-cov-names        [STR: ]             : Comma-separated covariates (column names in --cov) to treat as categorical, as REGENIE's --catCovarList. Their values may be any strings; each covariate with L levels among the analyzed samples is replaced by L-1 indicator columns (the first level in sorted order is the reference)
   --max-cat-levels       [INT: 10]           : Maximum number of levels of a categorical covariate, as REGENIE's --maxCatLevels (default: 10)
   --loco                 [STR: ]             : REGENIE step-1 LOCO predictions of a single trait (.loco file). For each chromosome, the covariate-adjusted phenotype is scaled to unit SD, that chromosome's prediction is subtracted, and the result is rescaled to the original SD for BETA/SE (the _RINT statistics use the RINT of the unscaled result, as REGENIE --apply-rerint)
   --pred                 [STR: ]             : REGENIE step-1 *_pred.list file for one or more traits: one '[pheno_id] [path_to_loco]' line per trait (relative paths are resolved against the list's directory). Each tested trait must be listed. Same LOCO adjustment as --loco
   --missing-str          [STR: NA]           : Comma-separated strings marking missing phenotype/covariate values (default: 'NA'). Samples with a missing phenotype are excluded from that trait's tests; samples with a missing covariate are excluded from all tests

== Output options ==
   --out                  [STR: ]             : Output prefix

== Auxiliary options ==
   --jump-thres-bp        [INT: 1000000]      : Jump threshold in base pairs for the variant index (default: 1000000)
   --max-chunk-vars       [INT: 1000]         : Maximum number of variants to store at once in memory (default: 1000)
   --icol-pivar-idx       [INT: 9]            : 1-based column index for the variant ID in the pvar file (default: 9)
   --colname-pheno-sample [STR: pheno]        : When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')
   --colname-geno-sample  [STR: geno]         : When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'geno')
   --colname-pair-trait   [STR: trait]        : Column name of the trait ID in the pait file
   --colname-pair-variant [STR: variant]      : Column name of the variant ID in the pair file
   --colname-pair-region  [STR: region]       : Column name of the region string in the pair file
   --skip-rint            [FLG: OFF]          : Skip tests based on rank-based inverse normal transformation (default: false)
   --mean-impute-geno     [FLG: OFF]          : Replace missing genotypes by the mean of the called genotypes and keep every sample with an observed phenotype (as rect-assoc does). By default, samples with a missing genotype are excluded from that variant's test
   --space-delim          [FLG: OFF]          : Input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) are delimited by single spaces instead of tabs (default: tab-delimited)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```