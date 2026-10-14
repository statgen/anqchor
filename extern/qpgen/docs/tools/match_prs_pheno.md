# qpgentools match-prs-pheno

## Summary 

`qpgentools match-prs-pheno` matches samples between PRS (Polygenic Risk Score) and phenotype matrices. It is typically used for sample identity verification or to detect sample swaps by correlating genetically predicted risks (PRS) with observed phenotypes across multiple traits. The tool identifies the best matching PRS profile for each phenotype sample.

A detailed description of the framework, the scores, the automatic shrinkage and the treatment of missing values is given on the [methods page](match_prs_pheno_methods.md).

A typical running example command is given below:

```bash
qpgentools match-prs-pheno --prs [prs_matrix] --pheno [pheno_matrix] --out [out_prefix]
```

## Required options

* `--prs` : Input PRS matrix file. The file should contain PRS values for multiple traits across samples.
* `--pheno` : Input phenotype matrix file. The file should contain observed phenotype values for the corresponding traits.
* `--out` : Output prefix to store the matching results.

## Additional Options

* `--cov-names` : Comma-separated covariates (column names in `--cov`) to use, as REGENIE's `--covarColList`. Ranges such as `PC{1:10}` are expanded. Covariates in `--cat-cov-names` are always used. Other columns are ignored and not parsed. Default: all columns. See [Categorical and selected covariates](../formats/phenotypes.md#categorical-and-selected-covariates).
* `--cat-cov-names` : Comma-separated covariates to treat as categorical, as REGENIE's `--catCovarList`. A covariate with L levels is replaced by L-1 indicator columns. Without `--cov-impute-mean`, levels are counted after samples with a missing covariate are dropped. With it, a missing level is imputed in every indicator by the level frequencies.
* `--max-cat-levels` : Maximum number of levels of a categorical covariate, as REGENIE's `--maxCatLevels` (default: 10).
* `--space-delim` : Read input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) as delimited by single spaces instead of tabs (default: tab-delimited). See [Phenotype file formats](../formats/phenotypes.md#delimiters).
* `--cov` : Input covariate matrix file (optional). If provided, phenotypes will be adjusted for covariates before matching.
* `--sample-tsv` : TSV file containing sample ID mapping between PRS and phenotype files. Format: `[PRS_SAMPLE_ID] [PHENO_SAMPLE_ID]`. 
* `--trait-tsv` : TSV file containing trait ID mapping between PRS and phenotype files. Format: `[PRS_TRAIT_ID] [PHENO_TRAIT_ID]`.
* `--weights` : Input weights file for each phenotype. Format: `[PHENO_ID] [WEIGHT]`. If not provided, weights are estimated based on the correlation between PRS and phenotype in overlapping samples.
* `--prs-format` : Format of the PRS file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'.
* `--pheno-format` : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'.
* `--cov-format` : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'.
* `--missing-str` : Comma-separated strings that represent missing values in the phenotype and covariate matrices (default: `NA`). The PRS matrix must not contain missing values; the program stops with an error if it does.
* `--missing-as-mean` : Impute missing phenotype values with the mean of the observed values of that trait and treat them as observed from then on (default: false). See [Handling missing values](#handling-missing-values).
* `--missing-as-half-min` : Impute missing phenotype values with half of the minimum observed value of that trait and treat them as observed, the usual convention when missing means below a detection limit, as in metabolomics (default: false). Assumes positive measurements; a warning is issued for traits whose minimum is not positive. Cannot be combined with `--missing-as-mean`.
* `--cov-impute-mean` : When covariates contain missing values, mean-impute them (using the mean of observed samples per covariate) instead of dropping samples with any missing covariate (default: false). See [Handling missing values](#handling-missing-values).
* `--rint` : Perform rank-based inverse normal transformation on phenotypes after covariate adjustment (default: false).
* `--mahalanobis` : Account for correlations between traits by whitening the profiles with a shrunk covariance matrix before computing the cosine similarity (default: false). With `--lambda 1` this is identical to the default independence score, computed by a slower route.
* `--weight-prs-mh` : Weight given to the PRS covariance when the PRS and phenotype covariances are blended into the Mahalanobis metric, as `w * COV_prs + (1 - w) * COV_pheno` (default: 0.5). Only used with `--mahalanobis`.
* `--lambda` : Fixed shrinkage parameter for the Mahalanobis covariance, applied as `(1 - lambda) * COV + lambda * I`. Must be between 0 and 1 (default: 0.0, i.e. no shrinkage). Values around 0.3 to 0.6 typically perform far better than 0 or 1; see `--auto-lambda`. Only used with `--mahalanobis`.
* `--no-norm` : Use the legacy independence score, the weighted inner product divided by the total absolute weight of the observed traits, instead of the weighted cosine similarity (default: false). Not applicable with `--mahalanobis`. See the [methods page](match_prs_pheno_methods.md#21-without-mahalanobis-weighted-cosine-similarity).
* `--exact-norm` : When phenotypes have missing values, compute each PRS norm over the traits observed for the phenotyped sample instead of over all traits (default: false). Without `--mahalanobis`, or with `--lambda 1`, this costs one extra matrix product; with other lambda values it costs one pass per distinct missingness pattern.
* `--min-weight` : Minimum weight (correlation) per trait to include in the analysis; traits below this are given zero weight (default: 0.0).
* `--z-threshold` : Z-score threshold for declaring a lenient match (default: 1.96).
* `--z-diff` : Z-score difference between the best and second-best match to declare a clear match (default: 2.0).
* `--threads` : Number of threads used by Eigen for the matrix operations (default: 1).

## Auto-lambda options

These options apply with `--mahalanobis` and control the automatic choice of the shrinkage parameter.

* `--auto-lambda` : Choose the Mahalanobis shrinkage parameter automatically (default: false). The phenotyped samples with a mapped PRS sample are scored against all PRS samples for candidate values of lambda, and the value maximizing the separation of the self match is selected by golden-section search over 0 to 1. Requires `--lambda 0` (the default) and a sample mapping or shared sample IDs. Costs roughly ten scoring passes. Only used with `--mahalanobis`.
* `--auto-lambda-metric` : Criterion maximized by `--auto-lambda` (default: `mean-log-softmax`, the mean log-probability of the self match under a softmax over the Z-scores of all candidates). Alternatives: `mean-z` (mean self Z-score) and `mrr` (mean reciprocal rank of the self match).
* `--auto-lambda-min-self` : Minimum number of phenotyped samples with a mapped PRS sample and at least one observed trait required for tuning (default: 100). Fewer is an error, since the criterion would be too noisy; check the sample mapping, set `--lambda` manually, or lower the threshold.
* `--auto-lambda-min-best` : Minimum number of mapped samples whose own PRS sample ranks first at the tuned lambda (default: 100; 0 disables the check). Fewer is an error. It usually indicates a mismatched sample mapping, for example completely mismatched IDs, which also invalidates the trait weights estimated from the same mapping. The log reports the self-best count and fraction at the tuned lambda.

## Handling missing values

Missing phenotype values (strings listed in `--missing-str`, `NA` by default) are ignored by default rather than imputed. Two alternatives impute them before any other processing and then treat them as observed: `--missing-as-mean` fills each missing cell with the trait's observed mean, and `--missing-as-half-min` fills it with half of the trait's observed minimum, the usual convention when a missing measurement indicates a value below the detection limit. With imputation, `N.Traits` counts all traits with non-zero weight. Under the default, concretely:

* Each trait is standardized using the mean and standard deviation of its observed values only. Missing cells contribute nothing to any downstream inner product.
* Per-trait weights (`Weight` in the weights file) are estimated from the matched samples in which the phenotype is observed; `N.Obs` reports that count.
* For each phenotyped sample, the score is computed over its observed traits only. `N.Traits` reports the number of observed traits with non-zero weight.
* The score is a weighted cosine similarity in both modes. Its numerator and the phenotype norm run over the observed traits only. The PRS norm is computed over all traits by default; `--exact-norm` restricts it to the observed traits of each phenotyped sample. With `--mahalanobis`, the trait covariance is estimated from the standardized matrix with missing cells at zero (mean imputation), and the cross term is restricted to the observed traits on both sides.
* With `--cov`, traits with missing values are adjusted by regressing on the covariates using only their observed samples. Traits without missing values are adjusted together in a single batch.
* With `--rint`, the rank-based inverse normal transformation is applied to the observed values of each trait.
* A phenotyped sample with no observed trait of non-zero weight is reported with `MatchStatus` `NO_OBS_TRAITS`, `N.Traits` 0 and `NA` in all numeric columns.

Missing covariate values are handled by dropping every sample with at least one missing covariate from the phenotype matrix before adjustment. The number of dropped samples is reported in the log, and these samples do not appear in the output. Use `--cov-impute-mean` to keep them with mean-imputed covariates instead.

The PRS matrix is assumed to be complete. If it contains a value matching `--missing-str`, the program stops with an error.

Matching on fewer than 10 traits is unreliable, so the command warns when fewer than 10 phenotypes have a non-zero weight, and when some phenotyped samples have fewer than 10 observed traits with non-zero weight. The warning gives their number and a few examples. Such samples are still matched and reported; `N.Traits` gives their count, so they can be filtered afterwards. There is no threshold on the number of samples per trait.

## Expected Output

The following files are expected to be generated:

* `[out_prefix].weights.tsv.gz` : (If `--weights` is not provided) Calculated weights for each trait.
    * Columns: `Trait`, `Weight`, `N.Obs` (number of matched samples with observed phenotype used to estimate the weight)
* `[out_prefix].match.assigned.tsv.gz` : Summary of the best match for each phenotype sample.
    * Columns:
        * `ID.Pheno`: Sample ID in the phenotype file.
        * `N.Traits`: Number of observed traits with non-zero weight used for this sample.
        * `MatchStatus`: Status of the match (e.g., `BEST_MATCH`, `NO_MATCH`, `LENIENT_MATCH`, or `NO_OBS_TRAITS` when no trait is observed).
        * `ID.self`: Corresponding sample ID in the PRS file (based on mapping or ID match).
        * `Z.self`: Z-score of the match with itself.
        * `COR.self`: Weighted correlation with itself.
        * `Rank.self`: Rank of the self-match among all PRS samples.
* `[out_prefix].match.all.tsv.gz` : Detailed output including the top 5 matches for each phenotype sample.
    * Columns: `ID.Pheno`, `N.Traits`, `MatchStatus`, `ID.self`, `Z.self`, `COR.self`, `Rank.self`, followed by `ID.1st`, `Z.1st`, `COR.1st` through `ID.5th`, `Z.5th`, `COR.5th` for the five best-matching PRS samples.
    * `MatchStatus` in this file is one of `SELF_BEST` (the mapped PRS sample ranks first), `SELF_LENIENT` (it does not rank first but its Z-score exceeds `--z-threshold`), `SINGLE_NEW_BEST` (the best match leads the second by more than `--z-diff`), `MULTI_NEW_BEST` (a clear gap appears further down the top five), `UNCLEAR`, `NO_SELF` (no mapped PRS sample and no clear new best), or `NO_OBS_TRAITS`. Phenotype samples without a mapped PRS sample (`ID.self` = `NA`) receive `NO_SELF`, or `SINGLE_NEW_BEST` / `MULTI_NEW_BEST` when the top matches are clearly separated, or `NO_OBS_TRAITS`.

## Full Usage 

The full usage of `qpgentools match-prs-pheno` can be viewed with the `--help` option:

```
$ ./qpgentools match-prs-pheno --help
[bin/qpgentools match-prs-pheno] -- Match PRS and phenotype matrices based on overlapping samples and compute weights for each phenotype

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input Options ==
   --prs                  [STR: ]             : Input PRS file
   --pheno                [STR: ]             : Input phenotype matrix
   --cov                  [STR: ]             : Input covariate matrix (optional)
   --sample-tsv           [STR: ]             : TSV file containing input sample IDs PRS and phenotype files in [PRS_SAMPLE_ID] [PHENO_SAMPLE_ID] format. If they use the same IDs, use only a single column if subsetting samples are needed
   --trait-tsv            [STR: ]             : TSV file containing input trait IDs PRS and phenotype files in [PRS_TRAIT_ID] [PHENO_TRAIT_ID] format. If they use the same IDs, use only a single column if subsetting traits are needed
   --weights              [STR: ]             : Input weights file for each phenotype in [PHENO_ID] [WEIGHT] format
   --prs-format           [STR: regenie]      : Format of the PRS file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --pheno-format         [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --cov-format           [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --cov-names            [STR: ]             : Comma-separated covariates (column names in --cov) to use, as REGENIE's --covarColList; ranges such as PC{1:10} are expanded. Covariates in --cat-cov-names are always used. Default: all columns
   --cat-cov-names        [STR: ]             : Comma-separated covariates (column names in --cov) to treat as categorical, as REGENIE's --catCovarList. Their values may be any strings; each covariate with L levels among the analyzed samples is replaced by L-1 indicator columns (the first level in sorted order is the reference)
   --max-cat-levels       [INT: 10]           : Maximum number of levels of a categorical covariate, as REGENIE's --maxCatLevels (default: 10)
   --missing-str          [STR: NA]           : Comma-separated strings representing missing values in the phenotype and covariate matrices (default: 'NA'). The PRS matrix must not contain missing values

== Output options ==
   --out                  [STR: ]             : Output prefix

== Analysis options ==
   --rint                 [FLG: OFF]          : Perform rank-based inverse normal transformation after covariate adjustment (default: false)
   --mahalanobis          [FLG: OFF]          : Use Mahalanobis distance for matching (default: false)
   --no-norm              [FLG: OFF]          : Without --mahalanobis, divide the weighted inner product by the total absolute weight of the observed traits instead of by the weighted norms of the two profiles (legacy score; default: false)
   --exact-norm           [FLG: OFF]          : When phenotypes have missing values, compute each PRS norm over the traits observed for the phenotyped sample instead of over all traits. One extra matrix product without --mahalanobis or with --lambda 1; one pass per missingness pattern otherwise (default: false)
   --weight-prs-mh        [FLT: 0.50]         : Weight for PRS distance when combining with weighted correlation (default: 0.5)
   --lambda               [FLT: 0.00]         : Regularization parameter for Mahalanobis distance, between 0 and 1 (default: 0.0)
   --min-weight           [FLT: 0.00]         : Minimum weight (in r) per trait to set to zero (default: 0.0)
   --z-threshold          [FLT: 1.96]         : Z-score threshold for lenient matching (default: 1.96)
   --z-diff               [FLT: 2.00]         : Z-score difference to declare a clear match (default: 2.0)
   --threads              [INT: 1]            : Number of threads to use (default: 1)

== Imputation options ==
   --cov-impute-mean      [FLG: OFF]          : Mean-impute missing covariate values instead of dropping samples with any missing covariate (default: false)
   --missing-as-mean      [FLG: OFF]          : Impute missing phenotype values with the mean of observed values for the trait, then treat them as observed (default: false, missing values are ignored)
   --missing-as-half-min  [FLG: OFF]          : Impute missing phenotype values with half of the minimum observed value for the trait, the usual convention for measurements below a detection limit, then treat them as observed (default: false, missing values are ignored)

== Auto-lambda options (with --mahalanobis) ==
   --auto-lambda          [FLG: OFF]          : Choose the shrinkage parameter lambda automatically by maximizing the separation of the mapped self matches (default: false)
   --auto-lambda-metric   [STR: mean-log-softmax] : Criterion maximized by --auto-lambda (default: 'mean-log-softmax'). Options: 'mean-log-softmax', 'mean-z', 'mrr'
   --auto-lambda-min-self [INT: 100]          : Minimum number of phenotyped samples with a mapped PRS sample and observed traits required for tuning; fewer is an error (default: 100)
   --auto-lambda-min-best [INT: 100]          : Minimum number of mapped samples whose own PRS ranks first at the tuned lambda; fewer is an error, indicating a mismatched sample mapping (default: 100; 0 disables)
   --space-delim          [FLG: OFF]          : Input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) are delimited by single spaces instead of tabs (default: tab-delimited)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```