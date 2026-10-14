# qpgentools region-assoc

## Summary 

`qpgentools region-assoc` loads **every variant in one genomic region** at once and tests it
against one or more traits. Unlike [`rect-assoc`](rect_assoc.md), which takes an explicit variant
list and streams it in chunks, `region-assoc` keeps the whole region in memory, which makes it
the entry point for region-based methods that need the full local LD structure --- in particular
**SuSiE fine-mapping** (`--susie`).

A typical running example command is given below:

```bash
qpgentools region-assoc --pgen-list [list] --pheno [pheno] --cov [cov] \
    --region chr1:1000000-2000000 --traits ENSG00000187634 --out [out_prefix]
```

To export the region's summary data for fine-mapping in R with
[susieR](https://github.com/stephenslab/susieR) instead of (or in addition to) the built-in
fine-mapper:

```bash
qpgentools region-assoc --pgen-list [list] --pheno [pheno] --cov [cov] \
    --region chr1:1000000-2000000 --traits ENSG00000187634 --out [out_prefix] \
    --out-suff --out-rss
```

To additionally fine-map the region with SuSiE:

```bash
qpgentools region-assoc --pgen-list [list] --pheno [pheno] --cov [cov] \
    --region chr1:1000000-2000000 --traits ENSG00000187634 --out [out_prefix] \
    --susie --susie-L 10 --susie-coverage 0.95 --output-lbf
```

## Required options

* Input genotype files : Either `--pgen-list` or `--pgen`, `--psam`, `--pivar` options are required. See [Genotype file formats](../formats/genotypes.md) for more details.
* `--pheno` : Input phenotype matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details.
* `--region` : Region to be analyzed, in `CHROM:BEG-END` format (e.g. `chr12:6534517-6538371`). The chromosome name must match the names used in the indexed pvar file.
* `--traits` or `--traitf` : Traits to be tested. `--traits` takes a comma-separated list of trait IDs; `--traitf` takes a file with one trait ID per line. Exactly one of the two must be provided.
* `--out` : Output prefix. Output file names are formed by appending the suffixes described in [Expected Output](#expected-output).

## Additional Options

* `--cov-names` : Comma-separated covariates (column names in `--cov`) to use, as REGENIE's `--covarColList`. Ranges such as `PC{1:10}` are expanded. Covariates in `--cat-cov-names` are always used. Other columns are ignored and not parsed, so they may hold text. Default: all columns. See [Categorical and selected covariates](../formats/phenotypes.md#categorical-and-selected-covariates).
* `--cat-cov-names` : Comma-separated covariates to treat as categorical, as REGENIE's `--catCovarList`. Their values may be any strings. A covariate with L levels among the analyzed samples is replaced by L-1 indicator columns.
* `--max-cat-levels` : Maximum number of levels of a categorical covariate, as REGENIE's `--maxCatLevels` (default: 10). More levels stop the run.
* `--space-delim` : Read input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) as delimited by single spaces instead of tabs (default: tab-delimited). See [Phenotype file formats](../formats/phenotypes.md#delimiters).
### Input and sample handling

* `--cov` : Input covariate matrix in Regenie or TSV format. When provided, **both** the phenotypes and the genotypes are residualized against the covariates (see [Covariate handling](#covariate-handling)).
* `--sample` : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files.
* `--loco` : REGENIE step-1 leave-one-chromosome-out (LOCO) predictions of a single trait (a `.loco` file). The prediction row for the chromosome of `--region` is subtracted from the phenotype (see [LOCO adjustment](#loco-adjustment)). A `*_pred.list` file is also accepted here for backward compatibility.
* `--pred` : REGENIE step-1 `*_pred.list` file for one or more traits, with one `[pheno_id] [path_to_loco]` line per trait. Every tested trait must be listed. It applies the same LOCO adjustment as `--loco`, and the two options cannot be combined.
* `--pheno-format` : Format of the phenotype file (default: `regenie`). Options: `regenie`, `tensorqtl`, `tsv-sample-col`, `tsv-sample-row`.
* `--cov-format` : Format of the covariate file (default: `regenie`). Options: `regenie`, `tsv-sample-col`, `tsv-sample-row`.
* `--missing-str` : Comma-separated strings marking missing phenotype or covariate values (default: `NA`). All tested traits must share one set of samples, so missing phenotype values are allowed only when a single trait is tested. In that case, samples with a missing value are excluded before sample matching, so covariate adjustment, genotype residualization, RINT, LOCO, and the exported statistics all use the observed samples, and `N`/`n` report that count. With several tested traits, a missing value in any of them stops the run. Missing covariate values are not supported.
* `--colname-pheno-sample` / `--colname-geno-sample` : When `--sample` is provided, the column names holding the phenotype-side and genotype-side sample IDs (defaults: `pheno` and `geno`).
* `--icol-pivar-idx` : 1-based column index for the variant index in the indexed pvar file (default: 9).

### Variant filters

* `--min-af` / `--max-af` : Minimum / maximum ALT allele frequency for a variant to be included (defaults: 0.0 and 1.0).
* `--min-ac` / `--max-ac` : Minimum / maximum ALT allele count for a variant to be included (defaults: 0.0 and 1e9).
* `--max-allowed-vars` : Maximum number of variants that may be held in memory for the region (default: 1000000). The run aborts if the region contains at least this many variants, since the whole region must fit in memory.
* `--jump-thres-bp` : Jump threshold in base pairs for the variant index (default: 1000000).

Monomorphic variants (AC = 0 or AC = AN) in the analyzed samples are always dropped, regardless of the filters above.

### Transformations

* `--rint-before-adj` : Perform rank-based inverse normal transformation before covariate adjustment (default: false).
* `--rint-after-adj` : Perform rank-based inverse normal transformation after covariate adjustment (default: false).

RINT uses Blom offsets, \(\Phi^{-1}\big((r - 3/8)/(n + 1/4)\big)\) with average ranks \(r\) for ties, as REGENIE's `--apply-rint` and `--apply-rerint` do.

### Test statistic

* `--score` : Use the REGENIE-style score test instead of the default Wald t-test (default: false). The association output then has `ZSTAT` columns instead of `TSTAT`.

The two tests give the same `BETA = x'y / x'x` but differ in how the residual variance is
estimated:

* **Wald t-test (default).** The residual variance is estimated after fitting the variant,
  `SE = sqrt((y'y - (x'y)^2 / x'x) / (N - 2 - k) / x'x)`, and `LOG10P` comes from a
  t-distribution with `df = N - 2 - k` degrees of freedom. `k` is the number of linearly
  independent covariates (the numerical rank of the intercept and covariates, minus 1; `k = 0`
  without `--cov`). With `--cov`, `x` and `y` are covariate-residualized, so this is exactly the
  OLS fit `y ~ 1 + covariates + genotype`. When the phenotype is transformed after the covariate
  adjustment (`--rint-after-adj`, or LOCO), it is no longer exactly orthogonal to the covariates,
  and `y'y` is replaced by its residual sum of squares on the intercept and covariates, which
  keeps the equivalence with that OLS fit.
* **Score test (`--score`).** The residual variance is estimated under the null model,
  `sigma0^2 = y'y / (N - 1 - k)`, where `1 + k` counts the intercept and the `k` linearly
  independent covariates. Then
  `SE = sqrt(sigma0^2 / x'x)`, `ZSTAT = BETA / SE`, and `LOG10P` comes from the standard normal.
  This matches REGENIE step 2 for quantitative traits.

For a variant that explains a fraction `R^2` of the phenotype variance, the Wald `SE` is smaller
than the score `SE` by about `sqrt(1 - R^2)`. The difference is negligible for most variants but
visible at strong hits.

### Output file names

* `--assoc-suffix` : Suffix for the marginal association output (default: `.assoc.tsv.gz`).
* `--susie-cs-suffix` : Suffix for the SuSiE credible-set output (default: `.susie.cs.tsv.gz`).
* `--susie-lbf-suffix` : Suffix for the SuSiE per-variant log Bayes factor output (default: `.susie.lbf.tsv.gz`).
* `--suff-suffix` : Suffix for the sufficient-statistics export written by `--out-suff` (default: `.suff.bin.gz`).
* `--rss-suffix` : Suffix for the RSS summary-statistics export written by `--out-rss` (default: `.rss.bin.gz`).

Suffixes ending in `.gz` produce bgzipped output; any other suffix produces plain text. The two
binary exports are always bgzipped regardless of suffix.

### Summary-data export for susieR

* `--out-suff` : Write the **sufficient statistics** `X'X`, `X'y`, `y'y` and `n` to `[out_prefix].suff.bin.gz` (default: off).
* `--out-rss` : Write the **RSS summary statistics** `z`, LD correlation matrix `R`, `n`, `bhat`, `shat` and `var_y` to `[out_prefix].rss.bin.gz` (default: off).

Both are computed from the same covariate-residualized, mean-centered genotype matrix `X` and
phenotype matrix `Y` used for the marginal association and for `--susie`, so fitting them with
`susieR` reproduces the built-in fine-mapping. See
the corresponding section under [Expected Output](#expected-output) for the file layout and how
to load them in R. Both files contain a `p x p` matrix of doubles, so their size grows as
`8 * p^2` bytes before compression: about 80 MB for 3,000 variants and 8 GB for 30,000. Restrict
the region (or the allele-frequency filters) accordingly.

## SuSiE fine-mapping

`region-assoc` implements SuSiE ("Sum of Single Effects") fine-mapping via the IBSS algorithm,
reimplemented in Eigen from Wang, Sarkar, Carbonetto & Stephens (2020), *JRSS-B* 82(5):1273-1300.
Field names in the output mirror `susieR`'s `susie()` fit object (`alpha`, `mu`, `mu2`,
`lbf_variable`, `lbf`, `pip`, `V`, `sigma2`), so results are directly comparable with `susieR`.

Fine-mapping runs on exactly the same data as the marginal association reported in the
`.assoc.tsv.gz` file: the same genotype matrix, the same covariate-adjusted phenotypes, and the
same set of samples. This is why the credible-set output can carry the marginal `BETA`, `SE`, and
`LOG10P` alongside the posterior quantities.

**Multiple traits.** Because all tested traits share the same samples, everything that does not
depend on the trait is computed once per region and reused by every trait:

* the centered (and standardized) genotype matrix;
* `X'X`, which is also reused by `--out-suff`/`--out-rss` and for credible-set purity;
* the eigendecomposition used by `inf`/`ash`;
* `X'Y` for all traits, from a single matrix product that also feeds the marginal test.

Each trait then only runs its own IBSS iterations. For example, at n = 12,000 and p = 3,000,
each additional trait costs about 0.1 s with either `none` or `inf`, whereas `inf` alone takes
about 3 s for the first trait.

### SuSiE options

* `--susie` : Run SuSiE fine-mapping for each tested trait in the region (default: off). Without this flag only the marginal association file is written.
* `--susie-L` : Maximum number of causal single effects (default: 10). Effectively capped at the number of variants in the region.
* `--susie-max-iter` : Maximum number of IBSS iterations (default: 100).
* `--susie-tol` : Convergence tolerance for the SuSiE objective (default: 1e-3). For `--unmappable-effects inf` and `ash`, which converge on PIPs (`inf`) or `alpha` (`ash`) rather than the ELBO, a default of 1e-4 is used instead, matching `susieR`.
* `--susie-coverage` : Target coverage of the credible sets (default: 0.95).
* `--susie-min-abs-corr` : Minimum purity, i.e. the minimum absolute pairwise correlation among the members of a credible set, required to report it (default: 0.5). Sets below this threshold are dropped, as are duplicate sets.
* `--susie-no-standardize` : Do not standardize genotype columns to unit variance before fitting. Columns are still mean-centered.
* `--output-lbf` : Also write the per-variant log Bayes factors (one column per single effect) to `[out_prefix].susie.lbf.tsv.gz`.

### Unmappable-effects models

`--unmappable-effects` selects how effects that no single variant can capture (e.g. untyped
causal variants, or polygenic background within the region) are modeled. It matches the
`--method` argument of `run_susie_v1.r`.

| Value | Model |
|---|---|
| `none` (default) | Standard SuSiE. |
| `inf` | **SuSiE-inf**: adds an infinitesimal effect, `theta_j ~ N(0, tau2)`, with `(sigma2, tau2)` estimated by method of moments and `theta` obtained as its BLUP. Convergence is assessed on PIPs. |
| `ash` | **SuSiE-ash**, following Algorithm 1 of the SuSiE 2.0 manuscript (McCreight et al., bioRxiv 2025.11.25.690514). `theta_j` gets a scale-mixture-of-normals prior, `sum_k pi_k * N(0, sa2_k * sigma2)`. Each iteration runs the single-effect regressions under `Omega = (tau2 XX' + sigma2 I)^-1`, as in `inf`, then fits Mr.ASH to convergence on `y - X b_bar` (see below). Convergence is assessed on `alpha`. |

For `inf` and `ash`, `tau2` is reported in the log and a per-variant `theta` column is added to
the LBF output.

**SuSiE-ash details.** After the single-effect regressions, each iteration:

1. Estimates provisional `(sigma2, tau2_0)` by the same method of moments as `inf`. If `tau2_0`
   is 0, it sets `theta = 0` and `tau2 = 0` and skips the remaining steps for that iteration.
2. Builds the variance grid from `t = tau2_0 / sigma2`: a point mass at 0, 10 log-spaced points on
   `[t/100, t/10]`, 6 on `[t/10, 3t]`, and 3 on `[3t, V_min/sigma2]`, where `V_min` is the
   smallest prior variance among active single effects. Points within 10% of the previous one are
   dropped.
3. Fits Mr.ASH to convergence on `y - X b_bar`, starting from the previous `theta` and a uniform
   `pi`. It updates `theta` by coordinate ascent, `pi_k` as the mean responsibility, and `sigma2`
   with the Mr.ASH update (manuscript Supplementary Notes S3).
4. Sets `tau2 = sigma2 * sum_k pi_k sa2_k` and refreshes `Omega`.

The grid band sizes and the zero-`tau2_0` shortcut are not specified in the manuscript. They follow
the `susieR` code from the time of the preprint (`initialize_mrash()`, December 2025). The
implementation reproduces an independent R version of Algorithm 1 that uses `susieR`'s `mr.ash`
to about 1e-11 in PIPs and `theta`.

Current `susieR` releases differ from the manuscript. Their ash mode adds an LD-based masking
layer, subtracts only "confident" effects before Mr.ASH, uses a fixed grid, and runs its
single-effect regressions without `Omega`. In addition, the individual-level `susie()` path
passes the unscaled, uncentered `X` to `mr.ash` together with standardized-scale effects, so it
disagrees with `susie_ss()` unless `X` is standardized beforehand. `region-assoc` follows the
manuscript, so results are not expected to match current `susieR` ash output.

`inf` and `ash` start from a thin eigendecomposition of the standardized genotype matrix,
taken on the smaller of `X'X` (p x p) or `XX'` (n x n), so the cost is
`O(n p min(n,p) + min(n,p)^3)` like `susieR`'s `svd(X)`. When CMake finds a system LAPACK/BLAS
(Accelerate on macOS; OpenBLAS/MKL on Linux, selectable with `-DBLA_VENDOR=...`), this step uses
the multithreaded library. Otherwise it falls back to Eigen's single-threaded routines, which
are several times slower for regions with thousands of variants.

* `--threads N` : Number of BLAS/LAPACK threads (default `0` = the library's own default, usually
  all cores). Supported for Accelerate, OpenBLAS, MKL, BLIS and FlexiBLAS. On a cluster, match it
  to the cores you requested, e.g. `--threads $SLURM_CPUS_PER_TASK`.

The backend actually loaded is detected at run time. It is shown in `region-assoc --help` (under
`--threads`) and logged as `SuSiE linear algebra: <backend> (threads: N)` when `inf`/`ash` runs.

As in `susieR`, single effects whose estimated prior variance is ~0 (`V <= 1e-9`) are treated
as inactive: they contribute neither to the reported PIPs nor to credible sets.

Two options exist for diagnostics, and are ignored unless `--unmappable-effects ash` is given:

* `--ash-fix-sa2` : Comma-separated prior-variance grid, in units of `sigma2` with a first entry of 0, held fixed instead of the data-driven grid.
* `--ash-fix-pi` : Comma-separated mixture weights held fixed (no `pi` update). Requires `--ash-fix-sa2` of the same length. With `pi = 1,0,...,0`, `theta` stays 0; because `sigma2` still comes from Mr.ASH rather than SuSiE's expected residual sum of squares, this is close to, but not identical to, `--unmappable-effects none`.

### Variants excluded from fine-mapping

Because a single effect shares one residual and one softmax across all variants, one bad column
would contaminate the whole fit rather than just its own coefficient. `region-assoc` therefore
zeroes out (and reports in the log) any variant column that is non-finite, monomorphic in the
analyzed samples, or fully explained by the covariates after residualization. Such variants get
zero prior weight, so their `alpha` and `pip` are 0, and the remaining variants are fit normally.
They still appear in the marginal association output.

## Covariate handling

When `--cov` is supplied, the phenotype matrix is adjusted for the covariates by linear
regression, and the genotype matrix is residualized against the *same* covariates
(Frisch-Waugh-Lovell). Residualizing only the phenotype would leave genotype variance that is
collinear with the covariates (e.g. genotype PCs) in the design, inflating `x'x`, attenuating the
correlation with the phenotype, and thus shrinking `BETA` and the test statistic. Residualizing
both sides makes the reported effect the *partial* effect of the variant, and matches the usual
`susieR` workflow where both `X` and `y` are residualized before fitting.

!!! note
    Without `--score`, `SE`, `TSTAT` and `LOG10P` use `df = N - 2 - k`, where `k` is the number of
    linearly independent covariates (see [Test statistic](#test-statistic)); the value is printed
    in the log. A redundant covariate (constant in the analyzed samples, or a linear combination of
    others) is adjusted for once and costs no degree of freedom. The numerical rank is computed on the intercept and covariates scaled to unit norm, counting the pivots of a column-pivoted QR above `1e-7` times the largest one, as R's `lm()` does. A covariate that is *nearly* a linear combination of others is therefore treated as redundant too. A variant whose genotype is (numerically) a linear combination of the covariates (residual variance below `1e-8` of its total, e.g. when the variant itself is a covariate) cannot be tested: its statistics are reported as NaN, and SuSiE
    gives it no prior weight. Versions before this change used
    `df = N - 2` with `--cov`, which made `SE` too small by a factor of
    `sqrt((N - 2 - k) / (N - 2))` (about 5% for 66 covariates and 729 samples), so `--cov` results
    differ from those versions in `SE`, `TSTAT` and `LOG10P` (and in the `--out-rss` `shat` and
    `z`), but not in `BETA`. `--score` results are unchanged for full-rank covariates.

!!! warning
    Covariate adjustment and the rank-based inverse normal transformations are currently supported
    only for phenotype and covariate matrices **without missing values**; the run aborts otherwise.

## LOCO adjustment

With `--loco` (one trait) or `--pred` (a `*_pred.list` for one or more traits), the phenotype is adjusted for the polygenic background estimated by REGENIE
step 1, following REGENIE step 2 for quantitative traits:

1. RINT the raw phenotype if `--rint-before-adj` is set.
2. Residualize the phenotype against the covariates (or just center it without `--cov`).
3. Scale the residual to unit SD, with `SD = ||r|| / sqrt(N - 1 - k)`, where `1 + k` counts the
   intercept and the `k` linearly independent covariates (REGENIE's `scale_Y`).
4. Subtract the LOCO prediction for the chromosome of `--region`.
5. RINT the result if `--rint-after-adj` is set. Otherwise, multiply it back by the SD from
   step 3.

Because of step 5, `BETA`, `SE`, the exported `--out-suff`/`--out-rss` statistics, and the SuSiE
effect sizes stay in units of the covariate-adjusted phenotype, as REGENIE reports them.
`TSTAT` and `LOG10P` do not depend on this rescaling. With `--rint-after-adj`, effects are on
the RINT scale, as with REGENIE's `--apply-rerint`.

A `.loco` file has a header `FID_IID <FID>_<IID> ...` and one row per chromosome (`1`-`22`,
and `23` for X). The region's chromosome is matched after dropping a `chr` prefix, with X
mapped to `23`. Sample IDs are matched to the `FID_IID` columns exactly, then as `<ID>_<ID>`,
and then as a unique `<FID>_<ID>` suffix. Samples without a non-missing LOCO prediction for
every tested trait are excluded, and a warning reports how many. In a `*_pred.list` file,
relative paths that do not exist as given are resolved against the list file's directory.

```bash
qpgentools region-assoc --pgen-list [list] --pheno [pheno] --cov [cov] \
    --traits T1,T2 --region chr2:1000000-2000000 \
    --pred [regenie_step1]_pred.list --susie --out [out_prefix]
```

With a single trait that has missing values, the samples with a missing phenotype are excluded
first. LOCO scaling (step 3) then uses only the observed samples, matching a run in which those
samples were removed from the phenotype file.

## Expected Output

### `[out_prefix].assoc.tsv.gz` --- marginal association

One row per variant in the region that passed the filters, in a wide format where the statistics
for each trait are appended horizontally.

* `#CHROM` : Chromosome
* `GENPOS` : Base position (1-based)
* `ID` : Variant ID in `[CHROM]:[GENPOS]:[REF]:[ALT]` format
* `ALLELE0` : Reference allele
* `ALLELE1` : Alternative allele
* `A1FREQ` : Frequency of the alternative allele (AC / AN)
* `N` : Number of samples with a called genotype (`N_RR + N_RA + N_AA`)
* `N_RR` : Count of Reference/Reference genotypes
* `N_RA` : Count of Reference/Alternative genotypes
* `N_AA` : Count of Alternative/Alternative genotypes
* `TEST` : Test type (always `Linear`)

Followed, for each trait `[TRAIT]`, by:

* `BETA.[TRAIT]` : Effect size
* `SE.[TRAIT]` : Standard error
* `TSTAT.[TRAIT]` : Wald t-statistic, or `ZSTAT.[TRAIT]` : score z-statistic with `--score`
* `LOG10P.[TRAIT]` : -log10 p-value of the selected test

### `[out_prefix].susie.cs.tsv.gz` --- credible sets

Written only with `--susie`. One row per (trait, credible set, variant) triple, i.e. only variants
that belong to a reported credible set appear here.

* `#trait` : Trait ID
* `variant` : Variant ID in `[CHROM]:[POS]:[REF]:[ALT]` format
* `pip` : Marginal posterior inclusion probability, `1 - prod_l (1 - alpha_lj)`
* `cs_id` : 1-based credible set index within this trait
* `alpha` : Posterior inclusion probability of this variant *within this single effect*
* `region` : Analyzed region, formatted as `CHROM_BEG_END`
* `n_region_vars` : Number of variants analyzed in this region (identical on every row of the file)
* `cs_size` : Number of variants in this credible set
* `cs_lbf` : Log Bayes factor (natural log) of the single effect defining this credible set, i.e. the prior-weighted average BF over all variants in the region. Constant across all rows of a credible set
* `var_lbf` : Log Bayes factor (natural log) of *this variant* under that single effect. Related to `alpha` by `var_lbf - cs_lbf = log(alpha) - log(prior)`
* `mu` : Posterior mean effect, conditional on inclusion
* `mu2` : Posterior second moment of the effect, conditional on inclusion
* `af` : Alternative allele frequency
* `n` : Sample size used in the marginal test
* `beta` : Marginal effect size (matches `BETA.[TRAIT]` in the association file)
* `se` : Marginal standard error (matches `SE.[TRAIT]`, so it follows `--score`)
* `log10p` : Marginal -log10 p-value (matches `LOG10P.[TRAIT]`, so it follows `--score`)

A trait with no credible set surviving the coverage and purity thresholds contributes no rows.

### `[out_prefix].susie.lbf.tsv.gz` --- per-variant log Bayes factors

Written only with `--susie --output-lbf`. One row per (trait, variant) pair for **every** variant
in the region, whether or not it belongs to a credible set.

* `#trait` : Trait ID
* `variant` : Variant ID
* `region` : Analyzed region, formatted as `CHROM_BEG_END`
* `n_region_vars` : Number of variants analyzed in this region (identical on every row of the file)
* `af` : Alternative allele frequency
* `pip` : Marginal posterior inclusion probability
* `theta` : Posterior mean unmappable effect on the standardized-genotype scale (equivalent to `susieR`'s `fit$theta`). **Present only with `--unmappable-effects inf` or `ash`.**
* `lbf.L1` ... `lbf.L[L]` : Log Bayes factor of this variant under each single effect, where `L` is `min(--susie-L, number of variants)`

### `[out_prefix].suff.bin.gz` and `[out_prefix].rss.bin.gz` --- summary-data export for susieR

Written only with `--out-suff` / `--out-rss`. These are gzip (BGZF) compressed binary files
that need no extra R packages to read: `gzfile()` + `readBin()` is enough. A ready-made reader is
provided in [`scripts/qpgen_susie_io.R`](https://github.com/hyunminkang/qpgen/blob/main/scripts/qpgen_susie_io.R).

**Quantities.** `X` is the `n x p` matrix of ALT-allele dosages, mean-imputed for missing
genotypes, mean-centered, and residualized against the covariates when `--cov` is given. `y` is
the covariate-residualized, mean-centered phenotype (one column per trait, `K` traits).
`n_cov` (stored in the file header) is `k`, the number of linearly independent covariates
*excluding* the intercept (0 without `--cov`), so the Wald df is `n - 2 - n_cov`.

After `--rint-after-adj` or LOCO, `y` is the transformed phenotype, which is no longer exactly
orthogonal to the covariates. `yty` and `var_y` are computed from its `y'y`, while the Wald test
in the association file uses its residual sum of squares on the intercept and covariates (see
[Test statistic](#test-statistic)). `BETA` and `X'y` are not affected, because `X` is
residualized. Only the residual variance differs: it is slightly larger in SuSiE and in fits
rebuilt from the exports than in the `.assoc` `SE`. The built-in `--susie` uses the same `y'y`,
so it agrees with the exports exactly. The log reports the gap as the largest share of `y'y` the
covariates explain. On the MAGE example (66 covariates, 729 samples, `--rint-after-adj`), the
share was at most 0.15%. Refitting `susie_ss()` with the residual sum of squares instead of
`y'y` then changed PIPs by at most 0.001 and posterior effects by at most 0.3%, and left every
credible set unchanged.

| File | Contents | susieR entry point |
|---|---|---|
| `.suff.bin.gz` | `XtX = X'X` (`p x p`), `Xty = X'y` (`p x K`), `yty = y'y` (`K`), `n` | `susie_ss(XtX, Xty, yty, n)` in susieR >= 0.14 (`susie_suff_stat()` in CRAN susieR <= 0.12) |
| `.rss.bin.gz` | `R = cov2cor(X'X)` (`p x p`), `z = bhat/shat`, `bhat`, `shat` (`p x K`), `var_y = y'y/(n-1)` (`K`), `n` | `susie_rss(z = z, R = R, n = n)` or `susie_rss(bhat, shat, R, n, var_y)` |

The sufficient-statistics file is the primary export. `X'X`, `X'y`, `y'y` and `n` are the exact
individual-level likelihood, and `susie_ss()` on them reproduces the built-in `--susie` fit
(same `sigma2`, credible sets and PIPs). Use this file whenever you fit susieR yourself.

The RSS file exists for tools that take marginal statistics plus an LD matrix. `bhat`, `shat` and
`z` in it are the same values as `BETA`, `SE` and `TSTAT` (or `ZSTAT` with `--score`) in the association file, and `R` is the
in-sample LD of the same residualized data. With `--cov`, `shat` and `z` therefore use
`df = n - 2 - n_cov` (see [Test statistic](#test-statistic)), while `susie_rss()` assumes
`df = n - 2` when it rebuilds `X'X`, `X'y` and `y'y` from `bhat` and `shat`. The helper
`susie_from_qpgen()` converts `shat` to that SE before calling `susie_rss()`, using the test
recorded in the file header:

* Wald (default): `shat * sqrt((n - 2 - n_cov) / (n - 2))`
* score (`--score`): `shat * sqrt((n - 1 - n_cov - z^2) / (n - 2))`, with `z = bhat / shat`

Both conversions are exact, so the fit matches `susie_ss()` to about 1e-9 in PIP. After
`--rint-after-adj` or LOCO, the Wald conversion is approximate, with the same small gap as
above. Apply the same factor when calling `susie_rss()` directly. Files from versions before
format version 2 do not record the test. For those, `susie_from_qpgen()` warns and assumes Wald;
pass `test = "score"` for files written with `--score`. Without the conversion, the MAGE example
`--score` files gave PIPs off by up to 0.21 and different credible sets.

The same applies to `susie_rss(z = z, R = R, n = n)`. susieR converts `z` assuming
`df = n - 2`, but the exported `z` is a t statistic with `df = n - 2 - n_cov`. Pass
`n = n - n_cov` so that the conversion uses the right df. With `--score` files this is only
approximate; use `bhat` and `shat` through `susie_from_qpgen()` instead. The difference is small when
`n_cov` is small relative to `n`. It is not negligible for eQTL-style designs with many
covariates: with 66 covariates and 729 samples, it is about 5% in `z`. Note that `susie_rss()` is **not** equivalent to the
sufficient-statistics fit under its defaults: it fixes the residual variance
(`estimate_residual_variance = FALSE`), which on the test data changed `sigma2` from 0.937 to
1.285 and moved PIPs by up to 0.05. With `estimate_residual_variance = TRUE` and `bhat`, `shat`,
`var_y` and `n` supplied, `susie_rss()` rebuilds `X'X`, `X'y` and `y'y` internally and matches
`susie_ss()` to about 1e-9. Passing `z` alone fixes `sigma2 = 1` and gives the RSS-model fit.

**Binary layout** (little-endian; strings are NUL-terminated; matrices are column-major, as R's
`matrix()` expects):

```
char[8]   magic       "QPGNSUFF" (suff) or "QPGN_RSS" (rss)
int32     version     2 (version 1 has no "test" field)
int32     n           number of samples
int32     p           number of variants
int32     K           number of traits
int32     n_cov       covariate degrees of freedom residualized out: the number of linearly
                      independent covariates, excluding the intercept (0 without --cov)
int32     test        test of the .assoc statistics and the RSS bhat/shat/z:
                      0 = Wald t-test, 1 = score test (--score)
str       region      CHROM:BEG-END
str[K]    trait_ids
str[p]    variant_ids CHROM:POS:REF:ALT (ALT is the dosage allele)
int32[p]  pos         1-based position
double[p] af          ALT allele frequency
--- suff body ---              --- rss body ---
double[p*p] XtX                double[p*p] R
double[p*K] Xty                double[p*K] z
double[K]   yty                double[p*K] bhat
                               double[p*K] shat
                               double[K]   var_y
```

**Loading in R.** With the helper script:

```r
library(susieR)
source("scripts/qpgen_susie_io.R")

ss  <- read_qpgen_suff("out.suff.bin.gz")   # list(XtX, Xty, yty, n, variants, traits, af, ...)
fit <- susie_from_qpgen(ss, trait = "ENSG00000187634", L = 10)   # = susie_ss(XtX, Xty, yty, n, L = 10)
fit$sets$cs                                  # credible sets (indices into ss$variants)
head(sort(susie_get_pip(fit), decreasing = TRUE))

rss  <- read_qpgen_rss("out.rss.bin.gz")    # list(R, z, bhat, shat, var_y, n, ...)
fit2 <- susie_from_qpgen(rss, trait = 1, L = 10,   # = susie_rss(bhat, shat, R, n, var_y, ...)
                         estimate_residual_variance = TRUE)   # needed to match fit
```

`susie_from_qpgen()` dispatches on the file type: `susie_ss()` for the sufficient statistics
(`susie_suff_stat()` on old susieR), `susie_rss()` for the RSS file. Any further arguments
(`L`, `coverage`, `min_abs_corr`, `estimate_residual_variance`, ...) pass through. Use
`fit$sets` for credible sets, or `susie_get_cs(fit, Xcorr = cov2cor(ss$XtX))`; calling
`susie_get_cs()` without `Xcorr` skips the purity filter.

Without the helper, the RSS file reads as:

```r
con <- gzfile("out.rss.bin.gz", "rb")
stopifnot(readChar(con, 8, useBytes = TRUE) == "QPGN_RSS")
version <- readBin(con, "integer", 1, size = 4, endian = "little")
d <- readBin(con, "integer", 4, size = 4, endian = "little"); n <- d[1]; p <- d[2]; K <- d[3]; n_cov <- d[4]
test <- if (version >= 2) readBin(con, "integer", 1, size = 4, endian = "little") else 0L  # 1 = --score
region   <- readBin(con, "character", 1)
traits   <- readBin(con, "character", K)
variants <- readBin(con, "character", p)
pos      <- readBin(con, "integer", p, size = 4, endian = "little")
af       <- readBin(con, "double", p)
R    <- matrix(readBin(con, "double", p * p), p, p, dimnames = list(variants, variants))
z    <- matrix(readBin(con, "double", p * K), p, K, dimnames = list(variants, traits))
bhat <- matrix(readBin(con, "double", p * K), p, K)
shat <- matrix(readBin(con, "double", p * K), p, K)
var_y <- readBin(con, "double", K)
close(con)
shat <- if (test == 0) shat * sqrt((n - 2 - n_cov) / (n - 2)) else      # Wald df n - 2 - n_cov -> n - 2
        shat * sqrt((n - 1 - n_cov - (bhat / shat)^2) / (n - 2))          # score SE -> Wald SE, df n - 2
fit <- susie_rss(bhat = bhat[, 1], shat = shat[, 1], R = R, n = n, var_y = var_y[1], L = 10)
```

The `.suff.bin.gz` file has the same header followed by `XtX` (`p*p` doubles), `Xty` (`p*K`)
and `yty` (`K`), to be passed as `susie_ss(XtX, Xty[, 1], yty[1], n, L = 10)`. Python users can read the same layout with `gzip.open()` and
`numpy.frombuffer(..., dtype="<f8").reshape((p, p), order="F")`.

!!! note
    The LD matrix `R` is the in-sample correlation of the *covariate-residualized* dosages, not
    of the raw genotypes. This is what makes `susie_rss` on these files equivalent to the
    individual-level fit; it is not a drop-in reference LD panel for other GWAS.

## Full Usage 

The full usage of `qpgentools region-assoc` can be viewed with the `--help` option:

```
$ ./qpgentools region-assoc --help
[bin/qpgentools region-assoc] -- Perform association analysis for a specific region and multiple phenotypes

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
   --sample               [STR: ]             : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
   --cov                  [STR: ]             : Input covariate matrix (optional)
   --loco                 [STR: ]             : REGENIE step-1 LOCO predictions of a single trait (.loco file). The covariate-adjusted phenotype is scaled to unit SD, the prediction for the region's chromosome is subtracted, and (unless --rint-after-adj) the result is rescaled to the original SD so BETA/SE stay in phenotype units. A *_pred.list file is also accepted here for backward compatibility
   --pred                 [STR: ]             : REGENIE step-1 *_pred.list file for one or more traits: one '[pheno_id] [path_to_loco]' line per trait (relative paths are resolved against the list's directory). Each tested trait must be listed. Same LOCO adjustment as --loco
   --pheno-format         [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --missing-str          [STR: NA]           : Comma-separated strings marking missing phenotype/covariate values (default: 'NA'). Missing phenotype values are allowed only when a single trait is tested; its samples with a missing value are excluded
   --cov-format           [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --cov-names            [STR: ]             : Comma-separated covariates (column names in --cov) to use, as REGENIE's --covarColList; ranges such as PC{1:10} are expanded. Covariates in --cat-cov-names are always used. Default: all columns
   --cat-cov-names        [STR: ]             : Comma-separated covariates (column names in --cov) to treat as categorical, as REGENIE's --catCovarList. Their values may be any strings; each covariate with L levels among the analyzed samples is replaced by L-1 indicator columns (the first level in sorted order is the reference)
   --max-cat-levels       [INT: 10]           : Maximum number of levels of a categorical covariate, as REGENIE's --maxCatLevels (default: 10)
   --traits               [STR: ]             : Trait IDs (comma-separated) to be tested (required)
   --traitf               [STR: ]             : Input file containing trait IDs to be tested (one per line; the first column is used)
   --region               [STR: ]             : Region string in the format CHROM:BEG-END (required)
   --min-af               [FLT: 0.00]         : Minimum allele frequency for variants to be tested (default: 0.0)
   --max-af               [FLT: 1.00]         : Maximum allele frequency for variants to be tested (default: 1.0)
   --min-ac               [FLT: 0.00]         : Minimum allele count for variants to be tested (default: 0.0)
   --max-ac               [FLT: 1000000000.00] : Maximum allele count for variants to be tested (default: 1e9)
   --threads              [INT: 0]            : Number of threads for the BLAS/LAPACK linear algebra of SuSiE-inf/ash (default: 0 = library default). Linked backend: Apple Accelerate

== Output options ==
   --out                  [STR: ]             : Output prefix
   --assoc-suffix         [STR: .assoc.tsv.gz] : Suffix for the association output file (default: '.assoc.tsv.gz')
   --susie-cs-suffix      [STR: .susie.cs.tsv.gz] : Suffix for the SuSiE credible set output file (default: '.susie.cs.tsv.gz')
   --susie-lbf-suffix     [STR: .susie.lbf.tsv.gz] : Suffix for the SuSiE log Bayes factor output file (default: '.susie.lbf.tsv.gz')
   --out-suff             [FLG: OFF]          : Write sufficient statistics (X'X, X'y, y'y, n) for susieR as a gzipped binary file [out]<suff-suffix>
   --out-rss              [FLG: OFF]          : Write RSS summary statistics (z, LD matrix R, n, bhat, shat, var_y) for susieR::susie_rss as a gzipped binary file [out]<rss-suffix>
   --suff-suffix          [STR: .suff.bin.gz] : Suffix for the sufficient statistics output file (default: '.suff.bin.gz')
   --rss-suffix           [STR: .rss.bin.gz]  : Suffix for the RSS summary statistics output file (default: '.rss.bin.gz')

== Auxiliary options ==
   --jump-thres-bp        [INT: 1000000]      : Jump threshold in base pairs for the variant index (default: 1000000)
   --max-allowed-vars     [INT: 1000000]      : Maximum number of allowed variants to store at once in memory (default: 1000000)
   --icol-pivar-idx       [INT: 9]            : 1-based column index for the variant ID in the pvar file (default: 9)
   --colname-pheno-sample [STR: pheno]        : When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')
   --colname-geno-sample  [STR: geno]         : When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'geno')
   --rint-before-adj      [FLG: OFF]          : Perform rank-based inverse normal transformation before covariate adjustment (default: false)
   --rint-after-adj       [FLG: OFF]          : Perform rank-based inverse normal transformation after covariate adjustment (default: false)
   --score                [FLG: OFF]          : Use the REGENIE-style score test (residual variance under the null with N - 1 - k df for k linearly independent covariates, normal p-value) and report ZSTAT instead of the Wald TSTAT

== SuSiE fine-mapping options ==
   --susie                [FLG: OFF]          : Run SuSiE fine-mapping for each tested phenotype in the region
   --susie-L              [INT: 10]           : Maximum number of causal single effects (default: 10)
   --susie-max-iter       [INT: 100]          : Maximum number of IBSS iterations (default: 100)
   --susie-coverage       [FLT: 0.95]         : Target coverage of credible sets (default: 0.95)
   --susie-min-abs-corr   [FLT: 0.50]         : Minimum absolute correlation (purity) required to report a credible set (default: 0.5)
   --susie-tol            [FLT: 1.0e-03]      : Convergence tolerance for the SuSiE objective (default: 1e-3)
   --susie-no-standardize [FLG: OFF]          : Do not standardize genotype columns to unit variance before SuSiE
   --unmappable-effects   [STR: none]         : Unmappable-effects model for SuSiE: 'none' (standard), 'inf' (SuSiE-inf, adds an infinitesimal effect), or 'ash' (SuSiE-ash as in the SuSiE 2.0 manuscript, Algorithm 1: adaptive-shrinkage background fit by Mr.ASH, SER under the implied Omega) (default: none)
   --ash-fix-pi           [STR: ]             : Diagnostics for --unmappable-effects ash: comma-separated mixture weights held fixed (no pi update). Requires --ash-fix-sa2 of the same length. pi=1,0,...,0 keeps theta=0 (sigma2 still comes from Mr.ASH, so this is close to, not identical to, --unmappable-effects none)
   --ash-fix-sa2          [STR: ]             : Diagnostics for --unmappable-effects ash: comma-separated prior-variance grid (in units of sigma2, first entry 0) held fixed instead of the data-driven grid rebuilt every iteration
   --output-lbf           [FLG: OFF]          : Also write per-variant log Bayes factors (one column per single effect) when running SuSiE
   --space-delim          [FLG: OFF]          : Input tables (phenotype/covariate matrices, pair files, sample/trait/variant lists) are delimited by single spaces instead of tabs (default: tab-delimited)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
