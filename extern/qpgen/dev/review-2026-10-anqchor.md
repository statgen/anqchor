# Review of qpgen for anqchor (2026-10): resolution log

Developer notes, not user documentation. This file records the review of qpgen (`dev`, from `938bf19`) made while building anqchor on the
MAGE chr21/22 example (2026-10-11). Item numbers follow the original review.

## Fixed

| # | Issue | Fix |
|---|---|---|
| 6 | `pair-assoc` without `--cov` did not center the phenotype or dosages, so `SE` was inflated (about 100-fold for a phenotype with mean 50) and dosage `BETA` was wrong. | `pair-assoc` now always fits `y ~ 1 + C + g`. Without `--cov`, this centers `y` and `g`. |
| 7 | With `--cov`, `pair-assoc` mean-imputed missing genotypes, so df counted all samples while `N` counted called genotypes. | Default is now complete-case: samples with a missing genotype (hard call or dosage) are dropped from that variant's test, and the covariates are refitted on the rest. `N = N_RR + N_RA + N_AA` and is the regression sample size. `--mean-impute-geno` keeps the old behaviour. |
| 8 | `pair-assoc` read `NA` as 0. | Added `--missing-str` (default `NA`). Missing phenotypes are dropped per trait, with the covariates refitted. Samples with a missing covariate are dropped from all tests. |
| 8+ | Every loader read non-numeric cells (`NA`, `XX`/`XY`) as 0 through `atof()`. The MAGE covariate file stores sex as `XX`/`XY`. | `PhenoMatrix::load_pheno_matrix()` now stops with an error on a value that is not a finite number and not a registered missing string. This affects all commands; `rect-assoc` now rejects `NA` instead of reading it as 0. |
| 9 | `pair-assoc` aborted when the pgen sample order was not lexicographic. | Overlapping IDs are put in genotype-file order, as `region-assoc` already did. |
| 9+ | `pair-assoc` ran on with 0 overlapping samples, for example with the wrong `--cov-format`. | Now an error. |
| 13 | `pair-assoc` computed a separate phenotype QR. | One covariate fit is shared by all traits without missing values. 20k traits × 1 variant runs in 18.1 s, against 18.9 s before, with identical output. |
| 14 | `pair-assoc` header started with `#TRAIT`. | Now `TRAIT`. |
| 15 | `pair-prs` progress lines did not add up. Monomorphic variants were never counted, and there was no final total. | Each variant is counted as scored, not found, or monomorphic, and a final summary line is printed. |
| 17 | `n_cov` meant two things in `region_assoc.md`. The `susie_export.h` header comment was stale. | The score test and LOCO text now use `1 + k`; `n_cov` always excludes the intercept. The header comment has been rewritten. |
| 4 | After `--rint-after-adj` or LOCO, `.suff`/`.rss` and `--susie` use `y'y`, while the `.assoc` Wald SE uses the residual sum of squares on `[1, C]`. | Documented, not changed (user decision). BETA and `X'y` are unaffected; only the residual variance differs. region-assoc now logs the gap. On MAGE (`--rint-after-adj`, 66 covariates) it was at most 0.15% of `y'y`: PIPs changed by at most 0.001, posterior effects by at most 0.3%, and the credible sets were identical. Only region-assoc has LOCO and the exports, so this is region-only. |
| 5 | `susie_from_qpgen()` applied the Wald `shat` rescale to `--score` files: PIPs were off by up to 0.21 and credible sets changed for 2 of 4 MAGE genes. The file did not record the test. | Export format version 2 adds a `test` header field. The R helper converts score SEs with `sqrt((n - 1 - n_cov - z^2) / (n - 2))`, which matches `susie_ss()` to about 1e-9. Version-1 files still read: the helper warns, assumes Wald, and accepts `test = "score"`. The docs also give the `n = n - n_cov` caveat for `susie_rss(z, R, n)`. |
| 5+ | The `.suff` X'X was not bit-exactly symmetric, so susieR warned. | The lower triangle is mirrored before writing. |
| 1 | The covariate rank used Eigen's unscaled default threshold. | `CovariateResidualizer` now computes the rank once, on unit-norm columns, with a relative tolerance of 1e-7 (as in R's `lm()`). The rank sets both df and the projection. pair-assoc's separate check was removed, so all three commands warn the same way. |
| 2 | Rank-deficient covariates gave unstable residuals. Exact constants and duplicates were already handled; the remaining failure was a near-duplicate (`c1` + 1e-9 noise), which shifted BETA by up to 73% with no warning. | Residuals are now always `(I - Q_r Q_r')v`. The old least-squares solve's huge, cancelling coefficients caused the error. Constant, duplicated and near-duplicate covariates now leave BETA unchanged, with a warning, in all three commands. |
| 3 | A variant explained by the covariates printed garbage in rect-/region-assoc (BETA -1.6e13). | `CovariateResidualizer::residualize_genotypes()` (shared by all three commands) zeroes such columns, so the tests report NaN and SuSiE excludes the variant. |
| FID | `load_psam`/`load_fam` keyed samples as `FID_IID` when FID differed from IID and was not `0`, but every command matches by IID. Subsetting then found no samples, and an empty subset made the reader load *all* samples, so genotype rows were misaligned with the phenotypes: **silently wrong results** (on synthetic data, a true signal went from LOG10P 4.1 to 0.05). It affected any command that subsets samples when the psam/fam has real family IDs. | Samples are keyed by IID in all three loaders, and duplicate IIDs are an error. An empty subset is an error. pair-, rect- and region-assoc check that the loaded genotype count equals the number of matched samples. With FID != IID, output is identical to the FID == IID run in all three commands. Note: `pgen2tsv`/`pair-prs` `--sample` lists now take IIDs, where `FID_IID` was needed before. |
| rect-MI | rect-assoc's mean imputation of missing genotypes was not explained. | Intentional (user decision). Documented in `rect_assoc.md` (Missing genotypes): `N` counts called genotypes while df uses all samples, and BETA is approximately unbiased under random missingness. |
| 10 | Input tables were split on any run of spaces or tabs, so a trait or sample ID containing a space was cut: an error in phenotype files, and silently truncated in pair files and `--sample` lists. | Input tables are tab-delimited by default; `--space-delim` reads single-space-delimited files instead, and the two are never mixed. This covers phenotype/covariate/PRS matrices, pair files and sample/trait/variant lists in pair-, rect- and region-assoc, pair-prs, match-prs-pheno and pgen2tsv. A header that splits into too few fields now gives an error naming the delimiter (it used to crash). psam/fam/pvar, `.loco` and `--pgen-list` files keep whitespace splitting. **Behaviour change:** space-delimited tables now need `--space-delim`. |
| 18 | No `--version` flag. | Version 0.2.0 (`project(qpgen VERSION 0.2.0)`). `qpgentools --version` prints `qpgentools 0.2.0 (git <commit>[-dirty])`, and every command logs it first. The commit comes from `cmake/qpgen_version.cmake`, which runs at every build. |
| obs | No warning when matching on few traits. | match-prs-pheno warns when fewer than 10 phenotypes have a non-zero weight, or when phenotyped samples have 1-9 observed traits (count and examples). Matching is unchanged; `N.Traits` is in the output. Per the user, sample counts per trait do not matter. |
| LOCO | `--loco`/`--pred` existed only in region-assoc. | Added to pair-assoc (per trait and chromosome, on the trait's samples) and rect-assoc (per chromosome; chunks split at chromosome changes, output in variant-list order). Both share region-assoc's arithmetic (`ind_assoc_input::apply_loco_chrom()`, public `loco_row_t`/`load_loco_row()`). For the same chromosome, results are identical to region-assoc, and pair-assoc complete-case output matches a NumPy reference. region-assoc output is unchanged. |
| 6b | rect-assoc without `--cov` did not center phenotypes, so `SE` was inflated (14.5x for a trait with mean 20; 0.08% on MAGE's inverse-normal expression). | The shared loader centers phenotypes when no covariates are given; rect-assoc now matches region-assoc. |
| cov | No REGENIE-style covariate selection or categorical covariates. Dead options `--icol-cov-id`, `--offset-pheno`, `--offset-cov` (and `--icol-pheno-id` outside region-assoc's `--traitf`). | Added `--cov-names` (REGENIE `--covarColList`, with `PC{1:10}` ranges; unselected columns are not parsed), `--cat-cov-names` (`--catCovarList`; L-1 indicators from levels counted on the analyzed samples) and `--max-cat-levels` (default 10) to pair-, rect- and region-assoc. Results equal hand-made indicator columns exactly, including with LOCO and SuSiE. match-prs-pheno takes the same three options; with `--cov-impute-mean`, indicators are imputed by level frequency. The dead options and `--icol-pheno-id` were removed (`--traitf` reads its first column). rect-assoc gained `--missing-str`, so `NA` is no longer read as a category level, and the shared loader stops clearly on missing phenotypes or covariates. |
| 19+ | `pair-assoc` cached the dosage buffer pointer once. With `--pgen-list`, each pgen has its own buffer, so dosage reads after the first pgen used a stale buffer. | The buffer is now fetched after every read. This is still untested on real dosage data; see item 19. |

Validation for 1-3: in every command, a constant, duplicated or near-duplicate covariate leaves
BETA unchanged, and a variant used as a covariate gives NaN. On MAGE, pair, rect and region
`.assoc` output and the SuSiE credible sets are identical at print precision. `.suff` matrices
differ by about 5e-15 of their scale.

Not changed: `bulk_adjust_for_covariates()` and `pheno_adj_cov_nxt_with_missing()` (used by
`match-prs-pheno`) still use their own least-squares solves.

Validation: `pair-assoc` output matches a NumPy OLS reference to print precision (about 5e-6).
This holds on MAGE (614 pairs, with and without `--cov`) and on a synthetic PLINK1 `.bed` set.
That set has missing genotypes, phenotypes and a covariate, a non-lexicographic sample order,
and an uncentered trait. It was run in pair and region mode with `--max-chunk-vars 7`, and
`--mean-impute-geno` was checked too. With `--cov` and no missing data, MAGE output is
byte-identical to `938bf19`. `rect-assoc` and `region-assoc` (`--out-suff`, `--susie`) are
byte-identical to `938bf19` on MAGE.

## Open

### Statistics and numerics

No open items. 4 and 5 are resolved; see the Fixed table.

### Missing data

No open items.

### Parsing and polish

11. **Allele matching** is exact `CHROM:POS:REF:ALT`, with no flips or normalisation. Not
    planned.
12. **`pair-assoc-v0`** is outdated (df `N - 2`, no genotype residualisation). Not planned;
    remove it when convenient.
16. **Region-pair chunking noise.** Different `--max-chunk-vars` values change floating-point
    summation order. The visible differences (e.g. 5e-17 vs 2e-17) are in values that are zero
    in exact arithmetic. Not a bug.
19. **Dosage pgen files are untested** for the new `pair-assoc` path: missing dosages flagged as
    -3 and complete-case analysis. The same goes for FWL and the df change. Test this when
    dosage data is available.
