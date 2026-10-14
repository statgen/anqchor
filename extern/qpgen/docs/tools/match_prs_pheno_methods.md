# match-prs-pheno: methods

This page describes how `qpgentools match-prs-pheno` works, from the input matrices to the match statuses in the output. It complements the [option reference](match_prs_pheno.md), which lists every parameter and the output columns. Section 1 gives the overall framework, section 2 the score with and without `--mahalanobis`, section 3 how the Z-scores turn into match calls, section 4 the automatic choice of the shrinkage parameter, and section 5 the treatment of missing values.

## 1. Framework

### The problem

Two matrices describe the same set of people. The PRS matrix holds, for every genotyped sample, a polygenic score for each of \(p\) traits. The phenotype matrix holds, for every phenotyped sample, the measured value of the same traits. If the sample identifiers are correct, the phenotype profile of a person should resemble their own PRS profile more than anyone else's, because each PRS predicts its trait with some correlation. The tool exploits this: for every phenotyped sample it scores every PRS sample as a candidate, ranks the candidates, and reports whether the sample's own PRS wins, how confidently, and who wins otherwise. It therefore detects sample swaps, mislabelled identifiers and contaminated samples, and can propose the correct assignment when the labelled one is wrong.

Because each individual trait is only weakly predicted by its PRS, no single trait can identify a person. The evidence is aggregated over many traits, weighted by how well each PRS predicts its phenotype, so that the aggregate score separates the true match from thousands of alternatives.

### Notation

| symbol | meaning |
|---|---|
| \(X\) | PRS matrix after standardization, \(n_x \times p\); row \(x_j\) is PRS sample \(j\) |
| \(Y\) | phenotype matrix after processing, \(n_y \times p\); row \(y_i\) is phenotyped sample \(i\) |
| \(O_i\) | set of traits observed (not missing) for phenotyped sample \(i\) |
| \(w_k\) | weight of trait \(k\), a non-negative estimate of how well the PRS predicts the phenotype |
| \(W\) | \(\mathrm{diag}(w_1,\dots,w_p)\) |
| \(\mathcal{P}\) | set of mapped pairs \((j,i)\): PRS sample \(j\) and phenotyped sample \(i\) are labelled as the same person |
| \(S(j,i)\) | raw similarity between PRS sample \(j\) and phenotyped sample \(i\) |
| \(Z(j,i)\) | standardized similarity, see section 3 |

### Pipeline

The command runs the following steps in this order.

1. **Load** the PRS, phenotype and optional covariate matrices. Cells matching one of the `--missing-str` strings (default `NA`) are flagged missing. The PRS matrix must be complete.
2. **Match traits** between the two matrices by identical names, or by the explicit `--trait-tsv` mapping. Only shared traits are kept, in the same order in both matrices.
3. **Impute missing phenotypes** if `--missing-as-mean` or `--missing-as-half-min` is set; see section 5. By default nothing is imputed.
4. **Adjust for covariates** if `--cov` is given: the phenotype and covariate matrices are restricted to shared samples, samples with a missing covariate are dropped (or mean-imputed with `--cov-impute-mean`), and each trait is replaced by its residual from a linear regression on an intercept and the covariates.
5. **Rank-based inverse normal transformation** (RINT) of each trait if `--rint` is set: observed values are replaced by \(\Phi^{-1}\big((r - 3/8)/(n + 1/4)\big)\) (Blom offsets, as in REGENIE), where \(r\) is the average rank among the \(n\) observed values of that trait.
6. **Map samples.** The mapped pairs \(\mathcal{P}\) are read from `--sample-tsv`, or formed from identical sample identifiers. Pairs are used to estimate the trait weights and to evaluate each phenotyped sample's own PRS. Phenotyped samples without a pair are still scored against all PRS samples but have no "self".
7. **Standardize** every trait of both matrices to mean zero and unit variance across samples (variance with \(n-1\)). Missing phenotype cells are set to exactly zero after standardization so that they contribute nothing to any sum below.
8. **Estimate trait weights**, unless supplied by `--weights`:
   $$ w_k = \frac{1}{|\mathcal{P}_k|}\sum_{(j,i)\in\mathcal{P}_k} x_{jk}\,y_{ik}, $$
   where \(\mathcal{P}_k\) is the subset of mapped pairs in which trait \(k\) is observed for the phenotyped sample. Since both matrices are standardized, \(w_k\) is essentially the correlation between the PRS and the phenotype for trait \(k\) among mapped samples. Weights below `--min-weight` (default 0, which also removes negative weights) are set to zero. The weights and the number of pairs behind each are written to `[out].weights.tsv.gz`. A trait with zero weight plays no further role.
9. **Score** every (PRS sample, phenotyped sample) pair, giving the \(n_x \times n_y\) matrix \(S\); see section 2.
10. **Standardize scores** within each phenotyped sample to Z-scores, rank the candidates, and write the match calls; see section 3.

Step 8 uses the labelled pairs as if they were correct. A modest fraction of wrong labels merely adds noise to the weights, so the method is robust to the sample swaps it is meant to detect, but the mapping must be mostly right. If it is mostly wrong, the weights are meaningless, and `--auto-lambda` will stop with an error (section 4).

## 2. Scores

Both scoring modes compute a weighted similarity between a PRS profile \(x_j\) and a phenotype profile \(y_i\) over the traits observed for \(i\). They differ in whether traits are treated as independent.

### 2.1 Without `--mahalanobis`: weighted cosine similarity

$$ S(j,i) = \frac{\sum_{k\in O_i} w_k\, x_{jk}\, y_{ik}}{\sqrt{\sum_{k} w_k\, x_{jk}^2}\;\sqrt{\sum_{k\in O_i} w_k\, y_{ik}^2}}. $$

Each trait contributes the product of the candidate's standardized PRS and the sample's standardized phenotype, weighted by how predictive the trait is. The two square roots are the weighted norms of the two profiles, so \(S\) is the cosine of the angle between them in a space where each trait axis is scaled by \(\sqrt{w_k}\), and it lies between minus one and one. The phenotype norm is constant within a phenotyped sample and cancels in the Z-scores of section 3; it only affects `COR`. The PRS norm differs between candidates and does affect the ranking: it prevents a candidate whose PRS profile is large in every trait, as happens with population structure, from scoring highly against everyone merely because of its magnitude.

By default the PRS norm runs over all traits. With missing phenotypes and `--exact-norm` it runs over the traits observed for sample \(i\), which costs one additional matrix product of the same size as the numerator.

**Legacy score (`--no-norm`).** Earlier versions used the weighted inner product divided by the total absolute weight of the observed traits, \(\sum_k w_k x_{jk} y_{ik} / \sum_{k\in O_i} |w_k|\). It has the same numerator, so the two differ only by the candidate norm in the denominator. Because that denominator is constant across candidates, the legacy score ranks candidates by the raw inner product and favours candidates with large PRS profiles. `--no-norm` restores it when earlier results must be reproduced.

This score assumes the traits carry independent evidence. When traits are correlated, as with related phenotypes such as height and weight, or through shared population structure in the PRS, a cluster of correlated traits is counted several times over.

### 2.2 With `--mahalanobis`: whitened cosine similarity

The Mahalanobis mode replaces the trait-wise product by a quadratic form that accounts for correlation between traits.

**Covariance.** Two \(p\times p\) covariance matrices are estimated from the standardized matrices, \(\Sigma_x = X^\top X/n_x\) and \(\Sigma_y = Y^\top Y/n_y\), and blended:

$$ \Sigma = \omega\,\Sigma_x + (1-\omega)\,\Sigma_y + 10^{-8} I, $$

with \(\omega\) given by `--weight-prs-mh` (default 0.5). Because both matrices are standardized, the diagonal of \(\Sigma\) is one and its trace over \(p\) is one, so the identity is the natural shrinkage target.

**Shrinkage.** The covariance is shrunk toward the identity with \(\lambda\) from `--lambda`, or chosen automatically by `--auto-lambda`:

$$ \Sigma_\lambda = (1-\lambda)\,\Sigma + \lambda\, I, \qquad 0 \le \lambda \le 1. $$

Shrinkage is essential in practice. With hundreds of traits, many eigenvalues of \(\Sigma\) are far below one, and inverting an unshrunk \(\Sigma\) amplifies those noisy directions. On real data, intermediate values (roughly 0.3 to 0.6) perform far better than either \(\lambda = 0\) (no shrinkage) or \(\lambda = 1\) (traits treated as uncorrelated).

**Metric.** The weighted precision matrix is

$$ M = W^{1/2}\,\Sigma_\lambda^{-1}\,W^{1/2}, $$

so that trait \(k\) enters with weight \(\sqrt{w_k}\) on each side of the quadratic form. Trait weights must be non-negative, which `--min-weight` at its default guarantees.

**Score.** For phenotyped sample \(i\) with observed traits \(O_i\), let \(M_{O_i}\) denote the sub-matrix of \(M\) restricted to rows and columns in \(O_i\). Then

$$ S(j,i) = \frac{x_{j,O_i}^\top\, M_{O_i}\, y_{i,O_i}}{\lVert x_j\rVert_M\;\lVert y_{i,O_i}\rVert_{M_{O_i}}}, \qquad \lVert v\rVert_A = \sqrt{v^\top A\, v}. $$

The numerator is a whitened inner product over the observed traits. The two norms turn it into a cosine similarity in the whitened space: a candidate whose PRS profile is large in every direction, as happens with population structure, is not favoured merely for its magnitude. The phenotype norm is constant within a phenotyped sample and cancels in the Z-scores, so it only affects `COR`. The PRS norm varies across candidates and does affect the ranking.

By default the PRS norm \(\lVert x_j\rVert_M\) is computed over all traits, once per candidate. With `--exact-norm` it is computed over the observed traits of each phenotyped sample, \(\lVert x_{j,O_i}\rVert_{M_{O_i}}\). When \(M\) is diagonal (\(\lambda = 1\)) that is a single matrix product; otherwise it requires a separate pass for every distinct missingness pattern. The two norms coincide when there are no missing values.

**Special cases.** With \(\lambda = 1\), \(\Sigma_\lambda = I\) and \(M = W\), so the score is exactly the weighted cosine similarity of section 2.1: `--mahalanobis --lambda 1` and the default independence mode produce identical output, the former by a slower route through a \(p \times p\) inversion. With \(\lambda = 0\), the score is the cosine similarity after full whitening by the blended covariance.

**Cost.** The covariance inversion is \(O(p^3)\) and the scoring is one \(n_x \times p \times n_y\) product, the same as the independence mode. Memory is dominated by the \(n_x\times n_y\) score matrix in both modes.

## 3. From scores to match calls

### Z-scores

Within each phenotyped sample \(i\), the scores against all \(n_x\) candidates are standardized:

$$ Z(j,i) = \frac{S(j,i) - \bar S_{\cdot i}}{\mathrm{sd}_j\, S(j,i)}. $$

\(Z(j,i)\) measures how far candidate \(j\) stands out from the crowd of all PRS samples for this phenotyped sample, in standard deviations. Because the standardization is per phenotyped sample, quantities that are constant across candidates (the denominators in section 2) have no influence on \(Z\), and Z-scores are comparable across phenotyped samples with different numbers of observed traits.

### Ranks and the two output files

For a mapped phenotyped sample, `Rank.self` is one plus the number of candidates with a strictly larger Z-score than its own PRS sample, and `Z.self` and `COR.self` are the corresponding values. The five candidates with the largest Z-scores are reported with their Z and COR.

`[out].match.assigned.tsv.gz` contains one row per mapped phenotyped sample and answers "does the labelled PRS sample look right?":

| MatchStatus | condition |
|---|---|
| `BEST_MATCH` | the labelled PRS ranks first |
| `LENIENT_MATCH` | it does not rank first, but \(Z_\text{self} \ge\) `--z-threshold` (default 1.96) |
| `NO_MATCH` | \(Z_\text{self} <\) `--z-threshold` |
| `NO_OBS_TRAITS` | the sample has no observed trait with non-zero weight |

`[out].match.all.tsv.gz` contains one row per phenotyped sample, mapped or not, and answers "who is the best match, and is it clear?". With \(Z_{(1)} \ge Z_{(2)} \ge \dots\) the top five Z-scores and \(\delta\) = `--z-diff` (default 2.0), the status is determined in this order:

| MatchStatus | condition |
|---|---|
| `SELF_BEST` | mapped, and the labelled PRS ranks first |
| `SINGLE_NEW_BEST` | \(Z_{(1)} > Z_{(2)} + \delta\): one candidate clearly leads |
| `MULTI_NEW_BEST` | \(Z_{(k)} > Z_{(k+1)} + \delta\) for some \(k\) in 2 to 4: a clear gap separates a small group of leading candidates from the rest, as with a contaminated or pooled sample |
| `SELF_LENIENT` | mapped, none of the above, and \(Z_\text{self} >\) `--z-threshold` |
| `UNCLEAR` | mapped, none of the above |
| `NO_SELF` | not mapped, and no clear leading candidate |
| `NO_OBS_TRAITS` | no observed trait with non-zero weight |

## 4. Automatic shrinkage: `--auto-lambda`

The shrinkage parameter matters more than any other setting in Mahalanobis mode, and its best value depends on the number of traits, the sample sizes and the correlation structure. `--auto-lambda` chooses it from the data using the mapped pairs as labels.

**Tuning set.** The phenotyped samples that have a mapped PRS sample and at least one observed trait with non-zero weight. Their number must be at least `--auto-lambda-min-self` (default 100); otherwise the criterion below would be too noisy and the command stops.

**Criterion.** For a candidate \(\lambda\), the tuning set is scored against all PRS samples exactly as in section 2 and standardized as in section 3. Writing \(Z_\text{self}(i)\) for the Z-score of the labelled PRS of tuning sample \(i\), and \(Z(j,i)\) for all candidates, `--auto-lambda-metric` selects one of:

| metric | definition | character |
|---|---|---|
| `mean-log-softmax` (default) | \(\frac{1}{m}\sum_i \Big[ Z_\text{self}(i) - \log \sum_j e^{Z(j,i)} \Big]\) | the mean log-probability of picking the labelled PRS under a softmax over all candidates. Smooth in \(\lambda\); saturates for easy samples, so borderline samples drive the optimum; bounded, so mislabelled pairs cannot dominate |
| `mean-z` | \(\frac{1}{m}\sum_i Z_\text{self}(i)\) | smooth, but linear in \(Z\), hence dominated by the easiest samples and more sensitive to mislabelled pairs |
| `mrr` | \(\frac{1}{m}\sum_i 1/\mathrm{rank}_\text{self}(i)\) | the outcome itself, robust, but piecewise constant in \(\lambda\), so the search may stall on a flat step |

**Search.** The criterion is maximized over \(\lambda \in [0,1]\) by golden-section search to a tolerance of 0.01, which takes 12 evaluations. Each evaluation costs one covariance inversion and one scoring pass over the tuning set, so the whole procedure costs roughly ten times a single run of the tuning set, with the phenotype norms omitted since they cancel in \(Z\). A final pass at the chosen \(\lambda\) reports how many tuning samples rank their own PRS first.

**Guard against mismatched labels.** If fewer than `--auto-lambda-min-best` (default 100) tuning samples rank their own PRS first at the chosen \(\lambda\), the command stops. Such a result almost always means the sample mapping is wrong, for example because identifiers do not correspond, and in that case the trait weights estimated from the same mapping are wrong as well. Setting the threshold to 0 disables the check.

`--auto-lambda` requires `--lambda` at its default of 0, and is only meaningful with `--mahalanobis`.

## 5. Missing values

Missing values are recognized by the strings in `--missing-str` (default `NA`, comma-separated for several). The PRS matrix must not contain any; the command stops if it does. Phenotypes and covariates may.

### 5.1 Default: missing phenotypes are ignored

A missing phenotype carries no information about the sample, so the trait is left out for that sample at every stage:

- **Covariate adjustment** regresses each trait on the covariates using only the samples in which the trait is observed. Traits without missing values are adjusted together in a single batch.
- **RINT** ranks only the observed values of each trait.
- **Standardization** computes the mean and standard deviation of each trait from its observed values, then sets missing cells to exactly zero. That zero is not an imputed value: every sum below is either unaffected by it or explicitly restricted to observed traits.
- **Weights** are estimated from the mapped pairs in which the trait is observed, and `N.Obs` in the weights file records how many.
- **Scores.** The numerator and the phenotype norm run over the observed traits only, because missing cells are zero; in Mahalanobis mode they use the sub-matrix \(M_{O_i}\), so unobserved traits do not enter on either side. The covariance \(\Sigma_y\) is estimated from the zero-filled standardized matrix, which amounts to mean imputation for that estimate only. The PRS norm uses all traits by default and the observed traits with `--exact-norm`; see section 2. With `--no-norm`, the legacy denominator is the total absolute weight of the observed traits.
- **`N.Traits`** reports, for every phenotyped sample, the number of observed traits with non-zero weight. A sample with none receives `NO_OBS_TRAITS`.

Two consequences are worth knowing. First, ignoring a trait and imputing its standardized mean give the same numerator, because an imputed zero contributes nothing to any candidate's score; with the default PRS norm over all traits, the Z-scores and ranks are then identical as well, and only `COR` and `N.Traits` differ. Second, the observed block \(M_{O_i}\) of the precision matrix is the precision of the observed traits conditional on the unobserved ones, not the inverse of their marginal covariance. When a missing trait is strongly correlated with an observed one, this gives the observed trait more weight than the marginal form would. A per-pattern marginal form is a possible future refinement.

### 5.2 `--missing-as-mean` and `--missing-as-half-min`

Both options fill missing phenotype cells immediately after trait matching (step 3 of the pipeline), before covariate adjustment and RINT, and then treat the filled cells as observed everywhere: in the regression, the ranks, the standardization, the weights and the scores. `N.Traits` then counts all traits with non-zero weight.

- `--missing-as-mean` fills each missing cell with the mean of the observed values of that trait. As noted above, this changes Z-scores little relative to ignoring, but the filled cells do participate in RINT ranks, in the standard deviation, in the weights, and in the Mahalanobis quadratic form.
- `--missing-as-half-min` fills each missing cell with half of the minimum observed value of that trait, the usual convention in metabolomics and other assay data when a missing measurement means "below the detection limit". It assumes positive measurements; a warning is issued for traits whose minimum is not positive, since half of it is then not below the observed values. With `--rint`, the filled cells form a tied block strictly below all observed values and are transformed to a low value, so the exact fraction of the minimum does not matter there. A quick diagnostic for informative missingness is whether the mapped samples' own standardized PRS for a trait is systematically negative when that trait is missing.

The two options cannot be combined. Traits with no observed value at all are left missing with a warning.

### 5.3 Missing covariates

By default, every sample with at least one missing covariate is dropped from the phenotype matrix before adjustment, the number dropped is logged, and those samples do not appear in the output at all. Entries of `--sample-tsv` that refer to dropped samples are skipped. With `--cov-impute-mean`, missing covariate values are instead replaced by the mean of the observed values of that covariate, and the samples are kept.

### 5.4 Choosing an option

| situation | recommendation |
|---|---|
| missing at random, or reasons unrelated to the trait value | default (ignore) |
| missing means below detection limit | `--missing-as-half-min` |
| reproducing an analysis that imputed means | `--missing-as-mean` |
| a few samples lack a covariate | default (drop) if they are few, `--cov-impute-mean` if they are many |

## 6. Practical notes

- Use `--rint` for phenotypes with skewed distributions; the scores rely on linear relationships between standardized values.
- Adjust for covariates that affect phenotypes but not PRS, such as age, sex and assay batch, so that the weights and scores reflect genetic prediction rather than shared covariate effects.
- In Mahalanobis mode, prefer `--auto-lambda`, or set `--lambda` in the 0.3 to 0.6 range; both extremes perform poorly on data with many correlated traits.
- Inspect `[out].weights.tsv.gz`: traits with very few observations behind their weight are estimated noisily, and traits with weights near zero contribute nothing. Thresholds on the number of observations per trait and per sample are a planned addition.
- `--z-threshold` and `--z-diff` control the calls, not the scores. Lower values call more matches and more clear winners; the defaults (1.96 and 2.0) are conservative.
