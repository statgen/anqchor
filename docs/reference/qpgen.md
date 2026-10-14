# qpgen integration

anqchor does no association statistics itself: [qpgen](https://github.com/hyunminkang/qpgen)
(`qpgentools`) computes polygenic scores, identity scores and association tests. anqchor
prepares the inputs, builds the command lines, runs them as steps and reads the outputs. The
qpgen flag names appear in one module, `anqchor.qpgen.commands`; the output layouts in one
module, `anqchor.qpgen.outputs`.

## Version requirement: qpgen 0.2.0 or later

qpgen is vendored as the git submodule `extern/qpgen` (branch `dev`), pinned at commit
`5b811e2`, which reports version 0.2.0:

```console
$ extern/qpgen/bin/qpgentools --version
qpgentools 0.2.0 (git 5b811e2)
```

Categorical covariates (`--cat-cov-names`, below) need a build from commit `5b811e2` or
later: `cb97eb0` added them to `pair-assoc`, `rect-assoc` and `region-assoc`, and `5b811e2`
to `match-prs-pheno`. Both report version 0.2.0, as does the earlier `76f2fd1`, so the
version string alone does not tell them apart; an older build does not accept
`--cat-cov-names`. `anqchor qpgen check` prints the commit of the checkout and fails with
`missing options: ...` when `match-prs-pheno`, `pair-assoc` or `rect-assoc` lacks
`--cat-cov-names` or `--max-cat-levels`.

Every qpgen command also logs this line first. Builds before 0.2.0 have no `--version`.

anqchor itself does not read the version: `anqchor qpgen check` verifies the subcommands and
prints the binary's sha256 and the git commit of the checkout it sits in, and `run.json`
records the same. Check `qpgentools --version` yourself when you use a binary built
elsewhere.

### What changed in qpgen 0.2.0

qpgen 0.2.0 collects the fixes made while anqchor was built and tested on the MAGE example.
The ones that change anqchor's results:

| Change | Effect on older builds |
|---|---|
| `pair-assoc` adjusts the **genotypes** for the covariates as well as the phenotype (Frisch-Waugh-Lovell), like `rect-assoc` | with covariates, `pair-assoc` effects and −log10 p were attenuated and differed from `rect-assoc` |
| With covariates, `pair-assoc`, `rect-assoc` and `region-assoc` use **df = N − 2 − k** (k = rank of the covariates) | `SE` too small by √((N−2−k)/(N−2)); −log10 p optimistic with many covariates |
| `pair-assoc` always fits an intercept; `rect-assoc` centres phenotypes without covariates | without covariates, `SE` was inflated for phenotypes with a large mean |
| Non-numeric values (such as `XX`/`XY`) are an error in every phenotype and covariate loader, except in covariates named with `--cat-cov-names` | they were read as 0 without a warning |
| `pair-assoc` reads `NA` as missing and drops that sample from the trait's tests (`--missing-str`); missing genotypes are dropped per variant (complete case) | `NA` was read as 0; missing genotypes were mean-imputed |
| Samples are keyed by IID in every loader, and duplicate IIDs are an error | a `.psam` with real family IDs (FID ≠ IID) could misalign genotypes and phenotypes silently |
| Input tables are tab-delimited by default (`--space-delim` for single-space files) | fields were split on any run of whitespace, so IDs with spaces were cut |
| Rank of the covariates computed on scaled columns; constant, duplicate and near-duplicate covariates leave `BETA` unchanged (with a warning) | near-duplicate covariates could shift `BETA` substantially without a warning |
| `pair-assoc` no longer requires lexicographic sample order, and stops when no samples overlap | it aborted, or ran on with zero samples |
| `pair-assoc` header starts with `TRAIT` (was `#TRAIT`) | anqchor accepts both |
| `match-prs-pheno` warns when fewer than 10 traits have a non-zero weight | no warning |
| `qpgentools --version` | none |
| Categorical covariates with `--cat-cov-names` (REGENIE `--catCovarList`) in `match-prs-pheno`, `pair-assoc`, `rect-assoc` and `region-assoc` (commits `cb97eb0`, `5b811e2`) | text covariates such as `XX`/`XY` had to be recoded to numbers first |

anqchor always writes tab-separated, regenie-layout matrices for qpgen (numeric except for
categorical covariates, which qpgen encodes) and keys samples by IID, so it does not depend on the delimiter and FID changes. The covariate
adjustment, degrees of freedom and missing-value changes do change the numbers: use 0.2.0.

## Build and discovery

```bash
git submodule update --init --recursive
anqchor qpgen build -j 8      # htslib -> qgenlib -> qpgen; produces extern/qpgen/bin/qpgentools
anqchor qpgen check           # path, sha256, git commit, required subcommands
```

Discovery order: `--qpgen` / `qpgen.binary` → `$ANQCHOR_QPGEN` → `qpgentools` on `PATH` →
`extern/qpgen/bin/qpgentools`. An explicit path or `$ANQCHOR_QPGEN` that does not name an
executable is an error (exit code 4).

## Commands and flags used

All genotype input is the prepared 4-column list with absolute prefixes,
`--pgen-list prepare/genotypes/pgen_list.tsv`. All phenotype and covariate input is in
regenie layout, so `--pheno-format` and `--cov-format` are never needed.

### Categorical covariates

Covariate files keep categorical covariates as text (`sex` is `XX`/`XY`). Every qpgen command
that takes `--cov` also gets `--cat-cov-names` with the categorical covariates in that file:
the names given with anqchor's `--cat-cov-names` (config `datasets.<ds>.cov_categorical`) plus
every column with a non-numeric value, detected when the Makefile is built. A command whose
covariates include no categorical one gets no `--cat-cov-names` (`--cov-include PC1-5,PEER1-10`).
qpgen follows REGENIE: a covariate with L levels becomes L − 1 indicator columns, the first
level in sorted order is the reference, and `--max-cat-levels` (default 10) caps the levels. anqchor
checks the level counts itself (see [covariates](../user-guide/inputs.md#covariates)), so when
a categorical covariate a command uses has N > 10 levels it adds `--max-cat-levels N`; with
10 or fewer levels the flag is omitted. The log shows the expansion:

```text
NOTICE [...] - Categorical covariate sex: 2 levels (XX, XY); reference XX; 1 indicator columns
NOTICE [...] - Covariate matrix has 6 columns after expanding categorical covariates
```

The results equal those of the earlier anqchor-side numeric recoding (`sex` → 1/2): the
covariates span the same column space.

| Command | Step | Flags anqchor passes |
|---|---|---|
| `pair-prs` | `prs.<ref>` | `--pgen-list`, `--pairs prepare/references/<ref>/pairs.tsv.gz`, `--out prs/<ref>/prs`, `--sample` (with `prs_samples`) |
| `match-prs-pheno` | `match.<ds>.qpgen` | `--prs`, `--pheno`, `--out`, `--cov` and `--cat-cov-names` (with `use_cov`), `--sample-tsv` (with a `sample_map`), `--trait-tsv` (when trait IDs differ), `--missing-str` (non-default), `--rint`, `--mahalanobis` + `--auto-lambda` [`--auto-lambda-min-self/-min-best`] or `--lambda`, `--min-weight` (non-zero), `--z-threshold`, `--z-diff`, `--threads` |
| `pair-assoc` | `cis.<ds>.<s>.assoc` | `--pgen-list`, `--pheno`, `--pairs prepare/datasets/<ds>/ref.<ref>/pairs.tsv.gz`, `--out <file>.tsv.gz`, `--cov` and `--cat-cov-names` (with covariates), `--skip-rint` (unless `--rint`) |
| `rect-assoc` | `cis.<ds>.<s>.rect.<t>` | `--pgen-list`, `--pheno`, `--var-list`, `--pheno-list`, `--cov`, `--cat-cov-names`, `--rint-after-adj` (with `--rint`), `--out` |
| `rect-assoc` | `trans.<ds>.<s>.rect.<t>` | as above plus `--min-maf`, `--min-mac`, `--max-chunk-vars`; `--pheno-list` unless `traits: all` and `trait_chroms: all` |

For example, the association step of `anqchor cis run ... --cov-include sex,PC1-5`:

```text
qpgentools pair-assoc --pgen-list prepare/genotypes/pgen_list.tsv \
  --pheno prepare/datasets/dataset/pheno.regenie.tsv.gz \
  --pairs prepare/datasets/dataset/ref.reference/pairs.tsv.gz \
  --out cis/dataset/strategies/default/pair_assoc.tsv.gz \
  --cov cis/dataset/strategies/default/cov.regenie.tsv.gz --cat-cov-names sex
```

and the identity step of `anqchor match all` with `--cov` (the `--out` prefix is a
temporary directory inside `match/dataset/`):

```text
qpgentools match-prs-pheno --prs prs/reference/prs.prs.tsv.gz \
  --pheno prepare/datasets/dataset/pheno.regenie.tsv.gz --out match/dataset/.qpgen.../qpgen \
  --cov prepare/datasets/dataset/cov.regenie.tsv.gz --cat-cov-names sex \
  --rint --z-threshold 1.96 --z-diff 2 --threads 1
```

Every step's exact command is in `steps/<step>.json` (`commands`) and its output in
`logs/<step>.log`. `match-prs-pheno` runs into a temporary prefix and its outputs are moved
into place, so a failed run leaves no partial files.

## The statistical model of the association tools

For each (variant, trait) with k covariates C (plus an intercept):

1. the phenotype y is residualised on C: ỹ = (I − P_C) y;
2. the genotype g is residualised on C: g̃ = (I − P_C) g;
3. β̂ = g̃ᵀỹ / g̃ᵀg̃, with the residual variance on **df = N − 2 − k** degrees of freedom,
   where k is the numerical rank of the covariates (columns of [1, C] minus 1); without
   covariates df = N − 2;
4. `rint` (`*_RINT` columns of `pair-assoc`, `--rint-after-adj` in `rect-assoc`) replaces ỹ
   by the rank-based inverse normal transform of ỹ and tests it the same way.

\[
\hat\beta = \frac{\tilde g^\top \tilde y}{\tilde g^\top \tilde g},\qquad
\widehat{\mathrm{SE}}^2 = \frac{\tilde y^\top \tilde y - \hat\beta^2\,\tilde g^\top \tilde g}{(N-2-k)\;\tilde g^\top \tilde g}
\]

Steps 1-3 are the Frisch-Waugh-Lovell form of the full regression `y ~ 1 + C + g`, so `BETA`,
`SE`, `TSTAT` and `LOG10P` equal an ordinary least-squares fit with the covariates in the
model. `pair-assoc` and `rect-assoc` therefore agree on shared pairs: on the example data a
cis run with sex, 5 PCs and 10 PEER factors replicates the same 231 (as given) or 230 (with
`--rint`) of 367 pairs at p < 10⁻⁵ with either engine.

Missing genotypes are handled differently: `pair-assoc` drops the sample from that variant's
test (complete case, `N` = called genotypes = regression sample size); `rect-assoc`
mean-imputes them (`N` counts called genotypes, df uses all samples). With few missing calls
the difference is negligible.

## Defences kept for older builds

`anqchor.qpgen.outputs` keeps two fallbacks so that older binaries still run:

- `pair-assoc` rows one field short of the header (raw `LOG10P` missing, before qpgen commit
  `89eab88`): raw `LOG10P` is recomputed from `TSTAT`, and the run's
  `cis/<ds>/summary.json` lists it under `recomputed_log10p`;
- `rect-assoc` `N` is always recomputed from `N_RR + N_RA + N_AA`.

The readers check the columns they rely on and fail with a descriptive error when qpgen's
layout changes. These fallbacks do not make an older build's statistics correct; use 0.2.0.

## Known qpgen behaviours

- Variants are looked up by exact `CHROM:POS:REF:ALT`, with no allele flips or normalisation;
  a variant not in the pvar is skipped. anqchor checks presence itself and reports
  `variant_absent`.
- `pair-assoc` and `rect-assoc` match samples by identical ID (no map).
- `rect-assoc` stops on a missing phenotype or covariate value (including a missing level of a
  categorical covariate); `pair-assoc` drops the sample. anqchor's prepared covariate files
  have no missing values (`cov_missing`).
- `match-prs-pheno` holds an `n_prs × n_pheno` score matrix in memory (729 × 729 in the
  example is trivial; 150,000 × 20,000 is about 24 GB). `prs_samples` limits the PRS samples.
- The `rect-assoc` output has `11 + 4 × n_traits` columns per variant; keep rectangles small
  (`trait_chroms`, `traits: reference`).

The remaining open qpgen items (allele matching, dosage pgen files not yet tested with the
new `pair-assoc` path) are listed in `extern/qpgen/dev/review-2026-10-anqchor.md`.
