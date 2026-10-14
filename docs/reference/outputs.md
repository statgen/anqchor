# Output files

Every file a fresh run writes, with its columns or keys. The lists below were taken from
runs of `match run` + `match report`, `cis run` (five runs) + `cis compare`, and `trans run`
(five runs) + `trans compare` on the example data. Paths are relative to the run (or comparison) directory; file names come
from one place, `anqchor.pipelines.layout.Layout`. Conventions:

- TSV files have a header (except the plain ID lists and maps noted below), use `NA` for
  missing values and `1`/`0` for booleans; `.gz` files are gzip/bgzip compressed.
- Every JSON file has `schema_version` (currently `1`; `3` for the cis `summary.json` of a
  run or comparison, which a comparison requires of every run; `2` for the trans
  `summary.json` of a run or comparison, which `trans compare` and `trans report` require).
- Paths inside summaries are relative to the run directory.
- In a flags-only run the dataset is called `dataset` and the reference `reference`, so
  `<ds>` = `dataset` and `<ref>` = `reference` below.

## Run-level files

| File | Content |
|---|---|
| `anqchor.resolved.yaml` | complete settings of the last command that ran steps in the directory (defaults filled in, paths absolute, flags included; a dry run or `anqchor steps` leaves it unchanged); a valid `--config` file |
| `run.json` | provenance; see [Running](../user-guide/running.md#provenance-runjson) |
| `Makefile` | generated; rewritten on every invocation |
| `.anqchor.lock` | run-directory lock held while a command runs (a second command on the same directory exits with code 2) |
| `steps/<step>.json` | step signature: `name kind pipeline inputs outputs commands params`, plus `params.qpgen_sha256` for steps that run qpgen and `anqchor_version` for Python steps |
| `steps/<step>.done`, `.seconds` | stamp and wall time of the last successful run |
| `logs/<step>.log` | step stdout and stderr |

## prepare

### `prepare/genotypes/`

| File | Content |
|---|---|
| `pgen_list.tsv` | qpgen pgen list, 4 columns `CHROM BEG END PREFIX`, absolute prefixes, no header |
| `links/` | symlinks giving non-standard file names a common prefix (only when needed) |
| `samples.txt` | genotype IIDs, one per line, psam order |
| `genotypes.json` | `pgen_list` or `pfile`, `n_samples`, `chrom_style` (`chr`/`plain`), `chromosomes`, `entries` (`chrom beg end pgen psam pivar`), `n_entries` |

### `prepare/references/<ref>/`

| File | Content |
|---|---|
| `pairs.tsv.gz` | normalised pairs: `trait variant beta se pval` (`pval` may be `NA`); variant IDs in the genotype chromosome style; invalid rows and duplicates removed |
| `variants.tsv` | one row per unique variant: `variant chrom pos ref alt is_indel present` |
| `variants.present.txt` | present variants, one CPRA per line, sorted, no header |
| `reference.json` | `reference kind n_pairs_input n_pairs n_invalid n_duplicates n_traits n_variants n_variants_list n_variants_list_invalid n_present n_absent n_absent_snv n_absent_indel absent_examples invalid_examples per_chrom chrom_renamed chrom_style presence_method` |

### `prepare/datasets/<ds>/`

| File | Content |
|---|---|
| `pheno.regenie.tsv.gz` | phenotypes, regenie layout `FID IID trait…`, the dataset's own sample IDs |
| `cov.regenie.tsv.gz` | covariates, regenie layout `FID IID sex PC1 …` (only with `cov`); categorical covariates keep their text levels (`sex` is `XX`/`XY`), constant columns are dropped and missing values are handled by `cov_missing`. There is no numeric copy: qpgen encodes the categorical columns (`--cat-cov-names`) |
| `cov.json` | `raw_columns` (the columns of `cov.regenie.tsv.gz`), `categorical_names` (the categorical covariates, passed to qpgen as `--cat-cov-names`), `categorical` (name → {level: code}, levels in sorted order), `names`, `columns` and `groups` (original name → columns of the equivalent numeric encoding: L − 1 per categorical covariate; used only to count covariates for degrees of freedom and to resolve `--cov-include`), `dropped_constant`, `dropped_all_missing`, `n_samples`, `missing_policy`, `n_missing_handled`, `dropped_samples` |
| `traits.tsv` | `trait chrom start end` for every phenotype trait (`NA` when unknown), genotype chromosome style |
| `samples.txt` | phenotype sample IDs |
| `sample_map.tsv` | `geno_id pheno_id`, no header, IDs present on both sides (identity pairs without a `sample_map`) |
| `dataset.json` | `dataset omics pheno pheno_format pheno_fast_path n_samples_pheno n_traits n_missing_values n_samples_geno n_samples_geno_overlap n_samples_mapped n_samples_analysable (phenotyped, genotyped and with covariates: usable by `cis`/`trans`) warnings (dropped or covariate-less samples, `sample_map` ID mismatch) n_pheno_unmapped pheno_unmapped_examples sample_map [n_sample_map_rows n_sample_map_rows_geno_missing n_sample_map_rows_pheno_missing n_geno_with_several_pheno] has_cov n_cov_columns cov_columns n_cov_samples n_pheno_without_cov cov_categorical (name → {level: code}) cov_dropped_constant cov_missing_policy cov_n_missing_handled has_loci n_traits_with_loci trait_annot_source trait_annot_rule trait_id_rule seconds_pheno seconds_total` |

### `prepare/datasets/<ds>/ref.<ref>/`

| File | Content |
|---|---|
| `trait_map.tsv` | `ref_trait pheno_trait`, no header, mapped traits only |
| `pairs.tsv.gz` | `trait variant beta se pval ref_trait present`: reference pairs rewritten to dataset trait IDs (`present` 0/1) |
| `variants.txt` | present variants of the mapped pairs, one CPRA per line, sorted |
| `traits.txt` | dataset trait IDs of the mapped reference traits |
| `mapping.json` | `rule trait_map n_ref_traits n_dataset_traits n_mapped_traits n_unmapped_traits unmapped_examples n_pairs n_pairs_mapped n_pairs_present n_variants_present n_traits_with_present_variant` |

## PRS

| File | Content |
|---|---|
| `prs/<ref>/prs.prs.tsv.gz` | `pair-prs` scores, regenie layout `FID IID trait…`, one row per genotyped sample |
| `prs/<ref>/prs.se.tsv.gz` | standard errors, same layout |

## match { #match }

Under `match/<ds>/`.

| File | Content |
|---|---|
| `qpgen.weights.tsv.gz` | qpgen: `Trait Weight N.Obs` |
| `qpgen.match.assigned.tsv.gz` | qpgen: mapped samples, `ID.Pheno N.Traits MatchStatus ID.self Z.self COR.self Rank.self`; `MatchStatus` is `BEST_MATCH`, `LENIENT_MATCH`, `NO_MATCH` or `NO_OBS_TRAITS` |
| `qpgen.match.all.tsv.gz` | qpgen: all phenotyped samples, the same seven columns plus `ID/Z/COR.{1st…5th}`; here `MatchStatus` is the detailed status (`SELF_BEST`, `SELF_LENIENT`, `SINGLE_NEW_BEST`, `MULTI_NEW_BEST`, `UNCLEAR`, `NO_SELF`, `NO_OBS_TRAITS`) that anqchor uses |
| `qpgen.params.json` | `argv qpgen dataset reference mahalanobis lambda_config lambda_mode (independent/auto/fallback/fixed) lambda_used rint use_cov z_threshold z_diff min_weight n_traits_used n_trait_map_rows n_dropped_* n_prs_samples n_samples_pheno n_samples_mapped n_samples_mapped_dataset n_sample_map_rows_dropped trait_tsv_used sample_tsv_used warnings log_excerpt seconds` |
| `trait_map.used.tsv` | trait map given to qpgen, no header |
| `samples.tsv` | one row per phenotyped sample, see below |
| `swaps.tsv` | samples whose best match is another individual: `kind (reciprocal/reciprocal-weak/one-way) pheno_id geno_id_assigned group best_id z_best z_self rank_self partner_pheno_id partner_geno_id_assigned partner_group partner_best_id` |
| `candidates.tsv` | rescue candidates (unmapped samples with a clear best match that is not another sample's confirmed match): `pheno_id candidate_geno_id z_best second_id z_second z_gap n_traits candidate_assigned_to` |
| `duplicates.tsv` | samples that match two genotype IDs about equally (rule in [match](../pipelines/match.md#status-rules)): `pheno_id geno_id_assigned group geno_id_1 z_1 geno_id_2 z_2 z_3 assigned_in_pair` |
| `related.tsv` | possible relatives (a block of two or more close candidates; samples in `duplicates.tsv` left out): `pheno_id geno_id_assigned group n_close close_ids close_z assigned_in_close` |
| `weights.tsv` | `trait ref_trait weight n_obs weight_le0` |
| `summary.json` | see below |

`samples.tsv`: `pheno_id geno_id_assigned qpgen_status group z_self cor_self rank_self n_traits
best_id z_best cor_best second_id z_second z_gap proposed_geno_id note id_1…id_5 z_1…z_5`
(meanings in [match](../pipelines/match.md#outputs)).

Also `report/match.html` after `match report` or `match all`.

`summary.json` keys: `pipeline dataset omics reference verdict verdict_status (good/warning/critical)
status_counts group_counts n_pheno n_mapped n_unmapped n_scored n_not_scored n_fail
n_fail_two_ids n_reciprocal_swaps n_rescue_candidates n_duplicates n_related n_status_disagreements failing (list of failing
samples) n_prs_samples n_traits_used n_traits_weight_le0 weights (min q05 q25 median mean q75
q95 max) mahalanobis lambda_mode lambda_used rint use_cov z_threshold z_diff
trait_tsv_used sample_tsv_used qpgen_argv qpgen_log_excerpt dataset_counts warnings` (`warnings` includes
the number of samples not scored for lack of covariates).

## cis run { #cis }

Under `cis/<ds>/`. A flags-only run has one strategy, `default` (so `<s>` = `default`), and
its rows are labelled with the **run label**, the name of the run directory (`run` and
`label` columns). A configuration file listing several strategies labels them by strategy
name. A run tests one phenotype transform: as given, or RINT with `--rint`
(`pipelines.cis.rint`).

| File | Content |
|---|---|
| `strategies/<s>/cov.regenie.tsv.gz` | the run's covariates, regenie layout, categorical covariates as text levels (absent without covariates) |
| `strategies/<s>/pheno.regenie.tsv.gz` | strategy phenotypes (only with a `pheno` override in a config file) |
| `strategies/<s>/strategy.json` | `pipeline dataset strategy cov_source cov_columns n_cov [cat_cov_names cov_include cov_empty cov_file cov_encoding] pheno`; `cov_columns` are the covariate names in the file, `n_cov` counts the numeric columns after encoding (degrees of freedom), `cat_cov_names` the categorical covariates passed to qpgen as `--cat-cov-names` (absent without covariates) |
| `strategies/<s>/pair_assoc.tsv.gz` | `pair-assoc` output (pair engine): `TRAIT CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N N_RR N_RA N_AA INFO BETA SE TSTAT LOG10P`, plus `BETA_RINT SE_RINT TSTAT_RINT LOG10P_RINT` with `--rint` |
| `strategies/<s>/rect_assoc.<t>.tsv.gz` | `rect-assoc` output (rect engine); `<t>` = `raw`, or `rint` with `--rint` |
| `pairs.long.tsv.gz` | one row per mapped pair, see below |
| `strategies.tsv` | one row per run, see below |
| `traits.tsv` | `trait ref_trait n_pairs n_variant_absent`, per run label: `best_log10p:<run>`, `level:<run>` (strongest level reached, in `level_order`), then `delta_log10p:<run>` (vs the first run; only when there is more than one), and `best_run` |
| `power_curve.tsv` | `run threshold_log10p n_tested n_replicated frac_replicated` (sign-concordant with p < 10^−threshold, as the levels; thresholds 0, 0.5, 1, …) |
| `summary.json` | see below |

Also `report/cis.html` after `cis report` or `cis all`.

`pairs.long.tsv.gz`, one row per (run, mapped pair):

| Column | Meaning |
|---|---|
| `run` | the run label |
| `trait`, `ref_trait`, `variant` | dataset trait ID, reference trait ID, CPRA |
| `beta_ref`, `se_ref` | reference effect (ALT allele) and standard error |
| `z_ref` | reference z-score: `beta_ref / se_ref`, or, without `se_ref`, the signed normal quantile of `pval_ref` |
| `pval_ref`, `log10p_ref` | the reference's `pval` column (or the column named by `columns.pval`) and its −log10, read from the text so that p-values below the double range (`1e-400`) keep their size; without a usable p-value, `log10p_ref` is the two-sided normal p-value of `z_ref` |
| `status` | `tested`, `variant_absent`, `not_tested` (no association row, or no finite statistic) |
| `a1freq`, `n`, `beta`, `se`, `tstat`, `log10p` | study statistics (ALT allele, the same allele as `beta_ref`); RINT statistics with `--rint` |
| `z` | signed study z-score `beta / se` (as `z_ref`); without a usable `se`, the normal quantile of the two-sided p-value with the sign of `beta` |
| `sign_concordant` | 1/0, `NA` when not tested |
| `rep_<level>` | 1/0 for each configured level (`sign nominal fdr bonferroni strict`) |

`strategies.tsv`, one row per run: `run label rint (1 with RINT) n_cov n_pairs n_tested
n_variant_absent n_not_tested n_rep_<level>… frac_rep_<level>… median_log10p mean_z2
beta_corr frac_sign_discordant pi1 bonferroni_log10p`; with the rect engine also
`n_control_cells lambda_gc_control n_control_p_lt_1e5 n_control_cross_chrom
lambda_gc_control_cross_chrom n_control_cross_chrom_p_lt_1e5`; last `mean_z2_common` (mean
z² over the pairs tested in every row's run) and `rel_neff` = (`mean_z2_common` − 1) / (first
run's `mean_z2_common` − 1), the effective sample size relative to the first run (1 for a
single run).

`summary.json` keys: `schema_version dataset reference engine rint levels (rule per level)
level_order (levels from weakest to strongest by their −log10 p cut-off for this run)
transform_semantics (as_given, rint, engines, degrees_of_freedom) runs (label strategy rint
n_cov per row) baseline n_ref_pairs n_pairs_mapped n_pairs_unmapped n_traits_mapped
n_traits_unmapped n_variant_absent n_absent_variants n_absent_variants_indel
n_absent_variants_snv absent_variant_examples unresolved_fraction unresolved_by_run
max_unresolved_fraction recomputed_log10p per_run (the strategies.tsv numbers by label)
ranking (labels per level, best first) recommended (run label key_level reason) warnings
sentinels` (reserved, `null`).

## cis compare { #cis-compare }

`anqchor cis compare RUN_DIR... --out DIR` writes into `DIR`. Each compared run is identified
by its label (`run` = `label` column); a run whose configuration listed several strategies
appears as `<label>:<strategy>`. Runs may differ in the phenotype transform. A run written
by an anqchor version with another cis summary schema is refused with a request to re-run it.

| File | Content |
|---|---|
| `runs.tsv` | one row per run: `label run_dir dataset genotypes n_geno_samples pheno n_pheno_samples n_samples_analysable cov cov_include (all, the names given, or NA) n_cov cov_columns rint reference reference_pairs engine n_pairs_mapped` |
| `strategies.tsv` | the columns of a run's `strategies.tsv` (`run` = `label` = run label), plus `n_common_tested` and `n_rep_<level>_common` (replication counted only over the pairs tested in every compared run); `mean_z2_common` and `rel_neff` are over those pairs, relative to the first run |
| `pairs.long.tsv.gz` | the runs' `pairs.long.tsv.gz` stacked, `run` = run label; levels re-derived with the first run's level settings |
| `traits.tsv` | as in a run, one column group per run label; `delta_log10p` is relative to the first run |
| `power_curve.tsv` | as in a run, per run label |
| `summary.json` | the keys of a run's `summary.json` (without `rint`, `n_traits_unmapped`, `recomputed_log10p`, `unresolved_by_run`), plus `kind` (`comparison`), `runs` (the `runs.tsv` rows) and `comparability` (`same_reference`, `reference_sha256` per label, `rint` per label); `baseline` is the first run's label; `warnings` include each run's warnings prefixed by its label |
| `compare.json` | manifest: `kind` (`cis`), `runs` (`label run_dir dataset settings`), `labels`, `argv`, `cwd` |
| `settings/<label>.anqchor.yaml` | a copy of each run's `anqchor.resolved.yaml` |
| `report/cis.html` | the comparison report (unless `--no-report`) |

## trans { #trans }

Under `trans/<ds>/` (`<ds>` is `dataset` in a flags-only run). `trans run` writes no HTML;
`trans report` or `trans all` writes `report/trans.html`. One run per directory, labelled
`run` (the run directory name). Summary schema 2.

| File | Content |
|---|---|
| `variants.txt` | rectangle variants, CPRA per line |
| `traits.txt` | rectangle genes (dataset IDs) |
| `lists.json` | `variants` {`source n_listed n_invalid n_renamed_chrom n_off_genome n_kept`}, `traits` {`source trait_chroms n_missing_from_dataset allowed_chroms n_candidates n_kept n_dropped_off_chrom n_dropped_no_locus examples_*`} |
| `strategies/<s>/cov.regenie.tsv.gz`, `strategy.json` | the run's covariates (`<s>` is `default` in a flags-only run) |
| `strategies/<s>/rect_assoc.<t>.tsv.gz` | `rect-assoc` output (`<t>` = `raw`, or `rint` with `--rint`): `#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N N_RR N_RA N_AA TEST` then `BETA.T SE.T TSTAT.T LOG10P.T` per gene |
| `reference_pairs.tsv` | one row per reference pair: `run trait variant ref_trait pair_class (cis/trans/no_locus) trans_kind (cross_chrom/same_chrom_distal) chrom pos trait_chrom trait_start distance beta_ref se_ref pval_ref status (tested/variant_absent/not_in_rectangle/not_tested) variant_tested a1freq n beta se z log10p sign_concordant rep_sign rep_nominal (p < 0.05) rep_strict (p < 1e-5)` |
| `hotspots.tsv` | one row per bin with a tested variant, in genome order: `chrom bin_start bin_end n_variants n_genes_<threshold>… top_trait top_variant top_log10p`; `<threshold>` is `p1e-4 p1e-5 p1e-6 p5e-8 p1e-8 p1e-10`, the configured levels (e.g. `p0.001`) and `bonferroni`; the top cell is the strongest trans cell of the bin among genes not in cis with it (NA without one at the loosest threshold) |
| `hotspot_genes.tsv` | one row per bin and gene with a trans cell below the loosest grid threshold (genes in cis with the bin excluded): `chrom bin_start bin_end trait trait_chrom trait_start variant beta se log10p` (the strongest variant of the bin for that gene) |
| `qq.tsv` | the cross-chromosome cells used for inflation, at most ~2,000 points (the 250 strongest, then log-spaced ranks): `n_cells rank expected_log10p observed_log10p band_lo_log10p band_hi_log10p` (95% range of the order statistic, Beta(rank, n − rank + 1)) |
| `variants.tsv` | one row per tested variant: `variant chrom pos bin_start a1freq n n_cis n_cross_chrom n_same_chrom_distal n_sig_trans (p < nominal) max_log10p_trans top_trans_trait max_log10p_cis top_cis_trait excluded_from_calibration` |
| `pairs.sig.tsv.gz` | every cell with p < `levels.nominal` (any class), strongest first, at most `max_sig_rows`: `variant chrom pos trait trait_chrom trait_start cell_class (cis/cross_chrom/same_chrom_distal/no_locus) reference_pair distance a1freq n beta se log10p` |
| `summary.json` | see below |

`summary.json` keys: `schema_version dataset reference run strategy rint n_cov cov_source
cov_columns n_variants n_traits (listed) n_variants_tested n_traits_tested n_samples n_cells
n_cells_trans cells (cis cross_chrom same_chrom_distal no_locus) variants traits trait_chroms
variant_selection trait_selection known_pairs n_reference_pairs`; `replication`: `levels
(nominal 0.05, strict 1e-5)` and per class `cis`, `trans`, `no_locus`: `n_pairs n_tested
n_rep_sign n_rep_nominal n_rep_strict median_log10p` (`trans` also `n_pairs_cross_chrom
n_pairs_same_chrom_distal`); `hotspots`: `bin_bp n_bins thresholds (key label alpha log10p
roles, loosest first) k counts (thresholds × k: bins with at least k genes) default (threshold
k n_hotspots) bonferroni_p`; `inflation`: `n_cells n_cells_excluded n_reference_trans_pairs
n_variants_excluded n_genes_excluded exclusion_window_bp median_log10p lambda_gc counts (alpha
label observed expected at p < 1e-3 and 1e-5)`; `n_sig_cells sig_rows_written params
definitions warnings memory_bound seconds`. Definitions: [trans](../pipelines/trans.md#definitions).

## trans compare { #trans-compare }

`anqchor trans compare RUN_DIR... --out DIR` writes into `DIR`. Each compared run is
identified by its label (the run directory name or `--label`; `<label>:<ds>` for a run with
several datasets). A run written by an anqchor version with another trans summary schema is
refused with a request to re-run it. Summary schema 2.

| File | Content |
|---|---|
| `runs.tsv` | one row per run: `label run_dir dataset genotypes n_geno_samples pheno n_samples cov cov_include (all, the names given, or NA) n_cov cov_columns rint reference reference_pairs`; per class `<c>` = `cis`, `trans`: `<c>_n_pairs <c>_n_tested <c>_n_rep_sign <c>_n_rep_nominal (p < 0.05) <c>_n_rep_strict (p < 1e-5)`, and over the reference pairs tested in every run `<c>_n_common <c>_n_rep_sign_common <c>_n_rep_nominal_common <c>_n_rep_strict_common`; then `n_variants_tested n_traits_tested n_cells n_cells_cis n_cells_cross_chrom n_cells_same_chrom_distal n_cells_trans cis_window_bp lambda_gc n_inflation_cells n_inflation_cells_excluded observed_p1e-3 expected_p1e-3 observed_p1e-5 expected_p1e-5 hotspot_bin_bp n_bins bonferroni_p hotspot_threshold hotspot_min_genes n_hotspots n_hotspots_shared` (`n_hotspots`: the hotspot count at the run's own opening threshold and minimum; `n_hotspots_shared`: at the first run's, the same setting for every run, empty when a run's grid lacks it) |
| `hotspot_grid.tsv` | one row per run, threshold and minimum genes: `run threshold (key, e.g. p1e-5, bonferroni) label alpha log10p k n_hotspots` |
| `hotspots.tsv` | the runs' `hotspots.tsv` stacked: `run` then the run's columns (`n_genes_<threshold>` missing in a run is `NA`) |
| `hotspot_genes.tsv` | the runs' `hotspot_genes.tsv` stacked, `run` first |
| `qq.tsv` | the runs' `qq.tsv` stacked, `run` first |
| `reference_pairs.tsv` | the runs' `reference_pairs.tsv` stacked (`run` = label), plus `tested_in_all` (1 when the pair, matched by `ref_trait` and `variant`, was tested in every run) |
| `summary.json` | `schema_version kind (comparison) pipeline (trans) dataset (comparison) labels reference runs (the input columns of runs.tsv) replication_levels n_common_pairs (cis, trans) thresholds (key label alpha log10p of every run's grid, loosest first) comparability (same_genotype_samples same_reference same_rectangle same_cells same_samples (same number and same phenotype samples) same_cis_window same_bin_width reference_sha256 rectangle_sha256 n_samples rint, per label where not a flag) per_run (per label: run dataset rint n_cov n_samples n_variants n_traits n_variants_tested n_traits_tested n_cells n_cells_trans cells n_reference_pairs replication hotspots inflation params warnings, as in a run's summary.json) definitions warnings` (comparability warnings, then each run's warnings prefixed by its label) |
| `compare.json` | manifest: `kind` (`trans`), `runs` (`label run_dir dataset settings`), `labels`, `argv`, `cwd` |
| `settings/<label>.anqchor.yaml` | a copy of each run's `anqchor.resolved.yaml` |
| `report/trans.html` | the comparison report (unless `--no-report`) |

## report

| File | Content |
|---|---|
| `report/match.html` | written by `match report` or `match all` |
| `report/cis.html` | written by `cis report` or `cis all` in a run directory, and by `cis compare` (or `cis report`) in a comparison directory |
| `report/trans.html` | written by `trans report` or `trans all` in a run directory, and by `trans compare` (or `trans report`) in a comparison directory |
| `report/index.html` | only from the top-level `anqchor run` with more than one pipeline: one tile per pipeline and dataset, with links |

Reports are self-contained interactive pages (plotly; every chart has a table view and TSV
download); with the default `--plotlyjs embed` the example reports are 5-10 MB. With
`--plotlyjs cdn` and no network, each chart area shows a notice suggesting
`--plotlyjs embed`. Reports are rendered only from the summaries above, never from qpgen output.
