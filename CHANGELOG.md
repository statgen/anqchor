# Changelog

All notable changes to anqchor are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased] - 0.1.0.dev0

First implementation of the design in `DESIGN.md` (milestones M0-M3 and most of M4), revised
after the first review (DESIGN.md section 11.2).

### Changed after the first review

- **`trans compare`** (DESIGN 11.2-R). Evaluate trans setups like cis: one `trans run` per
  setup, then `anqchor trans compare RUN_DIR... --out DIR [--label ...] [--no-report]`
  (about 1 s, no qpgen), and `trans report --out DIR` to re-render the page of a run or of a
  comparison. The comparison writes `runs.tsv` (inputs, cells, replicated reference pairs
  per class over each run's pairs and over the pairs tested in every run, λ_GC, observed and
  expected counts, hotspots), `hotspot_grid.tsv`, the runs' `hotspots.tsv`,
  `hotspot_genes.tsv`, `qq.tsv` and `reference_pairs.tsv` stacked with a `run` column,
  `summary.json` (`kind: comparison`, comparability checks, each run's summary),
  `compare.json` and `settings/<label>.anqchor.yaml`, and warns when the runs differ in
  genotype samples, reference, rectangle, cells, cis window or bin width. The comparison
  page (*Comparing N runs: trans-QTL evaluation with anQChor*) shows replicated pairs per
  run (cis and trans side by side), hotspot counts per run under shared threshold and
  minimum controls with the bins, grid and table of a selected run, and λ_GC per run, one
  QQ line per run and observed vs expected per run. `cis report` and `trans report` refuse
  each other's comparison directories. The label, settings-copy and manifest code is shared
  with `cis compare` (`steps/compare.py`); comparison file names are in `Layout`.
  Both comparisons refuse an `--out` that is a run directory or the other kind of
  comparison and a run directory given twice, and clean up a re-used comparison directory;
  `runs.tsv` of trans also has `n_hotspots_shared`.
- **trans: one run per setup, report redesigned** (DESIGN 11.2-Q). `trans run` now stops at
  the summary tables; new `trans report --out DIR` (re-renders `report/trans.html` from the
  tables in about a second, no qpgen) and `trans all` (both). A run tests one phenotype
  transform: as given, or with the new `--rint` (`pipelines.trans.rint`,
  `rect-assoc --rint-after-adj`); `transforms` under `pipelines.trans.strategies`, more than
  one trans strategy, and `pipelines.trans.levels.fdr` are now errors that say what to write
  instead. New `--hotspot-bin-bp` (`hotspot_bin_bp`, default 1 Mb); `hotspot_min_traits`
  now means genes per bin. Summaries (schema 2): `summary.json`, `reference_pairs.tsv`
  (every reference pair with its cis/trans class and replication at same sign, p < 0.05,
  p < 1e-5), `hotspots.tsv` and `hotspot_genes.tsv` (1-Mb bins: distinct genes with a trans
  signal, genes in cis with the bin not counted, over a grid of thresholds × minimum genes),
  `qq.tsv` (cross-chromosome cells with a 95% band), `variants.tsv`, `pairs.sig.tsv.gz`;
  `strategies.tsv`, `traits.tsv`, `known_pairs.tsv` and the recommendation are gone. The
  inflation cells leave out the regions of the reference's trans pairs. The report, *trans-QTL
  evaluation with anQChor*, shows replicated reference pairs (cis and trans), hotspots (with
  threshold and minimum-genes controls, bins along the genome and the grid) and inflation
  (QQ plot, λ_GC, observed vs expected). `examples/mage/trans.yaml` is the all-covariate run
  on the GTEx reference (`gtex_lcl`). On the example: λ_GC 1.74 without covariates (37
  hotspots at p < 1e-5, ≥ 5 genes), 1.00-1.02 with covariates (none).
- **Independent example reference** (DESIGN 11.2-P). The example's known cis-eQTLs are no
  longer the MAGE study's own top pairs (discovered in the same samples, which inflated
  identity matches and favoured the all-covariate cis setup by construction) but GTEx v10
  EBV-transformed lymphocyte eGenes with their sentinel variant:
  `examples/mage/data/reference/gtex_lcl_eqtl.chr21_22.tsv.gz` (replaces
  `mage_top_cis_eqtl.chr21_22.tsv.gz`), 367 pairs and 346 variants on chr21/22 (q ≤ 0.05,
  variant in the genotypes with MAF ≥ 5%, gene in the MAGE expression). The expression and
  gene loci now hold those 367 genes. The build script takes the GTEx table:
  `make_compact_data.sh FULL_DATA_DIR GTEX_EGENES [PLINK2]`. `examples/mage/match.yaml` and
  `cis.yaml` name the reference `gtex_lcl`. New example results: match 713 strict, 14
  lenient, 2 unclear (`NA11919` fails its own genotype clearly but has no clear
  alternative); cis at p < 1e-5 of 367 pairs: 178 (no covariates), 162 (sex, PC1-5), 231
  (+ PEER1-10), 251 (all 66). Documentation, `expected_mage.json` and the tests follow.
- **cis: one phenotype transform per run** (DESIGN 11.2-N). A cis run no longer evaluates
  raw and RINT automatically: it tests the phenotypes as given, or with the new flag
  `--rint` on `cis run` / `cis all` (config `pipelines.cis.rint: true`) the rank-based
  inverse normal transform after covariate adjustment (qpgen `pair-assoc` without
  `--skip-rint`, `rect-assoc --rint-after-adj`). Raw and RINT are compared as two runs with
  `cis compare`, which accepts runs that differ in the transform. cis strategies lose the
  `transforms` key; a configuration that still lists it is rejected with a message saying
  what to write instead (trans followed: see below).
- **cis summaries, schema 3** (2 had `z` as the normal quantile of the p-value; a new
  schema re-runs the summarize step). The transform dimension is gone: rows are labelled by
  `run` (= `label`): the run directory name for a single-strategy run, the run label in a
  comparison. `pairs.long.tsv.gz` gains `z_ref` (`beta_ref / se_ref`), `log10p_ref`
  (from the reference `pval`, read from the text so p-values below the double range keep
  their size; from `z_ref` when there is no p-value) and `z` (`beta / se`, aligned with `z_ref`; the signed normal
  quantile of the p-value only without `se`); `strategies.tsv` gains `rint`, `mean_z2`, `mean_z2_common` and `rel_neff`
  (effective sample size relative to the first run, (mean z² − 1) / (baseline mean z² − 1)
  over the pairs tested in every run); `traits.tsv` has `best_run`; `power_curve.tsv` has
  `run`; `summary.json` has `rint`, `runs`, `per_run`, `unresolved_by_run`, and
  `recommended.run`; `runs.tsv` of a comparison and its `comparability` record `rint` per
  run. `cis compare` refuses runs with another summary schema (re-run them) and re-derives
  `z` and `mean_z2` from each run's `beta` and `se`. `power_curve.tsv` counts p < 10^−t,
  as the levels (was p ≤). The MAGE
  example adds `peer60_rint`; on it the relative N_eff (vs nocov) is 0.83 (basic), 1.58
  (peer10), 2.27 (peer60) and 2.10 (peer60_rint) with the GTEx reference.
- **cis report redesigned** (DESIGN 11.2-O; single run and comparison, match report style):
  title "cis-QTL replication with anQChor" ("Comparing N runs: …"), then the power curve
  first (log(1 + t) axis labelled in −log10 p), signal strength (the effective sample
  size relative to the baseline as dots with 95% bootstrap intervals on a log axis, then
  |z| / −log10 p boxes per run over the pairs tested in every run, the reference hollow),
  significant hits per run (rows; same sign, p < 0.05, p < 1e-5), agreement with the
  reference (effect size, z, −log10 p with independent axes and a fit through 0, and z
  against the baseline with y = x, for a selected run; open diamonds for opposite signs),
  all pairs (100 a page, header and pair column pinned, opposite signs marked per run) and
  collapsed run details (a link to `strategies.tsv` instead of its raw columns). The "Raw vs RINT" section, the
  recommendation paragraph and the word "strategy" are gone. A summary warning about
  unmapped pairs says "the phenotypes" instead of the dataset's internal name.
- **Match report, round two: samples without an ID, duplicates, the table.**
  - One set of result names everywhere (tiles, cards, chart legend, filters, table):
    "Matches another sample", "ID already used", "Own ID too weak" (was "No match" /
    "Matches nothing"), "Matches several", "Matches two genotype IDs" (violet, its own
    glyph), "New ID found", "Lenient match", "Strict match", "No clear match", "Not scored"
    (with its own tile). The "Samples examined" tile is gone (the header has the count);
    tile numbers come first so they line up. Slider info icons are replaced by one-line
    hints. Own-ID-too-weak cards give the best genotype's lead and the lead needed.
  - A sample matching two genotype IDs, one of them its own (or with no ID), is not a
    mismatch: no proposed ID, not a rescue candidate, left out of the verdict's mismatches,
    `n_fail`, `failing` and `swaps.tsv` (new `n_fail_two_ids`); duplicates alone make the
    verdict a `warning`.
  - Table headers are short and on one line (own score, own rank, best other, its score,
    lead, suggested ID); rows without an assigned ID show the best genotype's lead in
    italics; lenient rows are no longer red.
  - A mismatch is never shown as a duplicate. "Possible duplicates" is replaced by
    **Duplicates**: samples that match two genotype IDs about equally (the top two stand
    clear of the third by `z_diff`, the second scores above `z_threshold` and at least 0.8 ×
    the first), most likely one person under two genotype IDs. The relatives fold is gone;
    an empty section is one line. Such a sample shows as "Matches two genotype IDs" (new
    tile when present) unless its own ID ranks first; qpgen's status is unchanged.
  - Samples without an assigned ID are counted ("No assigned ID" tile: new ID found, ID
    already used, unclear), get cards (one summary card when more than six), a table filter,
    and appear in the chart, which now plots every scored sample by its best and runner-up
    match score, with the clear-lead line and a grey "one genotype clearly best" region. Samples
    without an ID are hollow markers.
  - The sample table scrolls with the page (sticky header, no inner scroll box), has fixed
    column widths and right-aligned numbers, and plain headers (the lead column replaces
    "Margin"; no info icons). Phones show sample (best other genotype under it), result, own
    score and lead without sideways scrolling.
  - Summaries: new `duplicates.tsv` and `n_duplicates` (samples there are left out of
    `related.tsv`, and the verdict names them); `NOT_SCORED` samples keep their assigned ID
    (the summarize step now reads the prepared sample map) and are a `summary.json` warning
    ("N samples not scored: no covariates"), also named in the verdict; an unmapped sample
    whose clear best genotype is another sample's confirmed match is no longer a rescue
    candidate (report: "ID already used"); notes no longer say "duplicate or contaminated
    library" or "pooled or contaminated sample".
- **Match report redesigned** ("Sample identity matching with anQChor"). A dashboard of
  large tiles (examined, with an assigned ID, strict match, lenient match, matches a
  different sample, fits several samples when non-zero, matches nothing), two sliders
  directly below it that reclassify every sample in the browser, one card per sample to
  check, a single scatter chart, "Possible duplicates" (two samples with the same best
  genotype, or own and one other genotype scoring alike and both clearly ahead; the
  relative-level pairs of `related.tsv` fold into one collapsed line), a sortable sample
  table (500 rows per page, 100 on phones) and collapsed run details. Plain result names
  replace PASS/WARN/FAIL/NO_TEST and the qpgen codes, which remain in tooltips, the TSV
  download and a code table; qpgen's Z-score is called the "match score". The page shows the
  scoring mode and, for independent scoring, a tip that `--mahalanobis` matches slightly
  better. Summaries are unchanged. Shared report changes: reports link to the index page
  only in runs that write one (new `has_index` report-step parameter), larger type scale,
  compact theme toggle, tables default to 500 rows per page with a 100 | 500 | All selector.
- **Categorical covariates are encoded by qpgen.** qpgen (submodule at `5b811e2`, still
  `qpgentools 0.2.0`) accepts `--cat-cov-names` in `match-prs-pheno`, `pair-assoc`,
  `rect-assoc` and `region-assoc` (REGENIE `--catCovarList`: L levels → L − 1 indicators, the
  first sorted level is the reference, `--max-cat-levels` default 10, a missing level is a
  missing covariate). anqchor no longer recodes categorical covariates:
  `prepare/datasets/<ds>/cov.regenie.tsv.gz` and the strategies' `cov.regenie.tsv.gz` keep
  their text levels (`sex` stays `XX`/`XY`), and every qpgen call gets `--cat-cov-names` with
  the categorical covariates it uses (`--cov-include sex,PC1-5` passes `sex`,
  `--cov-include PC1-5,PEER1-10` none). Categorical = non-numeric columns (detected when the
  Makefile is built) plus the names given with the new flag `--cat-cov-names` on `match`,
  `cis` and `trans` (config `datasets.<ds>.cov_categorical`), for categories coded as
  numbers. `cov.json` gains `raw_columns` and `categorical_names`; `strategy.json` gains
  `cat_cov_names`. `anqchor prepare cov` keeps text levels by default and prints the
  `--cat-cov-names` value; `--encode` writes numeric indicator columns as before. A
  categorical covariate with more than 10 levels also gets `--max-cat-levels N` (qpgen's
  default cap is 10, anqchor's own check allows more). Results are
  identical to the numeric encoding (match outputs and cis statistics byte-identical on the
  example).
- **Separate run, report and all commands.** `match run` and `cis run` stop at the summary
  tables and write no HTML. New `match report --out RUN` and `cis report --out DIR` re-render
  the HTML from existing tables in about a second without running qpgen; `match all` and
  `cis all` run and then render.
- **cis evaluates one set of inputs per run.** `cis run` takes the genotypes, phenotypes,
  covariates and reference of one setup from flags, with `--cov-include` to use a subset of
  the covariate columns (comma-separated; ranges such as `PEER1-10` expand) and no strategy
  name. Runs are independent and can run in parallel.
- **Flags first.** All documented commands work from flags alone; `--config` is optional and
  flags override it.
- **No `report/index.html` for single-analysis runs**; it is written only when one invocation
  runs several pipelines (`anqchor run`).
- **qpgen 0.2.0 or later** (`qpgentools --version`; submodule first at `76f2fd1`, now at
  `5b811e2` for categorical covariates). It fixes, among
  others, `pair-assoc` without covariates (no centering), missing values read as 0, sample
  misalignment with real FIDs in the psam, and unstable residuals with near-duplicate
  covariates (see `extern/qpgen/dev/review-2026-10-anqchor.md`).
- **Example and tests use compact in-repo data.** `examples/mage/data/` (32 MB): 729 MAGE
  individuals × 256,645 variants with MAF ≥ 5% on chr21/22, 367 genes, 66 covariates and
  367 known cis-eQTL pairs from GTEx, rebuilt by `examples/mage/make_compact_data.sh`. Integration tests
  need no external data; `tests/integration/expected_mage.json` holds the numbers on this data.
- **Documentation** reorganised around three separate analyses, with flags-first commands,
  the new example numbers, a use-case overview and an example-data page; report screenshots
  removed until the report redesign; diagrams are short vertical mermaid charts.

### Added after the first review

- `anqchor cis compare RUN_DIR... --out DIR [--label L ...] [--no-report]`: compares finished
  cis runs with any inputs. Writes `runs.tsv`, `strategies.tsv` (one row per run × transform,
  with `n_rep_<level>` over each run's tested pairs and `n_rep_<level>_common` over the pairs
  tested in every run), `pairs.long.tsv.gz`, `traits.tsv`, `power_curve.tsv`, `summary.json`,
  `compare.json`, `settings/<label>.anqchor.yaml` (each run's exact settings) and
  `report/cis.html`.
- `--cov-include` for `cis run` / `cis all` and `trans run`; `--hotspot-min-traits` for `trans run`.
- Reproducing a run: `anqchor <analysis> run --config RUN/anqchor.resolved.yaml --out NEW`
  (the resolved file holds every setting, including flags).
- Optional per-analysis example configurations `examples/mage/match.yaml`, `cis.yaml`,
  `trans.yaml`.

### Fixed after the first review

- `--pairs` (and a configuration without `columns.pval`) reads a column named `pval` when the
  file has one, as `--help` says (`trait variant beta se [pval]`); flags-only runs reported
  `pval_ref` as `NA`, and comparing them with configuration runs of the same file warned
  about different reference pairs.
- `--dry-run` and `anqchor steps` no longer overwrite an existing `anqchor.resolved.yaml`, so
  the file keeps describing the settings of the outputs (and reproduces them).
- Top-level `anqchor report` on a configuration with one pipeline no longer fails asking for
  `report.index`.
- `match run` / `cis run` say when an existing report is older than the new results (instead
  of printing its path as current).
- `cis compare`: `n_common_tested` and `n_rep_<level>_common` are integers in `summary.json`.
- `--threads` help says that only identity matching uses it; the cis `transform_semantics`
  refer to qpgen 0.2.0 instead of development commits; the `trait_chroms` error no longer
  prints an empty `(e.g. )`.

### Removed after the first review

- `examples/mage/config.yaml` (the all-in-one example running three pipelines together) and
  the dependency of tests and examples on `tmp/mage_examples_data`.

### Not changed yet

- `trans run` still takes its covariate strategies from the configuration file and writes
  `report/trans.html` within the run; it will be restructured like cis.

### Added

- **CLI** (`anqchor`, Typer): `run`, `steps`, `report`, `match run`, `cis run`, `trans run`,
  `qpgen build|check`, `prepare index-pvar|pgen-list|pheno|cov`, and the hidden `_step`.
  Configuration file and/or flags (flags override the file); exit codes 2 (configuration),
  3 (step failed), 4 (qpgen missing).
- **Configuration** (pydantic v2, unknown keys rejected): genotypes, references, datasets,
  pipelines, report; resolved config written to every run directory.
- **Execution engine**: steps with declared inputs, outputs and parameters, written as a GNU
  make 3.81-compatible `Makefile`; step signatures re-run exactly the affected steps when a
  parameter changes; deleted outputs re-run their step and everything downstream;
  `--until`, `--only`, `--force`, `--dry-run`, `-j`, `--keep-going`; per-step logs;
  `run.json` provenance with qpgen sha256 and commit and per-step times.
- **qpgen integration**: submodule `extern/qpgen` (branch `dev`), `anqchor qpgen build`,
  binary discovery (`--qpgen`, `$ANQCHOR_QPGEN`, `PATH`, submodule), typed command builders and
  schema-checking output readers with fallbacks for older builds.
- **Input preparation**: pgen lists (4 or 6 columns, relative prefixes resolved near the list,
  written as 4-column lists with absolute prefixes, symlinks for non-standard names);
  chromosome style detection and normalisation; indexed-pvar creation; phenotype conversion
  from regenie, tensorQTL, tsv-sample-row and tsv-sample-col layouts; covariate encoding
  (categorical detection, constant columns dropped, missing-value policies; categorical
  covariates were first recoded to 1/2 codes and are now passed to qpgen as text, see above); trait
  annotation from BED-like files or tensorQTL matrices with automatic ID-rule matching;
  reference normalisation and presence check against the pvar with an SNV/indel split; trait-ID
  rules and trait maps; sample maps.
- **match pipeline**: `pair-prs` over all genotyped samples, `match-prs-pheno` with
  independence scoring by default (`mahalanobis: true` opt-in with auto/fallback/fixed λ),
  per-sample statuses and `PASS/WARN/FAIL/NO_TEST` groups, reciprocal and one-way swaps, rescue
  candidates, possible relatives (largest Z drop in the top five), trait weights; report with
  in-browser recomputation of statuses for user-chosen `z_threshold` / `z_diff` and
  "copy as config".
- **cis pipeline**: strategies (default, `include` subsets with ranges, separate files, none,
  phenotype overrides; raw and RINT), `pair` and `rect` engines, replication levels `sign`,
  `nominal`, `fdr`, `bonferroni`, `strict`, unresolved-pair check, per-trait table, power
  curves, π1, recommended strategy with a rationale; report.
- **trans pipeline**: variant × trait rectangle (`variants`, `traits`, `trait_chroms`), cell
  classes (known, other cis, same-chromosome distal, cross-chromosome, unknown), streaming
  summarizer with bounded memory, λ_GC, QQ, hotspots, significant pairs, known-pair replication,
  recommendation; report.
- **Reports**: self-contained HTML (plotly embedded or from CDN) with light/dark themes, table
  views and TSV downloads; an index page. Offline with the CDN, each chart area shows a notice;
  status dots carry a screen-reader status word.
- **Re-run safety**: step signatures include the `qpgentools` sha256 (qpgen steps) and the
  anqchor version (Python steps), so a rebuilt qpgen or an upgraded anqchor re-runs the affected
  steps; a step deletes its stamp before running, so a failed `--only`/`--force` re-run leaves
  no stale stamp; `--only` refuses when an upstream step is missing or out of date; an exclusive
  lock (`<run>/.anqchor.lock`) makes a second command on the same run directory exit with code 2;
  `--threads` reaches the step through `$ANQCHOR_THREADS` and is in no signature.
- **Selection flags**: with several pipelines, `--dataset` / `--strategy` skip the pipelines
  left with no matching dataset or strategy (an error only when none is left).
- **Input validation**: compression detected from file contents; IDs with whitespace or
  quotes rejected; infinite phenotype or covariate values rejected; a mostly numeric covariate
  with a few non-numeric cells, a non-numeric covariate with more than 20 levels (unless listed
  in `cov_categorical`) and colliding categorical level names are errors; empty cells are
  missing values; headerless trait annotations with numeric trait IDs need a header;
  nucleotide alleles upper-cased; samples dropped for missing covariates and `sample_map`
  ID mismatches reported in `dataset.json` (`n_samples_analysable`, `warnings`) and in the
  `cis`/`trans` summaries.
- **Statistics details**: one shared threshold rule for `cis` and `trans` (p < α for fixed
  levels and Bonferroni, BH q ≤ α for FDR); cis levels ranked by their actual −log10 p
  cut-off (`level_order`), with the recommendation keyed on the strongest level; rows without
  a finite statistic are `not_tested`; `reciprocal-weak` swaps (neither side a clear
  mismatch) listed but not counted; no `proposed_geno_id` when the best genotype is another
  sample's confirmed match.
- **Tests**: unit tests without qpgen; MAGE chr21/22 integration tests against
  `tests/integration/expected_mage.json`.
- **Documentation**: mkdocs-material site, README, example `examples/mage/`.

### Depends on qpgen `dev` changes

- `89eab88`: `pair-assoc` raw `LOG10P` column, `pair-assoc` lookup re-seek, `rect-assoc` `N`,
  tensorQTL header case.
- `37b2cb2`: 6-column pgen lists read as `PGEN PSAM PIVAR`.
- `f422acd`: `pair-assoc` residualises genotypes on the covariates (both cis engines agree).
- `938bf19`: residual df `N − 2 − k` with covariates in all association tools. On the full
  MAGE chr21/22 data (before the compact example): `peer60/raw` 511 known pairs at p < 1e-5
  (was 549 with df `N − 2`), 593 at Bonferroni; trans λ_GC ≈ 1.01.
- `b8df9c9`, `76f2fd1` (qpgen 0.2.0): the fixes from the anqchor review and
  `qpgentools --version`.
- `cb97eb0`, `5b811e2` (still 0.2.0): categorical covariates (`--cat-cov-names`) and named
  covariates (`--cov-names`) in the association tools and `match-prs-pheno`; anqchor uses
  `--cat-cov-names`.
