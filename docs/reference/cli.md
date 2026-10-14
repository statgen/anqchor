# Command line

Generated from `anqchor <command> --help` (version 0.1.0.dev0). Run any command with `--help`
for the same information. Commands are shown as run from the repository root.

```text
anqchor [--version] [-v/--verbose] [--qpgen PATH] COMMAND [ARGS]...
```

| Root option | Meaning |
|---|---|
| `--version` | show the version and exit |
| `-v`, `--verbose` | debug logging |
| `--qpgen PATH` | path to `qpgentools` (overrides discovery and `qpgen.binary`) |

| Command | Purpose |
|---|---|
| [`match run`](#anqchor-match-run), [`match report`](#anqchor-match-report), [`match all`](#anqchor-match-all) | sample identity matching |
| [`cis run`](#anqchor-cis-run), [`cis compare`](#anqchor-cis-compare), [`cis report`](#anqchor-cis-report), [`cis all`](#anqchor-cis-all) | replication of known cis pairs, and comparison of runs |
| [`trans run`](#anqchor-trans-run), [`trans compare`](#anqchor-trans-compare), [`trans report`](#anqchor-trans-report), [`trans all`](#anqchor-trans-all) | trans evaluation: replicated reference pairs, hotspots, inflation |
| [`steps`](#anqchor-steps) | list the steps of a run with their status |
| [`run`](#anqchor-run), [`report`](#anqchor-report) | several analyses from one configuration file |
| [`qpgen build`](#anqchor-qpgen-build), [`qpgen check`](#anqchor-qpgen-check) | build and check qpgen |
| [`prepare index-pvar`](#anqchor-prepare-index-pvar), [`pgen-list`](#anqchor-prepare-pgen-list), [`pheno`](#anqchor-prepare-pheno), [`cov`](#anqchor-prepare-cov) | standalone input utilities |

## Common run options

Shared by `match run|all`, `cis run|all`, `trans run|all` and the top-level `run`.

| Option | Meaning |
|---|---|
| `-c`, `--config FILE` | YAML configuration file (optional when the input flags are given) |
| `-o`, `--out PATH` | run directory (overrides `out_dir`) |
| `-j`, `--jobs INT` | steps to run in parallel (`make -j`); default 1 |
| `--threads INT` | threads for `qpgentools match-prs-pheno` (`qpgen.threads`), which only identity matching runs; cis and trans use single-threaded qpgen commands. Changing it re-runs nothing |
| `--dataset NAME` | restrict to this dataset of a multi-dataset configuration (repeatable) |
| `--strategy NAME` | restrict to this strategy of a configuration file (repeatable; `cis`, `trans` and `run` only) |
| `--dry-run` | write the Makefile, show what would run, run nothing; an existing `anqchor.resolved.yaml` is left unchanged (it is written only when the directory has none) |
| `--until STEP` | run up to and including this step |
| `--only STEP` | re-run only this step (its upstream steps must be up to date) |
| `--force` | re-run the selected steps and everything upstream |
| `--keep-going` | continue with independent steps after a failure |
| `--plotlyjs embed\|cdn` | report plotly.js: embedded (offline, default) or from a CDN (small files) |

## Input flags

Available on `match run|all`, `cis run|all` and `trans run|all`. Without `--config` they define
all inputs; with `--config` they override the file (see
[precedence](../user-guide/configuration.md#rules)). Paths are relative to the current
directory.

| Option | Meaning | Config key |
|---|---|---|
| `--geno-list FILE` | qpgen pgen list (`CHROM BEG END PREFIX`) | `genotypes.pgen_list` |
| `--pfile PREFIX` | single PLINK2 prefix with `.pgen`/`.psam`/`.pvar.idx.gz` | `genotypes.pfile` |
| `--pheno FILE` | phenotype matrix | `datasets.<name>.pheno` |
| `--pheno-format FMT` | `regenie` (default), `tensorqtl`, `tsv-sample-row`, `tsv-sample-col` | `datasets.<name>.pheno_format` |
| `--cov FILE` | covariate matrix | `datasets.<name>.cov` |
| `--cov-format FMT` | layout of `--cov` | `datasets.<name>.cov_format` |
| `--cat-cov-names NAMES` | comma-separated categorical covariates of `--cov` (e.g. `sex,batch`); their text levels are encoded by qpgen. Non-numeric covariate columns are treated as categorical automatically, so the flag is needed only for categories coded as numbers | `datasets.<name>.cov_categorical` |
| `--trait-annot FILE` | trait loci, BED-like or tensorQTL matrix (`cis`, `trans`) | `datasets.<name>.trait_annot` |
| `--trait-id-rule RULE` | `as_is`, `before_colon`, `strip_version`, `before_colon_strip_version` | `datasets.<name>.trait_id_rule` |
| `--sample-map FILE` | `geno_id pheno_id` map (`match` only) | `datasets.<name>.sample_map` |
| `--pairs FILE` | reference pairs, columns `trait variant beta se`, and `pval` when present | `references.<name>.pairs` |
| `--dataset-name NAME` | dataset name for these flags (default `dataset`) | |
| `--reference-name NAME` | reference name for these flags (default `reference`) | |

## anqchor match run

Run sample identity matching up to the summary tables (`pair-prs`, `match-prs-pheno`,
summaries). No HTML; render it with `anqchor match report --out DIR`.

```bash
anqchor match run --geno-list G --pheno P [--cov C] --pairs R [--sample-map M] --out DIR [options]
```

Common run options, input flags (no `--trait-annot`), and:

| Option | Config key | Meaning |
|---|---|---|
| `--mahalanobis` / `--independent` | `mahalanobis` | whiten correlated traits (slower) instead of scoring them as independent (default) |
| `--lambda auto\|X` | `lambda` | with `--mahalanobis`: `auto` or a value in [0, 1] |
| `--rint` / `--no-rint` | `rint` | rank-based inverse normal transform after covariate adjustment (default on) |
| `--use-cov` / `--no-use-cov` | `use_cov` | adjust phenotypes for covariates before matching (default on) |
| `--z-threshold X` | `z_threshold` | Z-score for a lenient self match (default 1.96) |
| `--z-diff X` | `z_diff` | Z-score gap that makes a new best match clear (default 2.0) |
| `--prs-samples FILE` | `prs_samples` | genotype IDs to compute PRS for (default: all) |

## anqchor match report

(Re-)render `report/match.html` from an existing run's summaries. Runs no analysis.

```bash
anqchor match report --out DIR [--plotlyjs embed|cdn]
```

| Option | Meaning |
|---|---|
| `-o`, `--out PATH` (required) | run directory written by `anqchor match run` |
| `--plotlyjs embed\|cdn` | plotly.js embedded or from a CDN |

## anqchor match all

`match run` followed by `match report`: analysis, summaries and the HTML report. Same options
as `match run`.

## anqchor cis run

Evaluate replication of known cis pairs for **one** set of inputs (genotypes, phenotypes,
covariates, phenotype transform, reference), up to the summary tables. No HTML: use
`cis report` or `cis compare`.

```bash
anqchor cis run --geno-list G --pheno P --pairs R [--cov C [--cov-include NAMES]] [--rint] [--trait-annot T] --out DIR [options]
```

Common run options, input flags (no `--sample-map`), and:

| Option | Config key | Meaning |
|---|---|---|
| `--cov-include NAMES` | `pipelines.cis.strategies` | use only these covariates of `--cov`: comma-separated names; ranges such as `PEER1-10` expand. Default: every covariate in `--cov`. Omit `--cov` for no covariates |
| `--rint` / `--no-rint` | `pipelines.cis.rint` | apply a rank-based inverse normal transform to the phenotypes after covariate adjustment. Default: test the phenotypes as given (one transform per run; compare runs with `cis compare`) |
| `--assoc-engine pair\|rect` | `assoc_engine` | `pair` (`pair-assoc`, default) or `rect` (`rect-assoc`) |

## anqchor cis compare

Compare finished cis runs (any inputs) and write comparison tables and the HTML report.

```bash
anqchor cis compare RUN_DIR [RUN_DIR...] --out DIR [--label L ...] [--no-report] [--plotlyjs embed|cdn]
```

| Argument / option | Meaning |
|---|---|
| `RUN_DIR...` (required) | run directories written by `anqchor cis run` |
| `-o`, `--out PATH` (required) | comparison directory to write |
| `--label STR` | label for each run, in order (repeat once per run); default: the run directory names (`<parent>_<name>` when two are the same) |
| `--report` / `--no-report` | also write `report/cis.html` (default `--report`) |
| `--plotlyjs embed\|cdn` | plotly.js embedded or from a CDN |

Prints `compared N run(s): <labels>`, the table directory and the report path. Labels must be
distinct; the number of `--label` values must equal the number of runs, a run directory may
appear only once, and `--out` must be a new directory or an earlier comparison of the same
kind (not a run directory); exit 2 otherwise.
Output files: [Output files](outputs.md#cis-compare).

## anqchor cis report

(Re-)render `report/cis.html` from existing tables. Runs no analysis.

```bash
anqchor cis report --out DIR [--plotlyjs embed|cdn]
```

| Option | Meaning |
|---|---|
| `-o`, `--out PATH` (required) | a cis run directory or a `cis compare` directory |
| `--plotlyjs embed\|cdn` | plotly.js embedded or from a CDN |

## anqchor cis all

`cis run` followed by `cis report` on the same run directory. Same options as `cis run`.

## anqchor trans run

Evaluate one set of inputs (genotypes, phenotypes, covariates, phenotype transform,
reference) over a variant x gene rectangle, up to the summary tables: replicated reference
pairs (cis and trans), hotspots and calibration. No HTML: use `trans compare` or
`trans report`.

```bash
anqchor trans run --geno-list G --pheno P --pairs R --trait-annot T [--cov C [--cov-include NAMES]] [--rint] --out DIR [options]
```

Common run options, input flags (no `--sample-map`), and:

| Option | Config key | Meaning |
|---|---|---|
| `--cov-include NAMES` | `pipelines.trans.strategies` | use only these covariates of `--cov`: comma-separated names; ranges such as `PEER1-10` expand. Default: every covariate in `--cov`. Omit `--cov` for no covariates |
| `--variants FILE` | `variants` (or a trans reference) | variant list (`CHROM:POS:REF:ALT` per line) instead of the reference's variants; without `--config` and `--pairs` it defines a variants-only reference |
| `--traits reference\|all\|FILE` | `traits` | rectangle genes (`file:` prefix optional); default `reference` |
| `--rint` / `--no-rint` | `pipelines.trans.rint` | apply a rank-based inverse normal transform to the phenotypes after covariate adjustment (`rect-assoc --rint-after-adj`). Default: test the phenotypes as given (one transform per run) |
| `--cis-window-bp INT` | `cis_window_bp` | a variant within this distance of a gene's start on the same chromosome is cis (default 1000000); every other pair is trans |
| `--min-maf X` | `min_maf` | minimum minor allele frequency |
| `--hotspot-bin-bp INT` | `hotspot_bin_bp` | bin width that thins trans signals for hotspots (default 1000000) |
| `--hotspot-min-traits INT` | `hotspot_min_traits` | genes with a trans signal that make a bin a hotspot (default 5; the report lets you change it) |

To evaluate several covariate sets or RINT, start one `trans run` per setup, each with its
own `--out`, then compare the runs with `trans compare`. A configuration file holds one
setup too (one entry under `pipelines.trans.strategies`).
`trait_chroms`, `levels` and `max_sig_rows` have no flag. Output files:
[Output files](outputs.md#trans).

## anqchor trans compare

Compare finished trans runs (any inputs) and write comparison tables and the HTML report.

```bash
anqchor trans compare RUN_DIR [RUN_DIR...] --out DIR [--label L ...] [--no-report] [--plotlyjs embed|cdn]
```

| Argument / option | Meaning |
|---|---|
| `RUN_DIR...` (required) | run directories written by `anqchor trans run` |
| `-o`, `--out PATH` (required) | comparison directory to write |
| `--label STR` | label for each run, in order (repeat once per run); default: the run directory names (`<parent>_<name>` when two are the same) |
| `--report` / `--no-report` | also write `report/trans.html` (default `--report`) |
| `--plotlyjs embed\|cdn` | plotly.js embedded or from a CDN |

Prints a `warning:` line per difference between the runs (genotype samples, reference,
samples, rectangle, tests, cis window, bin width, opening hotspot setting) and per run warning, then `compared N run(s):
<labels>`, the table directory and the report path. Labels must be distinct; the number of
`--label` values must equal the number of runs, a run directory may appear only once, and
`--out` must be a new directory or an earlier trans comparison (not a run directory or a cis
comparison); exit 2 otherwise. Output files:
[Output files](outputs.md#trans-compare).

## anqchor trans report

(Re-)render report/trans.html from existing tables. Runs no analysis.

```bash
anqchor trans report --out DIR [--plotlyjs embed|cdn]
```

| Option | Meaning |
|---|---|
| `-o`, `--out PATH` (required) | a trans run directory or a `trans compare` directory |
| `--plotlyjs embed\|cdn` | plotly.js embedded or from a CDN |

A `cis compare` directory is refused (exit 2) with a pointer to `cis report`, and
`cis report` refuses a `trans compare` directory the same way.

## anqchor trans all

`trans run` followed by `trans report` on the same run directory (one run, no comparison).
Same options as `trans run`.

## anqchor steps

List the steps of a run with their status (`done`, `stale`, `missing`) and descriptions.

```bash
anqchor steps --config FILE [--out PATH] [--pipeline P]... [--dataset D]... [--strategy S]...
```

`--config` is required; a run's own `anqchor.resolved.yaml` works. The command writes the
Makefile but runs nothing, and leaves the run's `anqchor.resolved.yaml` unchanged. It always
lists the report step, so after `match run`, `cis run` or `trans run` the report step shows as `missing` or
`stale`.

## anqchor run

Run every pipeline enabled in a configuration file, sharing common steps.

```bash
anqchor run --config FILE [--pipeline match|cis|trans]... [common run options]
```

`--config` is required and input flags are not accepted. `-p`, `--pipeline` (repeatable)
selects pipelines regardless of `enabled`. Reports are rendered as part of the run, and
`report/index.html` links them when more than one pipeline runs.

## anqchor report

Re-render HTML reports of an `anqchor run` directory from existing summaries; runs missing
upstream steps if needed.

```bash
anqchor report --config FILE [--out PATH] [--pipeline P]... [--plotlyjs embed|cdn] [--force]
```

`--force` re-renders even if up to date. Without `--pipeline` it renders every pipeline
enabled in the file. For a single analysis, `anqchor match report --out DIR`, `anqchor cis report --out DIR`
and `anqchor trans report --out DIR` are simpler.

## anqchor qpgen build

Build `qpgentools` from the submodule (htslib, qgenlib, then qpgen).

| Option | Meaning |
|---|---|
| `--src PATH` | qpgen checkout (default: the `extern/qpgen` submodule) |
| `-j`, `--jobs INT` | parallel compile jobs (default 4) |
| `--force` | rebuild htslib and qgenlib too |
| `--cmake-arg ARG` | extra argument for cmake (repeatable) |

Prints the binary path. Exit 4 on failure, with the tail of `extern/qpgen/anqchor-build.log`.

## anqchor qpgen check

Locate `qpgentools`, verify its subcommands (`pair-prs`, `match-prs-pheno`, `pair-assoc`,
`rect-assoc`), and print its path, discovery source, sha256 and git commit. It does not check
the qpgen version; run `qpgentools --version` (0.2.0 or later is required).

| Option | Meaning |
|---|---|
| `--qpgen PATH` | binary to check (default: discovery order) |

## anqchor prepare index-pvar

Build `PREFIX.pvar.idx.gz` (9 columns with a 1-based `INDEX`) and its tabix index.

| Option | Meaning |
|---|---|
| `--pfile PATH` (required) | PLINK2 prefix; reads `PREFIX.pvar`, `PREFIX.pvar.zst` or `PREFIX.pvar.gz` |
| `--force` | overwrite an existing `PREFIX.pvar.idx.gz` |

## anqchor prepare pgen-list

Write a qpgen pgen list (`CHROM BEG END PREFIX`) with absolute prefixes. With checking (the
default), each prefix must have `.pgen`, `.psam`, `.pvar.idx.gz` and `.tbi`, and `CHROM` is
taken from the first pvar record so it matches the pvar spelling.

| Option | Meaning |
|---|---|
| `--prefix-template STR` (required) | prefix with `{chrom}`, e.g. `geno/chr{chrom}` |
| `--chroms STR` (required) | chromosomes, e.g. `1-22,X` |
| `--out PATH` (required) | output pgen list |
| `--no-check` | do not check that the files exist |

## anqchor prepare pheno

Convert a phenotype matrix to regenie layout (`FID IID trait…`).

| Option | Meaning |
|---|---|
| `--in FILE` (required) | input matrix |
| `--format FMT` (required) | `regenie`, `tensorqtl`, `tsv-sample-row`, `tsv-sample-col` |
| `--out PATH` (required) | output regenie matrix (`.tsv.gz`) |
| `--missing-str STR` | string meaning "missing" (repeatable; default `NA`) |

## anqchor prepare cov

Convert covariates to regenie layout: validate values, drop constant columns, handle
missing values. Categorical covariates keep their text levels by default; the command
prints the `--cat-cov-names` value to pass to qpgen with the file.

| Option | Meaning |
|---|---|
| `--in FILE` (required) | input covariate matrix |
| `--format FMT` (required) | layout |
| `--out PATH` (required) | output regenie covariates (`.tsv.gz`) |
| `--categorical NAME` | covariate to treat as categorical (repeatable) |
| `--missing error\|drop_samples\|impute_mean` | missing-value policy (default `error`) |
| `--missing-str STR` | string meaning "missing" (repeatable; default `NA`) |
| `--json PATH` | also write the encoding summary as JSON |
| `--encode` / `--keep-categorical` | write categorical covariates as numeric indicator columns (`--encode`) instead of their text levels (default `--keep-categorical`; qpgen encodes them with `--cat-cov-names`) |

```bash
anqchor prepare cov --in examples/mage/data/pheno/mage_covariates.regenie.tsv.gz \
    --format regenie --out results/cov.regenie.tsv.gz
```

```text
wrote results/cov.regenie.tsv.gz: 729 samples x 66 covariates; categorical: sex ['XX', 'XY']
categorical covariates keep their text levels; use --cat-cov-names sex
```

With `--encode`, `sex` is written as one column coded 1 (`XX`) / 2 (`XY`), and a covariate
with k > 2 levels as k − 1 indicator columns `<name>_<level>` coded 2 (has the level) or 1;
the first sorted level is the reference. Both forms give the same association results.

## Exit codes

0 success · 2 configuration/validation error (or run directory locked by another command) ·
3 a step failed · 4 qpgen missing or unusable.
See [Running](../user-guide/running.md#exit-codes).
