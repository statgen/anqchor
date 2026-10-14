# Running

Each analysis runs on its own, into its own **run directory** (`--out`). This page explains the
commands, what a run directory holds, how to re-run, reproduce and compare runs, and what to
do when something fails.

## Commands

| Command | What it does |
|---|---|
| `anqchor match run` | identity matching up to the summary tables (no HTML) |
| `anqchor match report --out DIR` | renders `DIR/report/match.html` from the tables; runs no analysis |
| `anqchor match all` | `match run`, then `match report` |
| `anqchor cis run` | replication of known cis pairs for one set of inputs, up to the summary tables |
| `anqchor cis compare RUN_DIR... --out DIR` | puts finished cis runs side by side; writes comparison tables and `report/cis.html` |
| `anqchor cis report --out DIR` | renders `report/cis.html` of a cis run or of a comparison; runs no analysis |
| `anqchor cis all` | `cis run`, then `cis report` on that run |
| `anqchor trans run` | replicated reference pairs (cis and trans), hotspots and inflation for one set of inputs, up to the summary tables |
| `anqchor trans compare RUN_DIR... --out DIR` | puts finished trans runs side by side; writes comparison tables and `report/trans.html` |
| `anqchor trans report --out DIR` | renders `report/trans.html` of a trans run or of a comparison; runs no analysis |
| `anqchor trans all` | `trans run`, then `trans report` on that run |
| `anqchor steps` | lists the steps of a run with their status; runs nothing |

The top-level `anqchor run` and `anqchor report` (several analyses from one configuration
file) are described in the [command-line reference](../reference/cli.md#anqchor-run). This
guide always runs one analysis at a time.

### run, report, all

`run` does the work: it checks the inputs, runs qpgen and writes summary tables (TSV and JSON).
`report` reads only those tables and writes the HTML report. It takes about a second, never
runs qpgen, and can be repeated as often as you like, for example after upgrading anqchor or
with `--plotlyjs cdn` for a smaller file. `all` does both in one command.

```bash
D=examples/mage/data
IN=(--geno-list $D/geno/mage.pgen_list.tsv --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz
    --cov $D/pheno/mage_covariates.regenie.tsv.gz --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz)
anqchor match run "${IN[@]}" --out results/match
anqchor match report --out results/match
```

(`IN` is a shell array, so the commands work in both bash and zsh.)

`trans run`, `trans compare`, `trans report` and `trans all` work like the cis commands below.

### Flags or a configuration file

Every `run`/`all` command takes its inputs as flags (`--geno-list`, `--pheno`, `--cov`,
`--pairs`, ...), which is the usual way. A [configuration file](configuration.md)
(`--config FILE`) is optional; given together, flags override the file. Without `--config`,
the dataset is called `dataset` and the reference `reference` (`--dataset-name`,
`--reference-name` change that); these names appear in step names and paths such as
`match/dataset/`.

## Evaluating alternatives: one run per setup

A cis run evaluates **one** set of inputs. To find out which covariates, normalisation,
genotype set or reference works best, run each alternative into its own directory and compare
the finished runs:

```bash
D=examples/mage/data
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --out results/cis/nocov
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz --cov-include sex,PC1-5,PEER1-10 \
  --out results/cis/peer10
anqchor cis compare results/cis/nocov results/cis/peer10 --out results/cis/compare
```

```text
compared 2 run(s): nocov, peer10
tables: results/cis/compare
cis report: results/cis/compare/report/cis.html
```

Each run takes about 6 seconds on the example data; the comparison about a second. Runs are
labelled by their directory names; `--label` (once per run, in order) chooses other labels.
`--no-report` writes the tables only, and `anqchor cis report --out results/cis/compare`
renders the report later. The runs may differ in any input; the comparison also counts
replication over the pairs tested in every run (the `_common` columns of `strategies.tsv`), so
runs with different genotypes or references remain comparable. The
[cis use case](../use-cases/cis-strategies.md) walks through a full comparison.

trans works the same way: one `trans run` per setup (`--cov-include`, `--rint`), each into
its own directory
([the commands on the example](../use-cases/index.md#the-commands-on-the-example-data)),
then `trans compare RUN_DIR... --out DIR` puts them side by side and writes
`report/trans.html`; `trans report --out DIR` re-renders it. The
[trans use case](../use-cases/trans-strategies.md) walks through a full comparison.

### Parallel runs

Runs never share a directory, so independent runs can execute at the same time, for example
in separate terminals, as background jobs (`&` and `wait`) or as separate cluster jobs. Two
commands on the **same** run directory cannot: the second one stops at once with
`error: another anqchor command is running in DIR` (exit code 2).

Within a run, `-j N` lets make run up to N independent steps at once. `--threads N` gives
`match-prs-pheno` N threads; the other qpgen commands are single-threaded.

## How a run works

1. The configuration is built from the file and the flags, and validated; input files must
   exist (exit 2 otherwise). `qpgentools` is located (exit 4 if missing).
2. The complete settings are written to `<run>/anqchor.resolved.yaml` (not by a dry run or
   `anqchor steps` when the file exists).
3. anqchor builds the run's steps, writes one **signature file** per step
   (`steps/<name>.json`: inputs, outputs, commands and parameters) and a `Makefile`.
4. `make` runs what is out of date. Progress lines look like `[anqchor] running <step>`;
   each step's output goes to `logs/<step>.log`.
5. `run.json` is updated with provenance, and `ran: ...` lists the steps that ran.

The generated `Makefile` is plain GNU make (3.81 compatible). Each step's target is a stamp
`steps/<name>.done`, touched after the step succeeds and all its outputs exist.

## Re-runs and invalidation

Running the same command again does nothing (`ran: nothing (all steps up to date)`). A step
runs again when:

- **a parameter changed**: its signature changes, and anqchor removes the stamps of that step
  and everything downstream. Changing `--z-threshold` of `match run`, for example, re-runs
  `match.<ds>.qpgen` and `match.<ds>.summarize` and nothing else;
- **an input file is newer** than the step's stamp (make's usual rule);
- **one of its outputs was deleted**: the step and everything downstream run again;
- **qpgen was rebuilt or anqchor upgraded**: steps that run qpgen record the binary's sha256,
  and Python steps the anqchor version, so the affected steps run again.

`--threads` does not change results and changing it re-runs nothing.

To preview what a change would re-run, add `--dry-run` to the command with the changed
setting. On the identity-matching run of the [Quickstart](../getting-started/quickstart.md):

```bash
anqchor match run \
  --geno-list examples/mage/data/geno/mage.pgen_list.tsv \
  --pheno examples/mage/data/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --cov examples/mage/data/pheno/mage_covariates.regenie.tsv.gz \
  --pairs examples/mage/data/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --out results/match --z-threshold 2.5 --dry-run
# (first the recipes make would run, then)
# dry run: would run: match.dataset.qpgen, match.dataset.summarize
# dry run: Makefile written to /…/results/match/Makefile
```

## Partial runs

| Option | Effect |
|---|---|
| `--dry-run` | write the `Makefile`, print which steps would run, run nothing; existing signatures, stamps and `anqchor.resolved.yaml` are not changed |
| `--until STEP` | run `STEP` and everything it needs, nothing downstream |
| `--only STEP` | re-run exactly `STEP`; every upstream step must be up to date (`done`), otherwise exit 2 |
| `--force` | re-run the requested steps **and everything upstream** |
| `-j N` / `--jobs N` | run up to N steps in parallel |
| `--keep-going` | continue with independent steps after a failure |
| `--dataset NAME` | with a multi-dataset configuration file: restrict to these datasets (repeatable) |
| `--strategy NAME` | with a configuration file that lists several strategies (`cis`, `trans`): restrict to these (repeatable) |

`--until` and `--only` cannot be combined. Step names are those printed by `anqchor steps`.
It needs a configuration file; the run's own `anqchor.resolved.yaml` works:

```console
$ anqchor steps --config results/cis/peer10/anqchor.resolved.yaml --out results/cis/peer10
prepare.genotypes                  done     validate genotype files; write an absolute pgen list and the sample list
prepare.reference.reference        done     normalise reference reference and check which variants exist in the genotypes
prepare.dataset.dataset            done     convert dataset dataset to regenie layout, encode covariates, annotate traits
prepare.mapping.dataset.reference  done     map reference reference traits onto dataset dataset
cis.dataset.default.inputs         done     covariates/phenotypes for cis strategy default on dataset
cis.dataset.default.assoc          done     test the reference pairs on dataset with strategy default (qpgentools pair-assoc)
cis.dataset.summarize              done     replication of reference pairs on dataset by strategy
report.cis                         missing  cis report (report/cis.html)

run directory: /…/anqchor/results/cis/peer10   logs: /…/anqchor/results/cis/peer10/logs/<step>.log
```

Status is `done` (up to date), `stale` (an input or parameter changed) or `missing` (never run,
or an output was deleted). `steps` always lists the report step. After `cis run` it is
`missing`; after `cis report` it is `stale`, because `cis report` renders the page outside the
step graph. Both are expected.

## The run directory

A run of one analysis (here `cis run`) writes:

```text
results/cis/peer10/
├── anqchor.resolved.yaml          # complete settings of the run (defaults filled, paths absolute)
├── run.json                       # provenance (below)
├── Makefile                       # generated on every invocation
├── .anqchor.lock                  # held while a command runs in this directory
├── steps/<step>.json|.done|.seconds   # signatures, stamps, wall times
├── logs/<step>.log                # output of each step
├── prepare/
│   ├── genotypes/                 # pgen_list.tsv, samples.txt, genotypes.json
│   ├── references/<ref>/          # pairs.tsv.gz, variants.tsv, variants.present.txt, reference.json
│   └── datasets/<ds>/             # pheno.regenie.tsv.gz, cov.regenie.tsv.gz, cov.json, traits.tsv,
│       │                          # samples.txt, sample_map.tsv, dataset.json
│       └── ref.<ref>/             # trait_map.tsv, pairs.tsv.gz, variants.txt, traits.txt, mapping.json
├── cis/<ds>/                      # pairs.long.tsv.gz, strategies.tsv, traits.tsv, power_curve.tsv, summary.json
│   └── strategies/default/        # cov.regenie.tsv.gz, strategy.json, pair_assoc.tsv.gz
└── report/cis.html                # only after 'cis report' or 'cis all'
```

A `match` run has `prs/<ref>/` and `match/<ds>/` instead of `cis/`, and a `trans` run
`trans/<ds>/`. A run of a single analysis writes no `report/index.html`. Every file is
described in [Output files](../reference/outputs.md).

### The comparison directory

`cis compare` writes a separate directory:

```text
results/cis/compare/
├── runs.tsv                       # one row per compared run: label, directory, inputs, counts
├── strategies.tsv                 # one row per run: replication counts, relative N_eff
├── pairs.long.tsv.gz              # every pair of every run
├── traits.tsv  power_curve.tsv    # per-gene results, replication vs p-value threshold
├── summary.json                   # recommendation, warnings, comparability
├── compare.json                   # manifest: runs, labels, the compare command
├── settings/<label>.anqchor.yaml  # the exact settings of each compared run
└── report/cis.html
```

`settings/` keeps a copy of each run's `anqchor.resolved.yaml`, so a comparison records how
every run was made even after the run directories are deleted.

## Reproducing a run

`anqchor.resolved.yaml` holds every setting of the run, including those given as flags. To
repeat a run exactly, into a new directory:

```bash
anqchor cis run --config results/cis/peer10/anqchor.resolved.yaml --out results/cis/peer10_again
```

The same works with a copy from a comparison
(`--config results/cis/compare/settings/peer10.anqchor.yaml`) and for `match` and `trans`.
Change one setting with a flag on top (`--cov-include sex,PC1-5`) to make a variant of the run.
qpgen is found again by the usual discovery (`qpgen.binary: null`), so use the same qpgen
build; `run.json` records which one was used.

`anqchor.resolved.yaml` is rewritten by every command that runs steps in the directory, so
it always describes the settings of the current outputs. A `--dry-run` and `anqchor steps`
leave an existing file unchanged.

## Logs and failures

Each step writes `logs/<step>.log`. When a step fails, anqchor prints the step name, the last
lines of its log and the log path, then exits with code 3:

```text
[anqchor] FAILED cis.dataset.default.inputs (exit 1); log: logs/cis.dataset.default.inputs.log
error in step cis.dataset.default.inputs: unknown covariate name(s) 'FOO'; available: sex, PC1-5, PEER1-60
error: run failed in step cis.dataset.default.inputs; log: results/cis/x/logs/cis.dataset.default.inputs.log
```

A failed step leaves no stamp behind. Fix the cause and run the same command again; the steps
that already succeeded are not repeated.

`-v/--verbose` (a root option: `anqchor -v cis run ...`) turns on debug logging.

## Provenance: `run.json`

Rewritten after every invocation that is not a dry run:

| Key | Content |
|---|---|
| `schema_version` | `1` |
| `anqchor_version`, `python`, `platform`, `host` | environment |
| `config_sha256` | sha256 of `anqchor.resolved.yaml` |
| `qpgen` | `binary`, `sha256`, and `git_commit` of the git checkout the binary sits in |
| `invocations` | the last 20 invocations: `argv`, `started`, `finished`, `returncode` |
| `step_seconds` | wall time of the last successful run of each step |

## Exit codes

| Code | Meaning |
|---|---|
| 0 | success |
| 2 | configuration or validation error before any step ran (invalid YAML or keys, missing input file, unknown `--dataset`/`--strategy`/step name, `--only` with upstream steps not up to date, `report` on a directory without results), or the run directory is locked by another command |
| 3 | a step failed (the failing step and its log are printed) |
| 4 | `qpgentools` missing or unusable (`anqchor qpgen build`/`check` also use 4 for build failures and missing subcommands) |

`anqchor prepare ...` utilities exit with 2 on input errors.

## Environment variables

| Variable | Effect |
|---|---|
| `ANQCHOR_QPGEN` | path to `qpgentools` (after `--qpgen` / `qpgen.binary`) |
| `ANQCHOR_MAKE` | make executable (default `make`) |
| `ANQCHOR_TABIX`, `ANQCHOR_BGZIP` | htslib tools (default: `PATH`, then the qpgen submodule's htslib) |
