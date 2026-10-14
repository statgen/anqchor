# Configuration

Every analysis can be run from command-line flags alone, and that is the usual way to start.
A YAML configuration file is the alternative when you want the settings in a file, need
options that have no flag (for example `missing_str`, `trait_chroms` or the replication
`levels`), or want to repeat a run exactly.

Every run writes its complete settings, flags included, to `<run>/anqchor.resolved.yaml`
(defaults filled in, paths absolute). That file is a valid configuration file. After the
flags-only run

```bash
D=examples/mage/data
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz --cov-include sex,PC1-5,PEER1-10 \
  --out results/cis/peer10
```

the command

```bash
anqchor cis run --config results/cis/peer10/anqchor.resolved.yaml --out results/cis/peer10_again
```

reproduces it exactly (see
[Running](running.md#reproducing-a-run)).

The file is validated by pydantic models (`src/anqchor/config/model.py`). **Unknown keys are
errors**, so typos fail at load time with a message naming the dotted location (exit code 2):

```text
error: invalid configuration:
  pipelines.cis.levels.strcit: Extra inputs are not permitted
```

## Minimal example

The flags of the example data,

```bash
D=examples/mage/data
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --out results/cis/peer60
```

can be kept in a file instead. Saved as `results/mage.yaml` (relative paths are resolved
against the file's directory):

```yaml
genotypes:
  pgen_list: ../examples/mage/data/geno/mage.pgen_list.tsv
references:
  reference:
    pairs: ../examples/mage/data/reference/gtex_lcl_eqtl.chr21_22.tsv.gz
datasets:
  dataset:
    pheno: ../examples/mage/data/pheno/mage_expr.chr21_22.regenie.tsv.gz
    cov: ../examples/mage/data/pheno/mage_covariates.regenie.tsv.gz
```

```bash
anqchor cis run   --config results/mage.yaml --out results/cis/peer60_config
anqchor match run --config results/mage.yaml --out results/match_config
```

`anqchor <analysis> run --config FILE` runs that one analysis with the file's genotypes,
reference and dataset; each analysis is still a separate command with its own run
directory. With the names `dataset` and `reference` used by the flags, the outputs have the
same paths as a flags-only run. The `enabled` keys under `pipelines` only matter for the top-level
`anqchor run` (see [Command line](../reference/cli.md#anqchor-run)).

## Flags and configuration keys

| Flag | Configuration key |
|---|---|
| `--geno-list FILE` / `--pfile PREFIX` | `genotypes.pgen_list` / `genotypes.pfile` |
| `--pairs FILE` | `references.<name>.pairs` (`<name>` = `--reference-name`, default `reference`) |
| `--pheno FILE`, `--pheno-format` | `datasets.<name>.pheno`, `.pheno_format` (`<name>` = `--dataset-name`, default `dataset`) |
| `--cov FILE`, `--cov-format` | `datasets.<name>.cov`, `.cov_format` |
| `--cat-cov-names NAMES` | `datasets.<name>.cov_categorical` (comma-separated → list) |
| `--trait-annot FILE` | `datasets.<name>.trait_annot` |
| `--trait-id-rule RULE` | `datasets.<name>.trait_id_rule` |
| `--sample-map FILE` | `datasets.<name>.sample_map` |
| `--cov-include NAMES` (`cis`, `trans`) | `pipelines.<analysis>.strategies: [{name: default, cov: {include: [NAMES]}}]` |
| `--assoc-engine` (`cis`) | `pipelines.cis.assoc_engine` |
| `--rint` / `--no-rint` (`cis`, `trans`) | `pipelines.<analysis>.rint` |
| `--out DIR` | `out_dir` |
| `--threads N`, `--qpgen PATH` | `qpgen.threads`, `qpgen.binary` |
| `--plotlyjs` | `report.plotlyjs` |

The `match` and `trans` flags map to the keys of the same name under `pipelines.match` and
`pipelines.trans` (`--z-threshold` → `z_threshold`, `--mahalanobis` → `mahalanobis: true`,
`--cis-window-bp` → `cis_window_bp`, `--hotspot-min-traits` → `hotspot_min_traits`, ...). The [command-line reference](../reference/cli.md)
lists them.

!!! note "The reference p-value column"
    `--pairs` reads the `trait variant beta se` columns and, when the file has one, a column
    named `pval`, reported as `pval_ref`. A p-value column with another name needs
    `columns: {pval: NAME}` in a configuration file. Replication does not use it.

## Rules

**Paths.** Relative paths in the file are resolved against the **directory of the config
file**; paths given as command-line flags are resolved against the **current directory**. A
`file:` prefix is optional (`traits: file:my_traits.txt`). `~` is expanded.

**Precedence.** Command-line flags > YAML file > model defaults. Unset flags never override the
file. A genotype flag (`--geno-list`, `--pfile`) replaces the file's whole `genotypes`
section. Input flags such as `--pheno` describe a dataset named by `--dataset-name` (default
`dataset`); when the file defines exactly one dataset under another name, the flags override
that dataset instead of adding a second one (likewise `--pairs` and the single reference).
`--cov-include` replaces the file's whole `strategies` list for that analysis.

**Names.** Reference, dataset and strategy names use letters, digits, `_`, `.`, `-` and start
with a letter or digit. They appear in step and file names.

**Several references.** When more than one reference is defined, every *enabled* pipeline
must name its `reference`.

## `version`, `project`, `out_dir`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `version` | `1` | `1` | configuration schema version |
| `project` | string | `anqchor` | project name shown in reports |
| `out_dir` | path | `anqchor_out` | run directory; `--out` overrides it |

## `qpgen`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `binary` | path or `null` | `null` | explicit `qpgentools`; `null` = discovery (`--qpgen`, `$ANQCHOR_QPGEN`, `PATH`, `extern/qpgen/bin`) |
| `threads` | int ≥ 1 | `1` | threads passed to `match-prs-pheno --threads` (`--threads` flag). The other qpgen commands are single-threaded. Not part of any step signature: changing it re-runs nothing. |

## `genotypes`

Exactly one of:

| Key | Type | Meaning |
|---|---|---|
| `pgen_list` | path | qpgen pgen list, `CHROM BEG END PREFIX` (or 6 columns) |
| `pfile` | path | one PLINK2 prefix with `.pgen`, `.psam`, `.pvar.idx.gz` (+ `.tbi`) |

See [Inputs: genotypes](inputs.md#genotypes).

## `references.<name>`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `kind` | `cis` \| `trans` | `cis` | a `cis` reference needs `pairs`; a `trans` reference needs `pairs` or `variants` |
| `pairs` | path | — | pair file with a header (trait, variant, beta, se, optional pval) |
| `variants` | path | — | trans only: one `CHROM:POS:REF:ALT` per line |
| `columns.trait` | string | `trait` | column name of the trait ID |
| `columns.variant` | string | `variant` | column name of the CPRA variant ID |
| `columns.beta` | string | `beta` | effect of the ALT allele |
| `columns.se` | string | `se` | standard error |
| `columns.pval` | string or `null` | `null` | optional p-value column, carried into outputs; `null` reads a column named `pval` when present |

`match` and `cis` need a reference with `pairs` (effect sizes for the PRS, pairs to replicate).

## `datasets.<name>`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `omics` | string | `omics` | label shown in reports |
| `pheno` | path | **required** | phenotype matrix |
| `pheno_format` | `regenie` \| `tensorqtl` \| `tsv-sample-row` \| `tsv-sample-col` | `regenie` | [layout](inputs.md#phenotypes) |
| `trait_annot` | path | `null` | trait loci, BED-like or `trait chrom start end`; a `tensorqtl` phenotype provides loci itself |
| `cov` | path | `null` | covariate matrix (the dataset's **default** covariates) |
| `cov_format` | as `pheno_format` | `regenie` | covariate layout |
| `cov_categorical` | list of names | `[]` | force these covariates to be categorical, e.g. a batch number (non-numeric ones are detected automatically, see [covariates](inputs.md#covariates)); flag `--cat-cov-names a,b`. Categorical covariates keep their text levels in the prepared file and are passed to qpgen as `--cat-cov-names`, which encodes L levels as L − 1 indicators (anqchor adds `--max-cat-levels` when a covariate has more than qpgen's default of 10). A non-numeric column with more than 20 levels (and more than a tenth of the samples) must be listed here |
| `cov_missing` | `error` \| `drop_samples` \| `impute_mean` | `error` | missing covariate values |
| `missing_str` | list of strings | `["NA"]` | spellings of a missing value in phenotype and covariate files; an empty cell is always missing |
| `sample_map` | path | `null` | `match` only: `geno_id pheno_id` table; `null` = identical IDs |
| `trait_id_rule` | `as_is` \| `before_colon` \| `strip_version` \| `before_colon_strip_version` \| `file:<path>` | `as_is` | how reference trait IDs map to dataset trait IDs ([rules](concepts.md#trait-ids)); `file:<path>` = `trait_map: <path>` |
| `trait_map` | path | `null` | explicit `ref_trait pheno_trait` file; overrides the rule |

## `pipelines.match`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | bool | `true` | include in the top-level `anqchor run`; `anqchor match run` ignores it |
| `datasets` | list or `null` | `null` (all) | datasets to check |
| `reference` | name or `null` | `null` (the only reference) | reference for the PRS |
| `prs_samples` | `all` \| path | `all` | genotype IDs (one per line) to compute PRS for; `all` scores every genotyped sample, which improves calibration and enables rescue |
| `use_cov` | bool | `true` | adjust phenotypes for the dataset covariates |
| `rint` | bool | `true` | rank-based inverse normal transform after adjustment (`--rint`) |
| `mahalanobis` | bool | `false` | whiten correlated traits (`--mahalanobis`, slower); default scores traits as independent |
| `lambda` | `auto` \| number in [0, 1] | `auto` | with `mahalanobis`: `auto` → `--auto-lambda` when at least `min_self_for_auto` samples are mapped, else `--lambda 0.5` (with a warning); a number → `--lambda` |
| `min_self_for_auto` | int ≥ 1 | `100` | mapped-sample minimum for `--auto-lambda` |
| `z_threshold` | float | `1.96` | Z above which a non-best own genotype is a lenient match |
| `z_diff` | float ≥ 0 | `2.0` | Z gap that makes a leading genotype clear |
| `min_weight` | float | `0.0` | `--min-weight` (passed only when non-zero) |

## `pipelines.cis`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | bool | `true` | |
| `datasets` | list or `null` | `null` (all) | |
| `reference` | name or `null` | `null` | reference whose pairs are replicated |
| `assoc_engine` | `pair` \| `rect` | `pair` | `pair-assoc` over the pairs, or `rect-assoc` over the reference's variants × traits |
| `strategies` | list of [strategies](#strategies) | `[{name: default}]` | covariate/phenotype choices of this run; `--cov-include` sets a single one. First = baseline; names unique |
| `rint` | bool | `false` | rank-based inverse normal transform of the phenotypes after covariate adjustment (`--rint`); `false` tests the phenotypes as given. A cis run tests one transform; compare raw and RINT as two runs |
| `levels.nominal` | float in (0, 1) or `null` | `0.05` | sign and p < value |
| `levels.strict` | float or `null` | `1e-5` | sign and p < value |
| `levels.bonferroni` | bool | `true` | sign and p < 0.05 / n_tested |
| `levels.fdr` | float or `null` | `0.05` | sign and BH q ≤ value |
| `max_unresolved_fraction` | float in [0, 1] | `0.2` | fail when more of the mapped pairs are `variant_absent` or `not_tested` |

## `pipelines.trans`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | bool | `true` | |
| `datasets` | list or `null` | `null` (all) | |
| `reference` | name or `null` | `null` | |
| `variants` | `reference` \| path | `reference` | the reference's **present** variants, or a CPRA list file |
| `traits` | `reference` \| `all` \| path | `reference` | the reference's mapped traits (the initial test rectangle), every dataset trait, or a file of trait IDs (first column) |
| `trait_chroms` | `genotyped` \| `all` \| list | `genotyped` | keep traits whose locus is on these chromosomes: `genotyped` = the genotype set's chromosomes; a list such as `[21, 22]` (either style); `all` = no restriction. Traits without a locus are kept only with `all` |
| `known_pairs` | bool | `true` | replicate the reference pairs and leave the regions of its trans pairs out of the inflation cells (needs `pairs`) |
| `strategies` | list of one [strategy](#strategies) | `[{name: default}]` | the covariate/phenotype choice of this run (exactly one entry); `--cov-include` sets it |
| `rint` | bool | `false` | rank-based inverse normal transform of the phenotypes after covariate adjustment (`--rint`, `rect-assoc --rint-after-adj`); `false` tests the phenotypes as given. A trans run tests one transform |
| `cis_window_bp` | int ≥ 0 | `1000000` | a variant within this distance of the gene `start` on the same chromosome is cis; every other pair is trans (`--cis-window-bp`) |
| `min_maf` | float in [0, 0.5) | `0.01` | `rect-assoc --min-maf` |
| `min_mac` | float ≥ 0 | `1.0` | `rect-assoc --min-mac` |
| `levels.nominal` | float in (0, 1) | `1e-5` | a trans signal is p < value: the report's default hotspot threshold and the cut-off of `pairs.sig.tsv.gz` |
| `levels.strict` | float or `null` | `null` | an extra threshold for the hotspot grid (e.g. `1e-8`) |
| `levels.bonferroni` | bool | `true` | add p < 0.05 / trans cells tested to the hotspot grid |
| `hotspot_bin_bp` | int ≥ 1 | `1000000` | width of the bins that thin trans signals for hotspots (`--hotspot-bin-bp`) |
| `hotspot_min_traits` | int ≥ 1 | `5` | genes with a trans signal a bin needs to be a hotspot; the report's default (`--hotspot-min-traits`) |
| `max_chunk_vars` | int ≥ 1 | `100` | `rect-assoc --max-chunk-vars` |
| `max_sig_rows` | int ≥ 1 | `1000000` | cap on rows in `pairs.sig.tsv.gz` (strongest kept) |

Accepted spellings: `variants`/`traits: reference:<name>` name the reference and set
`reference`; `known_pairs: reference:<name>` means `true`, `known_pairs: null` means `false`.
All three must name the same reference.

!!! warning "Keep the rectangle small"
    `traits: all` with `trait_chroms: all` tests every variant against every trait. With a
    genome-wide expression matrix that is 20,000+ traits and an output of `11 + 4 × n_traits`
    columns per variant. Prefer `trait_chroms: genotyped` (the default) or `traits: reference`.

## Strategies

A *strategy* is a named covariate/phenotype choice inside one run. On the command line a run
has exactly one, called `default`: all covariates of `--cov`, the subset named by
`--cov-include`, or none when `--cov` is omitted. In a configuration file:

- **cis**: keep one strategy per run (as `examples/mage/cis.yaml` does) and compare runs with
  `anqchor cis compare`. A list of several strategies still works; their results then share
  one run directory and are labelled `<run label>:<strategy>` in a comparison.
- **trans**: exactly one strategy per run: start one `trans run` per covariate set
  (`--cov-include`), each into its own directory. A list of several strategies is an error
  that says so. Compare the runs with `anqchor trans compare`.

Each entry of `pipelines.cis.strategies` and `pipelines.trans.strategies`:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `name` | string | **required** | unique within the pipeline |
| `cov` | `default` \| `null` \| `{include: [...]}` \| `{file: path, format: ...}` \| path | `default` | `default` = all dataset covariates; `null` = none; `include` = a named subset (ranges `PEER1-15`; unknown names are errors); `file` = a separate covariate file, checked like the dataset's, categorical covariates kept as text (a bare path is shorthand). Each strategy passes qpgen `--cat-cov-names` with the categorical covariates it uses (recorded as `cat_cov_names` in `strategy.json`) |
| `pheno` | path | `null` | alternative phenotype matrix (e.g. another normalisation) for this strategy |
| `pheno_format` | layout | `regenie` | layout of `pheno` |

!!! note "One transform per run"
    Strategies have no `transforms` key: a cis or trans run tests the phenotypes as given,
    or with `pipelines.<analysis>.rint: true` (`--rint`) their RINT after covariate
    adjustment. Configuration files written for earlier versions with `transforms` under
    `pipelines.cis.strategies` or `pipelines.trans.strategies` are rejected with a message
    saying what to write instead: remove the key (`[raw]`), set `rint: true` (`[rint]`), or
    make two runs (`[raw, rint]`). `pipelines.trans.levels.fdr` is likewise rejected.

The forms of `cov`, written for a file next to `examples/mage/cis.yaml` (a run normally has
one of these entries):

```yaml
strategies:
  - {name: peer60, cov: default}                                    # every covariate (= --cov)
  - {name: peer10, cov: {include: [sex, PC1-5, PEER1-10]}}          # = --cov-include sex,PC1-5,PEER1-10
  - {name: nocov, cov: null}                                        # no covariates
  - {name: filecov, cov: {file: data/pheno/mage_covariates.regenie.tsv.gz, format: regenie}}
```

The last entry reads its covariates from a file of its own instead of the dataset's `cov`
(here the same file, so it gives the same results as `peer60`). A `pheno` entry works the
same way for a phenotype matrix with another normalisation.

## `report`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `plotlyjs` | `embed` \| `cdn` | `embed` | `embed` writes plotly.js into each page (works offline; the example reports are 5-10 MB each); `cdn` loads it from the CDN (small files; offline, each chart area shows a notice suggesting `--plotlyjs embed` and the tables still work) |

## Annotated examples

`examples/mage/cis.yaml`, the all-covariate cis run of the example data:

```yaml
--8<-- "examples/mage/cis.yaml"
```

`examples/mage/match.yaml` (identity matching) follows the same pattern. Each file's header
shows the equivalent flags-only command.
