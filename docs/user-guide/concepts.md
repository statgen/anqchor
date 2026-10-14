# Concepts

The terms below are used consistently in the commands, configuration, file names, logs and reports.

## Inputs

Genotype set
:   One PLINK2 dataset, possibly split by chromosome, described by a qpgen **pgen list**
    (`CHROM BEG END PREFIX`) or a single prefix (`pfile`). Each prefix needs `.pgen`, `.psam`
    and an indexed pvar `.pvar.idx.gz` + `.tbi`. All files of a set must share one sample list.
    See [Inputs](inputs.md#genotypes).

Reference
:   A set of known xQTL effects: a **pair file** with `trait variant beta se [pval]` per row.
    Variants are `CHROM:POS:REF:ALT` (CPRA), with `ALT` the effect allele. Kind `cis` (pairs
    required) or `trans` (pairs and/or a variant list). A reference feeds the PRS of `match`,
    the replication test of `cis`, and the variants and known positives of `trans`. For
    evaluating analysis choices, use a reference discovered in *other* samples (for example
    GTEx): pairs discovered in your own samples, with your own covariates, replicate best
    under those same choices by construction.

Dataset
:   One omics phenotype matrix (traits × samples) with optional covariates, trait annotation and
    sample map. A flags-only run has one dataset called `dataset`; a configuration file may
    define several (`rnaseq`, `proteomics`, ...). The `omics` label is free text shown in
    reports.

Trait
:   One molecular phenotype (gene, protein, metabolite, CpG), identified by a string ID.

Sample map
:   A table `geno_id pheno_id` linking genotype IDs to omics sample IDs, used **only by
    `match`**. Absent means identical IDs.

Trait map / trait-ID rule { #trait-ids }
:   How a reference's trait IDs are matched to a dataset's trait IDs. A **rule** is applied to
    both sides and IDs match when the transformed keys agree:

    | `trait_id_rule` | Transformation | `ENSG00000227232.5:WASH7P` becomes | `ENSG00000227232.5` becomes |
    |---|---|---|---|
    | `as_is` (default) | none | unchanged | unchanged |
    | `before_colon` | part before the first `:` | `ENSG00000227232.5` | unchanged |
    | `strip_version` | drop a trailing `.<digits>` | unchanged (no trailing version) | `ENSG00000227232` |
    | `before_colon_strip_version` | both | `ENSG00000227232` | `ENSG00000227232` |

    For anything else, give an explicit **trait map** file (`ref_trait pheno_trait`, two
    columns). When no trait maps, the error names the rule that would have worked best.
    Trait annotations (`trait_annot`) are matched to phenotype IDs with the same rules,
    automatically: the rule matching the most traits wins.

## Analyses

Analysis
:   One of the three questions anqchor answers, each run on its own with its own command:
    **identity matching** (`anqchor match`), **cis evaluation** (`anqchor cis`) and **trans
    evaluation** (`anqchor trans`). Internally each analysis is a *pipeline* of steps; the
    words are used interchangeably in file and step names (`match/`, `report.cis`).

PRS
:   Polygenic score of each genotyped sample for each reference trait, computed by
    `qpgentools pair-prs`. It depends only on the genotypes and the reference
    (`prs/<ref>/`).

Transform
:   How the phenotype is tested after covariate adjustment:

    - as given (`raw`): the phenotype is residualised on the covariates (with intercept) and
      regressed on the genotype, which is residualised on the same covariates.
    - `rint`: the covariate-adjusted phenotype is rank-based inverse-normal transformed (Blom
      offsets, average ranks for ties) and then tested the same way.

    A cis or trans run tests one transform: the phenotypes as given, or RINT with `--rint`;
    raw and RINT are compared as two runs (with `cis compare` or `trans compare`).

## Runs and comparisons

Run
:   One execution of an analysis on **one set of inputs**: genotypes, phenotypes, covariates
    (optionally a subset of their columns, `--cov-include`), for cis and trans the phenotype
    transform (`--rint`), and a reference. A run writes into
    its own **run directory** (`--out`). To evaluate an alternative, such as other covariates,
    another normalisation or another genotype set, start another run into another directory.
    Runs are independent and can run in parallel.

Comparison
:   `anqchor cis compare RUN_DIR... --out DIR` reads finished cis runs and puts them side by
    side in one set of tables and one report. The runs may differ in any input. Each run is
    identified by a **label** (its directory name unless `--label` is given). A comparison
    counts replication over each run's own tested pairs and, separately, over the pairs tested
    in *every* run (the `_common` columns), so runs with different genotypes or references stay
    comparable.

Run directory
:   The output directory of a run (`--out`). It holds the complete settings
    (`anqchor.resolved.yaml`), the provenance (`run.json`), the generated `Makefile`, step
    logs, every output and, once rendered, the report. See
    [Running](running.md#the-run-directory).

Run, report, all
:   Each analysis separates computing from rendering. `run` computes the summary tables;
    `report` renders the HTML from those tables in about a second and runs no analysis; `all`
    does both.

Strategy
:   A term of the configuration file, not of the command line. In a config file, a
    `strategies` list under `pipelines.cis` or `pipelines.trans` defines several named
    covariate/phenotype choices evaluated inside one run. Prefer one run per choice: a
    flags-only `cis run` or `trans run` has a single strategy called `default`. A cis run's
    and trans results are labelled with the run directory name (a trans configuration has
    exactly one strategy). cis runs are compared with `cis compare`, trans runs with
    `trans compare`.

## Execution

Step
:   One unit of work with declared input files, output files and commands, named with dots
    (`prepare.dataset.rnaseq`, `cis.rnaseq.default.assoc`). A step is either a qpgen command or
    a Python step run as `python -m anqchor _step <name> --run-dir .`. Steps form a graph and
    run through a generated `Makefile`, so a repeated command re-runs only what changed.

## Results

Match status and group { #match-groups }
:   qpgen's per-sample status is kept verbatim; anqchor adds a coarse group.

    | Group | qpgen statuses | Meaning |
    |---|---|---|
    | `PASS` | `SELF_BEST` | the labelled genotype is the best match |
    | `WARN` | `SELF_LENIENT`, `UNCLEAR` | the labelled genotype is not first, but its Z is above `z_threshold` (`SELF_LENIENT`), or nothing is clear (`UNCLEAR`) |
    | `FAIL` | `SINGLE_NEW_BEST`, `MULTI_NEW_BEST` (of a mapped sample) | one (or a block of) other genotype(s) matches clearly better |
    | `NO_TEST` | `NO_SELF`, `NO_OBS_TRAITS`, and every unmapped sample | there is no label to confirm |

    An *unmapped* sample with `SINGLE_NEW_BEST` is a **rescue candidate**. Rules:
    [match pipeline](../pipelines/match.md#status-rules).

Replication level { #replication-levels }
:   How strongly a known pair is reproduced. Every level requires **sign concordance**
    (`sign(beta_ref) == sign(beta)`). With the cis defaults and the 367 pairs of the example
    data, from weakest to strongest (the order is computed from each level's actual p
    cut-off, so it can change with the number of tested pairs):

    | Level | Rule (cis defaults) |
    |---|---|
    | `sign` | sign concordance alone (always reported) |
    | `nominal` | sign and p < 0.05 |
    | `fdr` | sign and Benjamini-Hochberg q ≤ 0.05 over the tested pairs of the run |
    | `bonferroni` | sign and p < 0.05 / n_tested |
    | `strict` | sign and p < 1e-5 |

    The cutoffs are configurable under `levels`; a level set to `null` (or `bonferroni: false`)
    is not reported. `trans` replicates its reference pairs at the same sign, p < 0.05 and
    p < 1e-5; its own `levels` (`nominal: 1e-5`, optional `strict`, Bonferroni over the trans
    cells) are thresholds of the hotspot grid.

Pair status
:   Every mapped reference pair is `tested`, `variant_absent` (its variant is not in the
    genotype pvar under that exact CPRA) or `not_tested` (variant present, but qpgen wrote no
    result or no finite statistic). Unresolved = absent + not tested; `cis` fails when more than
    `max_unresolved_fraction` (20%) of the pairs are unresolved.

cis and trans (trans)
:   A (variant, gene) cell of the trans rectangle, or a reference pair, is **cis** when the
    variant is on the gene's chromosome within `cis_window_bp` (1 Mb) of the gene start, and
    **trans** otherwise: on another chromosome (cross-chromosome) or farther on the same one.
    Genes without a locus are counted apart. Inflation (λ_GC, QQ) uses the cross-chromosome
    cells, leaving out the regions of the reference's trans pairs.

Hotspot (trans)
:   A bin of `hotspot_bin_bp` (1 Mb) whose variants have trans signals (p below a threshold)
    for at least k distinct genes, genes in cis with the bin not counted. The report shows
    the counts for a grid of thresholds and k.

λ_GC
:   Genomic control inflation: median χ²(1) of the cells divided by 0.4549. Values near 1 mean a
    calibrated null; anqchor warns above 1.1 or below 1/1.1.
