# Inputs

anqchor needs three kinds of input: a **genotype set**, a **reference** of known xQTLs, and
an omics **dataset** (phenotypes, optionally covariates and trait loci). On the command line
they are `--geno-list` (or `--pfile`), `--pairs`, and `--pheno` with `--cov`, `--trait-annot`
and `--sample-map`. Everything is validated and normalised by the `prepare.*` steps before qpgen
runs, so most input problems fail early with a clear message (exit code 3, with the step log
named).

The example data in `examples/mage/data/` has one file of each kind; the
[example data page](../use-cases/example-data.md) shows the first lines of each.

## Genotypes

PLINK2 files, one prefix per chromosome or one for everything. Each prefix needs:

| File | Notes |
|---|---|
| `PREFIX.pgen` | genotypes (hard calls or dosages) |
| `PREFIX.psam` | samples; header `#IID ...` or `#FID IID ...`; IIDs are the genotype sample IDs |
| `PREFIX.pvar.idx.gz` + `.tbi` | **indexed pvar**: a 9-column bgzipped, tabix-indexed VCF-like file (`#CHROM POS ID REF ALT QUAL FILTER INFO INDEX`) with a 1-based `INDEX` column in pgen row order |

The indexed pvar is what qpgen uses to look variants up by `CHROM:POS:REF:ALT`. Create it with:

```bash
anqchor prepare index-pvar --pfile geno/chr21
# wrote geno/chr21.pvar.idx.gz (127563 variants) and geno/chr21.pvar.idx.gz.tbi
```

It reads `PREFIX.pvar`, `PREFIX.pvar.zst` (needs `zstd`) or `PREFIX.pvar.gz`; a headerless
pvar is read in `.bim` column order. It needs `bgzip` and `tabix` and checks that the number of
records written equals the number read back. The pvar must be sorted for `tabix`.

### The pgen list

Describe a multi-file genotype set with a qpgen **pgen list**, one row per file, whitespace
separated, no header (`#` lines ignored):

```text
21  1  2000000000  /data/geno/chr21
22  1  2000000000  /data/geno/chr22
```

Columns are `CHROM BEG END PREFIX`. The 6-column form `CHROM BEG END PGEN PSAM PIVAR` is also
accepted; anqchor recognises the PSAM and PIVAR fields by name or content because qpgen's
parser and documentation disagree on their order. Write a list with:

```bash
anqchor prepare pgen-list --prefix-template '/data/geno/chr{chrom}' --chroms 1-22,X --out geno.list.tsv
```

`CHROM` is read from each file's first pvar record so it matches the pvar spelling; every
prefix must have `.pgen`, `.psam`, `.pvar.idx.gz` and `.tbi` (`--no-check` skips the checks).

For a single prefix use `--pfile PREFIX` (`genotypes.pfile` in a configuration file).

What `prepare.genotypes` does:

- Resolves **relative prefixes** against the list's directory and then up to three parent
  directories (lists are often written relative to a project directory); qpgen itself would
  resolve them against the current directory.
- Checks that every file exists, that all psam files list the same samples in the same order,
  and that sample IDs are unique.
- Detects the **chromosome style** from the first pvar record (`chr21` → `chr`, `21` → `plain`)
  and requires the list's `CHROM` column to use the same style.
- Writes `prepare/genotypes/pgen_list.tsv` in the **4-column form with absolute prefixes**.
  Files whose names are not `PREFIX.{pgen,psam,pvar.idx.gz}` are symlinked under
  `prepare/genotypes/links/` so that the 4-column form always works.

!!! note "Variant IDs"
    qpgen matches variants by exact `CHROM:POS:REF:ALT` string equality. anqchor rewrites the
    chromosome of reference variants to the genotype style (`chr22:…` ↔ `22:…`), but does not
    flip alleles or left-normalise indels. Variants that do not match are reported as
    `variant_absent` with an SNV/indel split.

## Phenotypes

A traits × samples matrix in one of four layouts (`pheno_format`):

| Format | Orientation | Header | Notes |
|---|---|---|---|
| `regenie` (default) | rows = samples | `FID IID trait1 trait2 …` | canonical internal format |
| `tensorqtl` | rows = traits | `#chr start end phenotype_id s1 s2 …` | first four column names matched case-insensitively with optional `#` (`#Chr` is fine; `chrom` accepted); also provides trait loci |
| `tsv-sample-row` | rows = samples | `ID trait1 …` | |
| `tsv-sample-col` | rows = traits | `ID s1 s2 …` | |

Files may be plain or gzip-compressed (detected from the file contents, not the extension),
tab or whitespace separated. Values must be numeric and finite (`inf` is an error naming the
trait and sample); missing values are spelled by `missing_str` (default `NA`).

Sample, trait, covariate and FID IDs must not contain whitespace or a double quote; anqchor
rejects such IDs. anqchor hands qpgen tab-separated files, and qpgen does not strip quotes, so
an ID such as `"A"` or `A B` would not match between the two (builds before qpgen 0.2.0
also split fields on any whitespace). A `tsv-sample-col` or `tensorqtl` header column with no
numeric value at all is reported as a probable annotation column, not a sample.

`prepare.dataset.<ds>` writes `prepare/datasets/<ds>/pheno.regenie.tsv.gz`. A tab-separated
regenie file with `FID == IID` and only `NA` as missing value is copied byte for byte;
anything else is parsed and written with round-trip precision, so values are never rounded.
The standalone equivalent is `anqchor prepare pheno --in F --format tensorqtl --out F.tsv.gz`.

!!! warning "Use identical sample IDs for `cis` and `trans`"
    `pair-assoc` and `rect-assoc` match phenotype and genotype samples **by identical ID**.
    The `sample_map` is used only by `match`. Run `match` first, then rename or drop samples in
    your phenotype and covariate files before running `cis` and `trans`.

!!! warning "Missing phenotype values"
    With qpgen 0.2.0, `pair-assoc` (the default cis engine) drops a sample with a missing
    phenotype from that trait's tests only. `rect-assoc` (the `rect` cis engine and `trans`)
    stops with an error on a missing value. anqchor does not check this before running qpgen:
    fill or remove missing phenotype values before running `trans` or `cis --assoc-engine
    rect`. `dataset.json` reports `n_missing_values`.

## Covariates

Optional (`--cov`, `--cov-format`), in any of the four layouts. The example covariates use the
regenie layout (one row per sample); a file with one row per covariate is `tsv-sample-col`.
`prepare.dataset` checks them and writes them in regenie layout to
`prepare/datasets/<ds>/cov.regenie.tsv.gz`. Categorical covariates keep their text levels
(`sex` stays `XX`/`XY`): qpgen encodes them itself when anqchor passes `--cat-cov-names`.

- **Empty cells** are missing values, never a category level.
- **Numeric covariates** must be numeric throughout: a column of mostly numbers with a few
  non-numeric cells (`nan`, `.`) is an error naming the first such cell; add the string to
  `missing_str` if it means "missing", or declare the covariate categorical. Infinite values
  are an error.
- **Categorical covariates** are every column with a non-numeric value (detected
  automatically) plus the ones you name with `--cat-cov-names` (config:
  `datasets.<ds>.cov_categorical`). Name a covariate when its categories are coded as
  numbers, for example a batch number `1`, `2`, `3`: without the flag it would be one numeric
  covariate with a linear effect. A non-numeric column with more than 20 levels (and more
  than a tenth of the samples) is an error unless named, since it would add one indicator
  column per level.
- **How qpgen encodes them** (REGENIE `--catCovarList` semantics): a covariate with L levels
  becomes L − 1 indicator columns; the first level in sorted order is the reference (`sex`:
  `XX` reference, one `XY` indicator). Which level is the reference does not change any
  result. qpgen rejects a categorical covariate with more than 10 levels unless told otherwise
  (`--max-cat-levels`), so anqchor adds `--max-cat-levels N` when a covariate it uses has
  N > 10 levels.
- **Constant columns** are dropped (they would be collinear with the intercept).
- **Missing values** follow `cov_missing`: `error` (default), `drop_samples` (remove samples
  with any missing covariate) or `impute_mean` (the most frequent level for a categorical
  covariate). qpgen therefore never sees a missing covariate. Samples dropped by
  `drop_samples` are logged as warnings, recorded in `dataset.json` (`n_samples_analysable`,
  `warnings`) and shown in the `cis` and `trans` summaries.

Each qpgen command gets `--cat-cov-names` with the categorical covariates it uses:
`--cov-include sex,PC1-5` passes `--cat-cov-names sex`, `--cov-include PC1-5,PEER1-10`
passes none. `prepare/datasets/<ds>/cov.json` records the columns as written
(`raw_columns`), the categorical covariates (`categorical_names`, `categorical`) and dropped
columns (see [outputs](../reference/outputs.md#prepare)). Standalone:

```bash
anqchor prepare cov --in examples/mage/data/pheno/mage_covariates.regenie.tsv.gz \
    --format regenie --out results/cov.regenie.tsv.gz --json results/cov.json
```

```text
wrote results/cov.regenie.tsv.gz: 729 samples x 66 covariates; categorical: sex ['XX', 'XY']
categorical covariates keep their text levels; use --cat-cov-names sex
```

Add `--encode` to write numeric indicator columns instead of text levels, for a tool without
categorical support.

To use only some covariates of the file, name them with `--cov-include` (or `cov: {include:
[...]}` in a configuration file). Ranges expand: `--cov-include sex,PC1-5,PEER1-10` means
`sex, PC1…PC5, PEER1…PEER10`. An unknown
name stops the run with an error listing the available names
(`unknown covariate name(s) 'FOO'; available: sex, PC1-5, PEER1-60`).

## Trait annotation

Trait loci are needed by `trans` (cell classes, `trait_chroms`) and by the `rect` cis engine's
cross-chromosome controls. Sources, in order:

1. `trait_annot`: a BED-like file (`chrom start end id` + anything, header optional, so a
   tensorQTL matrix works) or a `trait chrom start end` table with a header (column names such
   as `trait`/`gene_id`/`phenotype_id` and `chr`/`chrom`/`chromosome` are recognised). Without
   a header the layout is chosen by which column looks like a chromosome; when both or neither
   do (numeric trait IDs), anqchor asks for a header line;
2. otherwise the phenotype file itself when it is in `tensorqtl` format.

Annotation IDs are matched to phenotype trait IDs directly, or with whichever
[trait-ID rule](concepts.md#trait-ids) matches the most traits (a BED naming genes `ENSG…5`
matches expression IDs `ENSG…5:WASH7P` by `before_colon`). In the example data the IDs are
identical. Chromosome names are rewritten to the
genotype style. The cis window is measured from `start`. Result:
`prepare/datasets/<ds>/traits.tsv` (`trait chrom start end`, `NA` for unknown loci).

## References

A pair file with a header, tab or whitespace separated, plain or gzipped:

```text
trait                    variant          beta       se         pval
ENSG00000275993.3:SIK1B  21:6122794:C:T   -0.448959  0.0888711  8.054e-07
```

- Columns are mapped with `columns: {trait: …, variant: …, beta: …, se: …, pval: …}`; `pval` is
  optional and only carried into the outputs. Other columns are ignored.
- `variant` is `CHROM:POS:REF:ALT` and **`beta` is the effect of the `ALT` allele**: the PRS and
  the sign concordance test rely on it. Nucleotide alleles are upper-cased (`a` → `A`).
- Rows with a malformed variant, an empty trait or a non-finite `beta`/`se` are dropped
  (`n_invalid`); repeated (trait, variant) rows keep the first (`n_duplicates`).
- `prepare.reference.<ref>` checks every variant against the indexed pvar (`tabix -R`, or by
  streaming without tabix) and writes `variants.tsv` with a `present` flag. A reference with
  **no** present variant is an error (usually chromosome naming or genome build).

A `trans` reference may instead (or additionally) give `variants:`, a CPRA list with one
variant per line (first column; a non-CPRA first line is treated as a header).

## Sample maps (match only)

When omics sample IDs differ from genotype IDs, give `sample_map`: a table with a header
`geno_id pheno_id` (either order; the older `NWD_ID SAMPLE_ID` header is accepted too), one row
per omics sample. Rows naming IDs absent from the genotypes or the phenotypes are dropped and
counted in `dataset.json`. Unmapped phenotyped samples are still scored and can become rescue
candidates. `cis` and `trans` ignore the map; when some mapped samples do not have identical
genotype and phenotype IDs, `prepare` logs a warning and adds it to `dataset.json` `warnings`
(it does not fail, even with no identical IDs, because it does not know which pipelines will
run).

## Trait maps

When no trait-ID rule fits, give `trait_map`: a two-column file `ref_trait pheno_trait`
without header (a `ref_trait pheno_trait` header line is tolerated). `trait_id_rule:
file:<path>` is a shorthand for the same thing.
