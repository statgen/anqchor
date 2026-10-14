# Phenotype File Format for `qpgentools`

## Overview

`qpgentools` accepts phenotype or covariatefiles in the following formats:

1. Regenie format
2. TensorQTL format
3. Sample-column TSV format
4. Sample-row TSV format

## Specifying phenotype or covariate files

In `qpgentools`, many commands that accepted this modified PLINK2 format, including 
`qpgentools pair-assoc`, `qpgentools pair-prs`, and `qpgentools rect-assoc`, can accept the phenotype or covariate files in the above formats. The formats can be specified using the `--pheno-format` or `--cov-format` options.

## Delimiters

Phenotype, covariate and PRS matrices, and the other input tables (pair files, `--sample`,
trait and variant lists), are **tab-delimited**. A field may then contain spaces, e.g. a trait
name such as `IL6 (pg/ml)`. Each command takes `--space-delim` for files delimited by single
spaces instead; the two are never mixed. A header that splits into too few fields stops the run
with an error that names the delimiter in effect.

Genotype files (`.psam`, `.fam`, `.pvar`) and REGENIE files (`.loco`, `*_pred.list`) follow their
own formats and are not affected.

## Values and missing values

Every value must be a number, except in categorical covariates (below) and in covariate columns
not selected with `--cov-names`. A value that is not a number (for example sex coded as `XX`/`XY`
without `--cat-cov-names sex`) stops the run with an error, instead of being read as `0` as in
earlier versions.

`pair-assoc`, `rect-assoc`, `region-assoc` and `match-prs-pheno` take `--missing-str` (default `NA`):
a comma-separated list of strings that mark missing values. See each command's page for how
missing values are handled; `rect-assoc` stops on any missing value.

## Categorical and selected covariates

`pair-assoc`, `rect-assoc`, `region-assoc` and `match-prs-pheno` follow REGENIE's covariate options:

* `--cov-names` (REGENIE `--covarColList`) selects the covariate columns to use, e.g.
  `--cov-names "age,PC{1:10}"`; `PC{1:10}` expands to `PC1,...,PC10`. Columns in `--cat-cov-names`
  are always used. The other columns are ignored, are not parsed, and their missing values do not
  exclude samples. Without `--cov-names`, every column is used.
* `--cat-cov-names` (REGENIE `--catCovarList`) lists categorical covariates. Their values may be any
  strings, e.g. `XX`/`XY` or batch labels. After the samples are matched across the input files, a
  covariate with L observed levels is replaced by L-1 indicator columns named `<covariate>_<level>`.
  The levels are sorted as strings and the first is the reference. Which level is the reference
  does not change any result, because the model has an intercept. A level that only excluded
  samples have is not counted. A covariate with a single level is dropped with a warning, and
  `--max-cat-levels` (REGENIE `--maxCatLevels`, default 10) caps the number of levels.
* A missing value (`--missing-str`) in a used categorical covariate is a missing covariate:
  `pair-assoc` excludes the sample, and `rect-assoc` and `region-assoc` stop. `match-prs-pheno`
  excludes the sample, or with `--cov-impute-mean` imputes each indicator column with the level
  frequencies.

For example, with the MAGE covariate file, where sex is stored as `XX`/`XY`:

```
qpgentools pair-assoc ... --cov eQTL_covariates.tab.gz --cov-format tsv-sample-col \
    --cat-cov-names sex --cov-names "PC{1:5},PEER{1:60}"
```

## Regenie format

The Regenie phenotype format (with `--pheno-format regenie`) is a TSV file following convention:

* The first column is `FID`, representing family ID of the sample (PLINK2 convention)
* The second column is `IID`, representing individual ID of the sample (PLINK2 convention)
* From the third column, the remaining columns are the phenotype names. 

Each row represents a sample, and each column represents a phenotype.

Here is an example:

```
FID	IID	pheno1	pheno2	pheno3
1	1	1.2	3.4	5.6
2	2	2.3	4.5	6.7
...
```

(Columns are separated by tabs; use `--space-delim` for space-delimited files.)

Please refer to the [Regenie documentation](https://rgcgithub.github.io/regenie/options/#phenotype-file-format) for more details.

## TensorQTL format

The TensorQTL phenotype format is (with `--pheno-format tensorqtl`) a TSV file that represents the phenotype in a sample-column format.

* The first column is `#chr`, representing the chromosome of the molecular trait (e.g. gene)
* The second column is `start`, representing the base position (1-based) of the molecular trait (e.g. gene)
* The third column is `end`, representing the base position (1-based) of the molecular trait (e.g. gene)
* The fourth column is `phenotype_id`, representing the ID of the molecular trait (e.g. gene)
* From the fifth column, the remaining columns are the sample IDs. 

Each row represents a molecular trait, and each column represents a sample.

Here is an example:

```
#chr	start	end	phenotype_id	sample1	sample2	sample3
1	100000	100001	gene1	1.2	3.4	5.6
2	200000	200001	gene2	2.3	4.5	6.7
...
```

## Sample-row TSV format

The sample-row TSV format (with `--pheno-format tsv-sample-row`) is a TSV file that represents the phenotype in a sample-row format. This is very similar to the Regenie format, but instead of the `FID` and `IID` columns are expected for family and individual ID, only one sample ID column is expected. The rest of the columns are the phenotype names.

Here is an example:

```
SampleID	pheno1	pheno2	pheno3
sample1	1.2	3.4	5.6
sample2	2.3	4.5	6.7
...
```

## Sample-column TSV format

The sample-column TSV format (with `--pheno-format tsv-sample-column`) is a TSV file that represents the phenotype in a sample-column format. This is very similar to the TensorQTL format, but instead of the `#chr`, `start`, `end`, and `phenotype_id` columns are expected for chromosome, base position, base position, and molecular trait ID, only one phenotype ID column is expected. The rest of the columns are the sample IDs.

Here is an example:

```
TraitID	sample1	sample2	sample3
trait1	1.2	3.4	5.6
trait2	2.3	4.5	6.7
...
```

 


