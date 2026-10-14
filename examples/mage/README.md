# MAGE chr21/22 example

A compact, real xQTL dataset for trying anqchor: MAGE lymphoblastoid RNA-seq of 729
1000 Genomes individuals with NYGC whole-genome genotypes, restricted to chromosomes 21
and 22, and known cis-eQTLs from an independent study (GTEx v10, EBV-transformed
lymphocytes). It ships with the repository in `data/` (32 MB), so every command below works
right after installation and `anqchor qpgen build`.

| File (under `data/`) | Content |
|---|---|
| `geno/mage.pgen_list.tsv` (+ `mage_chr21.*`, `mage_chr22.*`) | 729 individuals × 256,645 variants with MAF ≥ 5%, unphased, indexed pvar |
| `pheno/mage_expr.chr21_22.regenie.tsv.gz` | expression, 729 samples × 367 genes (IDs `ENSG...:GENE`) |
| `pheno/mage_covariates.regenie.tsv.gz` | sex (`XX`/`XY`, categorical: anqchor passes `--cat-cov-names sex` and qpgen encodes it), PC1-5, PEER1-60: 66 covariates |
| `pheno/mage_genes.chr21_22.bed` | TSS loci of the 367 genes |
| `reference/gtex_lcl_eqtl.chr21_22.tsv.gz` | 367 known cis-eQTL pairs (`trait variant beta se pval`): GTEx eGenes (q ≤ 0.05) with their sentinel variant, 346 distinct variants |

Details, a `head` of each file and how the data was made:
[docs/use-cases/example-data.md](../../docs/use-cases/example-data.md).

The reference pairs were discovered in GTEx donors, not in these samples, so no analysis
setup is favoured by the way they were chosen and identity matches are not inflated by
the reference fitting these samples' own noise. Use an independent reference with your own
data too.

## Three separate analyses

Each analysis is its own command with its own run directory. Run from the repository
root. The examples use flags only; the configuration files in this directory are
optional (see below).

```bash
D=examples/mage/data
```

### Sample identity (`match`)

Does each expression sample come from the individual it is labelled with?

```bash
anqchor match all \
  --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --out results/match
```

About 8 s. `match all` = `match run` (analysis and summary tables) followed by
`match report` (writes `results/match/report/match.html`). To re-render the report only,
in about a second: `anqchor match report --out results/match`.

Expected result: 727 of 729 samples confirmed (14 only leniently); 2 unclear. `NA11919`'s
expression does not match its own genotypes (Z 0.41, rank 253 of 729); its best candidate,
`HG00237` (Z 3.60), leads the next by only 0.32, so no new ID is named. `NA12717` is unclear
because every candidate, its own included, is at background level. With `--mahalanobis`:
724 strict, 4 lenient, only `NA11919` unclear.

### cis evaluation (`cis`)

How many known cis pairs replicate with a given set of inputs? Run `cis run` once per
setup (each in its own directory; the runs are independent and can run in parallel),
then compare the runs:

```bash
# no covariates
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --out results/cis/nocov

# 6 covariates
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --cov-include sex,PC1-5 \
  --out results/cis/basic

# 16 covariates
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --cov-include sex,PC1-5,PEER1-10 \
  --out results/cis/peer10

# all 66 covariates
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --out results/cis/peer60

# all 66 covariates, RINT after covariate adjustment
anqchor cis run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --rint \
  --out results/cis/peer60_rint

# compare the five runs
anqchor cis compare results/cis/nocov results/cis/basic results/cis/peer10 results/cis/peer60 \
  results/cis/peer60_rint --out results/cis/compare
```

A run tests one phenotype transform: the expression as given, or with `--rint` its rank-based
inverse normal transform after covariate adjustment.

Each `cis run` takes about 6 s; `cis compare` about 1 s and writes
`results/cis/compare/report/cis.html` plus the comparison tables. Runs are labelled by
their directory names (`--label` overrides). `anqchor cis report --out DIR` re-renders
the report of a comparison or of a single run; `anqchor cis all` is `cis run` followed by
`cis report` on that run.

Expected replication (367 pairs tested) and effective sample size relative to nocov:

| Run | p < 1e-5 | Bonferroni | p < 0.05 | same sign | relative N_eff |
|---|---|---|---|---|---|
| nocov | 178 | 196 | 273 | 339 | 1 |
| basic | 162 | 189 | 266 | 343 | 0.83 |
| peer10 | 231 | 246 | 304 | 352 | 1.58 |
| peer60 | 251 | 269 | 318 | 353 | 2.27 |
| peer60_rint | 251 | 268 | 319 | 353 | 2.10 |

### trans evaluation (`trans`)

Do the known pairs replicate (cis and trans separately), do hotspots appear, and is a
trans scan calibrated? As for cis, run `trans run` once per covariate set, each into its
own directory. trans also needs the gene loci (`--trait-annot`) to tell cis from trans
pairs:

```bash
# no covariates
anqchor trans run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --trait-annot $D/pheno/mage_genes.chr21_22.bed \
  --out results/trans/nocov

# 6 covariates
anqchor trans run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --trait-annot $D/pheno/mage_genes.chr21_22.bed \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --cov-include sex,PC1-5 \
  --out results/trans/basic

# 16 covariates
anqchor trans run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --trait-annot $D/pheno/mage_genes.chr21_22.bed \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --cov-include sex,PC1-5,PEER1-10 \
  --out results/trans/peer10

# all 66 covariates
anqchor trans run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --trait-annot $D/pheno/mage_genes.chr21_22.bed \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --out results/trans/peer60

# all 66 covariates, RINT after covariate adjustment
anqchor trans run --geno-list $D/geno/mage.pgen_list.tsv \
  --pheno $D/pheno/mage_expr.chr21_22.regenie.tsv.gz \
  --pairs $D/reference/gtex_lcl_eqtl.chr21_22.tsv.gz \
  --trait-annot $D/pheno/mage_genes.chr21_22.bed \
  --cov $D/pheno/mage_covariates.regenie.tsv.gz \
  --rint \
  --out results/trans/peer60_rint
```

Each run takes about 8 s. It tests the reference's 346
variants × its 367 genes (126,982 cells: 6,186 cis, 53,508 on different chromosomes, 67,288
on the same chromosome more than 1 Mb apart) and writes its summary tables; `anqchor trans
report --out results/trans/peer10` renders the page of one run. Then put the runs side by
side, and re-render the page later if needed:

```bash
anqchor trans compare results/trans/nocov results/trans/basic results/trans/peer10 \
  results/trans/peer60 results/trans/peer60_rint --out results/trans/compare
anqchor trans report --out results/trans/compare   # re-render the page only
```

`trans compare` takes about 1 s and writes `results/trans/compare/report/trans.html`;
`trans report` rewrites that page from the comparison's tables.

Expected results (the GTEx reference has only cis pairs; inflation over the 53,508
cross-chromosome cells, of which 0.5 are expected at p < 1e-5):

| Run | cis pairs p < 1e-5 (of 367) | λ_GC | cells p < 1e-5 | hotspots (p < 1e-5, ≥ 5 genes) |
|---|---|---|---|---|
| nocov | 178 | 1.74 | 445 | 37 |
| basic | 162 | 1.02 | 1 | 0 |
| peer10 | 231 | 1.00 | 1 | 0 |
| peer60 | 251 | 1.00 | 0 | 0 |
| peer60_rint | 251 | 1.00 | 0 | 0 |

`--hotspot-min-traits` (default 5) sets the genes a 1-Mb bin (`--hotspot-bin-bp`) needs to
be a hotspot; the report shows every threshold and minimum.

## Optional configuration files

Flags are all you need. The YAML files hold the same settings with comments, if you prefer
a file (paths in them are relative to the file):

| File | Command |
|---|---|
| `match.yaml` | `anqchor match all --config examples/mage/match.yaml --out results/match` |
| `cis.yaml` | `anqchor cis run --config examples/mage/cis.yaml --out results/cis_config` |
| `trans.yaml` | `anqchor trans run --config examples/mage/trans.yaml --out results/trans_config` |

Each writes its own run directory, which `cis compare` / `trans compare` take like any
other. Flags given together with `--config` override the file. The trans analysis is described on
the [trans pipeline page](../../docs/pipelines/trans.md).

## Reproducing a run

Every run directory has `anqchor.resolved.yaml`, the complete settings including those
given as flags, and `run.json` with provenance (commands, anqchor version, qpgen binary
commit and checksum, step times). To repeat a run exactly:

```bash
anqchor cis run --config results/cis/peer10/anqchor.resolved.yaml --out results/cis/peer10_again
```

A comparison directory keeps a copy of each compared run's settings in
`settings/<label>.anqchor.yaml`.

Re-running a command on a finished run directory does nothing; changing a setting
re-runs only the affected steps.

## Rebuilding the data

`make_compact_data.sh` rebuilds `data/` from the full MAGE and NYGC files and the GTEx v10
eGene table for EBV-transformed lymphocytes
(`Cells_EBV-transformed_lymphocytes.v10.eGenes.txt.gz`, from the GTEx portal). It needs
`plink2`, `anqchor` and Python with pandas:

```bash
bash examples/mage/make_compact_data.sh FULL_DATA_DIR GTEX_EGENES [PLINK2]
```

Its header lists the input files it expects. It keeps the 729 individuals with expression
and every variant with MAF ≥ 5% among them; the GTEx eGenes on chr21/22 (q ≤ 0.05) whose
sentinel variant is in those genotypes and whose gene is in the MAGE expression data (367
of 516); and the expression, gene loci and covariates of those genes and samples.
