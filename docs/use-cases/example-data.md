# Example data

The repository ships a compact copy of a real xQTL dataset in `examples/mage/data/`
(32 MB). Every example in this documentation uses it, so you can run each command as
written right after [installation](../getting-started/installation.md).

## What it is

- **Expression:** RNA-seq of lymphoblastoid cell lines (LCLs) from the MAGE study
  (Multi-ancestry Analysis of Gene Expression), 729 individuals of the 1000 Genomes
  Project, normalised (TMM, then inverse normal transform per gene).
- **Genotypes:** the New York Genome Center (NYGC) high-coverage whole-genome sequencing
  of the 1000 Genomes samples.
- **Reference:** cis-eQTLs from an independent study, GTEx v10: for each gene with a
  significant eQTL (an eGene) in EBV-transformed lymphocytes, its sentinel (strongest)
  variant. GTEx donors are not MAGE individuals.
- **Region:** chromosomes 21 and 22 only.

| What | Size |
|---|---|
| Individuals (genotypes and expression) | 729 |
| Genotyped variants (MAF ≥ 5%) | 256,645 (127,563 on chr21, 129,082 on chr22) |
| Genes (expression and loci) | 367 |
| Covariates | 66: sex, PC1-5, PEER1-60 |
| Known cis-eQTL pairs | 367 (one per gene; 346 distinct variants) |

!!! note "Why the reference comes from another study"
    The known eQTLs were discovered in GTEx donors, not in these 729 individuals. That
    makes the comparisons meaningful: a reference discovered in the same samples would fit
    their own noise, make identity matches look stronger than they are and favour the
    covariate setup used for its discovery. Use an independent reference with your own
    data too.

## Files

```text
examples/mage/data/
├── geno/
│   ├── mage.pgen_list.tsv            pgen list: which file holds which chromosome
│   ├── mage_chr21.pgen               genotypes, chr21 (14.6 MB)
│   ├── mage_chr21.psam               sample IDs
│   ├── mage_chr21.pvar.idx.gz        indexed variant list (1.3 MB)
│   ├── mage_chr21.pvar.idx.gz.tbi    tabix index of the variant list
│   └── mage_chr22.*                  the same for chr22 (15.2 MB pgen)
├── pheno/
│   ├── mage_expr.chr21_22.regenie.tsv.gz      expression, 729 × 367 (0.6 MB)
│   ├── mage_covariates.regenie.tsv.gz         covariates, 729 × 66 (0.2 MB)
│   └── mage_genes.chr21_22.bed                gene TSS loci (17 KB)
└── reference/
    └── gtex_lcl_eqtl.chr21_22.tsv.gz          known cis-eQTL pairs from GTEx (12 KB)
```

### Genotypes

The **pgen list** tells qpgen which PLINK 2 file set holds which chromosome. Columns are
chromosome, first position, last position and the file prefix (relative prefixes are
looked up next to the list):

```console
$ cat examples/mage/data/geno/mage.pgen_list.tsv
21	1	1000000000	mage_chr21
22	1	1000000000	mage_chr22
```

The **psam** lists the samples, one ID per line:

```console
$ head -4 examples/mage/data/geno/mage_chr21.psam
#IID
HG00096
HG00100
HG00105
```

The **indexed pvar** (`.pvar.idx.gz` with its `.tbi`) lists the variants with an
`INDEX` column giving each variant's row in the `.pgen`. qpgen looks variants up through
the tabix index instead of reading the whole list. `anqchor prepare index-pvar` makes
one from an ordinary `.pvar`.

```console
$ gzip -dc examples/mage/data/geno/mage_chr21.pvar.idx.gz | head -5
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	INDEX
21	5031747	21:5031747:G:GA	G	GA	.	.	.	1
21	5031747	21:5031748:GA:G	GA	G	.	.	.	2
21	5083315	21:5083315:A:C	A	C	.	.	.	3
21	5083679	21:5083679:G:A	G	A	.	.	.	4
```

The `.pgen` files are binary. The genotypes are unphased hard calls.

### Expression

A samples × genes matrix in regenie layout: `FID`, `IID`, then one column per gene. Gene
IDs are `ENSEMBL_ID:GENE_NAME`.

```console
$ gzip -dc examples/mage/data/pheno/mage_expr.chr21_22.regenie.tsv.gz | head -3 | cut -f1-5
FID	IID	ENSG00000275993.3:SIK1B	ENSG00000278931.1:bP-2189O9.2	ENSG00000279720.1:SOWAHCP2
HG00096	HG00096	-0.238523	0.224454	2.54496
HG00100	HG00100	-0.217436	-0.309645	-0.0685407
```

The sample IDs are the genotype IDs, so no sample map is needed.

### Covariates

The same layout, one column per covariate. `sex` is text (`XX`/`XY`), so it is a
categorical covariate: anqchor keeps the text levels and passes `--cat-cov-names sex` to
qpgen, which replaces it by one indicator column (`XX` is the reference level).

```console
$ gzip -dc examples/mage/data/pheno/mage_covariates.regenie.tsv.gz | head -3 | cut -f1-9
FID	IID	sex	PC1	PC2	PC3	PC4	PC5	PEER1
HG00096	HG00096	XY	-0.0216243	0.0487524	0.00362148	0.0446186	0.00187148	-0.0565679
HG00100	HG00100	XX	-0.0215191	0.0503942	0.00281614	0.040869	0.000118772	0.0722641
```

The columns after `PEER1` are `PEER2` to `PEER60`. With `--cov-include` you can use a
subset, for example `--cov-include sex,PC1-5,PEER1-10` (ranges expand).

Non-numeric columns are categorical automatically. A covariate whose categories are coded as
numbers, for example a sequencing batch `1`, `2`, `3`, would be read as one numeric
covariate; name it with `--cat-cov-names batch` (config `cov_categorical: [batch]`) to treat
it as categorical. See [covariates](../user-guide/inputs.md#covariates).

### Gene loci

A BED-like file with each gene's transcription start site. Gene IDs are the same as in
the expression file. trans evaluation uses it to tell cis cells from trans cells; cis
evaluation needs it only with `--assoc-engine rect`.

```console
$ head -4 examples/mage/data/pheno/mage_genes.chr21_22.bed
#chrom	start	end	id
21	6111130	6111131	ENSG00000275993.3:SIK1B
21	8857259	8857260	ENSG00000278931.1:bP-2189O9.2
21	9083100	9083101	ENSG00000279720.1:SOWAHCP2
```

### Reference pairs

The known cis-eQTLs: trait, variant (`CHROM:POS:REF:ALT`), effect size, its standard
error and p-value. Each row is one GTEx eGene with its sentinel variant; `beta`, `se` and
`pval` are GTEx's `slope`, `slope_se` and `pval_nominal`, and `beta` is the effect of one
copy of the `ALT` allele. The trait is the gene's ID in the expression file (matched on the
Ensembl ID without its version).

```console
$ gzip -dc examples/mage/data/reference/gtex_lcl_eqtl.chr21_22.tsv.gz | head -4
trait	variant	beta	se	pval
ENSG00000275993.3:SIK1B	21:6122794:C:T	-0.448959	0.0888711	8.054e-07
ENSG00000278931.1:bP-2189O9.2	21:9370164:G:A	0.484109	0.10423	5.31901e-06
ENSG00000279720.1:SOWAHCP2	21:9088050:G:A	0.708996	0.180329	0.000107185
```

Every reference variant is present in the genotypes, so all 367 pairs are tested. The GTEx
p-values range from 1e-102 to 1e-4 (median 8e-11): some eGenes are strong, others barely
significant in GTEx, as in any real reference.

## What the example shows

- **Identity matching:** 727 of 729 samples are confirmed (14 only leniently). The
  expression sample labelled `NA11919` does not match its own genotypes (rank 253 of 729),
  and no other individual matches it clearly either, so the result is "unclear" rather than
  a named mismatch. Both `NA11919` and its best (but weak) candidate, `HG00237`, are in the
  data (genotypes and expression), so you can inspect the case yourself
  ([identity matching](identity-matching.md)).
- **cis evaluation:** the known pairs replicate much better with PEER factors: 178 of 367 at
  p < 10⁻⁵ without covariates, 251 with all 66 covariates
  ([cis evaluation](cis-strategies.md)).

## How it was made

`examples/mage/make_compact_data.sh` rebuilds `examples/mage/data/` from the full MAGE
and NYGC files and the GTEx v10 eGene table for EBV-transformed lymphocytes
(`Cells_EBV-transformed_lymphocytes.v10.eGenes.txt.gz`, part of the single-tissue cis-eQTL
data on the GTEx portal). It needs `plink2`, `anqchor` and Python with pandas:

```bash
bash examples/mage/make_compact_data.sh FULL_DATA_DIR GTEX_EGENES [PLINK2]
```

`FULL_DATA_DIR` holds the MAGE and NYGC files, `GTEX_EGENES` is the path of the GTEx eGene
table and `PLINK2` the plink2 binary (default: `plink2` on your `PATH`). For example:

```bash
bash examples/mage/make_compact_data.sh /data/mage \
    /data/gtex/Cells_EBV-transformed_lymphocytes.v10.eGenes.txt.gz
```

What it keeps, and why:

| Part | Kept | Why |
|---|---|---|
| Genotypes | the 729 individuals with expression; **every** variant with MAF ≥ 5% among them; unphased; psam with `#IID` only; pvar reduced to `CHROM POS ID REF ALT`, then indexed | any common variant can be queried, not just the reference variants, so trans rectangles and your own variant lists work |
| Reference | the GTEx eGenes on chr21/22 with q ≤ 0.05 whose sentinel variant is in those genotypes (same `CHROM:POS:REF:ALT`) and whose gene is in the MAGE expression data (same Ensembl ID, ignoring the version); columns `trait variant beta se pval` | an independent reference in the same cell type, and every pair can be tested |
| Expression | the 729 samples × the genes of those reference pairs; values rounded to 6 significant digits | the analyses only use genes with a known pair |
| Gene loci | the TSS of those genes, with the expression IDs and the genotypes' chromosome names (`21`, not `chr21`) | IDs and chromosomes match the other files without any mapping |
| Covariates | sex, PC1-5, PEER1-60 for the 729 samples, regenie layout, 6 significant digits | the covariates used by the MAGE study |

Unused columns (the other columns of the GTEx table; the extra pvar columns) are dropped.

!!! note "Provenance of the full files"
    The script starts from the NYGC phased chr21/22 genotypes of 3,202 1000 Genomes
    samples, the MAGE expression matrix of 729 samples × 20,154 genes with its tensorQTL
    BED, the MAGE covariate table (one row per covariate) and the GTEx v10 eGene table of
    EBV-transformed lymphocytes. GTEx tests 913 genes on chr21/22 and calls 516 of them
    eGenes at q ≤ 0.05. Of those 516, 367 are kept: 72 genes are not in the MAGE expression
    data, 67 have a sentinel variant with MAF below 5% in the 729 individuals, and 10 have
    a sentinel variant that the genotype release does not contain under the same
    `CHROM:POS:REF:ALT`.
