#!/usr/bin/env bash
# Rebuild the compact MAGE chr21/22 example data in examples/mage/data/ from the full
# example files. Requires plink2, anqchor (for 'anqchor prepare index-pvar') and python3
# with pandas.
#
#   bash examples/mage/make_compact_data.sh FULL_DATA_DIR GTEX_EGENES [PLINK2]
#
# GTEX_EGENES is the GTEx v10 eGene table of EBV-transformed lymphocytes
# (Cells_EBV-transformed_lymphocytes.v10.eGenes.txt.gz): one sentinel (top) variant per gene,
# from donors independent of MAGE. FULL_DATA_DIR holds the original MAGE files:
#   geno/nygc_tgp_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{21,22}.filtered.shapeit2-duohmm-phased.{pgen,psam,pvar.idx.gz}
#   exprs/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.729.gidgname.regenie.tsv.gz
#   exprs/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz
#   exprs/eQTL_covariates.tab.gz
#
# Output (examples/mage/data/):
#   genotypes  the 729 individuals with expression data; every variant with MAF >= 5% in
#              them; unphased; psam with #IID only; pvar reduced to CHROM POS ID REF ALT
#              and indexed for qpgen (pvar.idx.gz + tbi)
#   reference  GTEx LCL eGenes on chr21/22 (q <= 0.05) whose sentinel variant is in the
#              genotypes (same CHROM:POS:REF:ALT) and whose gene is in the MAGE expression
#              data; columns trait variant beta se pval (GTEx slope, slope_se, pval_nominal;
#              the slope refers to ALT); trait = the MAGE expression ID of the gene, matched
#              on the Ensembl ID without version
#   expression the 729 samples x the genes of those reference pairs, regenie layout,
#              values rounded to 6 significant digits
#   genes      TSS loci of those genes (IDs identical to the expression IDs)
#   covariates sex (XX/XY), PC1-5, PEER1-60 for the 729 samples, regenie layout,
#              6 significant digits
set -euo pipefail

FULL=${1:?usage: make_compact_data.sh FULL_DATA_DIR GTEX_EGENES [PLINK2]}
GTEX=${2:?usage: make_compact_data.sh FULL_DATA_DIR GTEX_EGENES [PLINK2]}
PLINK2=${3:-plink2}
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=$HERE/data
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
rm -rf "$OUT"
mkdir -p "$OUT/geno" "$OUT/pheno" "$OUT/reference"

# 1. the 729 individuals with expression data (IID column of the regenie expression file)
{ gzip -dc "$FULL/exprs/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.729.gidgname.regenie.tsv.gz" || true; } \
  | cut -f2 | tail -n +2 > "$TMP/samples.txt"

# 2. genotypes: those individuals, variants with MAF >= 5% among them
for c in 21 22; do
  P=$FULL/geno/nygc_tgp_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${c}.filtered.shapeit2-duohmm-phased
  # plink2 needs a plain pvar: take the first five columns of the indexed pvar
  { printf '#CHROM\tPOS\tID\tREF\tALT\n'; gzip -dc "$P.pvar.idx.gz" | grep -v '^#' | cut -f1-5; } > "$TMP/in.pvar"
  "$PLINK2" --pgen "$P.pgen" --pvar "$TMP/in.pvar" --psam "$P.psam" \
            --keep "$TMP/samples.txt" --maf 0.05 \
            --make-pgen erase-phase psam-cols=maybefid --out "$OUT/geno/mage_chr$c"
  rm -f "$OUT/geno/mage_chr$c.log"
  # CHROM:POS:REF:ALT of every kept variant, used to restrict the reference below
  awk -F'\t' '!/^#/ { print $1 ":" $2 ":" $4 ":" $5 }' "$OUT/geno/mage_chr$c.pvar" >> "$TMP/geno_variants.txt"
  anqchor prepare index-pvar --pfile "$OUT/geno/mage_chr$c" --force
  rm -f "$OUT/geno/mage_chr$c.pvar"
done
printf '21\t1\t1000000000\tmage_chr21\n22\t1\t1000000000\tmage_chr22\n' > "$OUT/geno/mage.pgen_list.tsv"

# 3. reference, expression, gene loci and covariates restricted to what is needed
python3 - "$FULL" "$OUT" "$TMP" "$GTEX" <<'PY'
import gzip
import sys
import pandas as pd
full, out, tmp, gtex = sys.argv[1] + "/", sys.argv[2] + "/", sys.argv[3] + "/", sys.argv[4]
gz = {"method": "gzip", "mtime": 0}
expr = full + "exprs/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.729.gidgname.regenie.tsv.gz"

geno = set(open(tmp + "geno_variants.txt").read().split())
with gzip.open(expr, "rt") as f:                                  # MAGE expression IDs ENSG...v:NAME
    mage_ids = {h.split(":")[0].split(".")[0]: h for h in f.readline().rstrip("\n").split("\t")[2:]}

eg = pd.read_csv(gtex, sep="\t")
eg = eg[eg["gene_chr"].isin(["chr21", "chr22"]) & (eg["qval"] <= 0.05)].copy()
n_all = len(eg)
eg["variant"] = (eg["chr"].str.replace("chr", "", regex=False) + ":" + eg["variant_pos"].astype(str)
                 + ":" + eg["ref"] + ":" + eg["alt"])
eg["trait"] = eg["gene_id"].str.split(".").str[0].map(mage_ids)
eg = eg[eg["variant"].isin(geno) & eg["trait"].notna()]
pairs = eg.rename(columns={"slope": "beta", "slope_se": "se", "pval_nominal": "pval"})
pairs = pairs[["trait", "variant", "beta", "se", "pval"]]
pairs.to_csv(out + "reference/gtex_lcl_eqtl.chr21_22.tsv.gz", sep="\t", index=False, compression=gz)
genes = list(dict.fromkeys(pairs["trait"]))                      # reference genes, file order

reg = pd.read_csv(expr, sep="\t", dtype={"FID": str, "IID": str}, usecols=lambda c: c in ("FID", "IID", *genes))
reg = reg[["FID", "IID", *genes]]
reg.to_csv(out + "pheno/mage_expr.chr21_22.regenie.tsv.gz", sep="\t", index=False,
           compression=gz, float_format="%.6g")

bed = pd.read_csv(full + "exprs/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz", sep="\t", usecols=[0, 1, 2, 3])
bed.columns = ["#chrom", "start", "end", "id"]
short = {g.split(":")[0]: g for g in genes}                      # ENSG...5 -> ENSG...5:GENE
bed = bed[bed["id"].isin(short)].copy()
bed["id"] = bed["id"].map(short)
bed["#chrom"] = bed["#chrom"].str.replace("chr", "", regex=False)  # genotype chromosome style
bed.to_csv(out + "pheno/mage_genes.chr21_22.bed", sep="\t", index=False)

cov = pd.read_csv(full + "exprs/eQTL_covariates.tab.gz", sep="\t", index_col=0, dtype=str).T
cov = cov.loc[list(reg["IID"])]
num = [c for c in cov.columns if c != "sex"]
cov[num] = cov[num].astype(float)
cov.insert(0, "IID", cov.index)
cov.insert(0, "FID", cov.index)
cov.to_csv(out + "pheno/mage_covariates.regenie.tsv.gz", sep="\t", index=False,
           compression=gz, float_format="%.6g")
print(f"genotype variants {len(geno)}; reference {len(pairs)} of {n_all} GTEx chr21/22 eGenes "
      f"({pairs['variant'].nunique()} variants); expression {len(reg)} samples x {len(genes)} genes "
      f"({len(bed)} with loci); covariates {cov.shape[1] - 2}")
PY
echo "done: $OUT"
