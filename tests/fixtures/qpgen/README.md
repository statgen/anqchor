# qpgen output fixtures

Excerpts of real `qpgentools` outputs (dev branch with the fixes of DESIGN.md 3.5),
produced on 2026-10-11 from the original, full-size MAGE chr21/22 files that
`examples/mage/make_compact_data.sh` compacts (pgen list rewritten to absolute paths, no
covariates). They are frozen inputs for the output readers; the live counterparts are
regenerated on the in-repo example data by `tests/unit/test_qpgen_live.py` and
`tests/integration/test_mage.py`.

| File | Command |
|---|---|
| `prs.prs.tsv.gz`, `prs.se.tsv.gz` | `pair-prs --pairs pair/MAGE_eQTL_top_cis.chr21_22.tsv.gz`; first 8 samples x first 6 traits |
| `match.*.tsv.gz` | `match-prs-pheno --rint` (independence); first rows plus `SELF_LENIENT`/`LENIENT_MATCH`, `UNCLEAR`/`NO_MATCH` rows |
| `matchauto.match.all.tsv.gz`, `matchauto.log` | `match-prs-pheno --rint --mahalanobis --auto-lambda`; the `SINGLE_NEW_BEST` row (NA11919 -> HG00237) and the lambda log lines |
| `pair_assoc.tsv.gz` | `pair-assoc --pairs pairs.sub.tsv` (10 pairs; `21:9357536:G:GTT` is absent from the pvar, so 9 rows of 20 fields) |
| `pair_assoc.short.tsv.gz` | the same with the raw `LOG10P` field deleted from every data row (19 fields; the unpatched dev-build layout) |
| `rect_assoc.tsv.gz` | `rect-assoc --pheno-list rect_traits.txt --var-list rect_variants.txt --min-maf 0.01` (6 present variants x 5 traits) |
| `rect_assoc.doubledN.tsv.gz` | the same with `N` doubled (the unpatched dev-build allele number) |
