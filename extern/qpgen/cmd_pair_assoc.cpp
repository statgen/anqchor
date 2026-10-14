#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "assoc_utils.h"
#include "qpgen_utils.h"
#include "qpgen.h"
#include "pheno.h"
#include "Eigen/Dense"
#include <cmath>
#include <algorithm>
#include <memory>

// pair-assoc model
// ----------------
// Each (trait, variant) pair is tested by least squares of y ~ 1 + C + g (C = --cov, empty
// without it), fitted by Frisch-Waugh-Lovell: y and g are both residualized on [1, C] and
// the slope of the residuals is tested with df = N - 2 - n_cov_df, where n_cov_df is the rank
// of [1, C] minus 1. Without --cov the intercept alone is fitted, i.e. y and g are centered.
//
// The analysis samples of a test are those with an observed phenotype and (by default) a
// called genotype (complete-case analysis), so N and the df always count the same samples.
// The covariate fit is therefore made on the samples of each test: once per trait for the
// samples with an observed phenotype, and once more for a variant with missing genotypes among
// them. Samples with a missing covariate are excluded from all tests. With --mean-impute-geno,
// missing genotypes are set to the mean of the called genotypes and all samples with an
// observed phenotype are kept (as rect-assoc does).
//
// The RINT phenotype is the rank-based inverse normal transform of the covariate-adjusted
// phenotype of the trait's samples (as rect-assoc --rint-after-adj); it is not orthogonal to
// the covariates, so its residual sum of squares on [1, C] replaces its y'y in the Wald test.
namespace {

// Rows of a matrix (or vector) selected by index
Eigen::MatrixXd select_rows(const Eigen::MatrixXd& m, const std::vector<int32_t>& rows) {
    Eigen::MatrixXd out(rows.size(), m.cols());
    for(size_t r = 0; r < rows.size(); ++r) out.row(r) = m.row(rows[r]);
    return out;
}

Eigen::VectorXd select_rows(const Eigen::VectorXd& v, const std::vector<int32_t>& rows) {
    Eigen::VectorXd out(rows.size());
    for(size_t r = 0; r < rows.size(); ++r) out(r) = v(rows[r]);
    return out;
}

void set_nan(slr_sumstat_t& ss, int32_t n_obs) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    ss.beta = ss.se = ss.tstat = ss.log10p = nan;
    ss.n_obs = n_obs;
}

// Samples, covariate fit and phenotypes of the trait under test
struct TraitData {
    int32_t phe_idx = -1;
    std::vector<int32_t> rows;          // overlapping-sample indices with an observed phenotype
    const Eigen::MatrixXd* cov = NULL;  // covariates on rows (n_t x k; k = 0 without --cov)
    std::shared_ptr<const CovariateResidualizer> resid; // fit of [1, C] on rows
    Eigen::MatrixXd own_cov;            // storage for cov when rows is a strict subset
    Eigen::VectorXd y;                  // phenotype residualized on [1, C]
    Eigen::VectorXd y_rint;             // RINT of y
    double rint_rss = -1.0;             // residual sum of squares of y_rint on [1, C]
    bool testable = false;              // enough samples for the Wald test

    // LOCO (--loco/--pred): y and y_rint depend on the chromosome; see apply_trait_loco()
    std::string loco_path;              // .loco file of the trait (empty without LOCO)
    double loco_sd = 0;                 // SD of the covariate residual (REGENIE's scale_Y)
    Eigen::VectorXd loco_base;          // covariate residual divided by loco_sd
    std::string loco_chrom;             // chromosome whose predictions y holds

    int32_t n() const { return (int32_t)rows.size(); }
};

// cov_all and resid_all are the covariates and their fit on all overlapping samples, shared
// by the traits without a missing phenotype
// loco_path: the trait's .loco file, or empty without LOCO. Samples absent from it are excluded.
void prepare_trait(TraitData& td, const PhenoMatrix& pheno_matrix, int32_t phe_idx,
                   const Eigen::MatrixXd& cov_all, const std::shared_ptr<const CovariateResidualizer>& resid_all,
                   bool skip_rint, const std::string& loco_path)
{
    td.phe_idx = phe_idx;
    td.rows.clear();
    td.loco_path = loco_path;
    td.loco_chrom.clear();
    const int32_t n_all = (int32_t)pheno_matrix.samp_ids.size();
    for(int32_t i = 0; i < n_all; ++i) {
        if ( pheno_matrix.pheno_mask(i, phe_idx) ) td.rows.push_back(i);
    }
    const char* phe_id = pheno_matrix.pheno_ids[phe_idx].c_str();
    if ( td.n() < n_all ) {
        notice("Trait %s: %d of %d samples have a missing phenotype and are excluded", phe_id, n_all - td.n(), n_all);
    }
    if ( !loco_path.empty() ) {
        loco_row_t ids;
        load_loco_row(loco_path, "", ids); // sample IDs only
        std::vector<int32_t> kept;
        for(int32_t r : td.rows) {
            if ( ids.has(pheno_matrix.samp_ids[r]) ) kept.push_back(r);
        }
        if ( kept.size() < td.rows.size() ) {
            warning("Trait %s: %zu of %d samples have no LOCO prediction in %s and are excluded",
                    phe_id, td.rows.size() - kept.size(), td.n(), loco_path.c_str());
        }
        td.rows.swap(kept);
    }
    if ( td.n() == n_all ) {
        td.cov = &cov_all;
        td.resid = resid_all;
    }
    else {
        td.own_cov = select_rows(cov_all, td.rows);
        td.cov = &td.own_cov;
        td.resid = std::make_shared<const CovariateResidualizer>(td.own_cov);
    }
    td.testable = ( td.n() - 2 - td.resid->n_cov_df() > 0 );
    if ( !td.testable ) {
        warning("Trait %s has %d samples with an observed phenotype, too few for the intercept, %d covariate degrees of freedom and the variant; its pairs are skipped",
                phe_id, td.n(), td.resid->n_cov_df());
        return;
    }
    td.y = select_rows(Eigen::VectorXd(pheno_matrix.pheno_mat.col(phe_idx)), td.rows);
    td.resid->residualize(td.y);
    if ( !loco_path.empty() ) {
        // REGENIE's scale_Y: SD = ||r|| / sqrt(n - n_cov), n_cov counting the intercept
        td.loco_sd = td.y.norm() / std::sqrt((double)(td.n() - 1 - td.resid->n_cov_df()));
        if ( !(td.loco_sd > 1e-12) ) {
            error("Covariate-adjusted phenotype %s has zero variance; cannot apply LOCO", phe_id);
        }
        td.loco_base = td.y / td.loco_sd; // y and y_rint are set per chromosome
        return;
    }
    if ( !skip_rint ) {
        td.y_rint = rint_without_missing(td.y);
        td.rint_rss = td.resid->residual_ss(td.y_rint)(0);
    }
}

// LOCO, as in region-assoc and REGENIE step 2: subtract the chromosome's prediction from the
// standardized covariate residual. y is the result in phenotype units (times loco_sd), and
// y_rint its RINT (REGENIE --apply-rerint). Neither is orthogonal to the covariates: y is
// residualized again (so that y'y is its residual sum of squares) and rint_rss is recomputed;
// genotypes are residualized, so x'y and BETA are unchanged by this. No-op without LOCO or
// when chrom is already applied.
void apply_trait_loco(TraitData& td, const PhenoMatrix& pheno_matrix, const std::string& chrom, bool skip_rint)
{
    if ( td.loco_path.empty() ) return;
    if ( !td.loco_chrom.empty() && loco_chrom_label(td.loco_chrom) == loco_chrom_label(chrom) ) return;
    loco_row_t row;
    load_loco_row(td.loco_path, chrom, row);
    Eigen::VectorXd t(td.n());
    for(int32_t r = 0; r < td.n(); ++r) {
        double v;
        if ( !row.get(pheno_matrix.samp_ids[td.rows[r]], v) ) {
            error("Missing LOCO prediction for sample %s and trait %s on chromosome %s in %s",
                  pheno_matrix.samp_ids[td.rows[r]].c_str(), pheno_matrix.pheno_ids[td.phe_idx].c_str(), chrom.c_str(), td.loco_path.c_str());
        }
        t(r) = td.loco_base(r) - v;
    }
    td.y = t * td.loco_sd;
    td.resid->residualize(td.y);
    if ( !skip_rint ) {
        td.y_rint = rint_without_missing(t);
        td.rint_rss = td.resid->residual_ss(td.y_rint)(0);
    }
    td.loco_chrom = chrom;
}

// Genotypes of the variants tested against one trait, on the trait's samples (TraitData::rows).
// Hard calls are centered by the mean of the called genotypes; missing genotypes (hard calls or
// dosages) are set to that mean (0 after centering) and flagged in the mask.
struct GenoChunk {
    Eigen::MatrixXd geno;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> mask;
    std::vector<cpra_t> cpras;
    std::vector<int32_t> ans;           // number of called alleles
    std::vector<double> acs;            // ALT allele counts (dosages)
    std::vector<int32_t> gc0s, gc1s, gc2s; // RR, RA, AA counts
    std::vector<double> infos;

    int32_t size() const { return (int32_t)cpras.size(); }

    void reset(int32_t n_rows) {
        geno.resize(n_rows, 16);
        mask.resize(n_rows, 16);
        cpras.clear(); ans.clear(); acs.clear();
        gc0s.clear(); gc1s.clear(); gc2s.clear(); infos.clear();
    }

    // Adds the variant whose genotypes mpr has just read. Returns false (and adds nothing) for
    // a variant that is monomorphic, or has no called genotype, among the trait's samples.
    bool add_variant(const MultiPgenIdxReader& mpr, const cpra_t& cpra, const std::vector<int32_t>& rows) {
        const int32_t n = (int32_t)rows.size();
        const int32_t j = size();
        if ( j >= geno.cols() ) {
            geno.conservativeResize(n, geno.cols() * 2);
            mask.conservativeResize(n, mask.cols() * 2);
        }
        int32_t an = 0;
        double ac = 0, sumsq = 0;
        int32_t gcs[3] = {0, 0, 0};
        const bool is_dosage = mpr.is_dosage_present();
        if ( is_dosage ) {
            const double* dbl_buf = mpr.get_dbl_buf(); // per-pgen buffer: fetch after each read
            for(int32_t r = 0; r < n; ++r) {
                const double d = dbl_buf[rows[r]];     // REF dosage; missing = -3
                if ( d < 0 ) { mask(r, j) = false; continue; }
                const double ds = 2.0 - d;             // ALT dosage
                geno(r, j) = ds;
                mask(r, j) = true;
                ac += ds;
                an += 2;
                sumsq += ds * ds;
                ++gcs[(ds < 0.5) ? 0 : ( (ds < 1.5) ? 1 : 2 )];
            }
        }
        else {
            const std::vector<int32_t>& int_buf = mpr.get_int_buf(); // REF count; missing = -9
            for(int32_t r = 0; r < n; ++r) {
                switch( int_buf[rows[r]] ) {
                case 0: geno(r, j) = 2.0; mask(r, j) = true; an += 2; ac += 2; ++gcs[2]; break;
                case 1: geno(r, j) = 1.0; mask(r, j) = true; an += 2; ac += 1; ++gcs[1]; break;
                case 2: geno(r, j) = 0.0; mask(r, j) = true; an += 2;          ++gcs[0]; break;
                default: mask(r, j) = false; break;
                }
            }
        }
        if ( an == 0 || ac == 0 || ac == an ) return false; // monomorphic or all missing

        const double mean = ac / (double)an * 2.0;
        double info;
        if ( is_dosage ) {
            // Var(g) / E(Var(g)) with E(Var(g)) = 2 af (1 - af)
            info = ( sumsq / 2.0 * an - 4.0 * ac * ac ) / ( 2.0 * ac * ( an - ac ) );
        }
        else {
            info = ((4.0 * gcs[2] + gcs[1]) / (an / 2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
        }
        for(int32_t r = 0; r < n; ++r) {
            geno(r, j) = mask(r, j) ? geno(r, j) - mean : 0.0;
        }
        cpras.push_back(cpra);
        ans.push_back(an);
        acs.push_back(ac);
        gc0s.push_back(gcs[0]);
        gc1s.push_back(gcs[1]);
        gc2s.push_back(gcs[2]);
        infos.push_back(info);
        return true;
    }
};

// Tests every variant of the chunk against the trait and writes one line per variant.
// Returns the number of variants whose genotype the covariates explain (reported as NaN).
int64_t test_and_write(htsFile* wf, const char* phe_id, const TraitData& td, GenoChunk& gc,
                       bool skip_rint, bool mean_impute_geno)
{
    const int32_t n = td.n();
    const int32_t m = gc.size();
    if ( m == 0 ) return 0;
    std::vector<slr_sumstat_t> ss(m), ss_rint(m);
    int64_t n_collinear = 0;

    // variants with a genotype for every sample of the trait share the trait's covariate fit
    std::vector<int32_t> complete, incomplete;
    for(int32_t j = 0; j < m; ++j) {
        if ( mean_impute_geno || gc.mask.col(j).head(n).all() ) complete.push_back(j);
        else incomplete.push_back(j);
    }

    if ( !complete.empty() ) {
        Eigen::MatrixXd G(n, complete.size());
        for(size_t c = 0; c < complete.size(); ++c) G.col(c) = gc.geno.col(complete[c]);
        n_collinear += td.resid->residualize_genotypes(G);
        std::vector<slr_sumstat_t> res, res_rint;
        // y is residualized, so its y'y is the null-model RSS
        simple_linear_regression_without_missing(td.y, G, res, td.resid->n_cov_df());
        if ( !skip_rint ) {
            simple_linear_regression_without_missing(td.y_rint, G, res_rint, td.resid->n_cov_df(), td.rint_rss);
        }
        for(size_t c = 0; c < complete.size(); ++c) {
            ss[complete[c]] = res[c];
            if ( !skip_rint ) ss_rint[complete[c]] = res_rint[c];
        }
    }

    // complete-case analysis: refit the covariates on the samples with a called genotype
    for(int32_t j : incomplete) {
        std::vector<int32_t> obs;
        for(int32_t r = 0; r < n; ++r) {
            if ( gc.mask(r, j) ) obs.push_back(r);
        }
        const int32_t n_obs = (int32_t)obs.size();
        CovariateResidualizer resid(select_rows(*td.cov, obs));
        const int32_t n_cov_df = resid.n_cov_df();
        if ( n_obs - 2 - n_cov_df <= 0 ) {
            set_nan(ss[j], n_obs);
            set_nan(ss_rint[j], n_obs);
            continue;
        }
        Eigen::MatrixXd x(n_obs, 1);
        Eigen::VectorXd y(n_obs), y_rint(skip_rint ? 0 : n_obs);
        for(int32_t k = 0; k < n_obs; ++k) {
            x(k, 0) = gc.geno(obs[k], j);
            y(k) = td.y(obs[k]);
            if ( !skip_rint ) y_rint(k) = td.y_rint(obs[k]);
        }
        n_collinear += resid.residualize_genotypes(x);
        // residualizing y again on the subset makes its y'y the null-model RSS of the subset
        resid.residualize(y);
        std::vector<slr_sumstat_t> res;
        simple_linear_regression_without_missing(y, x, res, n_cov_df);
        ss[j] = res[0];
        if ( !skip_rint ) {
            resid.residualize(y_rint);
            simple_linear_regression_without_missing(y_rint, x, res, n_cov_df);
            ss_rint[j] = res[0];
        }
    }

    for(int32_t j = 0; j < m; ++j) {
        const cpra_t& cpra = gc.cpras[j];
        hprintf(wf, "%s\t%s\t%d\t%s\t%s\t%s\t%.5g\t%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g",
            phe_id,                          // TRAIT
            cpra.chrom.c_str(),              // CHROM
            cpra.pos,                        // GENPOS
            cpra.to_string().c_str(),        // ID
            cpra.ref.c_str(),                // ALLELE0
            cpra.alts.c_str(),               // ALLELE1
            gc.acs[j] / (double)gc.ans[j],   // A1FREQ among the called genotypes
            ss[j].n_obs,                     // N = samples in the regression
            gc.gc0s[j], gc.gc1s[j], gc.gc2s[j], // N_RR, N_RA, N_AA
            gc.infos[j],                     // INFO
            ss[j].beta, ss[j].se, ss[j].tstat, ss[j].log10p);
        if ( !skip_rint ) {
            hprintf(wf, "\t%.6g\t%.6g\t%.6g\t%.6g", ss_rint[j].beta, ss_rint[j].se, ss_rint[j].tstat, ss_rint[j].log10p);
        }
        hprintf(wf, "\n");
    }
    return n_collinear;
}

} // anonymous namespace

int32_t cmd_pair_assoc(int32_t argc, char **argv)
{
    std::string pgenf;      // PLINK2 genotype file
    std::string psamf;      // PLINK2 sample file
    std::string pivarf;     // PLINK2 indexed pvar file, compatible with PgenIdxReader
    std::string pgenlistf;      // CHROM BEG END PGEN PSAM PIVAR containing the list of region-specific BED files
    std::string phef;
    std::string covf;
    std::string cov_names;       // covariates to use (--cov-names); empty = all
    std::string cat_cov_names;   // categorical covariates (--cat-cov-names)
    int32_t max_cat_levels = 10; // REGENIE's --maxCatLevels
    std::string pairf;
    std::string outf;
    std::string samplef;
    std::string missing_str("NA"); // comma-separated strings marking missing phenotype/covariate values
    std::string locof;      // REGENIE step-1 LOCO predictions of a single trait (.loco file)
    std::string predf;      // REGENIE step-1 *_pred.list: [pheno_id] [path_to_loco] per line
    int32_t jump_thres_bp = 1000000;
    int32_t max_chunk_vars = 1000; // Maximum number of variants to store at once in memory
    int32_t icol_pivar_idx = 9;
    std::string pheno_format("regenie");
    std::string cov_format("regenie");
    std::string colname_pheno_sample("pheno"); // column name for the sample IDs in the phenotype file
    std::string colname_geno_sample("geno"); // column name for the sample IDs in the genotype file
    std::string colname_pair_trait("trait"); // column name for the trait ID in the pair file
    std::string colname_pair_variant("variant"); // column name of the variant ID in the pair file
    std::string colname_pair_region("region"); // column name of the region string in the pair file
    bool skip_rint = false; // skip tests based on rank-based inverse normal transformation
    bool mean_impute_geno = false; // mean-impute missing genotypes instead of complete-case analysis

    bool space_delim = false; // space- instead of tab-delimited input tables

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("pgen", &pgenf, "Input PLINK2 genotype file")
    LONG_STRING_PARAM("psam", &psamf, "Input PLINK2 sample file")
    LONG_STRING_PARAM("pivar", &pivarf, "Input PLINK2 index pvar file (bgzipped and tabix)")
    LONG_STRING_PARAM("pgen-list", &pgenlistf, "Input file containing CHROM BEG END PGEN PSAM PIVAR")
    LONG_STRING_PARAM("pheno", &phef, "Input phenotype matrix")
    LONG_STRING_PARAM("pairs", &pairf, "Input file containing [trait_id] [variant_id] pairs to be tested")
    LONG_STRING_PARAM("sample", &samplef, "Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files")
    LONG_STRING_PARAM("cov", &covf, "Input covariate matrix (optional)")
    LONG_STRING_PARAM("pheno-format", &pheno_format, "Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("cov-format", &cov_format, "Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("cov-names", &cov_names, QPGEN_COV_NAMES_HELP)
    LONG_STRING_PARAM("cat-cov-names", &cat_cov_names, QPGEN_CAT_COV_NAMES_HELP)
    LONG_INT_PARAM("max-cat-levels", &max_cat_levels, QPGEN_MAX_CAT_LEVELS_HELP)
    LONG_STRING_PARAM("loco", &locof, "REGENIE step-1 LOCO predictions of a single trait (.loco file). For each chromosome, the covariate-adjusted phenotype is scaled to unit SD, that chromosome's prediction is subtracted, and the result is rescaled to the original SD for BETA/SE (the _RINT statistics use the RINT of the unscaled result, as REGENIE --apply-rerint)")
    LONG_STRING_PARAM("pred", &predf, "REGENIE step-1 *_pred.list file for one or more traits: one '[pheno_id] [path_to_loco]' line per trait (relative paths are resolved against the list's directory). Each tested trait must be listed. Same LOCO adjustment as --loco")
    LONG_STRING_PARAM("missing-str", &missing_str, "Comma-separated strings marking missing phenotype/covariate values (default: 'NA'). Samples with a missing phenotype are excluded from that trait's tests; samples with a missing covariate are excluded from all tests")


    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Auxiliary options", NULL)
    LONG_INT_PARAM("jump-thres-bp", &jump_thres_bp, "Jump threshold in base pairs for the variant index (default: 1000000)")
    LONG_INT_PARAM("max-chunk-vars", &max_chunk_vars, "Maximum number of variants to store at once in memory (default: 1000)")
    LONG_INT_PARAM("icol-pivar-idx", &icol_pivar_idx, "1-based column index for the variant ID in the pvar file (default: 9)")
    LONG_STRING_PARAM("colname-pheno-sample", &colname_pheno_sample, "When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')")
    LONG_STRING_PARAM("colname-geno-sample", &colname_geno_sample, "When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'geno')")
    LONG_STRING_PARAM("colname-pair-trait", &colname_pair_trait, "Column name of the trait ID in the pait file")
    LONG_STRING_PARAM("colname-pair-variant", &colname_pair_variant, "Column name of the variant ID in the pair file")
    LONG_STRING_PARAM("colname-pair-region", &colname_pair_region, "Column name of the region string in the pair file")
    LONG_PARAM("skip-rint", &skip_rint, "Skip tests based on rank-based inverse normal transformation (default: false)")
    LONG_PARAM("mean-impute-geno", &mean_impute_geno, "Replace missing genotypes by the mean of the called genotypes and keep every sample with an observed phenotype (as rect-assoc does). By default, samples with a missing genotype are excluded from that variant's test")
    LONG_PARAM("space-delim", &space_delim, QPGEN_SPACE_DELIM_HELP)
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();
    const char table_delim = table_delimiter(space_delim);
    if ( !cat_cov_names.empty() && covf.empty() ) {
        error("--cat-cov-names requires --cov");
    }
    if ( !cov_names.empty() && covf.empty() ) {
        error("--cov-names requires --cov");
    }

    // open the Plink file
    MultiPgenIdxReader mpr;
    if ( !pgenlistf.empty() ) { // list is provided
        if ( pgenf.empty() && pivarf.empty() && psamf.empty() ) {
            if ( !mpr.prep_pgen_list(pgenlistf.c_str()) ) {
                error("Failed to prepare pgen files with the following list file: %s", pgenlistf.c_str());
            }
        }
        else {
            error("When --pgen-list is provided, --pgen, --pivar, and --psam should not be provided");
        }
    }
    else {
        if ( pgenf.empty() || pivarf.empty() || psamf.empty() ) {
            error("When --pgen-list is not provided, --pgen, --pivar, and --psam should be provided");
        }
        if ( !mpr.set_single_chunk_pgen(pgenf.c_str(), pivarf.c_str(), psamf.c_str()) ) {
            error("Failed to add pgen files with the following files:\n%s\n%s\n%s", pgenf.c_str(), pivarf.c_str(), psamf.c_str());
        }
    }
    mpr.set_jump_thres_bp(jump_thres_bp);
    mpr.set_icol_pivar_idx(icol_pivar_idx - 1); // convert to 0-based index
    mpr.set_mean_impute_dosage(false); // missing dosages are flagged (-3) and handled below

    notice("Prepared to load PLINK2 genotype data with %d samples", mpr.get_all_sample_count());

    // load the sample IDs
    std::vector<std::string> subset_sample_ids;
    if ( !samplef.empty() ) {
        notice("Loading sample IDs from %s", samplef.c_str());
        tsv_reader tr_sample(samplef.c_str());
        tr_sample.delimiter = table_delim;
        while(tr_sample.read_line()) {
            subset_sample_ids.push_back(tr_sample.str_field_at(0));
        }
        notice("Loaded %d sample IDs from %s to subset the genotype and phenotype data", (int32_t)subset_sample_ids.size(), samplef.c_str());
    }

    // find overlapping sample IDs
    const std::vector<plink_samp_t>& geno_all_samps = mpr.get_all_samples();
    std::vector<std::string> geno_all_samp_ids;
    for(int32_t i=0; i < geno_all_samps.size(); ++i) {
        geno_all_samp_ids.push_back(geno_all_samps[i].indID);
    }

    notice("Loading phenotype matrix from %s", phef.c_str());

    // load the phenotype matrix
    PhenoMatrix pheno_matrix;
    pheno_matrix.add_missing_strs(missing_str);
    if ( !pheno_matrix.load_pheno_matrix(phef.c_str(), pheno_format.c_str(), table_delim) ) {
        error("Failed to load the phenotype matrix from file %s", phef.c_str());
    }

    notice("Loaded phenotype matrix with %d samples and %d phenotypes from %s", (int32_t)pheno_matrix.samp_ids.size(), (int32_t)pheno_matrix.pheno_ids.size(), phef.c_str());

    const bool have_cov = !covf.empty();
    PhenoMatrix cov_matrix;
    if ( have_cov ) {
        notice("Loading covariate matrix from %s", covf.c_str());
        cov_matrix.add_missing_strs(missing_str);
        cov_matrix.set_categorical_ids(split_comma_list(cat_cov_names));
        cov_matrix.set_selected_columns(expand_name_ranges(cov_names));
        if ( !cov_matrix.load_pheno_matrix(covf.c_str(), cov_format.c_str(), table_delim) ) {
            error("Failed to load the covariate matrix from file %s", covf.c_str());
        }
        notice("Loaded covariate matrix with %d samples and %d covariates from %s", (int32_t)cov_matrix.samp_ids.size(), (int32_t)cov_matrix.pheno_ids.size(), covf.c_str());
        cov_matrix.keep_covariate_columns(expand_name_ranges(cov_names), covf.c_str());
    }

    // overlap of the sample, genotype, covariate, and phenotype files
    std::vector<std::string> overlapping_sample_ids;
    {
        std::vector<std::string> ids = subset_sample_ids.empty() ? geno_all_samp_ids : std::vector<std::string>();
        if ( !subset_sample_ids.empty() ) identify_overlapping_ids(subset_sample_ids, geno_all_samp_ids, ids);
        if ( have_cov ) {
            std::vector<std::string> tmp;
            identify_overlapping_ids(ids, cov_matrix.samp_ids, tmp);
            ids.swap(tmp);
        }
        identify_overlapping_ids(ids, pheno_matrix.samp_ids, overlapping_sample_ids);
    }
    notice("%zu overlapping samples found among sample, genotype, phenotype, and covariate files", overlapping_sample_ids.size());

    // samples with a missing covariate cannot be adjusted and are excluded from all tests
    if ( have_cov && cov_matrix.has_missing ) {
        std::map<std::string, int32_t> cov_row;
        for(int32_t i = 0; i < (int32_t)cov_matrix.samp_ids.size(); ++i) cov_row[cov_matrix.samp_ids[i]] = i;
        std::vector<std::string> kept_ids;
        for(const std::string& id : overlapping_sample_ids) {
            if ( cov_matrix.pheno_mask.row(cov_row[id]).all() ) kept_ids.push_back(id);
        }
        notice("%zu of %zu overlapping samples have a missing covariate value and are excluded", overlapping_sample_ids.size() - kept_ids.size(), overlapping_sample_ids.size());
        overlapping_sample_ids.swap(kept_ids);
    }
    if ( overlapping_sample_ids.empty() ) {
        error("No samples are shared by the genotype, phenotype%s files%s. Check the sample IDs%s",
              have_cov ? ", covariate" : "", samplef.empty() ? "" : " and the --sample list",
              have_cov ? " and the covariate file layout (--cov-format)" : "");
    }

    // identify_overlapping_ids() returns IDs in lexicographic order, but the pgen reader requires
    // the sample subset in increasing genotype-file order. Put the overlapping IDs in genotype
    // order so that phenotype, covariate, and genotype rows all follow the same order.
    {
        std::map<std::string, int32_t> geno_rank;
        for(int32_t i = 0; i < (int32_t)geno_all_samp_ids.size(); ++i) geno_rank[geno_all_samp_ids[i]] = i;
        std::sort(overlapping_sample_ids.begin(), overlapping_sample_ids.end(),
                  [&geno_rank](const std::string& a, const std::string& b) { return geno_rank[a] < geno_rank[b]; });
    }

    if ( pheno_matrix.samp_ids != overlapping_sample_ids ) {
        notice("Subsetting the phenotype matrix to the overlapping samples");
        pheno_matrix.subset_sample_ids(overlapping_sample_ids);
    }
    if ( have_cov && cov_matrix.samp_ids != overlapping_sample_ids ) {
        notice("Subsetting the covariate matrix to the overlapping samples");
        cov_matrix.subset_sample_ids(overlapping_sample_ids);
    }
    // categorical covariates become indicator columns, with levels counted on the analyzed samples
    // (as in rect-assoc and region-assoc)
    if ( have_cov && cov_matrix.has_categorical() ) {
        cov_matrix.expand_categorical(max_cat_levels);
        notice("Covariate matrix has %zu columns after expanding categorical covariates", cov_matrix.pheno_ids.size());
    }
    if ( mpr.get_all_sample_count() != (int32_t)overlapping_sample_ids.size() ) { // same size => same order, as overlapping_sample_ids follows genotype order
        notice("Subsetting the genotype data to %zu overlapping samples", overlapping_sample_ids.size());
        mpr.subset_sample_ids(overlapping_sample_ids);
    }
    if ( mpr.get_loaded_sample_count() != (int32_t)overlapping_sample_ids.size() ) {
        error("%d genotype samples are loaded, but %zu samples were matched across the input files", mpr.get_loaded_sample_count(), overlapping_sample_ids.size());
    }
    const int32_t n_overlapping_samples = (int32_t)overlapping_sample_ids.size();

    // covariates of the analyzed samples (no columns without --cov: intercept only)
    const Eigen::MatrixXd cov_all = have_cov ? cov_matrix.pheno_mat : Eigen::MatrixXd(n_overlapping_samples, 0);
    const std::shared_ptr<const CovariateResidualizer> resid_all = std::make_shared<const CovariateResidualizer>(cov_all);
    if ( have_cov ) {
        // warn when [1, C] is (nearly) rank-deficient, e.g. a covariate that is constant in the
        // analyzed samples or (nearly) a linear combination of others
        const CovariateResidualizer& resid = *resid_all;
        if ( resid.rank() < resid.n_cols() ) {
            warning("The intercept and %d covariates have numerical rank %d among the %d analyzed samples: some covariates are constant or (nearly) linear combinations of others. They are adjusted for once, and the degrees of freedom use the rank; drop the redundant covariates.",
                    (int32_t)cov_all.cols(), resid.rank(), n_overlapping_samples);
        }
        notice("Covariate degrees of freedom: %d (rank of the intercept and %d covariates = %d); Wald tests use df = N - 2 - %d",
               resid.n_cov_df(), (int32_t)cov_all.cols(), resid.rank(), resid.n_cov_df());
    }
    int64_t n_collinear_tests = 0; // tests whose genotype is explained by the covariates

    // read the variant-trait pair lists
    tsv_reader tr_pair(pairf.c_str());
    tr_pair.delimiter = table_delim;
    std::map<int32_t, std::set<cpra_t> > pair_map;
    std::map<int32_t, std::set<cpra_t> >::iterator pair_map_it;
    std::map<int32_t, std::set<cbe_t> > region_map;
    std::map<int32_t, std::set<cbe_t> >::iterator region_map_it;

    notice("Reading the trait-variant/region pairs from %s", pairf.c_str());
    int32_t n_pairs = 0;

    // open the output file gz or plain based on the extension
    htsFile* wf = hts_open(outf.c_str(), outf.substr(outf.length() - 3).compare(".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s for writing", outf.c_str());
    }
    // write the header line
    hprintf(wf, "TRAIT\tCHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tN_RR\tN_RA\tN_AA\tINFO\tBETA\tSE\tTSTAT\tLOG10P");
    if ( !skip_rint ) {
        hprintf(wf, "\tBETA_RINT\tSE_RINT\tTSTAT_RINT\tLOG10P_RINT");
    }
    hprintf(wf, "\n");

    int32_t icol_pair_trait = -1;
    int32_t icol_pair_variant = -1;
    int32_t icol_pair_region = -1;
    bool region_mode = false;
    std::map<std::string, int32_t>::iterator phe_trait2idx_it;
    while( tr_pair.read_line() ) {
        if ( icol_pair_trait < 0 ) {
            for(int32_t i = 0; i < tr_pair.nfields; ++i) {
                if ( colname_pair_trait.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_trait = i;
                }
                else if ( colname_pair_variant.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_variant = i;
                    region_mode = false;
                }
                else if ( colname_pair_region.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_region = i;
                    region_mode = true;
                }
            }
            if ( icol_pair_trait < 0 ) {
                error("Cannot find column %s in the pair file %s", colname_pair_trait.c_str(), pairf.c_str());
            }
            if ( icol_pair_variant < 0 && icol_pair_region < 0 ) {
                error("Either %s or %s must exist in the pair file %s", colname_pair_variant.c_str(), colname_pair_region.c_str(), pairf.c_str());
            }
            if ( icol_pair_variant >= 0 && icol_pair_region >= 0 ) {
                error("Both %s and %s cannot exist in the pair file %s. Please provide either variant ID or region string", colname_pair_variant.c_str(), colname_pair_region.c_str(), pairf.c_str());
            }
        }
        else {
            if ( icol_pair_variant >= 0 ) {
                const char* var_id = tr_pair.str_field_at(icol_pair_variant);
                const char* phe_id = tr_pair.str_field_at(icol_pair_trait);
                phe_trait2idx_it = pheno_matrix.pheno_id2idx.find(phe_id);
                if ( phe_trait2idx_it == pheno_matrix.pheno_id2idx.end() ) {
                    notice("Skipping phenotype %s, which is not observed in %s", phe_id, phef.c_str());
                }
                else {
                    int32_t phe_idx = phe_trait2idx_it->second;
                    cpra_t cpra(var_id);
                    if ( pair_map[phe_idx].insert(cpra).second ) {
                        ++n_pairs;
                    }
                    else {
                        notice("Skipping duplicate pair %s\t%s", var_id, phe_id);
                    }
                }
            }
            else {
                const char* phe_id = tr_pair.str_field_at(icol_pair_trait);
                const char* region_str = tr_pair.str_field_at(icol_pair_region);
                phe_trait2idx_it = pheno_matrix.pheno_id2idx.find(phe_id);
                if ( phe_trait2idx_it == pheno_matrix.pheno_id2idx.end() ) {
                    notice("Skipping phenotype %s, which is not observed in %s", phe_id, phef.c_str());
                }
                else {
                    int32_t phe_idx = phe_trait2idx_it->second;
                    cbe_t cbe(region_str);
                    if ( region_map[phe_idx].insert(cbe).second ) {
                        ++n_pairs;
                    }
                    else {
                        notice("Skipping duplicate pair %s\t%s", region_str, phe_id);
                    }
                }
            }
        }
    }
    notice("Finished reading %d trait-%s pairs across %zu traits", n_pairs, region_mode ? "region" : "variant", region_mode ? region_map.size() : pair_map.size());

    // LOCO: the .loco file of each tested trait
    if ( !locof.empty() && !predf.empty() ) {
        error("--loco and --pred cannot be used together: --loco takes a single .loco file, --pred a *_pred.list file");
    }
    std::map<int32_t, std::string> loco_path_of; // phe_idx -> .loco file
    if ( !locof.empty() || !predf.empty() ) {
        std::vector<int32_t> tested;
        if ( region_mode ) { for(region_map_it = region_map.begin(); region_map_it != region_map.end(); ++region_map_it) tested.push_back(region_map_it->first); }
        else               { for(pair_map_it = pair_map.begin(); pair_map_it != pair_map.end(); ++pair_map_it) tested.push_back(pair_map_it->first); }
        std::vector<std::string> tested_ids;
        for(int32_t phe_idx : tested) tested_ids.push_back(pheno_matrix.pheno_ids[phe_idx]);
        const std::vector<std::string> paths = predf.empty()
            ? resolve_loco_files(locof, tested_ids, ind_assoc_input::LOCO_SINGLE)
            : resolve_loco_files(predf, tested_ids, ind_assoc_input::LOCO_PRED_LIST);
        for(size_t k = 0; k < tested.size(); ++k) loco_path_of[tested[k]] = paths[k];
        notice("LOCO predictions are subtracted per trait and chromosome for %zu traits", tested.size());
    }
    auto loco_path = [&loco_path_of](int32_t phe_idx) {
        std::map<int32_t, std::string>::const_iterator it = loco_path_of.find(phe_idx);
        return it == loco_path_of.end() ? std::string() : it->second;
    };

    TraitData td;
    GenoChunk gc;
    if ( region_mode ) {
        for(region_map_it = region_map.begin(); region_map_it != region_map.end(); ++region_map_it) {
            const int32_t phe_idx = region_map_it->first;
            const std::set<cbe_t>& regions = region_map_it->second;
            if ( regions.empty() ) continue;
            const char* phe_id = pheno_matrix.pheno_ids[phe_idx].c_str();

            notice("Processing %zu regions for phenotype %s", regions.size(), phe_id);
            prepare_trait(td, pheno_matrix, phe_idx, cov_all, resid_all, skip_rint, loco_path(phe_idx));
            if ( !td.testable ) continue;

            // process one region at a time. The variant may overlap, but let's not worry about it for now
            for(std::set<cbe_t>::const_iterator regions_it = regions.begin(); regions_it != regions.end(); ++regions_it) {
                const cbe_t& region = *regions_it;
                notice("Processing region %s for phenotype %s", region.to_string().c_str(), phe_id);
                if ( !mpr.read_pos(region.chrom.c_str(), region.beg1) ) continue;
                apply_trait_loco(td, pheno_matrix, region.chrom, skip_rint);
                gc.reset(td.n());
                do {
                    const plink_var_t& var = mpr.get_current_variant();
                    if ( var.pos > region.end0 ) {
                        break;
                    }
                    cpra_t cpra(var.to_string().c_str());
                    mpr.get_genos();
                    gc.add_variant(mpr, cpra, td.rows);
                    if ( gc.size() >= max_chunk_vars ) { // test the current chunk
                        notice("Processing %d variants for phenotype %s", gc.size(), phe_id);
                        n_collinear_tests += test_and_write(wf, phe_id, td, gc, skip_rint, mean_impute_geno);
                        gc.reset(td.n());
                    }
                }
                while ( mpr.read_pivar() );

                if ( gc.size() > 0 ) {
                    notice("Processing %d variants for phenotype %s", gc.size(), phe_id);
                    n_collinear_tests += test_and_write(wf, phe_id, td, gc, skip_rint, mean_impute_geno);
                }
            }
        }
    }
    else {
        for(pair_map_it = pair_map.begin(); pair_map_it != pair_map.end(); ++pair_map_it) {
            const int32_t phe_idx = pair_map_it->first;
            const std::set<cpra_t>& pairs = pair_map_it->second;
            if ( pairs.empty() ) continue;
            const char* phe_id = pheno_matrix.pheno_ids[phe_idx].c_str();

            notice("Processing %zu pairs for phenotype %s", pairs.size(), phe_id);
            prepare_trait(td, pheno_matrix, phe_idx, cov_all, resid_all, skip_rint, loco_path(phe_idx));
            if ( !td.testable ) continue;

            gc.reset(td.n());
            for(std::set<cpra_t>::iterator it = pairs.begin(); it != pairs.end(); ++it) {
                std::string cpra_s(it->to_string());
                if ( !mpr.read_pivar(cpra_s.c_str()) ) {
                    notice("Skipping pair %s\t%s, which is not found in the pvar file", cpra_s.c_str(), phe_id);
                    continue;
                }
                // with LOCO, the phenotype depends on the chromosome: test the pairs of the
                // previous chromosome before switching (pairs are sorted by chromosome)
                if ( !td.loco_path.empty() && ( td.loco_chrom.empty() || loco_chrom_label(td.loco_chrom) != loco_chrom_label(it->chrom) ) ) {
                    n_collinear_tests += test_and_write(wf, phe_id, td, gc, skip_rint, mean_impute_geno);
                    gc.reset(td.n());
                    apply_trait_loco(td, pheno_matrix, it->chrom, skip_rint);
                }
                mpr.get_genos();
                gc.add_variant(mpr, *it, td.rows);
            }
            n_collinear_tests += test_and_write(wf, phe_id, td, gc, skip_rint, mean_impute_geno);
        }
    }
    hts_close(wf); // close the output file
    if ( n_collinear_tests > 0 ) {
        notice("%lld tests were reported as NaN because the covariates explain the genotype (residual variance below 1e-8 of its total variance)", (long long)n_collinear_tests);
    }
    notice("Analysis finished");
    return 0;
}
