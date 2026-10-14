#ifndef __ASSOC_UTILS_H
#define __ASSOC_UTILS_H

#include "Eigen/Dense"
#include <memory>
#include <set>
#include <unordered_map>
#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"
#include "qpgen_utils.h"
#include "pheno.h"


class genotype_chunk {
public:
    Eigen::MatrixXd geno_mat; // genotype matrix (n_samples x n_variants)
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> geno_mask; // genotype missingness mask (n_samples x n_variants)
    std::vector<cpra_t> v_cpra; // vector of CPRA identifiers for the variants in the chunk
    std::vector<var_cnt_t> var_cnts; // vector of variant counts for the variants in the chunk
    std::vector<double> infos;
    int32_t n_variants; // number of variants in the chunk
    int32_t n_skipped;
    int32_t n_samples; // number of samples in the chunk
    uint64_t n_geno_missing; // number of missing genotypes in the chunk

    void clear() {
        geno_mat.resize(0, 0);
        geno_mask.resize(0, 0);
        v_cpra.clear();
        infos.clear();
        var_cnts.clear();
        n_variants = 0;
        n_skipped = 0;
        n_samples = 0;
        n_geno_missing = 0;
    }

    genotype_chunk() : n_variants(0), n_skipped(0), n_samples(0), n_geno_missing(0) {}
    ~genotype_chunk() {}
};

// Least-squares residualization against [1, C] (intercept + covariates, n x (p+1)) using a
// column-pivoting Householder QR that is computed once and reused for any number of columns.
//
// The QR is computed on [1, C] with unit-norm columns, which spans the same space (so the
// residuals are unchanged) but makes the rank independent of the covariates' scales. rank() is
// the number of pivots above kRankTol times the largest one, as R's lm() does with tol = 1e-7.
// The residuals are r = (I - Q_r Q_r') v with the first rank() Householder vectors, which stays
// accurate when covariates are nearly collinear (a least-squares solve does not: its huge,
// cancelling coefficients turn rounding errors into large residual errors).
//
// n_cov_df() = rank() - 1 is the number of residual degrees of freedom the covariates use
// beyond the intercept. After residualizing both the phenotype and the genotype (Frisch-Waugh-
// Lovell), the Wald test of y ~ 1 + C + g has df = n - 2 - n_cov_df().
class CovariateResidualizer {
public:
    static constexpr double kRankTol = 1e-7;
    // a genotype column whose residual sum of squares is below kCollinearTol of its centered
    // sum of squares is (numerically) a linear combination of the covariates
    static constexpr double kCollinearTol = 1e-8;

    explicit CovariateResidualizer(const Eigen::MatrixXd& covariates);
    Eigen::Index n_samples() const { return n_samples_; }
    int32_t rank() const { return rank_; }
    int32_t n_cov_df() const { return rank_ - 1; }
    // number of columns of [1, C] (the intercept plus the covariates)
    int32_t n_cols() const { return n_cols_; }
    // values (n x k) <- residuals of values on [1, C]; blocks of a matrix are accepted
    void residualize(Eigen::Ref<Eigen::MatrixXd> values) const;
    // Residualizes genotype columns and sets to zero those the covariates explain (see
    // kCollinearTol), so that their statistics are reported as NaN instead of noise. Returns
    // the number of zeroed columns. Columns are processed in blocks to bound peak memory.
    int32_t residualize_genotypes(Eigen::Ref<Eigen::MatrixXd> G) const;
    // Residual sum of squares of each column of Y on [1, C]. Columns are processed one at a
    // time, so a column's value does not depend on the other columns passed with it (pair-assoc,
    // one trait at a time, and rect-/region-assoc, all traits at once, get identical values).
    Eigen::RowVectorXd residual_ss(const Eigen::MatrixXd& Y) const;

private:
    Eigen::Index n_samples_;
    int32_t n_cols_;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
    int32_t rank_;
};

// ---- LOCO (REGENIE step-1 predictions) ---------------------------------------
// A REGENIE .loco file has a header "FID_IID <FID>_<IID> ..." and one row per chromosome
// ("1".."22", "23" for X) holding the leave-that-chromosome-out prediction for every sample.

// LOCO predictions of one trait for one chromosome, indexed by sample ID.
struct loco_row_t {
    std::unordered_map<std::string, int32_t> exact;  // FID_IID -> column
    std::unordered_map<std::string, int32_t> suffix; // IID candidate -> column (-1 if ambiguous)
    std::vector<double> values;                      // NaN when missing
    std::vector<bool> observed;

    // Match a phenotype/genotype sample ID (usually IID only) to a .loco column:
    // exact FID_IID match, then FID==IID, then a unique "<FID>_<id>" suffix match.
    int32_t find(const std::string& id) const {
        std::unordered_map<std::string, int32_t>::const_iterator it = exact.find(id);
        if ( it != exact.end() ) return it->second;
        it = exact.find(id + "_" + id);
        if ( it != exact.end() ) return it->second;
        it = suffix.find(id);
        if ( it != suffix.end() ) return it->second; // may be -1 when ambiguous
        return -1;
    }
    // true if the file has a column for sample `id` (works for a header-only row)
    bool has(const std::string& id) const { return find(id) >= 0; }
    // value for sample `id`; false if absent or missing
    bool get(const std::string& id, double& val) const {
        int32_t j = find(id);
        if ( j < 0 || j >= (int32_t)observed.size() || !observed[j] ) return false;
        val = values[j];
        return true;
    }
};

// Canonical chromosome label shared by region strings and .loco rows: strip a "chr" prefix and
// map X/Y/XY/M(T) to REGENIE's numeric codes.
std::string loco_chrom_label(const std::string& chrom);
// The .loco file of each trait. locof is a single .loco file (only when one trait is tested) or a
// REGENIE *_pred.list file ("TRAIT PATH" per line; relative paths are also tried against the
// list's directory). kind: ind_assoc_input::LOCO_ANY, LOCO_SINGLE or LOCO_PRED_LIST.
std::vector<std::string> resolve_loco_files(const std::string& locof, const std::vector<std::string>& pheno_ids, int32_t kind);
// Reads the predictions of chromosome chrom from a .loco file into row (aborts if the file has no
// such row). With an empty chrom, only the sample IDs of the header are read.
void load_loco_row(const std::string& path, const std::string& chrom, loco_row_t& row);

// class containing individual-level association analysis input data
class ind_assoc_input {
protected:
    bool load_pgenlist(const char* pgenlistf);
    bool load_pgen_files(const char* pgenf, const char* pivarf, const char* psamf);
    bool load_pheno_cov_matrices(const char* phef, const char* pheno_format, const char* covf, const char* cov_format);
    bool is_pheno_loaded() const { return !pheno_matrix.samp_ids.empty() && !pheno_matrix.pheno_ids.empty(); }
    //bool add_genotypes();

public:
    MultiPgenIdxReader mpr; // multi-pgen reader
    PhenoMatrix pheno_matrix; // phenotype matrix
    PhenoMatrix cov_matrix;   // covariate matrix
    std::vector<std::string> overlapping_sample_ids; // overlapping sample IDs among genotype, phenotype, and covariate files
    char table_delim;      // delimiter of the user-supplied tables (see table_delimiter())
    std::vector<std::string> cat_cov_ids; // categorical covariates (--cat-cov-names)
    std::vector<std::string> cov_col_ids; // covariates to use (--cov-names); empty = all
    int32_t max_cat_levels;               // maximum number of levels of a categorical covariate
    int32_t jump_thres_bp; // jump threshold in base pairs for the multi-pgen reader
    int32_t icol_pivar_idx; // column index for the variant ID in the pivar file
    std::vector<std::string> subset_sample_ids;
    std::vector<std::string> subset_pheno_ids;
    bool rint_before_adj; // Perform rank-based inverse normal transformation before covariate adjustment
    bool rint_after_adj;  // Perform rank-based inverse normal transformation after covariate adjustment
    bool geno_chunk_done; // true once the current region has been fully streamed
    double min_af; // minimum allele frequency for filtering variants
    double max_af; // maximum allele frequency for filtering variants
    double min_ac; // minimum allele count for filtering variants
    double max_ac; // maximum allele count for filtering variants
    genotype_chunk geno_chunk; // genotype chunk for the current region

    // LOCO adjustment using REGENIE step-1 predictions (optional)
    std::string loco_file;  // REGENIE .loco file (single trait) or *_pred.list file (TRAIT PATH per line)
    int32_t loco_kind;      // expected kind of loco_file: LOCO_ANY, LOCO_SINGLE or LOCO_PRED_LIST
    // chromosome whose leave-one-chromosome-out predictions are subtracted at load time; empty
    // for per-chromosome mode, where apply_loco_chrom() is called before each chromosome's tests
    std::string loco_chrom;
    // Per-trait SD of the covariate-residualized phenotype used to standardize it before
    // subtracting the LOCO predictions (REGENIE's scale_Y). Empty unless LOCO is applied.
    Eigen::RowVectorXd pheno_scale;
    Eigen::MatrixXd pheno_loco_base;   // covariate-residualized phenotype divided by pheno_scale
    std::vector<std::string> loco_paths; // .loco file of each trait (pheno_matrix.pheno_ids order)
    std::vector<loco_row_t> loco_rows;   // predictions of each trait for loco_rows_chrom
    std::string loco_rows_chrom;         // chromosome held in loco_rows ("" = sample IDs only)
    std::string loco_applied_chrom;      // chromosome whose predictions pheno_matrix holds
    std::set<std::string> loco_chroms_seen; // chromosomes applied so far (to report revisits)
    bool subset_pheno_drop_absent;       // subset_pheno_ids: drop IDs absent from the phenotype file

    // Covariate degrees of freedom for the Wald tests: rank of [1, C] minus 1 (the intercept is
    // already part of df = N - 2). 0 without --cov. Set by load_pheno_cov_matrices().
    int32_t n_cov_df;
    // Residualizer on [1, C] (QR computed once per run); NULL without --cov.
    std::shared_ptr<const CovariateResidualizer> cov_resid;
    // true when the phenotype was transformed after the covariate adjustment (--rint-after-adj
    // or LOCO), so it is no longer orthogonal to the covariates (see pheno_cov_rss()).
    bool pheno_transformed_after_cov;

    // Residual sums of squares of the columns of Y (analyzed phenotypes, in any trait order or
    // subset) on [1, C]. The Wald tests use them in place of y'y so that they equal the OLS fit
    // of the analyzed phenotype on [1, C, g]. Returns false, leaving rss empty, when Y is already
    // residualized (no --cov, or no transformation after the adjustment) and y'y is exact.
    bool pheno_cov_rss(const Eigen::MatrixXd& Y, Eigen::RowVectorXd& rss) const {
        rss.resize(0);
        if ( !pheno_transformed_after_cov ) return false;
        if ( cov_resid ) {
            rss = cov_resid->residual_ss(Y);
        }
        else { // without --cov, the fit is on the intercept only
            rss = ( Y.rowwise() - Y.colwise().mean() ).colwise().squaredNorm();
        }
        return true;
    }

    // When true, a SINGLE tested trait may have missing values: samples whose
    // phenotype is missing are excluded from the analyzed sample set (before the
    // genotype, covariate and LOCO matching), so every downstream step runs on the
    // observed samples. With several traits, missing values are an error, because
    // all traits must share one sample set. Off by default (other commands keep
    // their own missing-value handling).
    bool drop_missing_pheno_samples;

    ind_assoc_input() {
        table_delim = '\t';
        max_cat_levels = 10;
        jump_thres_bp = 1000000;
        icol_pivar_idx = 8;
        rint_before_adj = false;
        rint_after_adj = false;
        drop_missing_pheno_samples = false;
        loco_kind = LOCO_ANY;
        subset_pheno_drop_absent = false;
        geno_chunk_done = false;
        n_cov_df = 0;
        pheno_transformed_after_cov = false;
    }
    ~ind_assoc_input() {}

    void set_jump_thres_bp(int32_t thres) { 
        jump_thres_bp = thres; 
        mpr.set_jump_thres_bp(thres);
    }
    void set_icol_pivar_idx(int32_t idx) {
        icol_pivar_idx = idx;
        mpr.set_icol_pivar_idx(idx);
    }
    void set_rint_before_adj(bool val) { rint_before_adj = val; }
    void set_rint_after_adj(bool val) { rint_after_adj = val; }
    void set_drop_missing_pheno_samples(bool val) { drop_missing_pheno_samples = val; }

    // Enable LOCO adjustment. Must be called before the phenotype matrix is loaded.
    // locof: a REGENIE .loco file (header starting with FID_IID) when a single trait is tested,
    //        or a REGENIE *_pred.list file mapping each trait to its .loco file.
    // chrom: chromosome of the tested region (e.g. "chr1", "1", "chrX"); X maps to 23. Empty
    //        (or NULL) for per-chromosome mode: the predictions are then subtracted by
    //        apply_loco_chrom(), which must be called before testing each chromosome.
    // kind : LOCO_ANY detects the file type; LOCO_SINGLE / LOCO_PRED_LIST require it.
    enum { LOCO_ANY = 0, LOCO_SINGLE = 1, LOCO_PRED_LIST = 2 };
    void set_loco(const char* locof, const char* chrom, int32_t kind = LOCO_ANY) {
        if ( is_pheno_loaded() ) {
            error("Cannot set LOCO file after the phenotype matrix is loaded");
        }
        loco_file = locof;
        loco_chrom = ( chrom == NULL ) ? "" : chrom;
        loco_kind = kind;
    }
    bool use_loco() const { return !loco_file.empty(); }
    // Sets pheno_matrix to the LOCO-adjusted phenotype for chromosome chrom (and RINT after it
    // with --rint-after-adj), from the covariate-adjusted phenotype kept at load time. Reads the
    // predictions of chrom when they are not already loaded. No-op when chrom is already applied.
    void apply_loco_chrom(const std::string& chrom);

    void set_minmax_af(double _min_af, double _max_af) {
        min_af = _min_af;
        max_af = _max_af;
    }

    void set_minmax_ac(double _min_ac, double _max_ac) {
        min_ac = _min_ac;
        max_ac = _max_ac;
    }

    // delimiter of the phenotype/covariate matrices and the sample/trait lists (table_delimiter());
    // set it before the set_subset_*_file() calls
    inline void set_table_delimiter(char delim) { table_delim = delim; }
    // categorical covariates, expanded to indicator columns (see PhenoMatrix::expand_categorical())
    // covariates to use (plus the categorical ones); empty = all columns of --cov
    inline void set_covariate_columns(const std::vector<std::string>& ids) { cov_col_ids = ids; }
    inline void set_categorical_covariates(const std::vector<std::string>& ids, int32_t max_levels) {
        cat_cov_ids = ids;
        max_cat_levels = max_levels;
    }

    inline void set_subset_sample_file(const char* samplef) {
        if ( is_pheno_loaded() ) {
            error("Cannot set subset sample file after the phenotype matrix is loaded");
        }
        tsv_reader tr_sample(samplef);
        tr_sample.delimiter = table_delim;
        while(tr_sample.read_line()) {
            subset_sample_ids.push_back(tr_sample.str_field_at(0));
        }
    }

    inline void set_subset_sample_ids(const std::vector<std::string>& samp_ids) {
        if ( is_pheno_loaded() ) {
            error("Cannot set subset sample file after the phenotype matrix is loaded");
        }
        subset_sample_ids = samp_ids;
    }

    // trait IDs, one per line (first column)
    inline void set_subset_pheno_file(const char* phenof) {
        if ( is_pheno_loaded() ) {
            error("Cannot set subset phenotype file after the phenotype matrix is loaded");
        }
        tsv_reader tr_pheno(phenof);
        tr_pheno.delimiter = table_delim;
        while(tr_pheno.read_line()) {
            subset_pheno_ids.push_back(tr_pheno.str_field_at(0));
        }
    }

    // drop_absent: IDs absent from the phenotype file are dropped (with a notice) instead of
    // being an error, and the kept traits are in sorted order
    inline void set_subset_pheno_ids(const std::vector<std::string>& pheno_ids, bool drop_absent = false) {
        if ( is_pheno_loaded() ) {
            error("Cannot set subset phenotype file after the phenotype matrix is loaded");
        }
        subset_pheno_ids = pheno_ids;
        subset_pheno_drop_absent = drop_absent;
    }

    bool process_pgenlist(const char* pgenlistf, const char* phef, const char* pheno_format, const char* covf, const char* cov_format);
    bool process_single_pgen(const char* pgenf, const char* pivarf, const char* psamf, const char* phef, const char* pheno_format, const char* covf, const char* cov_format);
    //int32_t load_genotype_chunk(const char* chrom, int32_t beg, int32_t end, int32_t max_chunk_vars, Eigen::MatrixXd& geno_mat, Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& geno_mask);
    int32_t load_genotype_chunk(const char* chrom, int32_t beg, int32_t end, int32_t max_chunk_vars);
};

Eigen::MatrixXd bulk_adjust_for_covariates(const Eigen::MatrixXd& values, const Eigen::MatrixXd& covariates);

// Residuals of each column of values (n x k) on [1, covariates]; see CovariateResidualizer.
Eigen::MatrixXd pheno_adj_cov_nxt_without_missing(const Eigen::MatrixXd& values, const Eigen::MatrixXd& covariates);

// Covariate adjustment when some phenotype values are missing. Each trait is regressed on the covariates
// using only its observed rows; residuals are returned for observed cells and missing cells are set to 0.
Eigen::MatrixXd pheno_adj_cov_nxt_with_missing(const Eigen::MatrixXd& values,
                                               const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& mask,
                                               const Eigen::MatrixXd& covariates);

void center_rows(Eigen::MatrixXd &matrix);

// Structure to hold the results of the regression for one variable
struct _slr_sumstat_t {
    double beta;
    double se;
    double tstat;
    double log10p;
    int32_t n_obs;
};

typedef struct _slr_sumstat_t slr_sumstat_t;

double tstat2log10pval(double tstat, double df);
// Wald t-tests of a centered y on each centered column of X (no intercept term).
// n_cov_df: covariate degrees of freedom when y and X were residualized against covariates
// (CovariateResidualizer::n_cov_df(); 0 without covariates); df = n - 2 - n_cov_df.
// y_rss: residual sum of squares of y under the null model, used in place of y'y when y is not
// itself residualized (e.g. RINT after the covariate adjustment); negative = use y'y.
bool simple_linear_regression_without_missing(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, std::vector<slr_sumstat_t>& results,
                                              int32_t n_cov_df = 0, double y_rss = -1.0);

// Same with pairwise deletion of missing values; df = n_complete - 2 - n_cov_df.
bool simple_linear_regression_with_missing( const Eigen::VectorXd& y,
                                            const Eigen::Vector<bool, Eigen::Dynamic>& y_mask,
                                            const Eigen::MatrixXd& X,
                                            const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& X_mask,
                                            std::vector<slr_sumstat_t>& results,
                                            int32_t n_cov_df = 0);

// Wald t-tests for every (column of X, column of Y) pair; n_cov_df as above. y_rss, if not
// NULL, replaces diag(Y'Y) (see simple_linear_regression_without_missing).
bool simple_rect_regression_without_missing(
    const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& X,
    std::vector<std::vector<slr_sumstat_t> >& results,
    int32_t n_cov_df = 0,
    const Eigen::RowVectorXd* y_rss = NULL);

// Same Wald t-test from precomputed statistics, so callers that already hold
// them (e.g. region-assoc, which shares X'Y with SuSiE and the export) avoid a
// second pass over X: n samples, x_sq_norms = diag(X'X) (p), y_sq_norms =
// diag(Y'Y) (K) or the null-model residual sums of squares, xt_y = X'Y (p x K),
// and df = n - 2 - n_cov_df.
bool simple_rect_regression_from_stats(
    int n,
    const Eigen::VectorXd& x_sq_norms,
    const Eigen::RowVectorXd& y_sq_norms,
    const Eigen::MatrixXd& xt_y,
    std::vector<std::vector<slr_sumstat_t> >& results,
    int32_t n_cov_df = 0);

// -log10 of the two-sided standard normal p-value for a z statistic
double zstat2log10pval(double zstat);

// REGENIE-style score test (null-model residual variance with n - n_cov degrees
// of freedom, normal p-value). Y and X must be centered and covariate-residualized;
// n_cov counts the intercept (1 + CovariateResidualizer::n_cov_df(), i.e. the rank of
// [1, C]). The z statistic is stored in slr_sumstat_t::tstat.
bool simple_rect_score_test_without_missing(
    const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& X,
    int32_t n_cov,
    std::vector<std::vector<slr_sumstat_t> >& results);

// Same score test from precomputed statistics (see simple_rect_regression_from_stats).
bool simple_rect_score_test_from_stats(
    int n,
    const Eigen::VectorXd& x_sq_norms,
    const Eigen::RowVectorXd& y_sq_norms,
    const Eigen::MatrixXd& xt_y,
    int32_t n_cov,
    std::vector<std::vector<slr_sumstat_t> >& results);

Eigen::VectorXd rint_with_missing(
    const Eigen::VectorXd& values,
    const Eigen::Vector<bool, Eigen::Dynamic>& mask);

Eigen::VectorXd rint_without_missing(const Eigen::VectorXd& values);

Eigen::MatrixXd rint_matrix_without_missing(const Eigen::MatrixXd& matrix);

// Column-wise RINT using only observed cells; missing cells are set to 0 in the output.
Eigen::MatrixXd rint_matrix_with_missing(const Eigen::MatrixXd& matrix,
                                         const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& mask);

void standardize_matrix_columns_inplace(Eigen::MatrixXd& matrix);

// Column-wise standardization using mean/sd of observed cells only; missing cells are set to exactly 0
// so that they contribute nothing to downstream inner products.
void standardize_matrix_columns_inplace(Eigen::MatrixXd& matrix,
                                        const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& mask);

Eigen::VectorXd columnwise_dot(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);

int32_t assoc_single_trait( 
    htsFile* wf, // output file handle
    const char* pheno_id, // phenotype ID
    const Eigen::VectorXd& phe_vec,      // phenotype vector
    const Eigen::VectorXd& phe_rint_vec, // rinted phenotype vector 
    const Eigen::Vector<bool, Eigen::Dynamic>& phe_mask_vec, // phenotype mask vector
    bool phe_has_missing, // if the phenotype has missing values
    bool skip_rint, // if the rank-based inverse normal transformation should be skipped
    const Eigen::MatrixXd& geno_mat,    // genotype matrix
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& geno_mask,   // genotype mask matrix
    bool geno_has_missing, // if the genotype has missing values
    const std::vector<cpra_t>& v_cpra,      // variant pairs
    const std::vector<int32_t>& ans,         // allele counts
    const std::vector<double>& acs,         // allele counts
    const std::vector<double>& infos       // infor values
);

int32_t assoc_single_trait( 
    htsFile* wf, // output file handle
    const char* pheno_id, // phenotype ID
    const Eigen::VectorXd& phe_vec,      // phenotype vector
    const Eigen::VectorXd& phe_rint_vec, // rinted phenotype vector 
    const Eigen::Vector<bool, Eigen::Dynamic>& phe_mask_vec, // phenotype mask vector
    bool phe_has_missing, // if the phenotype has missing values
    bool skip_rint, // if the rank-based inverse normal transformation should be skipped
    const Eigen::MatrixXd& geno_mat,    // genotype matrix
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& geno_mask,   // genotype mask matrix
    bool geno_has_missing, // if the genotype has missing values
    const std::vector<cpra_t>& v_cpra,      // variant pairs
    const std::vector<int32_t>& ans,         // allele counts
    const std::vector<double>& acs,         // allele counts
    const std::vector<int32_t>& gc0s,        // genotype counts for genotype 0
    const std::vector<int32_t>& gc1s,        // genotype counts for genotype 1
    const std::vector<int32_t>& gc2s,        // genotype counts for genotype 2
    const std::vector<double>& infos,       // infor values
    bool n_from_geno_counts = false,  // report N as the number of called genotypes (gc0+gc1+gc2) instead of the regression sample size
    int32_t n_cov_df = 0,             // covariate degrees of freedom (df = n - 2 - n_cov_df); 0 without covariates
    double rint_y_rss = -1.0          // null-model residual sum of squares of phe_rint_vec (negative = use its y'y)
);

#endif // __ASSOC_UTILS_H