#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "assoc_utils.h"
#include "qpgen.h"
#include "pheno.h"
#include "qpgen_utils.h"
#include "Eigen/Dense"
#include <cmath>

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <algorithm>

// For a diagonal metric diag(w), the exact PRS norm of candidate j over the traits observed for phenotyped sample i is
// sqrt( sum_{k in O_i} w_k x_jk^2 ) = sqrt( (X.^2 diag(w) mask^T)(j,i) ). This divides S(j,i) by that norm times
// pheno_norms(i), processing phenotyped samples in blocks to bound memory.
static void divide_by_exact_diagonal_norms(Eigen::MatrixXd& S, const Eigen::MatrixXd& X, const Eigen::VectorXd& w,
                                           const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& mask,
                                           const Eigen::VectorXd& pheno_norms) {
    Eigen::MatrixXd X2W = X.cwiseAbs2() * w.asDiagonal();          // n_x x p
    const int32_t block = 512;
    for ( int32_t c0 = 0; c0 < S.cols(); c0 += block ) {
        const int32_t b = std::min(block, (int32_t)S.cols() - c0);
        Eigen::MatrixXd maskd = mask.middleRows(c0, b).cast<double>();   // b x p
        Eigen::MatrixXd nrm = ( X2W * maskd.transpose() ).cwiseSqrt();  // n_x x b
        for ( int32_t c = 0; c < b; ++c ) {
            S.col(c0 + c) = S.col(c0 + c).array() / ( nrm.col(c).array() * pheno_norms(c0 + c) + 1e-100 );
        }
    }
}

// Tunes the shrinkage parameter lambda of the Mahalanobis metric using the phenotyped samples that have a
// mapped PRS sample. For a candidate lambda, the mapped samples are scored against all PRS samples exactly as
// in the main analysis, and a separation criterion of the self match is computed:
//   mean-log-softmax : mean over mapped samples of log( exp(Z_self) / sum_j exp(Z_j) )  (default)
//   mean-z           : mean over mapped samples of Z_self
//   mrr              : mean over mapped samples of 1 / rank of the self match
// The criterion is maximized over [0,1] by golden-section search.
struct LambdaTuner {
    const Eigen::MatrixXd& total_cov;
    const Eigen::MatrixXd& X;          // standardized PRS matrix (n_prs x p)
    Eigen::MatrixXd Ymap;              // standardized phenotypes of mapped samples (n_map x p)
    Eigen::MatrixXd maskmap;           // 0/1 mask of Ymap (n_map x p), used only when has_missing
    bool has_missing;
    Eigen::VectorXd wsqrt;             // sqrt of trait weights
    std::vector<int32_t> self_rows;    // PRS row index of the self match for each column of Ymap
    std::string metric;
    int32_t n_eval;
    int32_t last_n_self_best;          // number of self matches ranked first in the last evaluation

    LambdaTuner(const Eigen::MatrixXd& total_cov_, const Eigen::MatrixXd& X_,
                const Eigen::MatrixXd& Y, const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& mask,
                bool has_missing_, const Eigen::VectorXd& weights,
                const std::vector<int32_t>& prs_idx, const std::vector<int32_t>& pheno_idx,
                const std::vector<int32_t>& n_traits_obs, const std::string& metric_, int32_t min_self)
        : total_cov(total_cov_), X(X_), has_missing(has_missing_), wsqrt(weights.cwiseSqrt()), metric(metric_), n_eval(0), last_n_self_best(0)
    {
        std::vector<int32_t> keep;
        for ( size_t k = 0; k < pheno_idx.size(); ++k ) {
            if ( n_traits_obs[pheno_idx[k]] > 0 ) keep.push_back((int32_t)k);
        }
        if ( (int32_t)keep.size() < std::max(min_self, 1) ) {
            error("--auto-lambda: only %d phenotyped samples have a mapped PRS sample and observed traits, fewer than --auto-lambda-min-self %d. Check the sample IDs or --sample-tsv mapping, set --lambda manually, or lower --auto-lambda-min-self",
                  (int32_t)keep.size(), min_self);
        }
        notice("--auto-lambda will tune on %d mapped samples", (int32_t)keep.size());
        Ymap.resize(keep.size(), Y.cols());
        maskmap.resize(has_missing ? (int32_t)keep.size() : 0, has_missing ? Y.cols() : 0);
        for ( size_t r = 0; r < keep.size(); ++r ) {
            Ymap.row(r) = Y.row(pheno_idx[keep[r]]);
            if ( has_missing ) maskmap.row(r) = mask.row(pheno_idx[keep[r]]).cast<double>();
            self_rows.push_back(prs_idx[keep[r]]);
        }
    }

    double evaluate(double lambda) {
        const int32_t p = total_cov.rows();
        Eigen::MatrixXd tot = (1.0 - lambda) * total_cov + lambda * Eigen::MatrixXd::Identity(p, p);
        Eigen::MatrixXd M = wsqrt.asDiagonal() * tot.inverse() * wsqrt.asDiagonal();
        Eigen::MatrixXd YM = Ymap * M;
        if ( has_missing ) YM = YM.cwiseProduct(maskmap);
        Eigen::MatrixXd S = X * YM.transpose();                                       // n_prs x n_map
        Eigen::VectorXd prs_norms = ( X * M ).cwiseProduct( X ).rowwise().sum().cwiseSqrt();
        S.array().colwise() /= ( prs_norms.array() + 1e-100 );                         // pheno norm cancels below
        const int32_t n = S.rows();
        double total = 0.0;
        int32_t n_self_best = 0;
        for ( int32_t c = 0; c < S.cols(); ++c ) {
            Eigen::VectorXd col = S.col(c);
            double mean = col.mean();
            double sd = std::sqrt( (col.array() - mean).square().sum() / (double)(n - 1) );
            if ( sd <= 0 ) continue;
            Eigen::ArrayXd z = (col.array() - mean) / sd;
            double z_self = z(self_rows[c]);
            if ( (z > z_self).count() == 0 ) ++n_self_best;
            if ( metric == "mean-z" ) {
                total += z_self;
            }
            else if ( metric == "mrr" ) {
                int32_t rank = 1 + (int32_t)(z > z_self).count();
                total += 1.0 / (double)rank;
            }
            else { // mean-log-softmax
                double zmax = z.maxCoeff();
                double lse = zmax + std::log( (z - zmax).exp().sum() );
                total += z_self - lse;
            }
        }
        ++n_eval;
        last_n_self_best = n_self_best;
        double value = total / (double)S.cols();
        notice("  auto-lambda evaluation %d: lambda = %.4f, %s = %.6f, self-best = %d / %d", n_eval, lambda, metric.c_str(), value, n_self_best, (int32_t)S.cols());
        return value;
    }

    // golden-section maximization on [a, b]
    double findOptimalLambda(double a = 0.0, double b = 1.0, double tol = 0.01) {
        const double phi = (std::sqrt(5.0) - 1.0) / 2.0;
        double x1 = b - phi * (b - a), x2 = a + phi * (b - a);
        double f1 = evaluate(x1), f2 = evaluate(x2);
        while ( (b - a) > tol ) {
            if ( f1 > f2 ) { b = x2; x2 = x1; f2 = f1; x1 = b - phi * (b - a); f1 = evaluate(x1); }
            else           { a = x1; x1 = x2; f1 = f2; x2 = a + phi * (b - a); f2 = evaluate(x2); }
        }
        return 0.5 * (a + b);
    }
};

int32_t cmd_match_prs_pheno(int32_t argc, char **argv)
{
    std::string prsf;       // PRS files
    std::string phef;
    std::string covf;
    std::string sample_mapf;    // TSV file containing input sample IDs PRS and phenotype files in [PRS_SAMPLE_ID] [PHENO_SAMPLE_ID] format. If they use the same IDs, use only a single column if subsetting samples are needed
    std::string trait_mapf;     // TSV file containing input trait IDs PRS and phenotype files in [PRS_TRAIT_ID] [PHENO_TRAIT_ID] format. If they use the same IDs, use only a single column if subsetting traits are needed
    std::string outf;
    std::string weightf;    // Weights file to use for each phenotype
    std::string prs_format("regenie");
    std::string pheno_format("regenie");
    std::string cov_format("regenie");
    std::string cov_names;      // covariates to use (--cov-names); empty = all
    std::string cat_cov_names;  // categorical covariates (--cat-cov-names)
    int32_t max_cat_levels = 10; // REGENIE's --maxCatLevels
    std::string missing_str("NA");  // comma-separated strings representing missing values
    bool cov_impute_mean = false;  // mean-impute missing covariates instead of dropping samples
    bool missing_as_mean = false;  // impute missing phenotypes with the trait mean instead of ignoring them
    bool missing_as_half_min = false; // impute missing phenotypes with half of the trait minimum (below detection limit)
    bool rint_after_adj = false;   // Perform rank-based inverse normal transformation after covariate adjustment
    bool use_mahalanobis = false;  // Use Mahalanobis distance for matching
    bool auto_lambda = false;     // tune the Mahalanobis shrinkage parameter on the mapped samples
    std::string auto_lambda_metric("mean-log-softmax"); // criterion maximized by --auto-lambda
    int32_t auto_lambda_min_self = 100;  // minimum number of mapped samples usable for tuning
    int32_t auto_lambda_min_best = 100;  // minimum number of self-best matches at the tuned lambda
    bool exact_norm = false;      // compute PRS norms over each phenotyped sample's observed traits instead of all traits
    bool no_norm = false;         // legacy independence score: inner product divided by total absolute weight, no profile norms
    double min_weight = 0.0;  // minimum weight (in r) per trait to set to zero
    double z_lenient_threshold = 1.96;  // Z-score threshold for lenient matching
    double z_diff_threshold = 2.0;      // Z-score difference to declare a clear match
    int32_t n_threads = 1;
    double lambda = 0.0; // regularization parameter for Mahalanobis distance
    double weight_prs_mh = 0.5; // weight for MH distance when combining with weighted correlation

    bool space_delim = false; // space- instead of tab-delimited input tables

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("prs", &prsf, "Input PRS file")
    LONG_STRING_PARAM("pheno", &phef, "Input phenotype matrix")
    LONG_STRING_PARAM("cov", &covf, "Input covariate matrix (optional)")
    LONG_STRING_PARAM("sample-tsv", &sample_mapf, "TSV file containing input sample IDs PRS and phenotype files in [PRS_SAMPLE_ID] [PHENO_SAMPLE_ID] format. If they use the same IDs, use only a single column if subsetting samples are needed")
    LONG_STRING_PARAM("trait-tsv", &trait_mapf, "TSV file containing input trait IDs PRS and phenotype files in [PRS_TRAIT_ID] [PHENO_TRAIT_ID] format. If they use the same IDs, use only a single column if subsetting traits are needed")
    LONG_STRING_PARAM("weights", &weightf, "Input weights file for each phenotype in [PHENO_ID] [WEIGHT] format")
    LONG_STRING_PARAM("prs-format", &prs_format, "Format of the PRS file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("pheno-format", &pheno_format, "Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("cov-format", &cov_format, "Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("cov-names", &cov_names, QPGEN_COV_NAMES_HELP)
    LONG_STRING_PARAM("cat-cov-names", &cat_cov_names, QPGEN_CAT_COV_NAMES_HELP)
    LONG_INT_PARAM("max-cat-levels", &max_cat_levels, QPGEN_MAX_CAT_LEVELS_HELP)
    LONG_STRING_PARAM("missing-str", &missing_str, "Comma-separated strings representing missing values in the phenotype and covariate matrices (default: 'NA'). The PRS matrix must not contain missing values")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Analysis options", NULL)
    LONG_PARAM("rint", &rint_after_adj, "Perform rank-based inverse normal transformation after covariate adjustment (default: false)")
    LONG_PARAM("mahalanobis", &use_mahalanobis, "Use Mahalanobis distance for matching (default: false)")
    LONG_PARAM("no-norm", &no_norm, "Without --mahalanobis, divide the weighted inner product by the total absolute weight of the observed traits instead of by the weighted norms of the two profiles (legacy score; default: false)")
    LONG_PARAM("exact-norm", &exact_norm, "When phenotypes have missing values, compute each PRS norm over the traits observed for the phenotyped sample instead of over all traits. One extra matrix product without --mahalanobis or with --lambda 1; one pass per missingness pattern otherwise (default: false)")
    LONG_DOUBLE_PARAM("weight-prs-mh", &weight_prs_mh, "Weight for PRS distance when combining with weighted correlation (default: 0.5)")
    LONG_DOUBLE_PARAM("lambda", &lambda, "Regularization parameter for Mahalanobis distance, between 0 and 1 (default: 0.0)")
    LONG_DOUBLE_PARAM("min-weight", &min_weight, "Minimum weight (in r) per trait to set to zero (default: 0.0)")
    LONG_DOUBLE_PARAM("z-threshold", &z_lenient_threshold, "Z-score threshold for lenient matching (default: 1.96)")
    LONG_DOUBLE_PARAM("z-diff", &z_diff_threshold, "Z-score difference to declare a clear match (default: 2.0)")
    LONG_INT_PARAM("threads", &n_threads, "Number of threads to use (default: 1)")

    LONG_PARAM_GROUP("Imputation options", NULL)
    LONG_PARAM("cov-impute-mean", &cov_impute_mean, "Mean-impute missing covariate values instead of dropping samples with any missing covariate (default: false)")
    LONG_PARAM("missing-as-mean", &missing_as_mean, "Impute missing phenotype values with the mean of observed values for the trait, then treat them as observed (default: false, missing values are ignored)")
    LONG_PARAM("missing-as-half-min", &missing_as_half_min, "Impute missing phenotype values with half of the minimum observed value for the trait, the usual convention for measurements below a detection limit, then treat them as observed (default: false, missing values are ignored)")

    LONG_PARAM_GROUP("Auto-lambda options (with --mahalanobis)", NULL)
    LONG_PARAM("auto-lambda", &auto_lambda, "Choose the shrinkage parameter lambda automatically by maximizing the separation of the mapped self matches (default: false)")
    LONG_STRING_PARAM("auto-lambda-metric", &auto_lambda_metric, "Criterion maximized by --auto-lambda (default: 'mean-log-softmax'). Options: 'mean-log-softmax', 'mean-z', 'mrr'")
    LONG_INT_PARAM("auto-lambda-min-self", &auto_lambda_min_self, "Minimum number of phenotyped samples with a mapped PRS sample and observed traits required for tuning; fewer is an error (default: 100)")
    LONG_INT_PARAM("auto-lambda-min-best", &auto_lambda_min_best, "Minimum number of mapped samples whose own PRS ranks first at the tuned lambda; fewer is an error, indicating a mismatched sample mapping (default: 100; 0 disables)")

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

    // validate option combinations before reading any input
    if ( prsf.empty() || phef.empty() || outf.empty() ) {
        error("--prs, --pheno and --out are required");
    }
    if ( missing_as_mean && missing_as_half_min ) {
        error("--missing-as-mean and --missing-as-half-min cannot be used together");
    }
    if ( no_norm && use_mahalanobis ) {
        error("--no-norm applies only to the independence score and cannot be combined with --mahalanobis");
    }
    if ( lambda < 0.0 || lambda > 1.0 ) {
        error("Invalid value for --lambda: %.4f. Must be between 0 and 1", lambda);
    }
    if ( auto_lambda ) {
        if ( !use_mahalanobis ) {
            error("--auto-lambda requires --mahalanobis; the independence score has no shrinkage parameter");
        }
        if ( lambda > 0 ) {
            error("Cannot use --auto-lambda together with a non-zero --lambda. Set --lambda 0 (default) to tune it automatically");
        }
        if ( auto_lambda_metric != "mean-log-softmax" && auto_lambda_metric != "mean-z" && auto_lambda_metric != "mrr" ) {
            error("Unknown --auto-lambda-metric '%s'. Options: mean-log-softmax, mean-z, mrr", auto_lambda_metric.c_str());
        }
    }
    if ( weight_prs_mh < 0.0 || weight_prs_mh > 1.0 ) {
        error("Invalid value for --weight-prs-mh: %.4f. Must be between 0 and 1", weight_prs_mh);
    }
    if ( n_threads < 1 ) {
        error("--threads must be at least 1");
    }

    if ( n_threads > 1 ) {
        Eigen::setNbThreads(n_threads);
        notice("Setting number of threads to %d", n_threads);
    }

    // Load the weight file
    std::map<std::string, double> phe_weights;
    if ( !weightf.empty() ) {
        notice("Loading phenotype weights from %s", weightf.c_str());
        tsv_reader tr_weight(weightf.c_str());
        tr_weight.delimiter = table_delim;
        while( tr_weight.read_line() ) {
            if ( tr_weight.nfields < 2 ) {
                error("Invalid format weights file %s in line %zu starting with %s. Must containing at least 2 fields", weightf.c_str(), (int32_t)weightf.size() + 1, tr_weight.str_field_at(0) );
            }
            std::string phe_id = tr_weight.str_field_at(0);
            double weight = tr_weight.double_field_at(1);
            phe_weights[phe_id] = weight;
        }
        if ( phe_weights.size() == 0 ) {
            error("No valid phenotype weights found in file %s", weightf.c_str());
        }
    }

    notice("Loading phenotype matrix from %s", phef.c_str());

    // load the phenotype matrix
    PhenoMatrix pheno_matrix;
    pheno_matrix.add_missing_strs(missing_str);
    if ( !pheno_matrix.load_pheno_matrix(phef.c_str(), pheno_format.c_str(), table_delim) ) {
        error("Failed to load the phenotype matrix from file %s", phef.c_str());
    }

    notice("Loaded phenotype matrix with %d samples and %d phenotypes from %s", (int32_t)pheno_matrix.samp_ids.size(), (int32_t)pheno_matrix.pheno_ids.size(), phef.c_str());
    if ( pheno_matrix.has_missing ) {
        int64_t n_missing = (int64_t)pheno_matrix.pheno_mask.size() - (int64_t)pheno_matrix.pheno_mask.cast<int64_t>().sum();
        notice("Phenotype matrix contains %lld missing values (%.3f%%); traits with missing values will be ignored for the corresponding individuals", (long long)n_missing, 100.0 * (double)n_missing / (double)pheno_matrix.pheno_mask.size());
    }

    notice("Loading PRS matrix from %s", prsf.c_str());
    // load the PRS matrix
    PhenoMatrix prs_matrix;
    prs_matrix.add_missing_strs(missing_str);
    if ( !prs_matrix.load_pheno_matrix(prsf.c_str(), prs_format.c_str(), table_delim) ) {
        error("Failed to load the PRS matrix from file %s", prsf.c_str());
    }
    if ( prs_matrix.has_missing ) {
        error("PRS matrix %s contains missing values, which is not supported. Please remove or impute missing PRS values before matching", prsf.c_str());
    }

    notice("Loaded PRS matrix with %d samples and %d traits from %s", (int32_t)prs_matrix.samp_ids.size(), (int32_t)prs_matrix.pheno_ids.size(), prsf.c_str());

    // load the covariance matrix
    PhenoMatrix cov_matrix;
    if ( !covf.empty() ) {
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

    // subset to overlapping phenotypes
    std::vector<int32_t> phe_pheno_indices;
    std::vector<int32_t> prs_pheno_indices;
    if ( trait_mapf.empty() ) { // identify overlapping phenotypes based on names
        std::vector<std::string> overlapping_phenos;
        identify_overlapping_ids(prs_matrix.pheno_ids, pheno_matrix.pheno_ids, overlapping_phenos);
        notice("Found %d overlapping phenotypes between PRS and phenotype matrices based on phenotype IDs", (int32_t)overlapping_phenos.size());
        if ( overlapping_phenos.size() == 0 ) {
            error("No overlapping phenotypes found between PRS and phenotype matrices. Please check if the phenotype IDs are consistent across the files or provide a trait ID mapping file using --trait-tsv option.");
        }
        pheno_matrix.subset_pheno_ids(overlapping_phenos);
        prs_matrix.subset_pheno_ids(overlapping_phenos);
    }
    else {
        notice("Loading trait ID mapping between PRS and phenotype files from %s", trait_mapf.c_str());
        tsv_reader tr_trait_map(trait_mapf.c_str());
        tr_trait_map.delimiter = table_delim;
        while( tr_trait_map.read_line() ) {
            if ( tr_trait_map.nfields > 2 ) {
                error("Invalid format trait ID mapping file %s in line %zu starting with %s. Must containing at least 2 fields", trait_mapf.c_str(), (int32_t)trait_mapf.size() + 1, tr_trait_map.str_field_at(0) );
            }
            std::string prs_id = tr_trait_map.str_field_at(0);
            std::string phe_id = tr_trait_map.str_field_at(tr_trait_map.nfields == 1 ? 0 : 1);
            std::map<std::string, int32_t>::const_iterator it_prs = prs_matrix.pheno_id2idx.find(prs_id);
            std::map<std::string, int32_t>::const_iterator it_phe = pheno_matrix.pheno_id2idx.find(phe_id);
            if ( it_prs == prs_matrix.pheno_id2idx.end() ) {
                error("PRS trait ID %s in mapping file %s not found in PRS matrix", prs_id.c_str(), trait_mapf.c_str());
            }
            if ( it_phe == pheno_matrix.pheno_id2idx.end() ) {
                error("Phenotype trait ID %s in mapping file %s not found in phenotype matrix", phe_id.c_str(), trait_mapf.c_str());
            }
            prs_pheno_indices.push_back( it_prs->second );
            phe_pheno_indices.push_back( it_phe->second );
        }
        pheno_matrix.subset_pheno_indices( phe_pheno_indices );
        prs_matrix.subset_pheno_indices( prs_pheno_indices );
    }

    // after the subsetting, the number of phenotypes in both matrices should be the same
    if ( pheno_matrix.pheno_ids.size() != prs_matrix.pheno_ids.size() ) {
        error("Number of phenotypes after subsetting do not match between phenotype and PRS matrices (%d vs %d)", (int32_t)pheno_matrix.pheno_ids.size(), (int32_t)prs_matrix.pheno_ids.size());
    }
    notice("Subsetted to %d overlapping phenotypes between PRS and phenotype matrices", (int32_t)pheno_matrix.pheno_ids.size());

    // optionally impute missing phenotypes so that they are treated as observed downstream
    if ( ( missing_as_mean || missing_as_half_min ) && pheno_matrix.has_missing ) {
        int64_t n_imputed = pheno_matrix.impute_missing(missing_as_half_min);
        notice("Imputed %lld missing phenotype values with %s of the trait (%s)", (long long)n_imputed,
               missing_as_half_min ? "half of the minimum" : "the mean", missing_as_half_min ? "--missing-as-half-min" : "--missing-as-mean");
    }

    // samples dropped from the phenotype matrix because of missing covariates
    std::set<std::string> dropped_samp_ids;

    // if covariates are provided, subset to overlapping samples
    if ( !covf.empty() ) {
        std::vector<std::string> overlapping_samples;
        identify_overlapping_ids(pheno_matrix.samp_ids, cov_matrix.samp_ids, overlapping_samples);
        notice("Found %d overlapping samples between phenotype and covariate matrices", (int32_t)overlapping_samples.size());
        if ( overlapping_samples.size() == 0 ) {
            error("No overlapping samples found between phenotype and covariate matrices. Please check if the sample IDs are consistent across the files.");
        }
        pheno_matrix.subset_sample_ids(overlapping_samples);
        cov_matrix.subset_sample_ids(overlapping_samples);
        notice("Subsetted phenotype and covariate matrices to %d overlapping samples", (int32_t)overlapping_samples.size());

        // Without --cov-impute-mean, samples with a missing covariate are dropped first, so that
        // categorical levels are counted on the retained samples. With it, the categorical
        // covariates are expanded first, so that their indicator columns are imputed (with the
        // level frequencies) rather than a missing value becoming a level.
        if ( cov_matrix.has_missing && !cov_impute_mean ) {
            // drop samples with any missing covariate from both matrices
            std::vector<std::string> keep_ids;
            for ( int32_t i = 0; i < cov_matrix.pheno_mat.rows(); ++i ) {
                if ( cov_matrix.pheno_mask.row(i).all() ) keep_ids.push_back( cov_matrix.samp_ids[i] );
                else dropped_samp_ids.insert( cov_matrix.samp_ids[i] );
            }
            if ( keep_ids.empty() ) {
                error("All %d samples have at least one missing covariate value. Consider using --cov-impute-mean", (int32_t)cov_matrix.pheno_mat.rows());
            }
            pheno_matrix.subset_sample_ids(keep_ids);
            cov_matrix.subset_sample_ids(keep_ids);
            notice("Dropped %d samples with at least one missing covariate value; %d samples retained. These samples will not appear in the output. Use --cov-impute-mean to mean-impute missing covariates instead",
                   (int32_t)dropped_samp_ids.size(), (int32_t)keep_ids.size());
        }
        // categorical covariates become indicator columns, with levels counted on the analyzed samples
        if ( cov_matrix.has_categorical() ) {
            cov_matrix.expand_categorical(max_cat_levels);
            notice("Covariate matrix has %zu columns after expanding categorical covariates", cov_matrix.pheno_ids.size());
        }
        if ( cov_matrix.has_missing && cov_impute_mean ) {
            // replace each missing covariate value with the mean of observed values in that column
            int32_t n_imputed = 0;
            for ( int32_t j = 0; j < cov_matrix.pheno_mat.cols(); ++j ) {
                double sum = 0.0;
                int32_t n_obs = 0;
                for ( int32_t i = 0; i < cov_matrix.pheno_mat.rows(); ++i ) {
                    if ( cov_matrix.pheno_mask(i, j) ) { sum += cov_matrix.pheno_mat(i, j); ++n_obs; }
                }
                double mean = ( n_obs > 0 ) ? sum / (double)n_obs : 0.0;
                for ( int32_t i = 0; i < cov_matrix.pheno_mat.rows(); ++i ) {
                    if ( !cov_matrix.pheno_mask(i, j) ) {
                        cov_matrix.pheno_mat(i, j) = mean;
                        cov_matrix.pheno_mask(i, j) = true;
                        ++n_imputed;
                    }
                }
            }
            cov_matrix.recompute_has_missing();
            notice("Mean-imputed %d missing covariate values (--cov-impute-mean)", n_imputed);
        }

        notice("Adjusting phenotypes by covariates using linear regression");
        if ( pheno_matrix.has_missing ) {
            pheno_matrix.pheno_mat = pheno_adj_cov_nxt_with_missing(pheno_matrix.pheno_mat, pheno_matrix.pheno_mask, cov_matrix.pheno_mat);
        }
        else {
            pheno_matrix.pheno_mat = pheno_adj_cov_nxt_without_missing(pheno_matrix.pheno_mat, cov_matrix.pheno_mat);
        }
    }

    if ( rint_after_adj ) {
        notice("Performing rank-based inverse normal transformation for all phenotypes after covariate adjustment");
        if ( pheno_matrix.has_missing ) {
            pheno_matrix.pheno_mat = rint_matrix_with_missing(pheno_matrix.pheno_mat, pheno_matrix.pheno_mask);
        }
        else {
            pheno_matrix.pheno_mat = rint_matrix_without_missing(pheno_matrix.pheno_mat);
        }
    }

    // construct matching ID maps
    std::vector<int32_t> matching_prs_samp_indices;
    std::vector<int32_t> matching_pheno_samp_indices;
    std::map<int32_t, int32_t> samp_idx_pheno2prs;
    bool has_sample_map = !sample_mapf.empty();
    bool has_trait_map = !trait_mapf.empty();

    // identify PRS samples to use
    if ( has_sample_map ) {
        notice("Loading sample ID mapping between PRS and phenotype files from %s", sample_mapf.c_str());
        tsv_reader tr_sample_map(sample_mapf.c_str());
        tr_sample_map.delimiter = table_delim;
        int32_t n_skipped_dropped = 0;
        while( tr_sample_map.read_line() ) {
            if ( tr_sample_map.nfields > 2 ) {
                error("Invalid format sample ID mapping file %s in line %zu starting with %s. Must containing at least 2 fields", sample_mapf.c_str(), (int32_t)sample_mapf.size() + 1, tr_sample_map.str_field_at(0) );
            }
            std::string prs_id = tr_sample_map.str_field_at(0);
            std::string phe_id = tr_sample_map.str_field_at(tr_sample_map.nfields == 1 ? 0 : 1);
            std::map<std::string, int32_t>::const_iterator it_prs = prs_matrix.samp_id2idx.find(prs_id);
            std::map<std::string, int32_t>::const_iterator it_phe = pheno_matrix.samp_id2idx.find(phe_id);
            if ( it_prs == prs_matrix.samp_id2idx.end() ) {
                error("PRS trait ID %s in mapping file %s not found in PRS matrix", prs_id.c_str(), sample_mapf.c_str());
            }
            if ( it_phe == pheno_matrix.samp_id2idx.end() ) {
                if ( dropped_samp_ids.find(phe_id) != dropped_samp_ids.end() ) {
                    ++n_skipped_dropped; // sample was dropped because of missing covariates
                    continue;
                }
                error("Phenotype trait ID %s in mapping file %s not found in phenotype matrix", phe_id.c_str(), sample_mapf.c_str());
            }
            matching_prs_samp_indices.push_back( it_prs->second );
            matching_pheno_samp_indices.push_back( it_phe->second );
            samp_idx_pheno2prs[it_phe->second] = it_prs->second;
        }
        if ( n_skipped_dropped > 0 ) {
            notice("Skipped %d entries in the sample ID mapping file whose phenotype samples were dropped due to missing covariates", n_skipped_dropped);
        }
        notice("Found %d overlapping samples between PRS and phenotype matrices based on sample ID mapping file", (int32_t)matching_pheno_samp_indices.size());
    }
    else { // identify overlapping samples based on names
        for(int32_t i=0; i < pheno_matrix.samp_ids.size(); ++i ) {
            const std::string& samp_id = pheno_matrix.samp_ids[i];
            std::map<std::string, int32_t>::const_iterator it = prs_matrix.samp_id2idx.find( samp_id );
            if ( it != prs_matrix.samp_id2idx.end() ) {
                matching_pheno_samp_indices.push_back( i );
                matching_prs_samp_indices.push_back( it->second );
                samp_idx_pheno2prs[i] = it->second;
            }
        }
        notice("Found %d overlapping samples between PRS and phenotype matrices based on sample IDs", (int32_t)matching_pheno_samp_indices.size());
    }

    const int32_t n_traits = (int32_t)pheno_matrix.pheno_ids.size();
    const int32_t n_pheno_samples = (int32_t)pheno_matrix.samp_ids.size();
    const int32_t n_prs_samples = (int32_t)prs_matrix.samp_ids.size();
    const bool pheno_has_missing = pheno_matrix.has_missing;

    // standardize each trait. For phenotypes, mean/sd are computed from observed values only and
    // missing cells are set to exactly 0 so that they contribute nothing to any inner product below.
    notice("Standardizing PRS and phenotype matrices");
    standardize_matrix_columns_inplace(prs_matrix.pheno_mat);
    if ( pheno_has_missing ) {
        standardize_matrix_columns_inplace(pheno_matrix.pheno_mat, pheno_matrix.pheno_mask);
    }
    else {
        standardize_matrix_columns_inplace(pheno_matrix.pheno_mat);
    }

    // for PRS X, and phenotype Y, compute column-wise correlation between X and Y
    Eigen::VectorXd weights( n_traits );
    if ( weightf.empty() ) {
        notice("Computing weights for each phenotype");
        std::vector<int32_t> n_obs_per_trait( n_traits, 0 );
        for(int32_t i=0; i < n_traits; ++i) {
            double r = 0;
            int32_t n_obs = 0;
            for(int32_t j=0; j < matching_pheno_samp_indices.size(); ++j) {
                int32_t pi = matching_pheno_samp_indices[j];
                if ( !pheno_matrix.pheno_mask( pi, i ) ) continue; // ignore samples with missing phenotype
                r += prs_matrix.pheno_mat( matching_prs_samp_indices[j], i ) * pheno_matrix.pheno_mat( pi, i );
                ++n_obs;
            }
            n_obs_per_trait[i] = n_obs;
            weights[i] = ( n_obs > 0 ) ? r / (double)n_obs : 0.0;
        }

        htsFile* wf = hts_open((outf + ".weights.tsv.gz").c_str(), "wz");
        if ( wf == NULL ) {
            error("Cannot open output file %s.weights.tsv for writing", outf.c_str());
        }
        hprintf(wf, "Trait\tWeight\tN.Obs\n");
        for ( int32_t i = 0; i < weights.size(); ++i) {
            hprintf(wf, "%s\t%.6f\t%d\n", pheno_matrix.pheno_ids[i].c_str(), weights[i], n_obs_per_trait[i]);
        }
        hts_close(wf);
    }
    else {
        notice("Using provided weights for each phenotype from %s", weightf.c_str());
        for(int32_t i=0; i < n_traits; ++i) {
            std::map<std::string, double>::const_iterator it = phe_weights.find( pheno_matrix.pheno_ids[i] );
            if ( it == phe_weights.end() ) {
                notice("No weight found for phenotype %s in weights file %s, setting weight to 0", pheno_matrix.pheno_ids[i].c_str(), weightf.c_str());
                weights[i] = 0.0;
            }
            else {
                weights[i] = it->second;
            }
        }
    }

    // set weights with abs value < min_weight to zero
    int32_t n_pass_weights = 0;
    for ( int32_t i = 0; i < weights.size(); ++i ) {
        if ( weights[i] < min_weight ) {
            weights[i] = 0.0;
        }
        else {
            n_pass_weights++;
        }
    }
    notice("%d / %zu phenotypes passed the minimum weight threshold of %.4f", n_pass_weights, weights.size(), min_weight);

    // Matching on fewer than kFewTraits traits is unreliable. Such individuals are still
    // matched and reported, with their trait count in the output, but a warning is printed.
    const int32_t kFewTraits = 10;
    if ( n_pass_weights < kFewTraits ) {
        warning("Only %d phenotypes have a non-zero weight; matching on fewer than %d traits is unreliable", n_pass_weights, kFewTraits);
    }

    // number of observed traits with non-zero weight per phenotyped individual
    std::vector<int32_t> n_traits_obs( n_pheno_samples, 0 );
    int32_t n_no_obs_traits = 0;
    int32_t n_few_obs_traits = 0;
    std::string few_obs_examples;
    for ( int32_t i = 0; i < n_pheno_samples; ++i ) {
        int32_t c = 0;
        for ( int32_t k = 0; k < n_traits; ++k ) {
            if ( weights[k] != 0.0 && pheno_matrix.pheno_mask(i, k) ) ++c;
        }
        n_traits_obs[i] = c;
        if ( c == 0 ) ++n_no_obs_traits;
        else if ( c < kFewTraits ) {
            if ( ++n_few_obs_traits <= 5 ) {
                few_obs_examples += ( few_obs_examples.empty() ? "" : ", " ) + pheno_matrix.samp_ids[i] + " (" + std::to_string(c) + ")";
            }
        }
    }
    if ( n_no_obs_traits > 0 ) {
        notice("%d phenotyped individuals have no observed traits with non-zero weight and will be reported as NO_OBS_TRAITS", n_no_obs_traits);
    }
    if ( n_few_obs_traits > 0 && n_pass_weights >= kFewTraits ) {
        warning("%d phenotyped individuals have fewer than %d observed traits with non-zero weight, so their matches are unreliable, e.g. %s%s. They are still matched; their trait count is in the output",
                n_few_obs_traits, kFewTraits, few_obs_examples.c_str(), n_few_obs_traits > 5 ? ", ..." : "");
    }

    if ( !no_norm && weights.minCoeff() < 0 ) {
        error("Negative trait weights are not allowed when profile norms are used (they arise from --weights or a negative --min-weight). Use --min-weight 0 or higher, or --no-norm");
    }

    // calculate all pair weighted correlations [n_prs x n_pheno] matrix
    Eigen::MatrixXd all_pair_wcor;
    if ( use_mahalanobis ) {
        notice("Using Mahalanobis distance for matching to account for correlation between traits (may take much longer time)");
        notice("Computing PRS covariance matrix");
        Eigen::MatrixXd prs_cov = prs_matrix.pheno_mat.transpose() * prs_matrix.pheno_mat / (double)prs_matrix.pheno_mat.rows();
        // With missing phenotypes, standardized missing cells are 0, so this is a mean-imputed covariance estimate
        notice("Computing phenotype covariance matrix%s", pheno_has_missing ? " (missing values mean-imputed)" : "");
        Eigen::MatrixXd pheno_cov = pheno_matrix.pheno_mat.transpose() * pheno_matrix.pheno_mat / (double)pheno_matrix.pheno_mat.rows();
        notice("Comptuing total covariance matrix");
        Eigen::MatrixXd total_cov = weight_prs_mh * prs_cov + (1.0 - weight_prs_mh) * pheno_cov + 1e-8 * Eigen::MatrixXd::Identity( prs_cov.rows(), prs_cov.cols() );
        if ( auto_lambda ) {
            notice("Tuning lambda on %d mapped samples by maximizing %s (--auto-lambda)", (int32_t)matching_pheno_samp_indices.size(), auto_lambda_metric.c_str());
            LambdaTuner tuner(total_cov, prs_matrix.pheno_mat, pheno_matrix.pheno_mat, pheno_matrix.pheno_mask, pheno_has_missing,
                              weights, matching_prs_samp_indices, matching_pheno_samp_indices, n_traits_obs, auto_lambda_metric, auto_lambda_min_self);
            lambda = tuner.findOptimalLambda(0.0, 1.0, 0.01);
            tuner.evaluate(lambda); // one more pass at the chosen value to report the self-best count
            int32_t n_tuned = (int32_t)tuner.self_rows.size();
            notice("Using automatically tuned lambda = %.4f after %d evaluations; %d / %d mapped samples (%.1f%%) rank their own PRS first",
                   lambda, tuner.n_eval, tuner.last_n_self_best, n_tuned, 100.0 * tuner.last_n_self_best / n_tuned);
            if ( auto_lambda_min_best > 0 && tuner.last_n_self_best < auto_lambda_min_best ) {
                error("--auto-lambda: only %d of %d mapped samples rank their own PRS first at the tuned lambda, fewer than --auto-lambda-min-best %d. This usually indicates a mismatched sample mapping (check --sample-tsv or the sample IDs), which also invalidates the trait weights estimated from it. Lower --auto-lambda-min-best or set it to 0 to proceed anyway",
                      tuner.last_n_self_best, n_tuned, auto_lambda_min_best);
            }
        }
        if ( lambda > 0 ) {
            notice("Adding regularization parameter %.4f to total covariance matrix", lambda);
            total_cov = (1-lambda) * total_cov + lambda * Eigen::MatrixXd::Identity( prs_cov.rows(), prs_cov.cols() );
        }
        notice("Inverting total covariance matrix");
        Eigen::MatrixXd total_cov_inv = total_cov.inverse();

        notice("Constructing weighted Mahalanobis matrix to multiply");
        Eigen::MatrixXd W_sqrt = weights.cwiseSqrt().asDiagonal();
        Eigen::MatrixXd M = W_sqrt * total_cov_inv * W_sqrt;

        // YM(i,.) = y_i^T M restricted to the traits observed for individual i. Because missing cells of Y are 0,
        // masking YM restricts the quadratic form to the observed block M_OO on both sides:
        //   numerator(j,i) = sum_{k,l in O_i} x_jk M_kl y_il
        notice("Computing numerator for all pair weighted correlations between PRS traits and phenotypes");
        Eigen::MatrixXd YM = pheno_matrix.pheno_mat * M;  // n_pheno x p
        if ( pheno_has_missing ) {
            YM = YM.cwiseProduct( pheno_matrix.pheno_mask.cast<double>() );
        }
        Eigen::MatrixXd numerator = prs_matrix.pheno_mat * YM.transpose();  // n_prs x n_pheno

        notice("Computing norms for PRS and phenotype matrices");
        // phenotype norms: sqrt( y_iO^T M_OO y_iO )
        Eigen::VectorXd pheno_norms = YM.cwiseProduct( pheno_matrix.pheno_mat ).rowwise().sum().cwiseSqrt();
        // PRS norms using all traits: sqrt( x_j^T M x_j ). This is the approximation used for all phenotyped
        // individuals unless --mh-exact-norm is set.
        Eigen::VectorXd prs_norms = ( prs_matrix.pheno_mat * M ).cwiseProduct( prs_matrix.pheno_mat ).rowwise().sum().cwiseSqrt();

        all_pair_wcor.resize(numerator.rows(), numerator.cols());
        notice("Computing all pair weighted correlations between PRS traits and phenotypes");
        if ( exact_norm && pheno_has_missing && lambda == 1.0 ) {
            // M is diagonal (= W) at lambda 1, so exact norms over observed traits are a single matrix product
            notice("Computing exact PRS norms over observed traits with the diagonal metric at lambda = 1 (--exact-norm)");
            all_pair_wcor = numerator;
            divide_by_exact_diagonal_norms(all_pair_wcor, prs_matrix.pheno_mat, weights, pheno_matrix.pheno_mask, pheno_norms);
        }
        else {
            for (int32_t i = 0; i < numerator.cols(); ++i) {
                all_pair_wcor.col(i) = numerator.col(i).array() / ( ( prs_norms * pheno_norms(i) ).array() + 1e-100 );
            }
        }

        if ( exact_norm && pheno_has_missing && lambda != 1.0 ) {
            // group phenotyped individuals by missingness pattern and recompute PRS norms as sqrt( x_jO^T M_OO x_jO )
            std::map<std::string, std::vector<int32_t> > pattern2indices;
            for ( int32_t i = 0; i < n_pheno_samples; ++i ) {
                if ( pheno_matrix.pheno_mask.row(i).all() ) continue; // full-trait norm is already exact
                std::string key( n_traits, '0' );
                for ( int32_t k = 0; k < n_traits; ++k ) if ( pheno_matrix.pheno_mask(i, k) ) key[k] = '1';
                pattern2indices[key].push_back(i);
            }
            notice("Computing exact PRS norms for %d unique missingness patterns among phenotyped individuals with missing values (--exact-norm)", (int32_t)pattern2indices.size());
            int32_t n_done = 0;
            for ( std::map<std::string, std::vector<int32_t> >::const_iterator it = pattern2indices.begin(); it != pattern2indices.end(); ++it ) {
                std::vector<int32_t> obs;
                for ( int32_t k = 0; k < n_traits; ++k ) if ( it->first[k] == '1' ) obs.push_back(k);
                if ( obs.empty() ) continue; // no observed traits: numerator is zero anyway
                Eigen::MatrixXd X_O = prs_matrix.pheno_mat( Eigen::placeholders::all, obs );
                Eigen::MatrixXd M_OO = M( obs, obs );
                Eigen::VectorXd norms_O = ( X_O * M_OO ).cwiseProduct( X_O ).rowwise().sum().cwiseSqrt();
                for ( size_t t = 0; t < it->second.size(); ++t ) {
                    int32_t i = it->second[t];
                    all_pair_wcor.col(i) = numerator.col(i).array() / ( ( norms_O * pheno_norms(i) ).array() + 1e-100 );
                }
                if ( ( ++n_done % 100 ) == 0 ) notice("Processed %d / %d missingness patterns", n_done, (int32_t)pattern2indices.size());
            }
        }
    }
    else {
        // missing phenotype cells are 0, so they drop out of the numerator sum_k w_k x_jk y_ik
        all_pair_wcor = prs_matrix.pheno_mat * ( weights.asDiagonal() * pheno_matrix.pheno_mat.transpose() );
        if ( no_norm ) {
            // legacy score: divide by the sum of absolute weights over the traits observed for each phenotyped individual
            notice("Computing all pair weighted inner products between PRS traits and phenotypes assuming independence (--no-norm)");
            Eigen::VectorXd abs_weights = weights.cwiseAbs();
            Eigen::VectorXd denom;
            if ( pheno_has_missing ) {
                denom = pheno_matrix.pheno_mask.cast<double>() * abs_weights;
            }
            else {
                denom = Eigen::VectorXd::Constant( n_pheno_samples, abs_weights.sum() );
            }
            for ( int32_t i = 0; i < n_pheno_samples; ++i ) {
                all_pair_wcor.col(i) /= ( denom(i) + 1e-100 );
            }
        }
        else {
            // weighted cosine similarity: identical to --mahalanobis --lambda 1, where M = diag(w)
            notice("Computing all pair weighted cosine similarities between PRS traits and phenotypes assuming independence");
            Eigen::VectorXd pheno_norms = ( pheno_matrix.pheno_mat.cwiseAbs2() * weights ).cwiseSqrt(); // sqrt( sum_{k in O_i} w_k y_ik^2 )
            if ( exact_norm && pheno_has_missing ) {
                notice("Computing exact PRS norms over observed traits (--exact-norm)");
                divide_by_exact_diagonal_norms(all_pair_wcor, prs_matrix.pheno_mat, weights, pheno_matrix.pheno_mask, pheno_norms);
            }
            else {
                Eigen::VectorXd prs_norms = ( prs_matrix.pheno_mat.cwiseAbs2() * weights ).cwiseSqrt();   // sqrt( sum_k w_k x_jk^2 )
                for ( int32_t i = 0; i < n_pheno_samples; ++i ) {
                    all_pair_wcor.col(i) = all_pair_wcor.col(i).array() / ( ( prs_norms * pheno_norms(i) ).array() + 1e-100 );
                }
            }
        }
    }

    notice("Standardizing all pair weighted correlation matrix");
    // copy the weighted correlation matrix
    Eigen::MatrixXd all_pair_z = all_pair_wcor;
    // Apply column-wise standardization
    standardize_matrix_columns_inplace(all_pair_z);

    // open the output file gz or plain based on the extension
    htsFile* wf1 = hts_open((outf + ".match.assigned.tsv.gz").c_str(), "wz");
    if ( wf1 == NULL ) {
        error("Cannot open output file %s.match.assigned.tsv.gz for writing", outf.c_str());
    }
    htsFile* wf2 = hts_open((outf + ".match.all.tsv.gz").c_str(), "wz");
    if ( wf2 == NULL ) {
        error("Cannot open output file %s.match.all.tsv.gz for writing", outf.c_str());
    }
 
    // write the header line
    hprintf(wf1, "ID.Pheno\tN.Traits\tMatchStatus\tID.self\tZ.self\tCOR.self\tRank.self\n");
    hprintf(wf2, "ID.Pheno\tN.Traits\tMatchStatus\tID.self\tZ.self\tCOR.self\tRank.self\tID.1st\tZ.1st\tCOR.1st\tID.2nd\tZ.2nd\tCOR.2nd\tID.3rd\tZ.3rd\tCOR.3rd\tID.4th\tZ.4th\tCOR.4th\tID.5th\tZ.5th\tCOR.5th\n");
    for ( int32_t i = 0; i < n_pheno_samples; ++i ) {
        int32_t self_idx = -1;
        if ( samp_idx_pheno2prs.find(i) != samp_idx_pheno2prs.end() ) {
            self_idx = samp_idx_pheno2prs[i];
        }

        // individuals without any observed trait cannot be scored
        if ( n_traits_obs[i] == 0 ) {
            if ( self_idx >= 0 ) {
                hprintf(wf1, "%s\t0\tNO_OBS_TRAITS\t%s\tNA\tNA\tNA\n",
                    pheno_matrix.samp_ids[i].c_str(), prs_matrix.samp_ids[self_idx].c_str());
            }
            hprintf(wf2, "%s\t0\tNO_OBS_TRAITS\t%s\tNA\tNA\tNA",
                pheno_matrix.samp_ids[i].c_str(), self_idx >= 0 ? prs_matrix.samp_ids[self_idx].c_str() : "NA");
            for ( int32_t k = 0; k < 5; ++k ) hprintf(wf2, "\tNA\tNA\tNA");
            hprintf(wf2, "\n");
            continue;
        }

        // find the best and second best matches
        int32_t self_rank = 1;
        double z_self = self_idx >=0 ? all_pair_z(self_idx, i) : -9999.0;
        double cor_self = self_idx >=0 ? all_pair_wcor(self_idx, i) : -9999.0;
        double idx_top[5] = { -1, -1, -1, -1, -1 };
        double z_top[5] = { -9999.0, -9999.0, -9999.0, -9999.0, -9999.0 };
        double cor_top[5] = { -9999.0, -9999.0, -9999.0, -9999.0, -9999.0 };
        for( int32_t j = 0; j < n_prs_samples; ++j ) {
            double z = all_pair_z(j, i);
            double cor = all_pair_wcor(j, i);
            if ( self_idx >= 0 && j != self_idx && z > z_self ) {
                self_rank++;
            }
            // check if this is in the top 5
            for ( int32_t k = 0; k < 5; ++k ) {
                if ( z > z_top[k] ) {
                    // shift down
                    for ( int32_t l = 4; l > k; --l ) {
                        z_top[l] = z_top[l-1];
                        cor_top[l] = cor_top[l-1];
                        idx_top[l] = idx_top[l-1];
                    }
                    z_top[k] = z;
                    cor_top[k] = cor;
                    idx_top[k] = j;
                    break;
                }
            }
        }
        if ( self_idx >= 0 ) {
            const char* match_status = ( self_rank == 1 ? "BEST_MATCH" : ( z_self < z_lenient_threshold ? "NO_MATCH" : "LENIENT_MATCH" ) );
            hprintf(wf1, "%s\t%d\t%s\t%s\t%.6f\t%.6f\t%d\n",
                pheno_matrix.samp_ids[i].c_str(),
                n_traits_obs[i],
                match_status,
                prs_matrix.samp_ids[self_idx].c_str(),
                z_self, cor_self,
                self_rank);
        }
        // self_rank is only meaningful when the individual has a mapped PRS sample; without one, the status is
        // NO_SELF unless the top matches are clearly separated
        const char* match_status = ( self_idx >= 0 ) ? "UNCLEAR" : "NO_SELF";
        if ( self_idx >= 0 && self_rank == 1 ) {
            match_status = "SELF_BEST";
        }
        else if ( z_top[0] > z_top[1] + z_diff_threshold ) {
            match_status = "SINGLE_NEW_BEST";
        }
        else {
            if ( self_idx >= 0 && z_self > z_lenient_threshold ) {
                match_status = "SELF_LENIENT";
            }
            for(int32_t k=2; k < 5; ++k ) {
                if ( z_top[k-1] > z_top[k] + z_diff_threshold ) {
                    match_status = "MULTI_NEW_BEST";
                    break;
                }
            }
        }
        if ( self_idx >= 0 ) {
            hprintf(wf2, "%s\t%d\t%s\t%s\t%.6f\t%.6f\t%d",
                pheno_matrix.samp_ids[i].c_str(),
                n_traits_obs[i],
                match_status,
                prs_matrix.samp_ids[self_idx].c_str(),
                z_self, cor_self,
                self_rank);
        }
        else {
            hprintf(wf2, "%s\t%d\t%s\tNA\tNA\tNA\tNA",
                pheno_matrix.samp_ids[i].c_str(),
                n_traits_obs[i],
                match_status);
        }
        for ( int32_t k = 0; k < 5; ++k ) {
            if ( idx_top[k] >=0 ) {
                hprintf(wf2, "\t%s\t%.6f\t%.6f",
                    prs_matrix.samp_ids[ (int32_t)idx_top[k] ].c_str(),
                    z_top[k],
                    cor_top[k]);
            }
            else {
                hprintf(wf2, "\tNA\tNA\tNA");
            }
        }
        hprintf(wf2, "\n");
    }
    hts_close(wf1); // close the output file
    hts_close(wf2); // close the output file

    notice("Analysis finished");
    return 0;
}
