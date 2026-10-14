#include "assoc_utils.h"
#include "qgenlib/qgen_error.h"
#include "qgenlib/tsv_reader.h"
#include <cmath>
#include <cstring>
#include <memory>
#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <sys/stat.h>


// calculate columnwise dot product between two matrices
Eigen::VectorXd columnwise_dot(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2) {
    const int32_t n_cols = mat1.cols();
    const int32_t n_rows = mat1.rows();

    if ( n_cols != mat2.cols() || n_rows != mat2.rows() ) {
        error("columnwise_dot(): Matrices must have the same dimensions for columnwise dot product. Got %d x %d and %d x %d.",
              n_rows, n_cols, mat2.rows(), mat2.cols());
    }
    Eigen::VectorXd result(n_cols);
    for ( int32_t j = 0; j < n_cols; ++j ) {
        result(j) = mat1.col(j).dot(mat2.col(j));
    }
    return result;
}
    

void standardize_matrix_columns_inplace(Eigen::MatrixXd& matrix) {
    const long n_cols = matrix.cols();
    const long n_rows = matrix.rows();

    notice("Standardizing columns of matrix with dimensions %d x %d.", n_rows, n_cols);

    for (long j = 0; j < n_cols; ++j) {
        Eigen::VectorXd col = matrix.col(j);
        double mean = col.mean();
        double stddev = std::sqrt((col.array() - mean).square().sum() / (n_rows - 1));

        if (stddev > 0) {
            matrix.col(j) = (col.array() - mean) / stddev;
        } else {
            // If stddev is zero, all values are the same. Set to zero vector.
            matrix.col(j).setZero();
        }
    }
}

Eigen::MatrixXd bulk_adjust_for_covariates(const Eigen::MatrixXd& values, const Eigen::MatrixXd& covariates) {
    const long g = values.rows();
    const long n = values.cols();
    const long p = covariates.rows();

    // --- Input Validation ---
    // Ensure the number of observations (columns) is the same in both matrices.
    if (n != covariates.cols()) {
        throw std::runtime_error("Value and covariate matrices must have the same number of columns (observations).");
    }

    // --- Transpose Inputs for Regression ---
    // The standard regression setup is (observations x variables).
    // We transpose our g x n and p x n inputs to n x g and n x p.
    Eigen::MatrixXd V_t = values.transpose(); // V_t is n x g
    Eigen::MatrixXd C_t = covariates.transpose(); // C_t is n x p

    // --- Augment Covariate Matrix with Intercept ---
    // We create an augmented covariate matrix `C_aug` of size n x (p+1).
    Eigen::MatrixXd C_aug(n, p + 1);
    C_aug.col(0).setOnes(); // First column is the intercept (all ones).
    C_aug.rightCols(p) = C_t; // The remaining columns are the transposed covariates.

    // --- Solve the Least Squares Problem ---
    // We solve for the coefficient matrix `beta` ((p+1) x g) where `V_t = C_aug * beta`.
    // Eigen's `colPivHouseholderQr()` provides a robust way to solve this.
    Eigen::MatrixXd beta = C_aug.colPivHouseholderQr().solve(V_t);

    // --- Calculate Residuals ---
    // The predicted values are `C_aug * beta`.
    // The residuals are the original (transposed) values minus the predicted values.
    // The resulting residual matrix is in n x g format.
    Eigen::MatrixXd residuals_t = V_t - (C_aug * beta);

    // --- Transpose Result Back to g x n Format ---
    return residuals_t.transpose();
}

CovariateResidualizer::CovariateResidualizer(const Eigen::MatrixXd& covariates)
    : n_samples_(covariates.rows()), n_cols_((int32_t)covariates.cols() + 1)
{
    // [1, C] with unit-norm columns: the first column is the intercept, the remaining columns
    // are the covariates. A column that is all zero stays zero and is not counted in the rank.
    Eigen::MatrixXd C_scaled(covariates.rows(), covariates.cols() + 1);
    C_scaled.col(0).setOnes();
    C_scaled.rightCols(covariates.cols()) = covariates;
    for (Eigen::Index j = 0; j < C_scaled.cols(); ++j) {
        const double nrm = C_scaled.col(j).norm();
        if ( nrm > 0 ) C_scaled.col(j) /= nrm;
    }
    qr.setThreshold(kRankTol);
    qr.compute(C_scaled);
    rank_ = (int32_t)qr.rank();
}

void CovariateResidualizer::residualize(Eigen::Ref<Eigen::MatrixXd> values) const {
    if ( values.rows() != n_samples_ ) {
        throw std::runtime_error("Value and covariate matrices must have the same number of rows (observations).");
    }
    if ( values.cols() == 0 ) return;
    // r = (I - Q_r Q_r') v with the first rank_ Householder vectors: the pivots beyond the rank
    // span rounding noise (e.g. a constant or duplicated covariate) and are not projected out
    Eigen::MatrixXd c = values;
    const auto Q = qr.householderQ().setLength(rank_);
    c.applyOnTheLeft(Q.adjoint());
    c.topRows(rank_).setZero();
    c.applyOnTheLeft(Q);
    values = c;
}

int32_t CovariateResidualizer::residualize_genotypes(Eigen::Ref<Eigen::MatrixXd> G) const {
    const Eigen::Index kBlockCols = 1024;
    const Eigen::Index n = G.rows();
    int32_t n_collinear = 0;
    for (Eigen::Index j0 = 0; j0 < G.cols(); j0 += kBlockCols) {
        const Eigen::Index w = std::min<Eigen::Index>(kBlockCols, G.cols() - j0);
        auto B = G.middleCols(j0, w);
        // centered sum of squares = sum(x^2) - n * mean(x)^2
        const Eigen::RowVectorXd col_means = B.colwise().mean();
        const Eigen::RowVectorXd ss_centered = B.colwise().squaredNorm() - ( (double)n * col_means.array().square() ).matrix();
        residualize(B);
        for (Eigen::Index j = 0; j < w; ++j) {
            if ( !( B.col(j).squaredNorm() > kCollinearTol * ss_centered(j) ) ) {
                B.col(j).setZero();
                ++n_collinear;
            }
        }
    }
    return n_collinear;
}

Eigen::RowVectorXd CovariateResidualizer::residual_ss(const Eigen::MatrixXd& Y) const {
    Eigen::RowVectorXd rss(Y.cols());
    Eigen::MatrixXd col(Y.rows(), 1);
    for (Eigen::Index j = 0; j < Y.cols(); ++j) {
        col = Y.col(j);
        residualize(col);
        rss(j) = col.squaredNorm();
    }
    return rss;
}

Eigen::MatrixXd pheno_adj_cov_nxt_without_missing(const Eigen::MatrixXd& values, const Eigen::MatrixXd& covariates) {
    // Ensure the number of observations (rows) is the same in both matrices.
    if (values.rows() != covariates.rows()) {
        throw std::runtime_error("Value and covariate matrices must have the same number of rows (observations).");
    }
    // Least squares on [1, C] via a column-pivoting Householder QR, computed once for all
    // columns of values; see CovariateResidualizer.
    const CovariateResidualizer resid(covariates);
    Eigen::MatrixXd residuals = values;
    resid.residualize(residuals);
    return residuals;
}

void standardize_matrix_columns_inplace(Eigen::MatrixXd& matrix,
                                        const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& mask) {
    const long n_cols = matrix.cols();
    const long n_rows = matrix.rows();
    if ( mask.rows() != n_rows || mask.cols() != n_cols ) {
        error("standardize_matrix_columns_inplace: mask dimensions (%ld x %ld) do not match matrix dimensions (%ld x %ld)",
              (long)mask.rows(), (long)mask.cols(), n_rows, n_cols);
    }

    notice("Standardizing columns of matrix with dimensions %ld x %ld using observed values only.", n_rows, n_cols);

    for (long j = 0; j < n_cols; ++j) {
        // two-pass mean / sd over observed cells
        double sum = 0.0;
        long n_obs = 0;
        for (long i = 0; i < n_rows; ++i) {
            if ( mask(i, j) ) { sum += matrix(i, j); ++n_obs; }
        }
        if ( n_obs < 2 ) {
            matrix.col(j).setZero();
            continue;
        }
        const double mean = sum / (double)n_obs;
        double ss = 0.0;
        for (long i = 0; i < n_rows; ++i) {
            if ( mask(i, j) ) { double d = matrix(i, j) - mean; ss += d * d; }
        }
        const double stddev = std::sqrt(ss / (double)(n_obs - 1));
        if ( stddev > 0 ) {
            for (long i = 0; i < n_rows; ++i) {
                matrix(i, j) = mask(i, j) ? (matrix(i, j) - mean) / stddev : 0.0;
            }
        }
        else {
            matrix.col(j).setZero();
        }
    }
}

Eigen::MatrixXd pheno_adj_cov_nxt_with_missing(const Eigen::MatrixXd& values,
                                               const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& mask,
                                               const Eigen::MatrixXd& covariates) {
    const long n = values.rows();
    const long g = values.cols();
    const long p = covariates.cols();

    if ( n != covariates.rows() ) {
        throw std::runtime_error("Value and covariate matrices must have the same number of rows (observations).");
    }
    if ( mask.rows() != n || mask.cols() != g ) {
        throw std::runtime_error("Value matrix and mask must have the same dimensions.");
    }

    Eigen::MatrixXd C_aug(n, p + 1);
    C_aug.col(0).setOnes();
    C_aug.rightCols(p) = covariates;

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, g);

    // split traits into complete and incomplete columns
    std::vector<long> complete_cols, incomplete_cols;
    for (long j = 0; j < g; ++j) {
        if ( mask.col(j).all() ) complete_cols.push_back(j);
        else incomplete_cols.push_back(j);
    }
    notice("Adjusting %ld complete traits in batch and %ld traits with missing values individually",
           (long)complete_cols.size(), (long)incomplete_cols.size());

    // complete traits: one QR decomposition shared across all of them
    if ( !complete_cols.empty() ) {
        Eigen::MatrixXd V(n, (long)complete_cols.size());
        for (size_t k = 0; k < complete_cols.size(); ++k) V.col(k) = values.col(complete_cols[k]);
        Eigen::MatrixXd beta = C_aug.colPivHouseholderQr().solve(V);
        Eigen::MatrixXd R = V - C_aug * beta;
        for (size_t k = 0; k < complete_cols.size(); ++k) result.col(complete_cols[k]) = R.col(k);
    }

    // incomplete traits: regress on observed rows only
    int32_t n_underdetermined = 0;
    for (size_t k = 0; k < incomplete_cols.size(); ++k) {
        const long j = incomplete_cols[k];
        std::vector<long> obs_rows;
        obs_rows.reserve(n);
        for (long i = 0; i < n; ++i) if ( mask(i, j) ) obs_rows.push_back(i);
        const long n_obs = (long)obs_rows.size();
        if ( n_obs <= p + 1 ) {
            // too few observations to fit the covariate model; leave residuals at zero so that the
            // trait becomes constant and receives zero weight downstream
            ++n_underdetermined;
            continue;
        }
        Eigen::MatrixXd C_obs(n_obs, p + 1);
        Eigen::VectorXd y_obs(n_obs);
        for (long r = 0; r < n_obs; ++r) {
            C_obs.row(r) = C_aug.row(obs_rows[r]);
            y_obs(r) = values(obs_rows[r], j);
        }
        Eigen::VectorXd beta = C_obs.colPivHouseholderQr().solve(y_obs);
        Eigen::VectorXd resid = y_obs - C_obs * beta;
        for (long r = 0; r < n_obs; ++r) result(obs_rows[r], j) = resid(r);
    }
    if ( n_underdetermined > 0 ) {
        warning("%d traits had no more observed values than covariates (plus intercept) and were set to zero residuals",
                n_underdetermined);
    }
    return result;
}

Eigen::MatrixXd rint_matrix_with_missing(const Eigen::MatrixXd& matrix,
                                         const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& mask) {
    const int32_t n_rows = matrix.rows();
    const int32_t n_cols = matrix.cols();
    if ( mask.rows() != n_rows || mask.cols() != n_cols ) {
        throw std::invalid_argument("rint_matrix_with_missing: matrix and mask must have the same dimensions.");
    }
    Eigen::MatrixXd result(n_rows, n_cols);
    for (int32_t j = 0; j < n_cols; ++j) {
        if ( mask.col(j).all() ) {
            result.col(j) = rint_without_missing(matrix.col(j));
        }
        else {
            Eigen::VectorXd r = rint_with_missing(matrix.col(j), mask.col(j));
            for (int32_t i = 0; i < n_rows; ++i) {
                result(i, j) = mask(i, j) ? r(i) : 0.0;
            }
        }
    }
    return result;
}

void center_rows(Eigen::MatrixXd &matrix) {
    // We iterate over each row of the matrix.
    // The .rowwise() method allows us to apply an operation to each row.
    // The .mean() method calculates the mean of the elements in a row.
    // This creates a column vector where each element is the mean of the corresponding row in the original matrix.
    Eigen::VectorXd row_means = matrix.rowwise().mean();

    notice("Centering rows of matrix with dimensions %d x %d.", matrix.rows(), matrix.cols());

    matrix.colwise() -= row_means;

    // // We then iterate over each row again to subtract the calculated mean.
    // for (int i = 0; i < matrix.rows(); ++i) {
    //     // .noalias() is a performance optimization in Eigen. It tells Eigen that
    //     // the result of the expression on the right-hand side does not alias
    //     // with the matrix on the left-hand side, allowing for a more efficient
    //     // subtraction operation without creating temporary objects.
    //     matrix.row(i).noalias() -= Eigen::VectorXd::Constant(matrix.cols(), row_means(i));
    // }
}

/**
 * @brief Computes the continued fraction part for the incomplete beta function.
 *
 * This is an implementation of Lentz's method as described in Numerical Recipes.
 * It's a key component for accurately calculating the incomplete beta function.
 *
 * @param a Parameter 'a' of the beta function.
 * @param b Parameter 'b' of the beta function.
 * @param x The upper limit of integration.
 * @return The result of the continued fraction evaluation.
 */
double betacf(double a, double b, double x) {
    const int MAXIT = 200;
    const double EPS = 1.0e-14;
    const double FPMIN = 1.0e-30;

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= MAXIT; m++) {
        int m2 = 2 * m;
        // Even step
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        // Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < EPS) {
            return h;
        }
    }
    // This warning can be enabled for debugging convergence issues.
    // std::cerr << "Warning: betacf did not converge for a=" << a << ", b=" << b << ", x=" << x << std::endl;
    return h;
}

/**
 * @brief Computes the natural logarithm of the regularized incomplete beta function I_x(a,b).
 *
 * This function computes log(I_x(a,b)) using a continued fraction evaluation (betacf).
 * It uses a symmetry property to maintain accuracy when x is large.
 * The formula used is:
 * log(I_x(a,b)) = a*log(x) + b*log(1-x) - log(a) - log(B(a,b)) + log(betacf(a,b,x))
 * where B(a,b) is the beta function.
 *
 * @param x The upper limit of integration (must be in [0, 1]).
 * @param a Parameter 'a' of the beta function (must be > 0).
 * @param b Parameter 'b' of the beta function (must be > 0).
 * @return The natural logarithm of I_x(a,b).
 */
double betaincln(double x, double a, double b) {
    if (x < 0.0 || x > 1.0) return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return -std::numeric_limits<double>::infinity();
    if (x == 1.0) return 0.0;

    // Use symmetry property for better accuracy when x is large
    if (x > (a + 1.0) / (a + b + 2.0)) {
        // Compute log(1 - I_{1-x}(b,a))
        // First compute log(I_{1-x}(b,a))
        double log_beta_func_sym = std::lgamma(b) + std::lgamma(a) - std::lgamma(b + a);
        double cf_sym = betacf(b, a, 1.0 - x);
        double log_power_term_sym = b * std::log(1.0 - x) + a * std::log(x);
        double log_i_sym = log_power_term_sym - std::log(b) - log_beta_func_sym + std::log(cf_sym);

        // log(1-y) is computed as log1p(-y) for precision when y is small.
        return std::log1p(-std::exp(log_i_sym));
    } else {
        // Direct computation for smaller x
        double log_beta_func = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
        double cf = betacf(a, b, x);
        double log_power_term = a * std::log(x) + b * std::log(1.0 - x);
        return log_power_term - std::log(a) - log_beta_func + std::log(cf);
    }
}

/**
 * @brief Converts a t-statistic to a 2-sided -log10 p-value.
 *
 * This function calculates the p-value for a given t-statistic and degrees
 * of freedom. It is designed to be numerically stable even for extremely
 * large t-statistic values by working in log-space to avoid underflow.
 * The p-value is derived from the student's t-distribution cumulative
 * distribution function, which is related to the regularized incomplete
 * beta function.
 *
 * @param tstat The t-statistic value.
 * @param df The degrees of freedom.
 * @return The 2-sided -log10 p-value. Returns NaN if df <= 0.
 */
double tstat2log10pval(double tstat, double df) {
    if (df <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double abs_t = std::abs(tstat);

    // The 2-sided p-value for the t-distribution is given by I_x(df/2, 1/2),
    // where x = df / (df + t^2).
    double x = df / (df + abs_t * abs_t);

    // To maintain numerical precision for very small p-values, we compute
    // the natural logarithm of the p-value directly.
    double log_p_value = betaincln(x, df / 2.0, 0.5);

    // Convert the natural log p-value to -log10 p-value.
    // -log10(p) = -ln(p) / ln(10)
    return -(log_p_value / std::log(10.0));
}
/**
 * @brief Performs simple linear regression for each column of a matrix against a vector.
 * * Assumes that the input vector y and each column of the matrix X are centered
 * (have a mean of zero).
 *
 * @param y An Eigen::VectorXd of size n (dependent variable).
 * @param X An Eigen::MatrixXd of size n x p (independent variables).
 * @return A std::vector of RegressionResult structs, one for each column of X.
 */
bool simple_linear_regression_without_missing(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, std::vector<slr_sumstat_t>& results,
                                              int32_t n_cov_df, double y_rss) {
    int n = y.size();
    int p = X.cols();

    // for(int32_t i = 0; i < 5; ++i) {
    //     notice("y[%d] = %.5g, X[%d,0] = %.5g", i, y[i], i, X(i, 0));
    // }

    if (n != X.rows()) {
        error("Vector and Matrix dimensions do not match : (%d x %d) vs (%d x %d)", n, 1, X.rows(), X.cols());
    }
    
    // Degrees of freedom for the t-distribution
    // For simple linear regression (one predictor), it's n - 2; covariates regressed out of
    // y and X (Frisch-Waugh-Lovell) use n_cov_df more
    int df = n - 2 - n_cov_df;
    if (df <= 0) {
        error("Not enough data points to perform regression (n - 2 - n_cov_df = %d must be > 0).", df);
    }

    results.clear();
    results.resize(p);

    // 1. Calculate squared norms for each column of X (x_i' * x_i)
    Eigen::VectorXd x_col_sq_norms = X.colwise().squaredNorm();

    // 2. Calculate X' * y
    Eigen::VectorXd xt_y = X.transpose() * y;

    // 3. Calculate betas for all columns
    // beta_i = (x_i' * y) / (x_i' * x_i)
    Eigen::ArrayXd betas = xt_y.array() / x_col_sq_norms.array();

    // 4. Calculate Sum of Squared Errors (SSE) for each regression
    // sse_i = y'y - (x_i'y)^2 / (x_i'x_i), with y'y replaced by the null-model RSS if given
    double y_sq_norm = ( y_rss >= 0 ) ? y_rss : y.squaredNorm();
    Eigen::ArrayXd sse = y_sq_norm - (xt_y.array().square() / x_col_sq_norms.array());

    // 5. Calculate standard errors for all betas
    // se(beta_i) = sqrt( (sse_i / df) / (x_i'x_i) )
    Eigen::ArrayXd std_errors = (sse / (df * x_col_sq_norms.array())).sqrt();

    // 6. Calculate t-statistics for all betas
    Eigen::ArrayXd t_values = betas / std_errors;

    for (int i = 0; i < p; ++i) {
        slr_sumstat_t& ss = results[i];
        ss.beta = betas[i];
        ss.se = std_errors[i];
        ss.tstat = t_values[i];
        ss.log10p = tstat2log10pval(ss.tstat, df);
        ss.n_obs = n; // Store the number of observations
    }

    //error("stop");

    return true;
}

bool simple_rect_regression_without_missing(
    const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& X,
    std::vector<std::vector<slr_sumstat_t>>& results,
    int32_t n_cov_df,
    const Eigen::RowVectorXd* y_rss) {
    const int n = Y.rows();

    if (n != X.rows()) {
        error("Matrix dimensions do not match : Y is (%d x %d) while X is (%d x %d)", Y.rows(), Y.cols(), X.rows(), X.cols());
    }

    if ( y_rss != NULL && y_rss->size() != Y.cols() ) {
        error("Null-model residual sums of squares are given for %d traits, but Y has %d columns", (int32_t)y_rss->size(), (int32_t)Y.cols());
    }

    const Eigen::VectorXd x_sq_norms = X.colwise().squaredNorm();
    const Eigen::RowVectorXd y_sq_norms = ( y_rss != NULL ) ? *y_rss : Eigen::RowVectorXd(Y.colwise().squaredNorm());
    const Eigen::MatrixXd xt_y = X.transpose() * Y;
    return simple_rect_regression_from_stats(n, x_sq_norms, y_sq_norms, xt_y, results, n_cov_df);
}

bool simple_rect_regression_from_stats(
    int n,
    const Eigen::VectorXd& x_sq_norms,
    const Eigen::RowVectorXd& y_sq_norms,
    const Eigen::MatrixXd& xt_y,
    std::vector<std::vector<slr_sumstat_t>>& results,
    int32_t n_cov_df) {
    const int num_x = (int)xt_y.rows();
    const int num_y = (int)xt_y.cols();

    // n - 2 for the intercept and the variant, minus the covariate degrees of freedom
    const int df = n - 2 - n_cov_df;
    if (df <= 0) {
        error("Not enough data points to perform regression (n - 2 - n_cov_df = %d must be > 0).", df);
    }

    results.clear();
    results.resize(num_x);

    if (num_x == 0 || num_y == 0) {
        return true;
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double df_double = static_cast<double>(df);

    for (int i = 0; i < num_x; ++i) {
        auto& row_results = results[i];
        row_results.resize(num_y);
        const double x_norm = x_sq_norms(i);

        if (x_norm <= 0.0) {
            for (int j = 0; j < num_y; ++j) {
                slr_sumstat_t& ss = row_results[j];
                ss.beta = nan;
                ss.se = nan;
                ss.tstat = nan;
                ss.log10p = nan;
                ss.n_obs = n;
            }
            continue;
        }

        const Eigen::ArrayXXd xty_row = xt_y.row(i).array();
        const Eigen::ArrayXXd betas = xty_row / x_norm;
        const Eigen::ArrayXXd sse = (y_sq_norms.array() - (xty_row.square() / x_norm)).max(0.0);
        const Eigen::ArrayXXd std_errors = ((sse / df_double) / x_norm).sqrt();

        for (int j = 0; j < num_y; ++j) {
            slr_sumstat_t& ss = row_results[j];
            const double beta = betas(0, j);
            const double se = std_errors(0, j);
            const bool valid_se = std::isfinite(se) && se > 0.0;
            const double tstat = valid_se ? beta / se : nan;
            const double log10p = valid_se ? tstat2log10pval(tstat, df) : nan;

            ss.beta = beta;
            ss.se = se;
            ss.tstat = tstat;
            ss.log10p = log10p;
            ss.n_obs = n;
        }
    }

    return true;
}

// -log10 of the two-sided standard normal p-value 2 * Phi(-|z|).
// Uses erfc() while it does not underflow, and the asymptotic Mills-ratio series
// for the far tail so that LOG10P stays finite for extreme statistics.
double zstat2log10pval(double zstat) {
    if ( !std::isfinite(zstat) ) {
        return std::isnan(zstat) ? std::numeric_limits<double>::quiet_NaN() : std::numeric_limits<double>::infinity();
    }
    const double z = std::fabs(zstat);
    if ( z < 30.0 ) {
        return std::max(0.0, -std::log10(std::erfc(z / std::sqrt(2.0)))); // erfc(z/sqrt2) = 2 * Phi(-z); max() avoids -0
    }
    // log Phi(-z) = -z^2/2 - log(z) - log(2*pi)/2 + log(1 - 1/z^2 + 3/z^4 - 15/z^6 + 105/z^8)
    const double z2 = z * z;
    const double series = 1.0 - 1.0 / z2 + 3.0 / (z2 * z2) - 15.0 / (z2 * z2 * z2) + 105.0 / (z2 * z2 * z2 * z2);
    const double log_p = std::log(2.0) - 0.5 * z2 - std::log(z) - 0.5 * std::log(2.0 * M_PI) + std::log(series);
    return -log_p / std::log(10.0);
}

// Score test for each (variant, trait) pair, matching REGENIE step 2 for
// quantitative traits. Y and X must be centered and covariate-residualized.
// The residual variance is estimated under the null model (no variant effect),
// sigma0^2 = y'y / (n - n_cov), where n_cov counts the intercept. Then
//   z    = x'y / sqrt(x'x * sigma0^2)
//   beta = x'y / x'x
//   se   = sqrt(sigma0^2 / x'x) = |beta / z|
// and log10p is from the two-sided standard normal. The z statistic is stored
// in the tstat field of slr_sumstat_t.
bool simple_rect_score_test_without_missing(
    const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& X,
    int32_t n_cov,
    std::vector<std::vector<slr_sumstat_t> >& results) {
    const int n = Y.rows();

    if (n != X.rows()) {
        error("Matrix dimensions do not match : Y is (%d x %d) while X is (%d x %d)", Y.rows(), Y.cols(), X.rows(), X.cols());
    }
    const Eigen::VectorXd x_sq_norms = X.colwise().squaredNorm();
    const Eigen::RowVectorXd y_sq_norms = Y.colwise().squaredNorm();
    const Eigen::MatrixXd xt_y = X.transpose() * Y;
    return simple_rect_score_test_from_stats(n, x_sq_norms, y_sq_norms, xt_y, n_cov, results);
}

bool simple_rect_score_test_from_stats(
    int n,
    const Eigen::VectorXd& x_sq_norms,
    const Eigen::RowVectorXd& y_sq_norms,
    const Eigen::MatrixXd& xt_y,
    int32_t n_cov,
    std::vector<std::vector<slr_sumstat_t> >& results) {
    const int num_x = (int)xt_y.rows();
    const int num_y = (int)xt_y.cols();
    const int df = n - n_cov;
    if (df <= 0) {
        error("Not enough data points for the score test (n - n_cov = %d must be > 0).", df);
    }

    results.clear();
    results.resize(num_x);
    if (num_x == 0 || num_y == 0) {
        return true;
    }

    const Eigen::RowVectorXd sigma0_sq = y_sq_norms / (double)df;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (int i = 0; i < num_x; ++i) {
        std::vector<slr_sumstat_t>& row_results = results[i];
        row_results.resize(num_y);
        const double x_norm = x_sq_norms(i);
        for (int j = 0; j < num_y; ++j) {
            slr_sumstat_t& ss = row_results[j];
            ss.n_obs = n;
            const double se = (x_norm > 0.0) ? std::sqrt(sigma0_sq(j) / x_norm) : nan;
            if ( !(std::isfinite(se) && se > 0.0) ) {
                ss.beta = ss.se = ss.tstat = ss.log10p = nan;
                continue;
            }
            ss.beta = xt_y(i, j) / x_norm;
            ss.se = se;
            ss.tstat = ss.beta / se; // z statistic
            ss.log10p = zstat2log10pval(ss.tstat);
        }
    }
    return true;
}

bool simple_linear_regression_with_missing( const Eigen::VectorXd& y,
                                            const Eigen::Vector<bool, Eigen::Dynamic>& y_mask,
                                            const Eigen::MatrixXd& X,
                                            const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& X_mask,
                                            std::vector<slr_sumstat_t>& results,
                                            int32_t n_cov_df) {
    const int n_total = y.size();
    const int p = X.cols();

    // --- Input Validation ---
    if (n_total != X.rows() || n_total != y_mask.size() ||
        X.rows() != X_mask.rows() || X.cols() != X_mask.cols()) {
        error("Data and mask dimensions do not match.");
        return false;
    }

    results.clear();
    results.resize(p);

    // --- Main loop to iterate over each column of X ---
    for (int i = 0; i < p; ++i) {
        std::vector<double> y_complete;
        std::vector<double> x_complete;
        y_complete.reserve(n_total);
        x_complete.reserve(n_total);

        // 1. Pairwise deletion: Collect all pairs where mask is true
        for (int j = 0; j < n_total; ++j) {
            // THE CORE CHANGE IS HERE: Check the boolean mask instead of isnan()
            if (y_mask(j) && X_mask(j, i)) {
                y_complete.push_back(y(j));
                x_complete.push_back(X(j, i));
            }
        }

        const int n_complete = y_complete.size();
        slr_sumstat_t& ss = results[i];
        ss.n_obs = n_complete;

        // 2. Check for sufficient data to perform regression
        const int df = n_complete - 2 - n_cov_df;
        if (df <= 0) {
            ss.beta = std::numeric_limits<double>::quiet_NaN();
            ss.se = std::numeric_limits<double>::quiet_NaN();
            ss.tstat = std::numeric_limits<double>::quiet_NaN();
            ss.log10p = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // 3. Map the std::vectors to Eigen vectors for calculation (no copy)
        Eigen::Map<Eigen::VectorXd> y_vec(y_complete.data(), n_complete);
        Eigen::Map<Eigen::VectorXd> x_vec(x_complete.data(), n_complete);

        // 4. Perform regression calculations on the complete data
        const double x_sq_norm = x_vec.squaredNorm();
        const double y_sq_norm = y_vec.squaredNorm();
        const double xt_y = x_vec.dot(y_vec);

        if (x_sq_norm == 0) {
            ss.beta = std::numeric_limits<double>::quiet_NaN();
            ss.se = std::numeric_limits<double>::quiet_NaN();
            ss.tstat = std::numeric_limits<double>::quiet_NaN();
            ss.log10p = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        
        ss.beta = xt_y / x_sq_norm;
        
        const double sse = y_sq_norm - (xt_y * xt_y) / x_sq_norm;
        ss.se = std::sqrt((sse / df) / x_sq_norm);
        
        if (ss.se > 0) {
            ss.tstat = ss.beta / ss.se;
            ss.log10p = tstat2log10pval(ss.tstat, df);
        } else {
            ss.tstat = std::numeric_limits<double>::quiet_NaN();
            ss.log10p = std::numeric_limits<double>::quiet_NaN();
        }
    }

    return true;
}

// Inverse of the standard normal CDF, accurate to about 1e-16 relative error.
// Algorithm AS 241 (PPND16), Wichura M.J. (1988) Applied Statistics 37:477-484,
// the same algorithm R uses for qnorm().
double inverseNormalCDF(double p) {
    if (p <= 0.0 || p >= 1.0) {
        // Return infinity or NaN for out-of-range probabilities
        if (p == 0.0) return -std::numeric_limits<double>::infinity();
        if (p == 1.0) return std::numeric_limits<double>::infinity();
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double q = p - 0.5;
    if (std::fabs(q) <= 0.425) { // central region
        const double r = 0.180625 - q * q;
        return q * (((((((r * 2509.0809287301226727 +
                   33430.575583588128105) * r + 67265.770927008700853) * r +
                 45921.953931549871457) * r + 13731.693765509461125) * r +
               1971.5909503065514427) * r + 133.14166789178437745) * r +
             3.387132872796366608)
          / (((((((r * 5226.495278852545925 +
                   28729.085735721942674) * r + 39307.89580009271061) * r +
                 21213.794301586595867) * r + 5394.1960214247511077) * r +
               687.1870074920579083) * r + 42.313330701600911252) * r + 1.0);
    }

    // tails: r = sqrt(-log(min(p, 1-p)))
    double r = std::sqrt(-std::log(q < 0.0 ? p : 1.0 - p));
    double val;
    if (r <= 5.0) {
        r -= 1.6;
        val = (((((((r * 7.7454501427834140764e-4 +
                   0.0227238449892691845833) * r + 0.24178072517745061177) *
                 r + 1.27045825245236838258) * r +
                3.64784832476320460504) * r + 5.7694972214606914055) *
              r + 4.6303378461565452959) * r +
             1.42343711074968357734)
            / (((((((r *
                     1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + 0.0151986665636164571966) * r +
                   0.14810397642748007459) * r + 0.68976733498510000455) *
                 r + 1.6763848301838038494) * r +
                2.05319162663775882187) * r + 1.0);
    }
    else {
        r -= 5.0;
        val = (((((((r * 2.01033439929228813265e-7 +
                   2.71155556874348757815e-5) * r +
                  0.0012426609473880784386) * r + 0.026532189526576123093) *
                r + 0.29656057182850489123) * r +
               1.7848265399172913358) * r + 5.4637849111641143699) *
             r + 6.6579046435011037772)
            / (((((((r *
                     2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                    r + 1.8463183175100546818e-5) * r +
                   7.868691311456132591e-4) * r + 0.0148753612908506148525)
                 * r + 0.13692988092273580531) * r +
                0.59983220655588793769) * r + 1.0);
    }
    return q < 0.0 ? -val : val;
}

// Blom's rank-based inverse normal score, qnorm((rank - 3/8) / (n + 1/4)),
// matching REGENIE's rint_pheno() (--apply-rint / --apply-rerint).
// `rank` is 1-based; ties should receive the average rank.
static inline double blom_rint_score(double rank, int32_t n) {
    return inverseNormalCDF((rank - 0.375) / ((double)n + 0.25));
}

// Struct to hold value and its original index for sorting purposes
struct ValueIndex {
    double value;
    int32_t original_index;
};


/**
 * @brief Performs rank-based inverse normal transformation (unmasked version).
 *
 * This function assumes all values are valid and non-missing. It is more
 * efficient than the masked version as it avoids the filtering step.
 *
 * @param values An Eigen::VectorXd containing the numerical data.
 * @return An Eigen::VectorXd of the same size with the transformed values.
 */
Eigen::VectorXd rint_without_missing(const Eigen::VectorXd& values) {
    int32_t n = values.size();
    if (n == 0) {
        return Eigen::VectorXd();
    }

    // --- 1. Indexing ---
    // Create a vector of structs to hold value and original index.
    std::vector<ValueIndex> indexed_values(n);
    for (int32_t i = 0; i < n; ++i) {
        indexed_values[i] = {values(i), i};
    }

    // --- 2. Sorting ---
    std::sort(indexed_values.begin(), indexed_values.end(),
                     [](const ValueIndex& a, const ValueIndex& b) {
                         return a.value < b.value;
                     });

    // --- 3. Ranking with Tie Handling (Average Rank) ---
    std::vector<double> ranks(n);
    for (int32_t i = 0; i < n; ) {
        int32_t j = i;
        while (j < n && indexed_values[j].value == indexed_values[i].value) {
            j++;
        }
        double sum_of_ranks = (double)(j - i) / 2.0 * ((i + 1) + j);
        double average_rank = sum_of_ranks / (j - i);
        for (int k = i; k < j; ++k) {
            ranks[k] = average_rank;
        }
        i = j;
    }

    // for(int32_t i=0; i < 5; ++i) {
    //     notice("%.5g\t%.5g\t%d\t%.5g", values(i), indexed_values[i].value, indexed_values[i].original_index, ranks[i]);
    // }


    // --- 4. Inverse Normal Transformation ---
    Eigen::VectorXd result(n);
    for (int i = 0; i < n; ++i) {
        result(indexed_values[i].original_index) = blom_rint_score(ranks[i], n);
    }

    // for(int32_t i=0; i < 5; ++i) {
    //     notice("%.5g -> %.5g", values(i), result(i));
    // }

    return result;
}


/**
 * @brief Performs rank-based inverse normal transformation on a numeric vector (masked version).
 *
 * This function takes a vector of values and a boolean mask. It ranks the
 * non-masked values, handles ties by assigning the average rank, and then
 * applies the inverse normal transformation to these ranks.
 *
 * @param values An Eigen::VectorXd containing the numerical data.
 * @param mask An Eigen::Vector<bool, Eigen::Dynamic> of the same size as `values`.
 * `true` indicates a value to be included in the transformation,
 * `false` indicates a missing value to be ignored.
 * @return An Eigen::VectorXd of the same size as the input. Transformed values
 * are placed in their original positions. Positions corresponding to
 * `false` in the mask are set to NaN.
 */
Eigen::VectorXd rint_with_missing(
    const Eigen::VectorXd& values,
    const Eigen::Vector<bool, Eigen::Dynamic>& mask) {

    // --- 1. Pre-computation and Filtering ---
    if (values.size() != mask.size()) {
        throw std::invalid_argument("Input 'values' and 'mask' vectors must have the same size.");
    }

    std::vector<ValueIndex> filtered_values;
    filtered_values.reserve(values.size()); // Reserve capacity
    for (int32_t i = 0; i < values.size(); ++i) {
        if (mask(i)) {
            filtered_values.push_back({values(i), i});
        }
    }
    
    int32_t n_unmasked = filtered_values.size();
    if (n_unmasked == 0) {
        return Eigen::VectorXd::Constant(values.size(), NAN);
    }

    // --- 2. Sorting ---
    std::sort(filtered_values.begin(), filtered_values.end(),
                     [](const ValueIndex& a, const ValueIndex& b) {
                         return a.value < b.value;
                     });

    // --- 3. Ranking with Tie Handling (Average Rank) ---
    std::vector<double> ranks(n_unmasked);
    for (int32_t i = 0; i < n_unmasked; ) {
        int32_t j = i;
        while (j < n_unmasked && filtered_values[j].value == filtered_values[i].value) {
            j++;
        }
        double sum_of_ranks = (double)(j - i) / 2.0 * ((i + 1) + j);
        double average_rank = sum_of_ranks / (j - i);
        for (int k = i; k < j; ++k) {
            ranks[k] = average_rank;
        }
        i = j;
    }

    // --- 4. Inverse Normal Transformation ---
    Eigen::VectorXd result = Eigen::VectorXd::Constant(values.size(), NAN);
    for (int i = 0; i < n_unmasked; ++i) {
        result(filtered_values[i].original_index) = blom_rint_score(ranks[i], n_unmasked);
    }

    return result;
}

/**
 * @brief Performs column-wise rank-based inverse normal transformation on a matrix (unmasked version).
 *
 * This function assumes all values are valid and non-missing. It performs RINT on each column independently.
 *
 * @param matrix An Eigen::MatrixXd containing the numerical data.
 * @return An Eigen::MatrixXd of the same size with the transformed values.
 */
Eigen::MatrixXd rint_matrix_without_missing(const Eigen::MatrixXd& matrix) {
    int32_t n_rows = matrix.rows();
    int32_t n_cols = matrix.cols();
    
    if (n_rows == 0 || n_cols == 0) {
        return Eigen::MatrixXd(n_rows, n_cols);
    }

    Eigen::MatrixXd result(n_rows, n_cols);
    
    // Reuse the vector for sorting to avoid repeated allocations
    std::vector<ValueIndex> indexed_values(n_rows);

    for (int32_t j = 0; j < n_cols; ++j) {
        // 1. Indexing
        for (int32_t i = 0; i < n_rows; ++i) {
            indexed_values[i].value = matrix(i, j);
            indexed_values[i].original_index = i;
        }

        // 2. Sorting
        std::sort(indexed_values.begin(), indexed_values.end(),
                     [](const ValueIndex& a, const ValueIndex& b) {
                         return a.value < b.value;
                     });

        // 3. Ranking and Transformation
        for (int32_t i = 0; i < n_rows; ) {
            int32_t k = i + 1;
            while (k < n_rows && indexed_values[k].value == indexed_values[i].value) {
                k++;
            }
            
            // Average rank calculation
            // The ranks for the tied group are i+1, i+2, ..., k
            // Average rank = ( (i+1) + k ) / 2.0
            double average_rank = (i + 1 + k) / 2.0;
            double transformed_value = blom_rint_score(average_rank, n_rows);
            
            for (int l = i; l < k; ++l) {
                result(indexed_values[l].original_index, j) = transformed_value;
            }
            i = k;
        }
    }
    return result;
}

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
) 
{
    std::vector<slr_sumstat_t> sumstats;
    std::vector<slr_sumstat_t> sumstats_rint;
    if ( !geno_has_missing && !phe_has_missing ) {
        simple_linear_regression_without_missing(phe_vec, geno_mat, sumstats);
        if ( ! skip_rint ) {
            simple_linear_regression_without_missing(phe_rint_vec, geno_mat, sumstats_rint);
        }
    }
    else {
        simple_linear_regression_with_missing(phe_vec, phe_mask_vec, geno_mat, geno_mask, sumstats);
        if ( ! skip_rint ) {
            simple_linear_regression_with_missing(phe_rint_vec, phe_mask_vec, geno_mat, geno_mask, sumstats_rint);
        }
    }

    // print the results
    for(int32_t i=0; i < (int32_t)v_cpra.size(); ++i) {
        const slr_sumstat_t& ss = sumstats[i];
        hprintf(wf, "%s\t%s\t%d\t%s\t%s\t%s\t%.5g\t%.5g\t%d\t%.6g\t%.6g\t%.6g\t%.6g",
            pheno_id, // TRAIT
            v_cpra[i].chrom.c_str(), // CHROM
            v_cpra[i].pos,           // POS
            v_cpra[i].to_string().c_str(), // ID
            v_cpra[i].ref.c_str(),  // REF
            v_cpra[i].alts.c_str(), // ALT
            (double)acs[i] / (double)ans[i], // AF
            infos[i],  // INFO - placeholder, not calculated
            ss.n_obs, // N - number of samples
            ss.beta, // BETA
            ss.se,   // SE
            ss.tstat, // TSTAT
            ss.log10p); // LOG10P
        if ( ! skip_rint ) {
            const slr_sumstat_t& ss_rint = sumstats_rint[i];
            hprintf(wf, "\t%.6g\t%.6g\t%.6g\t%.6g", // BETA_RINT, SE_RINT, TSTAT_RINT, LOG10P_RINT
                ss_rint.beta, ss_rint.se, ss_rint.tstat, ss_rint.log10p);
        }
        hprintf(wf, "\n");
    }
    return (int32_t)v_cpra.size(); // return the number of variants processed
}

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
    const std::vector<int32_t>& gc0s,     // genotype counts of 0
    const std::vector<int32_t>& gc1s,     // genotype counts of 1
    const std::vector<int32_t>& gc2s,     // genotype counts of
    const std::vector<double>& infos,       // infor values
    bool n_from_geno_counts, // report N as the number of called genotypes instead of the regression sample size
    int32_t n_cov_df,        // covariate degrees of freedom (df = n - 2 - n_cov_df)
    double rint_y_rss        // null-model residual sum of squares of phe_rint_vec (negative = its y'y)
) 
{
    std::vector<slr_sumstat_t> sumstats;
    std::vector<slr_sumstat_t> sumstats_rint;
    if ( !geno_has_missing && !phe_has_missing ) {
        simple_linear_regression_without_missing(phe_vec, geno_mat, sumstats, n_cov_df);
        if ( ! skip_rint ) {
            simple_linear_regression_without_missing(phe_rint_vec, geno_mat, sumstats_rint, n_cov_df, rint_y_rss);
        }
    }
    else {
        simple_linear_regression_with_missing(phe_vec, phe_mask_vec, geno_mat, geno_mask, sumstats, n_cov_df);
        if ( ! skip_rint ) {
            simple_linear_regression_with_missing(phe_rint_vec, phe_mask_vec, geno_mat, geno_mask, sumstats_rint, n_cov_df);
        }
    }

    // print the results
    for(int32_t i=0; i < (int32_t)v_cpra.size(); ++i) {
        const slr_sumstat_t& ss = sumstats[i];
        hprintf(wf, "%s\t%s\t%d\t%s\t%s\t%s\t%.5g\t%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g",
            pheno_id, // TRAIT
            v_cpra[i].chrom.c_str(), // CHROM
            v_cpra[i].pos,           // POS
            v_cpra[i].to_string().c_str(), // ID
            v_cpra[i].ref.c_str(),  // REF
            v_cpra[i].alts.c_str(), // ALT
            (double)acs[i] / (double)ans[i], // AF
            n_from_geno_counts ? ( gc0s[i] + gc1s[i] + gc2s[i] ) : ss.n_obs, // N - number of samples
            gc0s[i], // GC0 - genotype count 0
            gc1s[i], // GC1 - genotype count 1
            gc2s[i], // GC2 - genotype count 2
            infos[i],  // INFO - placeholder, not calculated
            ss.beta, // BETA
            ss.se,   // SE
            ss.tstat, // TSTAT
            ss.log10p); // LOG10P
        if ( ! skip_rint ) {
            const slr_sumstat_t& ss_rint = sumstats_rint[i];
            hprintf(wf, "\t%.6g\t%.6g\t%.6g\t%.6g", // BETA_RINT, SE_RINT, TSTAT_RINT, LOG10P_RINT
                ss_rint.beta, ss_rint.se, ss_rint.tstat, ss_rint.log10p);
        }
        hprintf(wf, "\n");
    }
    return (int32_t)v_cpra.size(); // return the number of variants processed
}


// ---- LOCO (REGENIE step-1 predictions) helpers -------------------------------
// (loco_row_t and the file format are described in assoc_utils.h)
namespace {

bool loco_file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

bool is_loco_missing_str(const char* s) {
    return strcmp(s, "NA") == 0 || strcmp(s, "nan") == 0 || strcmp(s, "NaN") == 0 || strcmp(s, ".") == 0;
}

} // namespace

std::string loco_chrom_label(const std::string& chrom) {
    std::string c = chrom;
    if ( c.size() > 3 && strncasecmp(c.c_str(), "chr", 3) == 0 ) c = c.substr(3);
    if ( c == "X" || c == "x" ) return "23";
    if ( c == "Y" || c == "y" ) return "24";
    if ( c == "XY" || c == "xy" ) return "25";
    if ( c == "M" || c == "MT" || c == "m" || c == "mt" ) return "26";
    return c;
}

// Determine the .loco file for each trait. `locof` is either a single .loco file
// (only valid when exactly one trait is tested) or a REGENIE *_pred.list file
// with lines "TRAIT PATH". Relative paths in the list that do not exist as given
// are resolved against the directory of the list file.
// kind: ind_assoc_input::LOCO_ANY (detect), LOCO_SINGLE (.loco required) or
// LOCO_PRED_LIST (*_pred.list required).
std::vector<std::string> resolve_loco_files(const std::string& locof, const std::vector<std::string>& pheno_ids, int32_t kind) {
    std::vector<std::string> paths;
    tsv_reader tr(locof.c_str());
    if ( !tr.read_line() ) {
        error("LOCO file %s is empty", locof.c_str());
    }
    if ( strcmp(tr.str_field_at(0), "FID_IID") == 0 ) { // a single .loco file
        if ( kind == ind_assoc_input::LOCO_PRED_LIST ) {
            error("--pred %s is a single .loco file (header FID_IID), not a *_pred.list file. Use --loco for a single .loco file",
                  locof.c_str());
        }
        if ( pheno_ids.size() != 1 ) {
            error("--loco %s is a single .loco file, but %zu traits are tested. Provide a REGENIE *_pred.list file (TRAIT PATH per line) with --pred instead",
                  locof.c_str(), pheno_ids.size());
        }
        paths.push_back(locof);
        return paths;
    }
    if ( kind == ind_assoc_input::LOCO_SINGLE ) {
        error("--loco %s is not a .loco file (its header does not start with FID_IID). Use --pred for a *_pred.list file", locof.c_str());
    }

    // *_pred.list file
    std::string dir;
    size_t slash = locof.rfind('/');
    if ( slash != std::string::npos ) dir = locof.substr(0, slash + 1);
    std::map<std::string, std::string> trait2path;
    do {
        if ( tr.nfields == 0 ) continue;
        if ( tr.nfields < 2 ) {
            error("Invalid line in LOCO list file %s: expected 'TRAIT PATH'", locof.c_str());
        }
        std::string path(tr.str_field_at(1));
        if ( !loco_file_exists(path) && path[0] != '/' && !dir.empty() && loco_file_exists(dir + path) ) {
            path = dir + path;
        }
        trait2path[tr.str_field_at(0)] = path;
    } while ( tr.read_line() );

    for(size_t k = 0; k < pheno_ids.size(); ++k) {
        std::map<std::string, std::string>::const_iterator it = trait2path.find(pheno_ids[k]);
        if ( it == trait2path.end() ) {
            error("Trait %s is not listed in the LOCO list file %s", pheno_ids[k].c_str(), locof.c_str());
        }
        paths.push_back(it->second);
    }
    return paths;
}

void load_loco_row(const std::string& path, const std::string& chrom, loco_row_t& row) {
    const std::string target = loco_chrom_label(chrom);
    row = loco_row_t();
    tsv_reader tr(path.c_str());
    if ( !tr.read_line() || strcmp(tr.str_field_at(0), "FID_IID") != 0 ) {
        error("LOCO file %s does not start with a 'FID_IID' header", path.c_str());
    }
    const int32_t n = tr.nfields - 1;
    std::vector<std::string> ids(n);
    for(int32_t j = 0; j < n; ++j) {
        ids[j] = tr.str_field_at(j + 1);
        row.exact[ids[j]] = j;
    }
    for(int32_t j = 0; j < n; ++j) {
        const std::string& s = ids[j];
        for(size_t p = s.find('_'); p != std::string::npos; p = s.find('_', p + 1)) {
            std::string suf = s.substr(p + 1);
            if ( suf.empty() ) continue;
            std::unordered_map<std::string, int32_t>::iterator it = row.suffix.find(suf);
            if ( it == row.suffix.end() ) row.suffix[suf] = j;
            else if ( it->second != j ) it->second = -1; // ambiguous
        }
    }
    if ( chrom.empty() ) return; // sample IDs only
    while ( tr.read_line() ) {
        if ( tr.nfields == 0 ) continue;
        if ( loco_chrom_label(tr.str_field_at(0)) != target ) continue;
        if ( tr.nfields != n + 1 ) {
            error("LOCO file %s: chromosome %s row has %d values, but the header has %d samples",
                  path.c_str(), tr.str_field_at(0), tr.nfields - 1, n);
        }
        row.values.resize(n);
        row.observed.resize(n);
        for(int32_t j = 0; j < n; ++j) {
            const char* s = tr.str_field_at(j + 1);
            row.observed[j] = !is_loco_missing_str(s);
            row.values[j] = row.observed[j] ? tr.double_field_at(j + 1) : std::numeric_limits<double>::quiet_NaN();
        }
        return;
    }
    error("LOCO file %s has no row for chromosome %s", path.c_str(), chrom.c_str());
}


bool ind_assoc_input::load_pheno_cov_matrices(const char* phef, const char* pheno_format, const char* covf, const char* cov_format) {
    // if the non-necessary arguments are empty, set them to NULL
    if ( covf != NULL && strlen(covf) == 0 ) covf = NULL;
    if ( phef != NULL && strlen(phef) == 0 ) {
        error("Phenotype file is required for association analysis");
    }

    // get the list of all genotyped samples
    const std::vector<plink_samp_t>& geno_all_samps = mpr.get_all_samples();
    std::vector<std::string> geno_all_samp_ids;
    for(int32_t i=0; i < geno_all_samps.size(); ++i) {
        geno_all_samp_ids.push_back(geno_all_samps[i].indID);
    }

    notice("Loading phenotype matrix from %s", phef);
    if ( !pheno_matrix.load_pheno_matrix(phef, pheno_format, table_delim) ) {
        error("Failed to load the phenotype matrix from file %s", phef);
    }
    if ( subset_pheno_ids.size() > 0 ) {
        if ( subset_pheno_drop_absent ) {
            // keep the listed traits found in the phenotype file, in sorted order
            std::vector<std::string> found_ids;
            identify_overlapping_ids(subset_pheno_ids, pheno_matrix.pheno_ids, found_ids);
            notice("%zu of %zu listed phenotypes are found in the phenotype matrix of %zu phenotypes",
                   found_ids.size(), subset_pheno_ids.size(), pheno_matrix.pheno_ids.size());
            if ( found_ids.empty() ) {
                error("None of the listed phenotypes are in the phenotype file %s", phef);
            }
            subset_pheno_ids.swap(found_ids);
        }
        notice("Subsetting the phenotype matrix to the specified phenotypes");
        pheno_matrix.subset_pheno_ids(subset_pheno_ids);
    }
    
    notice("Loaded phenotype matrix with %d samples and %d phenotypes from %s", (int32_t)pheno_matrix.samp_ids.size(), (int32_t)pheno_matrix.pheno_ids.size(), phef);

    if ( covf != NULL && strlen(covf) > 0 ) {
        notice("Loading covariate matrix from %s", covf);
        cov_matrix.set_categorical_ids(cat_cov_ids);
        cov_matrix.set_selected_columns(cov_col_ids);
        if ( !cov_matrix.load_pheno_matrix(covf, cov_format, table_delim) ) {
            error("Failed to load the covariate matrix from file %s", covf);
        }
        notice("Loaded covariate matrix with %d samples and %d covariates from %s", (int32_t)cov_matrix.samp_ids.size(), (int32_t)cov_matrix.pheno_ids.size(), covf);
        cov_matrix.keep_covariate_columns(cov_col_ids, covf);
    }

    // load the covariate matrix
    // NB: populate the class member (not a local) so load_genotype_chunk can use it later
    overlapping_sample_ids.clear();
    if ( !subset_sample_ids.empty() ) {
        std::vector<std::string> temp1_ids;
        std::vector<std::string> temp2_ids;
        identify_overlapping_ids(subset_sample_ids, geno_all_samp_ids, temp1_ids);
        if ( covf == NULL || std::string(covf).length() == 0 ) {
            identify_overlapping_ids(temp1_ids, pheno_matrix.samp_ids, overlapping_sample_ids);
        }
        else {
            identify_overlapping_ids(temp1_ids, cov_matrix.samp_ids, temp2_ids);
            identify_overlapping_ids(temp2_ids, pheno_matrix.samp_ids, overlapping_sample_ids);
        }
    }
    else {
        if ( covf == NULL )  {
            identify_overlapping_ids(geno_all_samp_ids, pheno_matrix.samp_ids, overlapping_sample_ids);
        }
        else {
            std::vector<std::string> temp1_ids;
            identify_overlapping_ids(geno_all_samp_ids, cov_matrix.samp_ids, temp1_ids);
            identify_overlapping_ids(temp1_ids, pheno_matrix.samp_ids, overlapping_sample_ids);
        }
    }
    notice("%zu overlapping samples found among sample, genotype, phenotype, and covariate files", (int32_t)overlapping_sample_ids.size());

    // Single-trait missing values: analyze only the samples with an observed
    // phenotype. Removing them here, before any subsetting, means the phenotype,
    // covariate and genotype matrices (and the covariate residualization, RINT and
    // LOCO scaling) are all built on the observed samples, so nothing downstream
    // sees a missing value.
    if ( drop_missing_pheno_samples && pheno_matrix.has_missing ) {
        const int32_t n_traits = (int32_t)pheno_matrix.pheno_ids.size();
        if ( n_traits > 1 ) {
            error("The phenotype matrix has missing values across %d traits. Missing phenotype values are allowed only when a single trait is tested (all traits must share the same samples); run each trait separately", n_traits);
        }
        std::set<std::string> missing_ids;
        for(int32_t i = 0; i < (int32_t)pheno_matrix.samp_ids.size(); ++i) {
            if ( !pheno_matrix.pheno_mask(i, 0) ) missing_ids.insert(pheno_matrix.samp_ids[i]);
        }
        std::vector<std::string> kept_ids;
        kept_ids.reserve(overlapping_sample_ids.size());
        for(const std::string& id : overlapping_sample_ids) {
            if ( missing_ids.find(id) == missing_ids.end() ) kept_ids.push_back(id);
        }
        notice("Trait %s: %zu of %zu overlapping samples have a missing phenotype and are excluded; analyzing %zu samples",
               pheno_matrix.pheno_ids[0].c_str(), overlapping_sample_ids.size() - kept_ids.size(),
               overlapping_sample_ids.size(), kept_ids.size());
        if ( kept_ids.empty() ) {
            error("No overlapping samples have an observed value for trait %s", pheno_matrix.pheno_ids[0].c_str());
        }
        overlapping_sample_ids.swap(kept_ids);
    }

    // identify_overlapping_ids() returns IDs in lexicographic order, but the pgen
    // reader requires the sample subset in increasing genotype-file order. Put the
    // overlapping IDs in genotype order so that phenotype, covariate, and genotype
    // rows all follow the same (genotype-file) order.
    {
        std::map<std::string, int32_t> geno_rank;
        for(int32_t i = 0; i < (int32_t)geno_all_samp_ids.size(); ++i) geno_rank[geno_all_samp_ids[i]] = i;
        std::sort(overlapping_sample_ids.begin(), overlapping_sample_ids.end(),
                  [&geno_rank](const std::string& a, const std::string& b) { return geno_rank[a] < geno_rank[b]; });
    }

    // LOCO: keep only samples that have a LOCO prediction for every tested trait. With a fixed
    // chromosome (region-assoc), its predictions are loaded now and a missing value also excludes
    // the sample. In per-chromosome mode only the sample IDs are read here; a missing value in a
    // chromosome's row is an error when that chromosome is applied.
    loco_rows.clear();
    loco_rows_chrom.clear();
    loco_applied_chrom.clear();
    loco_chroms_seen.clear();
    if ( use_loco() ) {
        loco_paths = resolve_loco_files(loco_file, pheno_matrix.pheno_ids, loco_kind);
        loco_rows.resize(loco_paths.size());
        for(size_t k = 0; k < loco_paths.size(); ++k) {
            if ( !loco_chrom.empty() ) {
                notice("Loading LOCO predictions for trait %s on chromosome %s from %s",
                       pheno_matrix.pheno_ids[k].c_str(), loco_chrom.c_str(), loco_paths[k].c_str());
            }
            load_loco_row(loco_paths[k], loco_chrom, loco_rows[k]);
        }
        loco_rows_chrom = loco_chrom;
        std::vector<std::string> kept_ids;
        kept_ids.reserve(overlapping_sample_ids.size());
        for(size_t i = 0; i < overlapping_sample_ids.size(); ++i) {
            bool ok = true;
            double v;
            for(size_t k = 0; ok && k < loco_rows.size(); ++k) {
                ok = loco_chrom.empty() ? loco_rows[k].has(overlapping_sample_ids[i])
                                        : loco_rows[k].get(overlapping_sample_ids[i], v);
            }
            if ( ok ) kept_ids.push_back(overlapping_sample_ids[i]);
        }
        if ( kept_ids.size() < overlapping_sample_ids.size() ) {
            warning("%zu of %zu overlapping samples have no LOCO prediction for at least one trait and are excluded",
                    overlapping_sample_ids.size() - kept_ids.size(), overlapping_sample_ids.size());
        }
        if ( kept_ids.empty() ) {
            error("No overlapping samples have LOCO predictions. Check that the .loco FID_IID IDs match the sample IDs");
        }
        overlapping_sample_ids.swap(kept_ids);
        notice("%zu samples retained after matching LOCO predictions", overlapping_sample_ids.size());
    }

    if ( pheno_matrix.samp_ids != overlapping_sample_ids ) {
        notice("Subsetting the phenotype matrix to the overlapping samples");
        pheno_matrix.subset_sample_ids(overlapping_sample_ids);
    }
    else {
        notice("No need to subset the phenotype matrix, as the sample IDs already match the overlapping samples");
    }
    if ( covf != NULL ) {
        if ( cov_matrix.samp_ids != overlapping_sample_ids ) {
            notice("Subsetting the covariate matrix to the overlapping samples");
            cov_matrix.subset_sample_ids(overlapping_sample_ids);
        }
        else {
            notice("No need to subset the covariate matrix, as the sample IDs already match the overlapping samples");
        }
        // categorical covariates become indicator columns, with levels counted on the analyzed samples
        if ( cov_matrix.has_categorical() ) {
            cov_matrix.expand_categorical(max_cat_levels);
            notice("Covariate matrix has %zu columns after expanding categorical covariates", cov_matrix.pheno_ids.size());
        }
    }
    if ( mpr.get_all_sample_count() != (int32_t)overlapping_sample_ids.size() ) { // same size => same order, as overlapping_sample_ids follows genotype order
        notice("Subsetting the genotype data to %zu overlapping samples", (int32_t)overlapping_sample_ids.size());
        mpr.subset_sample_ids(overlapping_sample_ids);
    }
    else {
        notice("No need to subset the genotype data, as the sample IDs already match the overlapping samples"); 
    }
    if ( mpr.get_loaded_sample_count() != (int32_t)overlapping_sample_ids.size() ) {
        error("%d genotype samples are loaded, but %zu samples were matched across the input files", mpr.get_loaded_sample_count(), overlapping_sample_ids.size());
    }
    int32_t n_overlapping_samples = (int32_t)overlapping_sample_ids.size();

    // adjust phenotype matrix by covariates
    if ( rint_before_adj ) {
        notice("Performing rank-based inverse normal transformation for all phenotypes before covariate adjustment");
        if ( !pheno_matrix.has_missing ) {
            pheno_matrix.pheno_mat = rint_matrix_without_missing(pheno_matrix.pheno_mat);
        }
        else {
            error("Rank-based inverse normal transformation adjustment is currently only supported for phenotype matrices with missing values");
        }
    }

    // all traits share the analyzed samples here (region-assoc drops the missing samples of a
    // single trait above), so a missing phenotype value cannot be handled
    if ( pheno_matrix.has_missing ) {
        error("The phenotype matrix has missing values. This command tests all traits on the same samples and needs complete phenotypes (pair-assoc handles missing values per trait)");
    }
    if ( covf != NULL && cov_matrix.has_missing ) {
        error("The covariate matrix has missing values, which this command does not support. Remove the samples with missing covariates (pair-assoc excludes them)");
    }

    // perform covariate adjustment
    n_cov_df = 0;
    cov_resid.reset();
    pheno_transformed_after_cov = false;
    if ( covf != NULL ) {
        notice("Adjusting phenotypes by covariates using linear regression");
        if ( pheno_matrix.has_missing || cov_matrix.has_missing ) {
            error("Covariate adjustment is currently only supported for phenotype and covariate matrices without missing values");
        }
        else {
            std::shared_ptr<CovariateResidualizer> resid(new CovariateResidualizer(cov_matrix.pheno_mat));
            resid->residualize(pheno_matrix.pheno_mat);
            cov_resid = resid;
            // the Wald tests lose one degree of freedom per (linearly independent) covariate
            n_cov_df = cov_resid->n_cov_df();
            const int32_t n_cov_cols = (int32_t)cov_matrix.pheno_mat.cols();
            if ( n_cov_df < n_cov_cols ) {
                warning("The intercept and %d covariates have numerical rank %d among the %d analyzed samples: some covariates are constant or linear combinations of others. They are adjusted for once, and the degrees of freedom use the rank; drop the redundant covariates.",
                        n_cov_cols, cov_resid->rank(), n_overlapping_samples);
            }
            notice("Covariate degrees of freedom: %d (rank of the intercept and %d covariates = %d); Wald tests use df = N - 2 - %d",
                   n_cov_df, n_cov_cols, cov_resid->rank(), n_cov_df);
        }
    }
    else if ( !pheno_matrix.has_missing ) {
        // Without --cov the model still has an intercept: center the phenotypes. The tests assume
        // centered phenotypes (y'y is their residual sum of squares), so an uncentered phenotype
        // would inflate SE.
        Eigen::RowVectorXd means = pheno_matrix.pheno_mat.colwise().mean();
        pheno_matrix.pheno_mat.rowwise() -= means;
    }

    // LOCO adjustment, following REGENIE step 2 (Pheno.cpp residualize_phenotypes,
    // Data.cpp compute_res):
    //   (2) scale the covariate residual to unit SD, SD = ||r|| / sqrt(n - n_cov)
    //       where n_cov counts the intercept (REGENIE's scale_Y);
    //   (3) subtract the LOCO prediction of the tested chromosome.
    // Unless --rint-after-adj follows, the result is multiplied back by scale_Y so the
    // phenotype, and hence BETA/SE, the exported sufficient/RSS statistics, and the
    // SuSiE effect sizes, stay in the units of the covariate-adjusted phenotype
    // (REGENIE reports BETA as stat * scale_Y). T-statistics and p-values do not
    // depend on this rescaling. With --rint-after-adj, the RINT output defines the
    // scale, as with REGENIE's --apply-rerint.
    if ( use_loco() ) {
        if ( pheno_matrix.has_missing ) {
            error("--loco is currently only supported for phenotype matrices without missing values");
        }
        const int32_t n = (int32_t)pheno_matrix.pheno_mat.rows();
        const int32_t K = (int32_t)pheno_matrix.pheno_mat.cols();
        // (without covariates, the phenotypes were centered above)
        const int32_t n_cov = 1 + n_cov_df; // rank of [1, C]
        if ( n - n_cov <= 0 ) {
            error("Not enough samples (%d) for %d covariates (including the intercept) to scale phenotypes for LOCO", n, n_cov);
        }
        pheno_scale = pheno_matrix.pheno_mat.colwise().norm() / std::sqrt((double)(n - n_cov));
        for(int32_t k = 0; k < K; ++k) {
            if ( !(pheno_scale(k) > 1e-12) ) {
                error("Covariate-adjusted phenotype %s has zero variance; cannot apply LOCO", pheno_matrix.pheno_ids[k].c_str());
            }
        }
        // standardized covariate residual; apply_loco_chrom() subtracts each chromosome's
        // predictions from it
        pheno_loco_base.resize(n, K);
        for(int32_t k = 0; k < K; ++k) {
            pheno_loco_base.col(k) = pheno_matrix.pheno_mat.col(k) / pheno_scale(k);
        }
        if ( !loco_chrom.empty() ) {
            apply_loco_chrom(loco_chrom);
        }
        else {
            notice("LOCO predictions are subtracted per chromosome, before each chromosome's tests");
        }
    }
    else if ( rint_after_adj ) {
        notice("Performing rank-based inverse normal transformation for all phenotypes after covariate adjustment");
        if ( !pheno_matrix.has_missing ) {
            pheno_matrix.pheno_mat = rint_matrix_without_missing(pheno_matrix.pheno_mat);
        }
        else {
            error("Rank-based inverse normal transformation adjustment is currently only supported for phenotype matrices with missing values");
        }
    }

    // After RINT or LOCO, the phenotype is no longer orthogonal to the covariates; the Wald
    // tests then use its residual sum of squares on [1, C] (pheno_cov_rss()). The phenotype
    // itself (used by --score, SuSiE and the exports) is left unchanged.
    // Without --cov this is the intercept-only fit: subtracting LOCO predictions can leave the
    // phenotype uncentered.
    pheno_transformed_after_cov = ( cov_resid && rint_after_adj ) || use_loco();
    return true;
}

void ind_assoc_input::apply_loco_chrom(const std::string& chrom) {
    if ( !use_loco() ) error("apply_loco_chrom() called without --loco/--pred");
    if ( pheno_loco_base.size() == 0 && pheno_matrix.pheno_ids.size() > 0 ) error("apply_loco_chrom() called before the phenotypes are loaded");
    const std::string label = loco_chrom_label(chrom);
    if ( !loco_applied_chrom.empty() && loco_chrom_label(loco_applied_chrom) == label ) return;
    if ( loco_chroms_seen.count(label) ) {
        notice("Chromosome %s is tested again after another chromosome: its LOCO predictions are read again (group the variants by chromosome to avoid this)", chrom.c_str());
    }
    const int32_t n = (int32_t)pheno_loco_base.rows();
    const int32_t K = (int32_t)pheno_loco_base.cols();
    if ( loco_rows_chrom.empty() || loco_chrom_label(loco_rows_chrom) != label ) {
        for(int32_t k = 0; k < K; ++k) load_loco_row(loco_paths[k], chrom, loco_rows[k]);
        loco_rows_chrom = chrom;
    }

    // LOCO predictions aligned with the (subsetted) sample and trait order
    Eigen::MatrixXd loco_mat(n, K);
    for(int32_t i = 0; i < n; ++i) {
        for(int32_t k = 0; k < K; ++k) {
            double v;
            if ( !loco_rows[k].get(pheno_matrix.samp_ids[i], v) ) {
                error("Missing LOCO prediction for sample %s and trait %s on chromosome %s",
                      pheno_matrix.samp_ids[i].c_str(), pheno_matrix.pheno_ids[k].c_str(), chrom.c_str());
            }
            loco_mat(i, k) = v;
        }
    }

    pheno_matrix.pheno_mat.resize(n, K);
    for(int32_t k = 0; k < K; ++k) {
        const double sd = pheno_scale(k);
        pheno_matrix.pheno_mat.col(k) = pheno_loco_base.col(k) - loco_mat.col(k);
        if ( !rint_after_adj ) pheno_matrix.pheno_mat.col(k) *= sd;
        if ( K <= 20 ) {
            notice("LOCO-adjusted trait %s with chromosome %s predictions (residual SD before LOCO = %.6g%s)",
                   pheno_matrix.pheno_ids[k].c_str(), chrom.c_str(), sd,
                   rint_after_adj ? "; RINT follows, so the output is on the RINT scale" : "");
        }
    }
    if ( K > 20 ) {
        notice("LOCO-adjusted %d traits with chromosome %s predictions%s", K, chrom.c_str(),
               rint_after_adj ? "; RINT follows, so the output is on the RINT scale" : "");
    }
    if ( rint_after_adj ) {
        notice("Performing rank-based inverse normal transformation for all phenotypes after the LOCO adjustment");
        pheno_matrix.pheno_mat = rint_matrix_without_missing(pheno_matrix.pheno_mat);
    }
    loco_applied_chrom = chrom;
    loco_chroms_seen.insert(label);
}

bool ind_assoc_input::process_pgenlist(const char* pgenlistf, const char* phef, const char* pheno_format, const char* covf, const char* cov_format) {
    if ( !load_pgenlist(pgenlistf) ) {
        return false;
    }
    mpr.set_icol_pivar_idx(icol_pivar_idx);
    if ( !load_pheno_cov_matrices(phef, pheno_format, covf, cov_format) ) {
        return false;
    }
    return true;
}

bool ind_assoc_input::process_single_pgen(const char* pgenf, const char* pivarf, const char* psamf, const char* phef, const char* pheno_format, const char* covf, const char* cov_format) {
    if ( !load_pgen_files(pgenf, pivarf, psamf) ) {
        return false;
    }
    mpr.set_icol_pivar_idx(icol_pivar_idx);
    if ( !load_pheno_cov_matrices(phef, pheno_format, covf, cov_format) ) {
        return false;
    }
    return true;
}

bool ind_assoc_input::load_pgenlist(const char* pgenlistf) {
    if ( !mpr.prep_pgen_list(pgenlistf) ) {
        error("Failed to prepare pgen files with the following list file: %s", pgenlistf);
        return false;
    }
    return true;
}

bool ind_assoc_input::load_pgen_files(const char* pgenf, const char* pivarf, const char* psamf) {
    if ( !mpr.set_single_chunk_pgen(pgenf, pivarf, psamf) ) {
        error("Failed to add pgen files with the following files:\n%s\n%s\n%s", pgenf, pivarf, psamf);
        return false;
    }
    return true;
}

// bool ind_assoc_input::add_genotypes() {
//     const plink_var_t& var = mpr.get_current_variant();

//     gcs[0] = gcs[1] = gcs[2] = 0; // initialize genotype counts
//     std::string cpra_s(var.to_string());
//     mpr.get_genos();

//     if ( icol >= n_col_est ) {
//         geno_mat.conservativeResize(n_overlapping_samples, n_col_est * 2);
//         geno_mask.conservativeResize(n_overlapping_samples, n_col_est * 2);
//         n_col_est *= 2;
//     }
//     int32_t an = 0;
//     double ac = 0;
//     bool skip  = false;
//     if ( mpr.is_dosage_present() ) {
//         if ( dbl_buf == NULL) {
//             dbl_buf = mpr.get_dbl_buf();
//         }
//         double sumsq = 0;
//         for(int32_t i =0; i < n_overlapping_samples; ++i) {
//             double ds = 2.0 - dbl_buf[i];
//             geno_mat(i, icol) = ds;
//             geno_mask(i, icol) = true; // not missing
//             ac += ds;
//             an += 2;
//             ++gcs[(ds < 0.5) ? 0 : ( (ds < 1.5) ? 1 : 2 )];
//             sumsq += (ds * ds);
//         }
//         if ( ac == 0 || an == ac ) {
//             skip = true; // skip monomorphic variants
//         }
//         else if ( ac < min_ac || ac > max_ac || ac < min_af * an || ac > max_af * an ) {
//             skip = true; // skip variants outside the specified allele count/frequency range
//         }
//         else {
//             double info = ( sumsq / 2.0 * an - 4.0 * ac * ac ) / ( 2.0 * ac * ( an - ac ) );
//             infos.push_back(info);
//         }
//     }
//     else {
//         for(int32_t i = 0; i < n_overlapping_samples; ++i) {
//             switch(int_buf[i]) { // make sure to convert 1-based index to 0-based
//             case 0:
//                 an += 2;
//                 ac += 2;
//                 ++gcs[2];
//                 break;
//             case 1:
//                 an += 2;
//                 ++ac;
//                 ++gcs[1];
//                 break;
//             case 2:
//                 an += 2;
//                 ++gcs[0];
//                 break;
//             }
//         }
//         if ( ac == 0 || an == ac ) {
//             skip = true; // skip monomorphic variants
//         }
//         else if ( ac < min_ac || ac > max_ac || ac < min_af * an || ac > max_af * an ) {
//             skip = true; // skip variants outside the specified allele count/frequency range
//         }
//         else {
//             double mean = (double)ac / (double)an * 2.0;
//             double info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
//             infos.push_back(info);
//             for(int32_t i = 0; i < n_overlapping_samples; ++i) {
//                 switch(int_buf[i]) {
//                 case 0:
//                     geno_mat(i, icol) = 2.0 - mean; // homalt
//                     geno_mask(i, icol) = true; // not missing
//                     break;
//                 case 1:
//                     geno_mat(i, icol) = 1.0 - mean; // het
//                     geno_mask(i, icol) = true; // not missing
//                     break;
//                 case 2:
//                     geno_mat(i, icol) = 0.0 - mean; // homref
//                     geno_mask(i, icol) = true; // not missing
//                     break;
//                 default:
//                     geno_mat(i, icol) = 0; // missing - mean imputation
//                     geno_mask(i, icol) = false; // missing
//                     ++n_geno_missing;
//                     break;
//                 }
//             }
//         }
//     }

//     if ( !skip ) {
//         v_cpra.push_back(cpra_t(cpra_s.c_str()));
//         acs.push_back(ac);
//         ans.push_back(an);
//         gc0s.push_back(gcs[0]);
//         gc1s.push_back(gcs[1]);
//         gc2s.push_back(gcs[2]);
//         ++icol;
//     }
//     else {
//         ++nskip;
//     }    
// }

// int32_t ind_assoc_input::load_genotype_chunk(const char* chrom, int32_t beg, int32_t end, int32_t max_chunk_vars, Eigen::MatrixXd& geno_mat, Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& geno_mask) {
//     // read genotype data from the pgen files
//     notice("Loading genotype data from pgen files with chromosome %s, position %d to %d, and maximum chunk size of %d variants", chrom == NULL ? "NULL" : chrom, beg, end, max_chunk_vars);
//     int32_t n_overlapping_samples = (int32_t)overlapping_sample_ids.size();
//     if ( chrom != NULL ) { // start reading from a specific position
//         geno_chunk_done = false;
//         if ( !mpr.read_pos(chrom, beg) ) { // no variant at/after the requested start
//             geno_chunk_done = true;
//         }
//         notice("Reading genotype data from %s:%d to %s:%d", chrom, beg, chrom, end);
//     }
//     else { // continue streaming from the current position left by the previous chunk
//         notice("Reading genotype data from the current position to %s:%d", geno_chunk_done ? "(end)" : mpr.get_current_variant().schrom.c_str(), end);
//     }

//     if ( geno_chunk_done ) { // the region has already been fully streamed
//         geno_mat.resize(n_overlapping_samples, 0);
//         geno_mask.resize(n_overlapping_samples, 0);
//         return false;
//     }

//     std::vector<cpra_t> v_cpra;
//     std::vector<int32_t> ans;
//     std::vector<double> acs;
//     std::vector<int32_t> gc0s, gc1s, gc2s;
//     std::vector<double> infos;
//     const std::vector<int32_t>& int_buf = mpr.get_int_buf();
//     const double* dbl_buf = mpr.get_dbl_buf();
//     int32_t n_col_est = 10;
//     geno_mat.resize(n_overlapping_samples, n_col_est);
//     geno_mask.resize(n_overlapping_samples, n_col_est);
//     uint32_t n_geno_missing = 0;
//     int32_t icol = 0;
//     int32_t nskip = 0;

//     notice("Reading genotype data for %d overlapping samples", n_overlapping_samples);

//     int32_t gcs[3] = {0, 0, 0};
//     while ( true ) {
//         const plink_var_t& var = mpr.get_current_variant();
//         if ( var.pos > end ) { // moved past the requested region
//             geno_chunk_done = true;
//             break;
//         }

//         gcs[0] = gcs[1] = gcs[2] = 0; // initialize genotype counts
//         std::string cpra_s(var.to_string());
//         mpr.get_genos();

//         if ( icol >= n_col_est ) {
//             geno_mat.conservativeResize(n_overlapping_samples, n_col_est * 2);
//             geno_mask.conservativeResize(n_overlapping_samples, n_col_est * 2);
//             n_col_est *= 2;
//         }
//         int32_t an = 0;
//         double ac = 0;
//         bool skip  = false;
//         if ( mpr.is_dosage_present() ) {
//             if ( dbl_buf == NULL) {
//                 dbl_buf = mpr.get_dbl_buf();
//             }
//             double sumsq = 0;
//             for(int32_t i =0; i < n_overlapping_samples; ++i) {
//                 double ds = 2.0 - dbl_buf[i];
//                 geno_mat(i, icol) = ds;
//                 geno_mask(i, icol) = true; // not missing
//                 ac += ds;
//                 an += 2;
//                 ++gcs[(ds < 0.5) ? 0 : ( (ds < 1.5) ? 1 : 2 )];
//                 sumsq += (ds * ds);
//             }
//             if ( ac == 0 || an == ac ) {
//                 skip = true; // skip monomorphic variants
//             }
//             else if ( ac < min_ac || ac > max_ac || ac < min_af * an || ac > max_af * an ) {
//                 skip = true; // skip variants outside the specified allele count/frequency range
//             }
//             else {
//                 double info = ( sumsq / 2.0 * an - 4.0 * ac * ac ) / ( 2.0 * ac * ( an - ac ) );
//                 infos.push_back(info);
//             }
//         }
//         else {
//             for(int32_t i = 0; i < n_overlapping_samples; ++i) {
//                 switch(int_buf[i]) { // make sure to convert 1-based index to 0-based
//                 case 0:
//                     an += 2;
//                     ac += 2;
//                     ++gcs[2];
//                     break;
//                 case 1:
//                     an += 2;
//                     ++ac;
//                     ++gcs[1];
//                     break;
//                 case 2:
//                     an += 2;
//                     ++gcs[0];
//                     break;
//                 }
//             }
//             if ( ac == 0 || an == ac ) {
//                 skip = true; // skip monomorphic variants
//             }
//             else if ( ac < min_ac || ac > max_ac || ac < min_af * an || ac > max_af * an ) {
//                 skip = true; // skip variants outside the specified allele count/frequency range
//             }
//             else {
//                 double mean = (double)ac / (double)an * 2.0;
//                 double info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
//                 infos.push_back(info);
//                 for(int32_t i = 0; i < n_overlapping_samples; ++i) {
//                     switch(int_buf[i]) {
//                     case 0:
//                         geno_mat(i, icol) = 2.0 - mean; // homalt
//                         geno_mask(i, icol) = true; // not missing
//                         break;
//                     case 1:
//                         geno_mat(i, icol) = 1.0 - mean; // het
//                         geno_mask(i, icol) = true; // not missing
//                         break;
//                     case 2:
//                         geno_mat(i, icol) = 0.0 - mean; // homref
//                         geno_mask(i, icol) = true; // not missing
//                         break;
//                     default:
//                         geno_mat(i, icol) = 0; // missing - mean imputation
//                         geno_mask(i, icol) = false; // missing
//                         ++n_geno_missing;
//                         break;
//                     }
//                 }
//             }
//         }

//         if ( !skip ) {
//             v_cpra.push_back(cpra_t(cpra_s.c_str()));
//             acs.push_back(ac);
//             ans.push_back(an);
//             gc0s.push_back(gcs[0]);
//             gc1s.push_back(gcs[1]);
//             gc2s.push_back(gcs[2]);
//             ++icol;
//         }
//         else {
//             ++nskip;
//         }

//         // advance to the next variant right away so the next chunk never re-reads this one
//         if ( !mpr.read_pivar() ) {
//             geno_chunk_done = true; // reached the end of available variants
//             break;
//         }

//         if ( icol >= max_chunk_vars ) { // chunk is full; cur_var already points to the next variant
//             break;
//         }
//     }

//     // shrink to the number of variants actually loaded so callers never see garbage columns
//     geno_mat.conservativeResize(n_overlapping_samples, icol);
//     geno_mask.conservativeResize(n_overlapping_samples, icol);

//     notice("Loaded %d variants in the current chunk until %s, skipped %d", icol, mpr.get_current_variant().to_string().c_str(), nskip);

//     return icol;
// }

int32_t ind_assoc_input::load_genotype_chunk(const char* chrom, int32_t beg, int32_t end, int32_t max_chunk_vars) {
    // read genotype data from the pgen files
    notice("Loading genotype data from pgen files with chromosome %s, position %d to %d, and maximum chunk size of %d variants", chrom == NULL ? "NULL" : chrom, beg, end, max_chunk_vars);
    int32_t n_overlapping_samples = (int32_t)overlapping_sample_ids.size();
    if ( chrom != NULL ) { // start reading from a specific position
        geno_chunk_done = false;
        if ( !mpr.read_pos(chrom, beg) ) { // no variant at/after the requested start
            geno_chunk_done = true;
        }
        notice("Reading genotype data from %s:%d to %s:%d", chrom, beg, chrom, end);
    }
    else { // continue streaming from the current position left by the previous chunk
        notice("Reading genotype data from the current position to %s:%d", geno_chunk_done ? "(end)" : mpr.get_current_variant().schrom.c_str(), end);
    }

    if ( geno_chunk_done ) { // the region has already been fully streamed
        geno_chunk.geno_mat.resize(n_overlapping_samples, 0);
        geno_chunk.geno_mask.resize(n_overlapping_samples, 0);
        return false;
    }

    geno_chunk.clear();
    const std::vector<int32_t>& int_buf = mpr.get_int_buf();
    const double* dbl_buf = mpr.get_dbl_buf();
    int32_t n_col_est = 10;
    geno_chunk.geno_mat.resize(n_overlapping_samples, n_col_est);
    geno_chunk.geno_mask.resize(n_overlapping_samples, n_col_est);
    geno_chunk.n_geno_missing = 0;

    geno_chunk.n_samples = n_overlapping_samples;
    int32_t icol = 0;

    notice("Reading genotype data for %d overlapping samples", n_overlapping_samples);

    int32_t gcs[3] = {0, 0, 0};
    while ( true ) {
        const plink_var_t& var = mpr.get_current_variant();
        if ( var.pos > end ) { // moved past the requested region
            geno_chunk_done = true;
            break;
        }

        gcs[0] = gcs[1] = gcs[2] = 0; // initialize genotype counts
        std::string cpra_s(var.to_string());
        mpr.get_genos();

        if ( icol >= n_col_est ) {
            geno_chunk.geno_mat.conservativeResize(n_overlapping_samples, n_col_est * 2);
            geno_chunk.geno_mask.conservativeResize(n_overlapping_samples, n_col_est * 2);
            n_col_est *= 2;
        }
        int32_t an = 0;
        double ac = 0;
        bool skip  = false;
        if ( mpr.is_dosage_present() ) {
            if ( dbl_buf == NULL) {
                dbl_buf = mpr.get_dbl_buf();
            }
            double sumsq = 0;
            for(int32_t i =0; i < n_overlapping_samples; ++i) {
                double ds = 2.0 - dbl_buf[i];
                geno_chunk.geno_mat(i, icol) = ds;
                geno_chunk.geno_mask(i, icol) = true; // not missing
                ac += ds;
                an += 2;
                ++gcs[(ds < 0.5) ? 0 : ( (ds < 1.5) ? 1 : 2 )];
                sumsq += (ds * ds);
            }
            if ( ac == 0 || an == ac ) {
                skip = true; // skip monomorphic variants
            }
            else if ( ac < min_ac || ac > max_ac || ac < min_af * an || ac > max_af * an ) {
                skip = true; // skip variants outside the specified allele count/frequency range
            }
            else {
                double info = ( sumsq / 2.0 * an - 4.0 * ac * ac ) / ( 2.0 * ac * ( an - ac ) );
                geno_chunk.infos.push_back(info);
            }
        }
        else {
            for(int32_t i = 0; i < n_overlapping_samples; ++i) {
                switch(int_buf[i]) { // make sure to convert 1-based index to 0-based
                case 0:
                    an += 2;
                    ac += 2;
                    ++gcs[2];
                    break;
                case 1:
                    an += 2;
                    ++ac;
                    ++gcs[1];
                    break;
                case 2:
                    an += 2;
                    ++gcs[0];
                    break;
                }
            }
            if ( ac == 0 || an == ac ) {
                skip = true; // skip monomorphic variants
            }
            else if ( ac < min_ac || ac > max_ac || ac < min_af * an || ac > max_af * an ) {
                skip = true; // skip variants outside the specified allele count/frequency range
            }
            else {
                double mean = (double)ac / (double)an * 2.0;
                double info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
                geno_chunk.infos.push_back(info);
                for(int32_t i = 0; i < n_overlapping_samples; ++i) {
                    switch(int_buf[i]) {
                    case 0:
                        geno_chunk.geno_mat(i, icol) = 2.0 - mean; // homalt
                        geno_chunk.geno_mask(i, icol) = true; // not missing
                        break;
                    case 1:
                        geno_chunk.geno_mat(i, icol) = 1.0 - mean; // het
                        geno_chunk.geno_mask(i, icol) = true; // not missing
                        break;
                    case 2:
                        geno_chunk.geno_mat(i, icol) = 0.0 - mean; // homref
                        geno_chunk.geno_mask(i, icol) = true; // not missing
                        break;
                    default:
                        geno_chunk.geno_mat(i, icol) = 0; // missing - mean imputation
                        geno_chunk.geno_mask(i, icol) = false; // missing
                        ++geno_chunk.n_geno_missing;
                        break;
                    }
                }
            }
        }

        if ( !skip ) {
            geno_chunk.v_cpra.push_back(cpra_t(cpra_s.c_str())); 
            geno_chunk.var_cnts.push_back(var_cnt_t(an, ac, gcs[0], gcs[1], gcs[2]));
            // geno_chunk.acs.push_back(ac);
            // geno_chunk.ans.push_back(an);
            // geno_chunk.gc0s.push_back(gcs[0]);
            // geno_chunk.gc1s.push_back(gcs[1]);
            // geno_chunk.gc2s.push_back(gcs[2]);
            ++icol;
        }
        else {
            //++nskip;
            ++geno_chunk.n_skipped;
        }

        // advance to the next variant right away so the next chunk never re-reads this one
        if ( !mpr.read_pivar() ) {
            geno_chunk_done = true; // reached the end of available variants
            break;
        }

        if ( icol >= max_chunk_vars ) { // chunk is full; cur_var already points to the next variant
            break;
        }
    }

    // shrink to the number of variants actually loaded so callers never see garbage columns
    geno_chunk.geno_mat.conservativeResize(n_overlapping_samples, icol);
    geno_chunk.geno_mask.conservativeResize(n_overlapping_samples, icol);

    // Residualize the genotypes against the same covariates used to adjust the
    // phenotype. The phenotype matrix is covariate-adjusted in load_pheno_cov_matrices(),
    // so for the association test (and SuSiE) to estimate the *partial* effect of each
    // variant, the genotypes must be projected onto the orthogonal complement of the
    // covariate space as well (Frisch-Waugh-Lovell). Adjusting only the phenotype leaves
    // genotype variance that is collinear with the covariates (e.g. genotype PCs) in the
    // design, which inflates x'x, attenuates the correlation with the phenotype, and thus
    // shrinks beta and the test statistic. This mirrors susieR usage where both X and y
    // are residualized against the covariates before fitting.
    if ( icol > 0 && cov_matrix.pheno_mat.cols() > 0 ) {
        if ( cov_matrix.pheno_mat.rows() != n_overlapping_samples ) {
            error("Covariate matrix has %d rows but %d overlapping samples were expected",
                  (int32_t)cov_matrix.pheno_mat.rows(), n_overlapping_samples);
        }
        // same fit as the phenotypes, with the QR computed once per run. Variants the covariates
        // explain are zeroed, so the tests report NaN and SuSiE gives them no prior weight.
        std::shared_ptr<const CovariateResidualizer> resid = cov_resid;
        if ( !resid ) resid = std::make_shared<const CovariateResidualizer>(cov_matrix.pheno_mat);
        const int32_t n_collinear = resid->residualize_genotypes(geno_chunk.geno_mat);
        if ( n_collinear > 0 ) {
            notice("%d variants in this chunk are explained by the covariates (residual variance below %g of the total) and are reported as NaN",
                   n_collinear, CovariateResidualizer::kCollinearTol);
        }
    }


    notice("Loaded %d variants in the current chunk until %s, skipped %d", icol, mpr.get_current_variant().to_string().c_str(), geno_chunk.n_skipped);

    geno_chunk.n_variants = icol;
    return icol;
}
