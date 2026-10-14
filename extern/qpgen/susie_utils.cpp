#include "susie_utils.h"
#include "qgenlib/qgen_error.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>
#include <string>

namespace susie {

// Effects whose estimated prior variance is at most this are inactive: their
// alpha stays uniform (1/p) and carries no information. Same as susieR's
// prior_tol default.
static const double SUSIE_PRIOR_TOL = 1e-9;

// Reported PIP = 1 - prod_l (1 - alpha_lj) over ACTIVE effects only, as in
// susieR's susie_get_pip(prior_tol = 1e-9). Counting inactive effects would add
// a spurious floor of 1 - (1 - 1/p)^{#inactive} to every variant (0.077 at p =
// 100 with 8 inactive effects).
static VectorXd reported_pip(const MatrixXd& alpha, const VectorXd& V) {
    const int L = (int)alpha.rows();
    const int p = (int)alpha.cols();
    VectorXd pip(p);
    for (int j = 0; j < p; ++j) {
        double prod = 1.0;
        for (int l = 0; l < L; ++l)
            if (V(l) > SUSIE_PRIOR_TOL) prod *= (1.0 - alpha(l, j));
        pip(j) = 1.0 - prod;
    }
    return pip;
}

// ---------------------------------------------------------------------------
// Dense kernels for the one-time cross-products / decompositions
// ---------------------------------------------------------------------------
// With QPGEN_USE_LAPACK (set by CMake when a system LAPACK/BLAS is found) these
// call the (multithreaded) BLAS/LAPACK directly; otherwise Eigen's built-in,
// single-threaded routines. At n = 10K, p = 4000 on Accelerate: cross-product
// 0.27s vs 3.3s, eigendecomposition 5.7s vs 37s.
#ifdef QPGEN_USE_LAPACK
extern "C" {
void dsyrk_(const char* uplo, const char* trans, const int* n, const int* k,
            const double* alpha, const double* a, const int* lda,
            const double* beta, double* c, const int* ldc);
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
            const double* alpha, const double* a, const int* lda,
            const double* b, const int* ldb, const double* beta, double* c, const int* ldc);
void dsyevd_(const char* jobz, const char* uplo, const int* n, double* a, const int* lda,
             double* w, double* work, const int* lwork, int* iwork, const int* liwork, int* info);
void dsyevr_(const char* jobz, const char* range, const char* uplo, const int* n,
             double* a, const int* lda, const double* vl, const double* vu,
             const int* il, const int* iu, const double* abstol, int* m, double* w,
             double* z, const int* ldz, int* isuppz, double* work, const int* lwork,
             int* iwork, const int* liwork, int* info);
}
#endif

// Lower triangle of C = X'X (trans) or XX' (!trans).
static void sym_crossprod_lower(const MatrixXd& X, bool trans, MatrixXd& C) {
    const int m = trans ? (int)X.cols() : (int)X.rows();
    C.setZero(m, m);
#ifdef QPGEN_USE_LAPACK
    const int k = trans ? (int)X.rows() : (int)X.cols();
    const int lda = (int)X.rows();
    const double one = 1.0, zero = 0.0;
    dsyrk_("L", trans ? "T" : "N", &m, &k, &one, X.data(), &lda, &zero, C.data(), &m);
#else
    if (trans) C.selfadjointView<Eigen::Lower>().rankUpdate(X.transpose());
    else       C.selfadjointView<Eigen::Lower>().rankUpdate(X);
#endif
}

// Symmetric eigendecomposition from the lower triangle of C. On return w holds
// the eigenvalues in ascending order and C the corresponding eigenvectors.
static void sym_eigen_inplace(MatrixXd& C, VectorXd& w) {
    const int m = (int)C.rows();
    w.resize(m);
    if (m == 0) return;
#ifdef QPGEN_USE_LAPACK
    int info = 0;
    if (m <= 16384) {
        // divide-and-conquer: fastest, but needs ~2m^2 doubles of workspace
        int lwork = -1, liwork = -1, iwq = 0;
        double wq = 0.0;
        dsyevd_("V", "L", &m, C.data(), &m, w.data(), &wq, &lwork, &iwq, &liwork, &info);
        if (info == 0) {
            lwork = (int)wq; liwork = iwq;
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);
            dsyevd_("V", "L", &m, C.data(), &m, w.data(), work.data(), &lwork,
                    iwork.data(), &liwork, &info);
        }
    } else {
        // MRRR: O(m) workspace (plus the m x m output), for very large m
        MatrixXd Z(m, m);
        std::vector<int> isuppz(2 * (size_t)m);
        int nfound = 0, lwork = -1, liwork = -1, iwq = 0, il = 0, iu = 0;
        double wq = 0.0, vl = 0.0, vu = 0.0, abstol = 0.0;
        dsyevr_("V", "A", "L", &m, C.data(), &m, &vl, &vu, &il, &iu, &abstol, &nfound,
                w.data(), Z.data(), &m, isuppz.data(), &wq, &lwork, &iwq, &liwork, &info);
        if (info == 0) {
            lwork = (int)wq; liwork = iwq;
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);
            dsyevr_("V", "A", "L", &m, C.data(), &m, &vl, &vu, &il, &iu, &abstol, &nfound,
                    w.data(), Z.data(), &m, isuppz.data(), work.data(), &lwork,
                    iwork.data(), &liwork, &info);
        }
        if (info == 0) C = std::move(Z);
    }
    if (info != 0)
        error("SuSiE: LAPACK eigendecomposition of the %d x %d genotype cross-product failed (info=%d)", m, m, info);
#else
    Eigen::SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(C.selfadjointView<Eigen::Lower>());   // reads the lower triangle only
    if (es.info() != Eigen::Success)
        error("SuSiE: eigendecomposition of the %d x %d genotype cross-product failed", m, m);
    w = es.eigenvalues();
    C = es.eigenvectors();
#endif
}

// out = A' B
static void gemm_tn(const MatrixXd& A, const MatrixXd& B, MatrixXd& out) {
#ifdef QPGEN_USE_LAPACK
    const int m = (int)A.cols(), n = (int)B.cols(), k = (int)A.rows();
    out.resize(m, n);
    if (m == 0 || n == 0) return;
    const int lda = (int)A.rows(), ldb = (int)B.rows();
    const double one = 1.0, zero = 0.0;
    dgemm_("T", "N", &m, &n, &k, &one, A.data(), &lda, B.data(), &ldb, &zero, out.data(), &m);
#else
    out.noalias() = A.transpose() * B;
#endif
}

// ---------------------------------------------------------------------------
// Single-effect regression (SER) and prior-variance optimization
// ---------------------------------------------------------------------------

struct SER {
    VectorXd alpha;   // p
    VectorXd mu;      // p
    VectorXd mu2;     // p
    VectorXd lbf_var; // p
    double   lbf;     // scalar log Bayes factor for the single effect
    double   V;
};

static double ser_loglik(double v, const ArrayXd& betahat2,
                         const ArrayXd& s2, const ArrayXd& logpi) {
    ArrayXd lbf = 0.5 * (s2 / (s2 + v)).log()
                + 0.5 * (betahat2 / s2) * (v / (s2 + v));
    ArrayXd w = lbf + logpi;
    double m = w.maxCoeff();
    return m + std::log((w - m).exp().sum());
}

static double optimize_prior_variance(const ArrayXd& betahat2,
                                      const ArrayXd& s2,
                                      const ArrayXd& logpi,
                                      double lo, double hi) {
    // golden-section search on log(v), then compare against the null (v = 0)
    const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
    double a = std::log(lo), b = std::log(hi);
    double c = b - gr * (b - a), d = a + gr * (b - a);
    double fc = ser_loglik(std::exp(c), betahat2, s2, logpi);
    double fd = ser_loglik(std::exp(d), betahat2, s2, logpi);
    for (int i = 0; i < 100 && (b - a) > 1e-6; ++i) {
        if (fc > fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = ser_loglik(std::exp(c), betahat2, s2, logpi); }
        else         { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = ser_loglik(std::exp(d), betahat2, s2, logpi); }
    }
    double v_opt   = std::exp((a + b) / 2.0);
    double ll_opt  = ser_loglik(v_opt, betahat2, s2, logpi);
    double ll_null = ser_loglik(0.0,   betahat2, s2, logpi);
    return (ll_null > ll_opt) ? 0.0 : v_opt;
}

// xtr = X' r, with r the residual excluding the effect being updated.
static SER single_effect_regression(const ArrayXd& xtr,
                                     const ArrayXd& xtx, double sigma2,
                                     double V_init, const ArrayXd& logpi,
                                     bool estimate_v, double v_lo, double v_hi) {
    const int p = (int)xtr.size();
    ArrayXd betahat = xtr / xtx;
    ArrayXd s2 = sigma2 / xtx;
    ArrayXd betahat2 = betahat.square();

    double v = estimate_v ? optimize_prior_variance(betahat2, s2, logpi, v_lo, v_hi)
                          : V_init;

    SER out; out.V = v;
    if (v <= 0.0) {
        out.alpha   = VectorXd::Constant(p, 1.0 / p);
        out.mu      = VectorXd::Zero(p);
        out.mu2     = VectorXd::Zero(p);
        out.lbf_var = VectorXd::Zero(p);
        out.lbf     = 0.0;
        return out;
    }

    ArrayXd lbf = 0.5 * (s2 / (s2 + v)).log()
                + 0.5 * (betahat2 / s2) * (v / (s2 + v));
    ArrayXd w = lbf + logpi;
    double m = w.maxCoeff();
    ArrayXd a = (w - m).exp();
    a /= a.sum();

    double lbf_model = m + std::log((w - m).exp().sum()); // log average BF under prior

    ArrayXd post_var = 1.0 / (1.0 / s2 + 1.0 / v);
    ArrayXd mu  = (post_var / s2) * betahat;
    ArrayXd mu2 = mu.square() + post_var;

    out.alpha   = a.matrix();
    out.mu      = mu.matrix();
    out.mu2     = mu2.matrix();
    out.lbf_var = lbf.matrix();
    out.lbf     = lbf_model;
    return out;
}

// ---------------------------------------------------------------------------
// IBSS fit
// ---------------------------------------------------------------------------

SusieFit fit_susie(SusieDesign& d, const VectorXd& y, const VectorXd& Xty, const SusieOptions& opt) {
    const MatrixXd& X = d.X();
    const int n = d.n();
    const int p = d.p();
    const int L = std::max(1, std::min(opt.L, p));

    ArrayXd xtx = d.xtx();
    // Columns with no variance (monomorphic, or fully explained by the covariates
    // after the Frisch-Waugh-Lovell residualization) carry no information. They must
    // be given zero prior weight rather than merely a guarded x'x: dividing by a
    // sentinel x'x yields betahat2/s2 = 0/0 = NaN, and because every variable shares
    // the same softmax and residual, a single NaN propagates to the whole fit
    // (all-NaN alpha/lbf, sigma2 frozen at its initial value, no credible sets).
    int n_active = 0;
    for (int j = 0; j < p; ++j) {
        if (std::isfinite(xtx(j)) && xtx(j) > 0.0) ++n_active;
    }
    if (n_active == 0) {
        error("fit_susie(): all %d genotype columns have zero variance; nothing to fine-map", p);
    }
    ArrayXd logpi(p);
    const double logpi_active = -std::log((double)n_active);
    for (int j = 0; j < p; ++j) {
        if (std::isfinite(xtx(j)) && xtx(j) > 0.0) {
            logpi(j) = logpi_active;
        } else {
            xtx(j)   = 1.0;                                        // keep the arithmetic finite
            logpi(j) = -std::numeric_limits<double>::infinity();    // prior weight 0 -> alpha 0
        }
    }
    double var_y = (y.array() - y.mean()).square().sum() / std::max(1, n);

    MatrixXd alpha = MatrixXd::Zero(L, p);
    MatrixXd mu    = MatrixXd::Zero(L, p);
    MatrixXd mu2   = MatrixXd::Zero(L, p);
    MatrixXd lbf_variable = MatrixXd::Zero(L, p);
    VectorXd lbf   = VectorXd::Zero(L);
    VectorXd V     = VectorXd::Constant(L, opt.scaled_prior_variance * var_y);
    double sigma2  = var_y > 0.0 ? var_y : 1.0;

    VectorXd b_bar = VectorXd::Zero(p);

    // Two equivalent ways to form X'r_l for each single effect, picked by shape:
    //   p <= n : sufficient statistics (as susieR's susie_ss). X'X is the design's
    //            cached one (formed once per region, multithreaded BLAS when
    //            available), and each effect then costs ONE p x p product:
    //            X'r_l = X'y - X'X b_bar + X'X b_l.
    //   n <  p : individual-level, with X b_l cached per effect so each effect
    //            costs two n x p products (X'r_l, then X b_l for the update).
    // Either way the ERSS needs no further passes over X. (Recomputing the
    // residual and every ||X b_l||^2 from scratch took ~3L+1 passes over the n x p
    // matrix per iteration, far slower than SuSiE-inf's eigenspace updates.)
    const bool use_suff = (p <= n);
    static const MatrixXd empty;
    const MatrixXd& XtX = use_suff ? d.XtX() : empty;   // p x p (use_suff)
    MatrixXd Q = MatrixXd::Zero(use_suff ? p : n, L);  // use_suff ? X'X b_l : X b_l
    const double yty = y.squaredNorm();

    SusieFit fit; fit.converged = false; fit.niter = 0;
    double obj_prev = -std::numeric_limits<double>::infinity();
    MatrixXd alpha_prev = alpha;

    for (int iter = 0; iter < opt.max_iter; ++iter) {
        // Q_sum = X'X b_bar (use_suff) or X b_bar; rebuilt each sweep to avoid drift
        VectorXd Q_sum = Q.rowwise().sum();
        for (int l = 0; l < L; ++l) {
            VectorXd b_l = (alpha.row(l).array() * mu.row(l).array()).matrix();
            VectorXd Q_minus = Q_sum - Q.col(l);          // effect l removed
            ArrayXd xtr;
            if (use_suff) xtr = (Xty - Q_minus).array();
            else          xtr = (X.transpose() * (y - Q_minus)).array();
            SER ser = single_effect_regression(xtr, xtx, sigma2, V(l), logpi,
                                               opt.estimate_prior_variance,
                                               opt.prior_v_min, opt.prior_v_max);
            alpha.row(l)        = ser.alpha.transpose();
            mu.row(l)           = ser.mu.transpose();
            mu2.row(l)          = ser.mu2.transpose();
            lbf_variable.row(l) = ser.lbf_var.transpose();
            lbf(l)              = ser.lbf;
            V(l)                = ser.V;
            VectorXd b_l_new = (alpha.row(l).array() * mu.row(l).array()).matrix();
            b_bar = b_bar - b_l + b_l_new;
            if (use_suff) Q.col(l).noalias() = XtX * b_l_new;
            else          Q.col(l).noalias() = X * b_l_new;
            Q_sum = Q_minus + Q.col(l);
        }

        // Expected residual sum of squares, matching susieR's get_ER2():
        //   ERSS = ||y - X b_bar||^2 - sum_l ||X b_l||^2 + sum_l sum_j d_j alpha_lj mu2_lj
        // The exact ||X b_l||^2 (a full quadratic form) is required so that LD
        // cross-terms are retained; approximating it by sum_j d_j (alpha_lj mu_lj)^2
        // inflates sigma2 for correlated X and over-shrinks the single effects.
        // Both terms come from the cached Q: with sufficient statistics
        // ||y - X b||^2 = y'y - 2 b'X'y + b'X'X b and ||X b_l||^2 = b_l'(X'X b_l).
        Q_sum = Q.rowwise().sum();
        double erss = use_suff ? yty - 2.0 * b_bar.dot(Xty) + b_bar.dot(Q_sum)
                               : (y - Q_sum).squaredNorm();
        for (int l = 0; l < L; ++l) {
            ArrayXd a  = alpha.row(l).array();
            ArrayXd m2 = mu2.row(l).array();
            if (use_suff) {
                VectorXd b_l = (a * mu.row(l).transpose().array()).matrix();
                erss -= b_l.dot(Q.col(l));
            } else {
                erss -= Q.col(l).squaredNorm();
            }
            erss += (xtx * (a * m2)).sum();
        }
        if (opt.estimate_residual_variance && std::isfinite(erss) && erss > 0.0) sigma2 = erss / n;

        double elbo = -0.5 * n * std::log(2.0 * M_PI * sigma2) - erss / (2.0 * sigma2);
        for (int l = 0; l < L; ++l) elbo += lbf(l);
        fit.elbo.push_back(elbo);

        fit.niter = iter + 1;
        bool done = false;
        if (opt.convergence_method == SusieOptions::ELBO) {
            if (iter > 0 && std::abs(elbo - obj_prev) < opt.tol) done = true;
            obj_prev = elbo;
        } else {
            double da = (alpha - alpha_prev).cwiseAbs().maxCoeff();
            if (iter > 0 && da < opt.tol) done = true;
            alpha_prev = alpha;
        }
        if (done) { fit.converged = true; break; }
    }

    VectorXd pip = reported_pip(alpha, V);

    VectorXd Xb = X * b_bar;
    fit.alpha = alpha; fit.mu = mu; fit.mu2 = mu2;
    fit.lbf_variable = lbf_variable; fit.lbf = lbf;
    fit.pip = pip; fit.Xr = Xb; fit.fitted = Xb;
    fit.V = V; fit.sigma2 = sigma2; fit.intercept = 0.0;
    return fit;
}

// ---------------------------------------------------------------------------
// SuSiE-inf (unmappable infinitesimal effects), matching susieR's
// unmappable_effects = "inf" with estimate_residual_method = "MoM" and
// PIP-based convergence.
//
// Model: y = sum_l X b_l + X theta + e, theta_j ~ N(0, tau2), e_i ~ N(0, sigma2).
// The infinitesimal effect makes the residual covariance tau2 XX' + sigma2 I.
// All per-effect single-effect regressions are performed in the eigenspace of
// the (standardized) design so that Omega = (tau2 XX' + sigma2 I)^{-1} enters
// through Omega-weighted sufficient statistics X'Omega y and diag(X'Omega X).
// ---------------------------------------------------------------------------

// Negative SER log-likelihood as a function of the prior variance V, on the
// Omega-whitened scale (predictor weights pw = diag(X'Omega X), residual r =
// X'Omega (y - X b_{-l})). Mirrors susieR's neg_loglik on the inf path with the
// shat2 inflation factor equal to 1 (no finite-reference-R correction).
static double susie_inf_negll(double V, const ArrayXd& pw, const ArrayXd& res,
                              const ArrayXd& logpi) {
    ArrayXd denom = 1.0 + V * pw;
    ArrayXd lbf = -0.5 * denom.log() + 0.5 * V * res.square() / denom;
    ArrayXd w = lbf + logpi;
    double m = w.maxCoeff();
    return -(m + std::log((w - m).exp().sum()));
}

// Optimize V over [0,1] (golden section) then keep whichever of {optimum, V_init}
// has the higher likelihood, and finally snap to 0 when the null (V=0, loglik 0)
// is at least as good. Mirrors susieR's optimize_scalar_prior_variance for inf.
static double susie_inf_optimize_V(const ArrayXd& pw, const ArrayXd& res,
                                   const ArrayXd& logpi, double V_init) {
    const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
    double a = 0.0, b = 1.0;
    double c = b - gr * (b - a), d = a + gr * (b - a);
    double fc = susie_inf_negll(c, pw, res, logpi);
    double fd = susie_inf_negll(d, pw, res, logpi);
    for (int i = 0; i < 100 && (b - a) > 1e-6; ++i) {
        if (fc < fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = susie_inf_negll(c, pw, res, logpi); }
        else         { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = susie_inf_negll(d, pw, res, logpi); }
    }
    double v_opt = (a + b) / 2.0;
    double f_opt = susie_inf_negll(v_opt, pw, res, logpi);
    double f_init = susie_inf_negll(V_init, pw, res, logpi);
    double v_best = (f_init < f_opt) ? V_init : v_opt;
    double f_best = (f_init < f_opt) ? f_init : f_opt;
    // null (V=0) has loglik 0; keep it if it is at least as good
    if (f_best >= 0.0) return 0.0;
    return v_best;
}

// ---------------------------------------------------------------------------
// SusieDesign: trait-independent, region-level work
// ---------------------------------------------------------------------------

SusieDesign::SusieDesign(MatrixXd&& X, bool standardize) : X_(std::move(X)) {
    const int n = (int)X_.rows();
    const int p = (int)X_.cols();
    xtx_.setZero(p);
    col_scale_.setZero(p);

    // Every single-effect regression shares one residual and one softmax across all p
    // variables, so a single bad column contaminates the entire fit rather than just
    // its own coefficient (unlike the marginal per-variant regressions, which stay
    // correct for every other variant). Two kinds of bad column are screened here:
    //   - non-finite entries (e.g. a NaN dosage that survived mean-imputation);
    //   - zero variance (monomorphic in the analyzed samples, an imputed variant with
    //     an identical dosage for everybody, or a variant fully explained by the
    //     covariates after the FWL residualization).
    // Both are zeroed out and reported; the fits then give them zero prior weight
    // so their PIP/alpha are 0 and the remaining variants are fit normally.
    int n_nonfinite_cols = 0, n_novar_cols = 0;
    std::vector<int> bad_cols;
    for (int j = 0; j < p; ++j) {
        bool finite = X_.col(j).allFinite();
        double mean = finite ? X_.col(j).mean() : 0.0;
        double ss   = finite ? (X_.col(j).array() - mean).square().sum() : 0.0;
        double sd   = std::sqrt(ss / std::max(1, n - 1));
        if (!finite || !(ss > 0.0) || !(sd > 0.0)) {
            if (!finite) ++n_nonfinite_cols; else ++n_novar_cols;
            if ((int)bad_cols.size() < 10) bad_cols.push_back(j);
            X_.col(j).setZero();
            continue;
        }
        if (standardize) { X_.col(j) = (X_.col(j).array() - mean) / sd; col_scale_(j) = sd; }
        else             { X_.col(j) = X_.col(j).array() - mean;        col_scale_(j) = 1.0; }
        xtx_(j) = X_.col(j).squaredNorm();
    }
    n_unusable_ = n_nonfinite_cols + n_novar_cols;
    if (n_unusable_ > 0) {
        std::string idx;
        for (size_t i = 0; i < bad_cols.size(); ++i) {
            if (i) idx += ",";
            idx += std::to_string(bad_cols[i]);
        }
        if ((int)bad_cols.size() < n_unusable_) idx += ",...";
        notice("SuSiE: excluding %d of %d variants with no usable genotype variance "
               "(%d with non-finite values, %d monomorphic/collinear-with-covariates); "
               "0-based column indices: %s",
               n_unusable_, p, n_nonfinite_cols, n_novar_cols, idx.c_str());
    }
}

MatrixXd SusieDesign::crossprod(const MatrixXd& Y) const {
    MatrixXd XtY;
    gemm_tn(X_, Y, XtY);
    return XtY;
}

const MatrixXd& SusieDesign::XtX() {
    if (!has_XtX_) {
        sym_crossprod_lower(X_, true, XtX_);
        XtX_.triangularView<Eigen::StrictlyUpper>() = XtX_.transpose(); // full symmetric for GEMV
        has_XtX_ = true;
    }
    return XtX_;
}

// Thin eigendecomposition of the design, X = U D V', shared by SuSiE-inf and
// SuSiE-ash (and by every trait on this design). Equivalent to susieR's svd(X)
// path: only the variant-space eigenvectors V and the eigenvalues d^2 are needed.
//
// The decomposition is taken on whichever cross-product is SMALLER:
//   p <= n : X'X (p x p, the cached one), whose eigenvectors are V directly;
//   n <  p : XX' (n x n), with V = X'U D^{-1} formed by a single GEMM.
// This keeps the cost at O(n p min(n,p) + min(n,p)^3), like a thin SVD. (Always
// using the n x n Gram matrix made the cost O(n^3) and the memory O(n^2) even
// for a handful of variants, e.g. ~8 min at n = 10K, p = 50.) Components with
// ~0 singular value contribute ~0 everywhere and are dropped.
const ThinEigen& SusieDesign::thin_eigen() {
    if (has_eigen_) return eigen_;
    const bool variant_space = (p() <= n());
    const int m = variant_space ? p() : n();

    MatrixXd C;
    if (variant_space) C = XtX();                    // copy: the solver overwrites it
    else               sym_crossprod_lower(X_, false, C);
    VectorXd evals;                                  // ascending
    sym_eigen_inplace(C, evals);                     // C <- eigenvectors

    const double ev_tol = std::max(evals.maxCoeff(), 0.0) * 1e-8;
    std::vector<int> keep;
    for (int k = 0; k < m; ++k) if (evals(k) > ev_tol) keep.push_back(k);
    const int r = (int)keep.size();

    eigen_.eigval.resize(r);
    MatrixXd E(m, r);                                // kept eigenvectors of C
    for (int idx = 0; idx < r; ++idx) {
        eigen_.eigval(idx) = evals(keep[idx]);
        E.col(idx) = C.col(keep[idx]);
    }
    C.resize(0, 0);
    if (variant_space) {
        eigen_.Vmat = std::move(E);
    } else {
        // V = X'U D^{-1}
        const ArrayXd dk = eigen_.eigval.array().sqrt();
        gemm_tn(X_, E, eigen_.Vmat);
        eigen_.Vmat.array().rowwise() /= dk.transpose();
    }
    has_eigen_ = true;
    return eigen_;
}

SusieFit fit_susie_inf(SusieDesign& d, const VectorXd& y, const VectorXd& Xty, const SusieOptions& opt) {
    const MatrixXd& X = d.X();
    const int n = d.n();
    const int p = d.p();
    const int L = std::max(1, std::min(opt.L, p));
    const double var_y = y.squaredNorm() / std::max(1, n - 1); // y is pre-centered

    // --- thin eigendecomposition of the design (cached, shared by all traits) ---
    const ThinEigen& te = d.thin_eigen();
    const VectorXd& eigval = te.eigval;             // d_k^2
    const MatrixXd& Vmat   = te.Vmat;               // p x r
    const VectorXd VtXty   = Vmat.transpose() * Xty; // r (trait-specific)
    MatrixXd Vsq = Vmat.array().square();           // p x r
    const double yty = y.squaredNorm();
    ArrayXd logpi = ArrayXd::Constant(p, -std::log((double)p));

    // --- initialization (matches susieR init: uniform alpha, V = spv*var_y) ---
    MatrixXd alpha = MatrixXd::Constant(L, p, 1.0 / p);
    MatrixXd mu    = MatrixXd::Zero(L, p);
    MatrixXd mu2   = MatrixXd::Zero(L, p);
    MatrixXd lbf_variable = MatrixXd::Zero(L, p);
    VectorXd lbf   = VectorXd::Zero(L);
    VectorXd V     = VectorXd::Constant(L, opt.scaled_prior_variance * var_y);
    double sigma2  = var_y;
    double tau2    = 0.0;
    VectorXd theta = VectorXd::Zero(p);

    // Omega-derived caches, refreshed whenever (tau2, sigma2) change.
    ArrayXd omega_var = tau2 * eigval.array() + sigma2;                  // r
    VectorXd pw = Vsq * (eigval.array() / omega_var).matrix();          // p = diag(X'Omega X)
    VectorXd XtOmegay = Vmat * (VtXty.array() / omega_var).matrix();    // p

    VectorXd pip_prev = VectorXd::Zero(p);
    MatrixXd alpha_prev = alpha;

    // VtB.col(l) = V' b_l with b_l = alpha_l * mu_l, kept in sync with (alpha, mu)
    // so each SER needs only two p x r matrix-vector products and the MoM step
    // reuses them instead of recomputing V' b_l for every effect.
    MatrixXd VtB = MatrixXd::Zero(eigval.size(), L);                             // r x L

    SusieFit fit; fit.converged = false; fit.niter = 0;

    for (int iter = 0; iter < opt.max_iter; ++iter) {
        VectorXd Vtb_sum = VtB.rowwise().sum();                                     // r = V' b_bar
        // ---- single-effect regressions (Omega-weighted) ----
        for (int l = 0; l < L; ++l) {
            VectorXd Vtb     = Vtb_sum - VtB.col(l);                                // r = V' b_{-l}
            VectorXd XtOmegaXb = Vmat * (Vtb.array() * eigval.array() / omega_var).matrix();
            ArrayXd res = (XtOmegay - XtOmegaXb).array();                           // p (= X'Omega r)

            double Vl = susie_inf_optimize_V(pw.array(), res, logpi, V(l));
            V(l) = Vl;

            if (Vl <= 0.0) {
                alpha.row(l).setConstant(1.0 / p);
                mu.row(l).setZero();
                mu2.row(l).setZero();
                lbf_variable.row(l).setZero();
                lbf(l) = 0.0;
                VtB.col(l).setZero();
                Vtb_sum = Vtb;
                continue;
            }
            ArrayXd denom = 1.0 + Vl * pw.array();
            ArrayXd lbfj  = -0.5 * denom.log() + 0.5 * Vl * res.square() / denom;
            ArrayXd w = lbfj + logpi;
            double m = w.maxCoeff();
            ArrayXd a = (w - m).exp();
            double s = a.sum();
            a /= s;
            ArrayXd post_var = Vl / denom;            // V*shat2/(V+shat2), shat2 = 1/pw
            ArrayXd post_mean = post_var * res;       // post_var/shat2 * betahat
            alpha.row(l) = a.transpose();
            mu.row(l)    = post_mean.transpose();
            mu2.row(l)   = (post_var + post_mean.square()).transpose();
            lbf_variable.row(l) = lbfj.transpose();
            lbf(l) = m + std::log(s);
            VtB.col(l).noalias() = Vmat.transpose() * (a * post_mean).matrix();
            Vtb_sum = Vtb + VtB.col(l);
        }
        fit.niter = iter + 1;

        // ---- PIP-based convergence (susieR forces this for inf) ----
        VectorXd pip(p);
        for (int j = 0; j < p; ++j) {
            double prod = 1.0;
            for (int l = 0; l < L; ++l) prod *= (1.0 - alpha(l, j));
            pip(j) = 1.0 - prod;
        }
        if (iter > 0) {
            double da = (alpha - alpha_prev).cwiseAbs().maxCoeff();
            double dp = (pip - pip_prev).cwiseAbs().maxCoeff();
            if (std::max(da, dp) < opt.tol) { fit.converged = true; break; }
        }
        alpha_prev = alpha;
        pip_prev   = pip;

        // ---- update variance components (method of moments) + theta BLUP ----
        VectorXd b   = (alpha.array() * mu.array()).colwise().sum().transpose();
        const VectorXd& Vtb_all = Vtb_sum;                                          // r = V' b
        // theta uses the current (pre-update) Omega caches and tau2, as in susieR.
        VectorXd XtOmegaXb_all = Vmat * (Vtb_all.array() * eigval.array() / omega_var).matrix();
        theta = tau2 * (XtOmegay - XtOmegaXb_all);

        // diag(V' M V), M = bbar bbar' - sum_l bl bl' + diag(sum_l alpha_l (mu_l^2 + 1/omega_l))
        ArrayXd diagVtMV = Vtb_all.array().square();
        ArrayXd tmpD = ArrayXd::Zero(p);
        for (int l = 0; l < L; ++l) {
            diagVtMV -= VtB.col(l).array().square();
            ArrayXd omega_l = pw.array() + 1.0 / V(l);   // V(l)==0 -> +inf -> 1/omega_l == 0
            tmpD += alpha.row(l).transpose().array() *
                    (mu.row(l).transpose().array().square() + 1.0 / omega_l);
        }
        diagVtMV += (Vsq.transpose() * tmpD.matrix()).array();

        const double sumEv  = eigval.sum();
        const double sumEv2 = eigval.array().square().sum();
        Eigen::Matrix2d A;
        A << (double)n, sumEv, sumEv, sumEv2;
        double x1 = yty - 2.0 * b.dot(Xty) + (eigval.array() * diagVtMV).sum();
        double x2 = Xty.squaredNorm()
                    - 2.0 * (Vtb_all.array() * VtXty.array() * eigval.array()).sum()
                    + (eigval.array().square() * diagVtMV).sum();
        Eigen::Vector2d sol = A.colPivHouseholderQr().solve(Eigen::Vector2d(x1, x2));
        if (sol(0) > 0.0 && sol(1) > 0.0) { sigma2 = sol(0); tau2 = sol(1); }
        else                              { sigma2 = x1 / n; tau2 = 0.0; }
        if (sigma2 <= 0.0) sigma2 = var_y * 1e-8;

        // refresh Omega caches with the new (tau2, sigma2)
        omega_var = tau2 * eigval.array() + sigma2;
        pw        = Vsq * (eigval.array() / omega_var).matrix();
        XtOmegay  = Vmat * (VtXty.array() / omega_var).matrix();
    }

    VectorXd pip = reported_pip(alpha, V);

    VectorXd b_final = (alpha.array() * mu.array()).colwise().sum().transpose();
    VectorXd Xr = X * (b_final + theta);
    fit.alpha = alpha; fit.mu = mu; fit.mu2 = mu2;
    fit.lbf_variable = lbf_variable; fit.lbf = lbf;
    fit.pip = pip; fit.Xr = Xr; fit.fitted = Xr;
    fit.V = V; fit.sigma2 = sigma2; fit.intercept = 0.0;
    fit.tau2 = tau2; fit.theta = theta;
    return fit;
}

// ---------------------------------------------------------------------------
// SuSiE-ash, following Algorithm 1 and Supplementary Notes S1-S3 of the SuSiE 2.0
// manuscript (McCreight et al., bioRxiv 2025.11.25.690514v2).
//
// Model: y = X beta + X theta + e,  e ~ N(0, sigma2 I),  beta = sum_l b_l gamma_l,
//        theta_j ~ sum_k pi_k N(0, sigma2 * sa2_k)       (sa2_1 = 0, point mass)
// theta enters the SER only through its marginal variance tau2 = sigma2 *
// sum_k pi_k sa2_k, i.e. through Omega = (tau2 XX' + sigma2 I)^{-1}. Each outer
// iteration:
//   Step 1  Omega-weighted SER for every effect l (the SuSiE-inf machinery),
//           on the leave-one-effect-out residual y - X sum_{l' != l} b_l'.
//   Step 2  Converge on max |alpha^(t) - alpha^(t-1)| < tol.
//   Step 3  Provisional (sigma2, tau2_0) by method of moments (as SuSiE-inf);
//           build the sa2 grid from tau2_0 / sigma2 (S3); fit Mr.ASH to
//           convergence on r = y - X b_bar: coordinate-ascent normal-means
//           updates of theta (S19-S26), pi_k = mean_j phi_jk, and sigma2 (S3).
//   Step 4  tau2 = sigma2 * sum_k pi_k sa2_k; refresh the Omega caches.
//
// Where the manuscript is silent, choices follow the susieR code of that period
// (initialize_mrash(), Dec 2025): the grid band sizes (10/6/3 points), the 10%
// relative thinning of grid points, the active-effect threshold for the grid's
// upper bound, starting Mr.ASH at the MoM sigma2 with uniform pi, and skipping
// Mr.ASH (theta = 0, tau2 = 0) when the MoM tau2_0 is 0.
//
// Deliberately NOT ported from later susieR versions (neither is in the
// manuscript): the LD-based masking / "confident effect" subtraction layer, and
// the Beta-Binomial slot-activity weights. Note that susieR's individual-level
// ash path passes the raw (unscaled, uncentered) X to mr.ash together with
// effects on the standardized scale, so it disagrees with its own sufficient-
// statistics path unless X is standardized beforehand; the design here is
// always centered/standardized, so that issue cannot arise.
// ---------------------------------------------------------------------------

// Solve the SuSiE-inf 2x2 MoM system for (sigma2, tau2) and return the
// non-negative pair (or (x1/n, 0) when the system is not proper).
static void mom_variance_components(
    int n, double yty, const VectorXd& b, const VectorXd& Xty,
    const VectorXd& eigval, const VectorXd& Vtb_all, const VectorXd& VtXty,
    const ArrayXd& diagVtMV, double var_y_floor,
    double& sigma2, double& tau2) {
    const double sumEv  = eigval.sum();
    const double sumEv2 = eigval.array().square().sum();
    Eigen::Matrix2d A;
    A << (double)n, sumEv, sumEv, sumEv2;
    double x1 = yty - 2.0 * b.dot(Xty) + (eigval.array() * diagVtMV).sum();
    double x2 = Xty.squaredNorm()
                - 2.0 * (Vtb_all.array() * VtXty.array() * eigval.array()).sum()
                + (eigval.array().square() * diagVtMV).sum();
    Eigen::Vector2d sol = A.colPivHouseholderQr().solve(Eigen::Vector2d(x1, x2));
    if (sol(0) > 0.0 && sol(1) > 0.0) { sigma2 = sol(0); tau2 = sol(1); }
    else                              { sigma2 = std::max(x1 / n, var_y_floor); tau2 = 0.0; }
}

// Data-driven sa2 grid (Supplementary Notes S3), in units of sigma2 (theta_j ~
// N(0, sigma2 * sa2_k)). t = tau2_0 / sigma2 from the MoM step; upper = smallest
// active sparse-effect prior variance / sigma2. Three log-spaced bands:
// [t/100, t/10] dense (10 points), [t/10, 3t] moderate (6), [3t, upper] coarse
// (3, only if upper > 3t); points within 10% of the previous kept point are
// dropped, and the point mass at zero is prepended.
static ArrayXd ash_variance_grid(double t, double upper) {
    std::vector<double> g;
    auto logspace = [&g](double a, double b, int m) {
        for (int i = 0; i < m; ++i)
            g.push_back(std::exp(std::log(a) + (std::log(b) - std::log(a)) * i / (double)(m - 1)));
    };
    logspace(t / 100.0, t / 10.0, 10);
    logspace(t / 10.0, 3.0 * t, 6);
    if (upper > 3.0 * t) logspace(3.0 * t, upper, 3);
    std::sort(g.begin(), g.end());
    std::vector<double> kept;
    for (double v : g) {
        if (!(v > 0.0) || !std::isfinite(v)) continue;
        if (kept.empty() || (v - kept.back()) / kept.back() > 0.1) kept.push_back(v);
    }
    ArrayXd sa2(kept.size() + 1);
    sa2(0) = 0.0;
    for (size_t k = 0; k < kept.size(); ++k) sa2(k + 1) = kept[k];
    return sa2;
}

// Mr.ASH (Kim et al. 2024; manuscript S19-S26 and S3) fit to convergence on the
// residual r0 = y - X b_bar, updating theta (warm-started), pi (unless fixed)
// and sigma2 (if requested). When the design caches X'X (p <= n) the
// coordinate updates run on X'r (O(p) per coordinate); otherwise on the n-vector
// residual (O(n) per coordinate). Returns the number of sweeps.
static int mr_ash_fit(SusieDesign& d, const VectorXd& Xtr0, double rtr0, const VectorXd& r0,
                      const ArrayXd& sa2, ArrayXd& mix_pi, bool update_pi,
                      double& sigma2, bool update_sigma2,
                      VectorXd& theta, int max_iter, double convtol) {
    const int n = d.n();
    const int p = d.p();
    const int K = (int)sa2.size();
    const ArrayXd& w = d.xtx();
    const bool use_xtx = (p <= n);
    const MatrixXd& X = d.X();
    static const MatrixXd empty;
    const MatrixXd& XtX = use_xtx ? d.XtX() : empty;
    int p_eff = 0;
    for (int j = 0; j < p; ++j) if (w(j) > 0.0) ++p_eff;
    if (p_eff == 0) return 0;

    // residual state for the current theta: X'r (use_xtx) or r
    VectorXd Xtr, r;
    if (use_xtx) Xtr = Xtr0 - XtX * theta;
    else         r = r0 - X * theta;

    ArrayXd logw(K), phi(K), mu(K);
    int iter = 0;
    for (; iter < max_iter; ++iter) {
        const VectorXd theta_old = theta;
        const ArrayXd log_pi = mix_pi.max(1e-300).log();
        ArrayXd pi_new = ArrayXd::Zero(K);
        double a1 = 0.0;                          // sum_j (x_j' r_{-j}) * theta_j
        for (int j = 0; j < p; ++j) {
            const double wj = w(j);
            if (!(wj > 0.0)) { theta(j) = 0.0; continue; }
            const double tj_old = theta(j);
            const double bjwj = (use_xtx ? Xtr(j) : X.col(j).dot(r)) + wj * tj_old;
            // NM posterior per component (S21-S26): mu_k = bjwj / (w_j + 1/sa2_k),
            // log BF_k = -0.5 log(1 + sa2_k w_j) + bjwj mu_k / (2 sigma2)
            double mmax = -std::numeric_limits<double>::infinity();
            for (int k = 0; k < K; ++k) {
                if (sa2(k) > 0.0) {
                    mu(k) = bjwj / (wj + 1.0 / sa2(k));
                    logw(k) = log_pi(k) - 0.5 * std::log1p(sa2(k) * wj) + mu(k) * bjwj / (2.0 * sigma2);
                } else {
                    mu(k) = 0.0;
                    logw(k) = log_pi(k);
                }
                if (logw(k) > mmax) mmax = logw(k);
            }
            phi = (logw - mmax).exp();
            phi /= phi.sum();
            const double tj = (phi * mu).sum();
            theta(j) = tj;
            pi_new += phi;
            a1 += bjwj * tj;
            const double delta = tj - tj_old;
            if (delta != 0.0) {
                if (use_xtx) Xtr.noalias() -= XtX.col(j) * delta;
                else         r.noalias()   -= X.col(j) * delta;
            }
        }
        pi_new /= (double)p_eff;
        if (update_pi) mix_pi = pi_new;

        // sigma2 (S3): [||r - X theta||^2 + sum_j sum_{k>=2} phi (w_j + 1/sa2_k)(mu^2 + s^2)
        //               - sum_j w_j theta_j^2] / (n + p (1 - pi_1))
        // which, with s^2 = sigma2 / (w_j + 1/sa2_k), equals
        //   (varobj + p (1 - pi_1) sigma2) / (n + p (1 - pi_1)),
        //   varobj = ||r - X theta||^2 - sum_j w_j theta_j^2 + sum_j bjwj theta_j.
        if (update_sigma2) {
            double rss = use_xtx ? rtr0 - theta.dot(Xtr0) - theta.dot(Xtr) : r.squaredNorm();
            double varobj = rss - (w * theta.array().square()).sum() + a1;
            double nonnull = (double)p_eff * (1.0 - pi_new(0));
            double s2 = (varobj + nonnull * sigma2) / ((double)n + nonnull);
            if (std::isfinite(s2) && s2 > 0.0) sigma2 = s2;
        }

        // mr.ash convergence rule: relative L2 change of theta
        if ((theta - theta_old).norm() < convtol * std::max(1.0, theta.norm())) { ++iter; break; }
    }
    return iter;
}

SusieFit fit_susie_ash(SusieDesign& d, const VectorXd& y, const VectorXd& Xty, const SusieOptions& opt) {
    const MatrixXd& X = d.X();
    const int n = d.n();
    const int p = d.p();
    const int L = std::max(1, std::min(opt.L, p));
    const bool fix_pi  = !opt.ash_fix_pi.empty();
    const bool fix_sa2 = !opt.ash_fix_sa2.empty();
    const double var_y = y.squaredNorm() / std::max(1, n - 1);

    // --- thin eigendecomposition of the design (cached, same as fit_susie_inf) ---
    const ThinEigen& te = d.thin_eigen();
    const VectorXd& eigval = te.eigval;
    const MatrixXd& Vmat   = te.Vmat;
    const VectorXd VtXty   = Vmat.transpose() * Xty;
    MatrixXd Vsq = Vmat.array().square();
    const double yty = y.squaredNorm();
    ArrayXd logpi_ser = ArrayXd::Constant(p, -std::log((double)p));

    // --- SER state (initialization as SuSiE-inf: uniform alpha, V = spv*var_y) ---
    MatrixXd alpha = MatrixXd::Constant(L, p, 1.0 / p);
    MatrixXd mu    = MatrixXd::Zero(L, p);
    MatrixXd mu2   = MatrixXd::Zero(L, p);
    MatrixXd lbf_variable = MatrixXd::Zero(L, p);
    VectorXd lbf   = VectorXd::Zero(L);
    VectorXd V     = VectorXd::Constant(L, opt.scaled_prior_variance * var_y);
    double sigma2  = var_y;
    double tau2    = 0.0;

    // --- Mr.ASH state ---
    VectorXd theta = VectorXd::Zero(p);
    ArrayXd sa2(1);  sa2(0) = 0.0;
    ArrayXd mix_pi(1); mix_pi(0) = 1.0;
    if (fix_sa2) {
        sa2.resize((int)opt.ash_fix_sa2.size());
        for (int k = 0; k < sa2.size(); ++k) sa2(k) = opt.ash_fix_sa2[k];
    }
    if (fix_pi) {
        // the data-driven grid changes size every iteration, so a fixed pi
        // needs a fixed grid of the same length
        if (!fix_sa2 || opt.ash_fix_sa2.size() != opt.ash_fix_pi.size())
            error("SuSiE-ash: --ash-fix-pi requires --ash-fix-sa2 with the same number of components");
        mix_pi.resize((int)opt.ash_fix_pi.size());
        for (int k = 0; k < mix_pi.size(); ++k) mix_pi(k) = opt.ash_fix_pi[k];
        if (mix_pi.sum() > 0.0) mix_pi /= mix_pi.sum();
    }

    // Omega caches (theta enters through Omega, not by subtracting X*theta from y).
    ArrayXd omega_var = tau2 * eigval.array() + sigma2;
    VectorXd pw = Vsq * (eigval.array() / omega_var).matrix();
    VectorXd XtOmegay = Vmat * (VtXty.array() / omega_var).matrix();

    // VtB.col(l) = V' b_l, kept in sync with (alpha, mu) as in fit_susie_inf
    MatrixXd VtB = MatrixXd::Zero(eigval.size(), L);
    MatrixXd alpha_prev = alpha;
    SusieFit fit; fit.converged = false; fit.niter = 0;

    for (int iter = 0; iter < opt.max_iter; ++iter) {
        // ---- Step 1: Omega-weighted single-effect regressions ----
        VectorXd Vtb_sum = VtB.rowwise().sum();
        for (int l = 0; l < L; ++l) {
            VectorXd Vtb = Vtb_sum - VtB.col(l);
            VectorXd XtOmegaXb = Vmat * (Vtb.array() * eigval.array() / omega_var).matrix();
            ArrayXd res = (XtOmegay - XtOmegaXb).array();

            double Vl = susie_inf_optimize_V(pw.array(), res, logpi_ser, V(l));
            V(l) = Vl;
            if (Vl <= 0.0) {
                alpha.row(l).setConstant(1.0 / p);
                mu.row(l).setZero();
                mu2.row(l).setZero();
                lbf_variable.row(l).setZero();
                lbf(l) = 0.0;
                VtB.col(l).setZero();
                Vtb_sum = Vtb;
                continue;
            }
            ArrayXd denom = 1.0 + Vl * pw.array();
            ArrayXd lbfj  = -0.5 * denom.log() + 0.5 * Vl * res.square() / denom;
            ArrayXd w = lbfj + logpi_ser;
            double m = w.maxCoeff();
            ArrayXd a = (w - m).exp();
            double s = a.sum();
            a /= s;
            ArrayXd post_var = Vl / denom;
            ArrayXd post_mean = post_var * res;
            alpha.row(l) = a.transpose();
            mu.row(l)    = post_mean.transpose();
            mu2.row(l)   = (post_var + post_mean.square()).transpose();
            lbf_variable.row(l) = lbfj.transpose();
            lbf(l) = m + std::log(s);
            VtB.col(l).noalias() = Vmat.transpose() * (a * post_mean).matrix();
            Vtb_sum = Vtb + VtB.col(l);
        }
        fit.niter = iter + 1;

        // ---- Step 2: convergence on alpha (Algorithm 1) ----
        if (iter > 0 && (alpha - alpha_prev).cwiseAbs().maxCoeff() < opt.tol) {
            fit.converged = true;
            break;
        }
        alpha_prev = alpha;

        // ---- Step 3a: provisional (sigma2, tau2_0) by MoM, as SuSiE-inf ----
        VectorXd b = (alpha.array() * mu.array()).colwise().sum().transpose();
        const VectorXd& Vtb_all = Vtb_sum;
        ArrayXd diagVtMV = Vtb_all.array().square();
        ArrayXd tmpD = ArrayXd::Zero(p);
        for (int l = 0; l < L; ++l) {
            diagVtMV -= VtB.col(l).array().square();
            ArrayXd omega_l = pw.array() + 1.0 / V(l);   // V(l)==0 -> 1/omega_l == 0
            tmpD += alpha.row(l).transpose().array() *
                    (mu.row(l).transpose().array().square() + 1.0 / omega_l);
        }
        diagVtMV += (Vsq.transpose() * tmpD.matrix()).array();
        double sigma2_mom = sigma2, tau2_mom = tau2;
        mom_variance_components(n, yty, b, Xty, eigval, Vtb_all, VtXty, diagVtMV,
                                var_y * 1e-8, sigma2_mom, tau2_mom);

        // ---- Step 3b: sa2 grid from tau2_0 / sigma2 (S3) ----
        if (!fix_sa2) {
            if (!(tau2_mom > 0.0)) {
                // no polygenic signal by MoM: skip Mr.ASH (theta = 0, tau2 = 0)
                sa2.resize(1); sa2(0) = 0.0;
                mix_pi.resize(1); mix_pi(0) = 1.0;
                theta.setZero();
                if (opt.estimate_residual_variance) sigma2 = sigma2_mom;
                tau2 = 0.0;
                omega_var = tau2 * eigval.array() + sigma2;
                pw        = Vsq * (eigval.array() / omega_var).matrix();
                XtOmegay  = Vmat * (VtXty.array() / omega_var).matrix();
                continue;
            }
            // upper bound: smallest prior variance among active sparse effects
            // (V > 0.01 * var(y); susieR used V > 0.01 for standardized y)
            double min_V = std::numeric_limits<double>::infinity();
            for (int l = 0; l < L; ++l) if (V(l) > 0.01 * var_y) min_V = std::min(min_V, V(l));
            if (!std::isfinite(min_V)) min_V = 0.01 * var_y;
            sa2 = ash_variance_grid(tau2_mom / sigma2_mom, min_V / sigma2_mom);
            if (!fix_pi) mix_pi = ArrayXd::Constant(sa2.size(), 1.0 / (double)sa2.size());
        }

        // ---- Step 3c: Mr.ASH to convergence on r = y - X b_bar ----
        double sigma2_ash = opt.estimate_residual_variance ? sigma2_mom : sigma2;
        const VectorXd Xtr0 = Xty - (p <= n ? (VectorXd)(d.XtX() * b) : (VectorXd)(X.transpose() * (X * b)));
        double rtr0 = 0.0;
        VectorXd r0;
        if (p <= n) rtr0 = yty - 2.0 * b.dot(Xty) + b.dot(Xty - Xtr0);   // b'X'X b = b'(X'y - X'r0)
        else        r0 = y - X * b;
        mr_ash_fit(d, Xtr0, rtr0, r0, sa2, mix_pi, !fix_pi,
                   sigma2_ash, opt.estimate_residual_variance,
                   theta, opt.ash_inner_max, opt.ash_inner_tol);

        // ---- Step 4: variance components from Mr.ASH; refresh Omega ----
        sigma2 = sigma2_ash;
        tau2 = sigma2 * (mix_pi * sa2).sum();
        omega_var = tau2 * eigval.array() + sigma2;
        pw        = Vsq * (eigval.array() / omega_var).matrix();
        XtOmegay  = Vmat * (VtXty.array() / omega_var).matrix();
    }

    VectorXd pip = reported_pip(alpha, V);

    VectorXd b_final = (alpha.array() * mu.array()).colwise().sum().transpose();
    VectorXd Xr = X * (b_final + theta);
    fit.alpha = alpha; fit.mu = mu; fit.mu2 = mu2;
    fit.lbf_variable = lbf_variable; fit.lbf = lbf;
    fit.pip = pip; fit.Xr = Xr; fit.fitted = Xr;
    fit.V = V; fit.sigma2 = sigma2; fit.intercept = 0.0;
    fit.tau2 = tau2; fit.theta = theta;
    fit.ash_sa2 = sa2.matrix();
    fit.ash_pi  = mix_pi.matrix();
    return fit;
}

// ---------------------------------------------------------------------------
// Credible sets
// ---------------------------------------------------------------------------

// Minimum absolute correlation among a set of design columns. For large sets the
// estimate is based on a random subsample of up to 100 columns (mirrors
// susieR's get_purity()).
static double set_purity(const std::vector<int32_t>& vars, SusieDesign& d) {
    const int k = (int)vars.size();
    if (k <= 1) return 1.0;

    std::vector<int32_t> use = vars;
    if (k > 100) {
        std::mt19937 rng(12345);
        std::shuffle(use.begin(), use.end(), rng);
        use.resize(100);
    }
    const int m = (int)use.size();
    // Pairwise inner products: looked up in the cached X'X when the design has
    // one, otherwise one symmetric rank-n update over the sampled columns (BLAS
    // dsyrk when available) rather than m(m-1)/2 separate length-n dot products.
    MatrixXd G(m, m);
    if (d.has_XtX()) {
        const MatrixXd& XtX = d.XtX();
        for (int j = 0; j < m; ++j)
            for (int i = j; i < m; ++i) G(i, j) = XtX(use[i], use[j]);
    } else {
        const MatrixXd& X = d.X();
        MatrixXd Xs(X.rows(), m);
        for (int i = 0; i < m; ++i) Xs.col(i) = X.col(use[i]);
        sym_crossprod_lower(Xs, true, G);            // m x m, lower triangle
    }
    double min_abs = 1.0;
    for (int j = 0; j < m; ++j) {
        for (int i = j + 1; i < m; ++i) {
            double denom_ij = std::sqrt(G(i, i) * G(j, j));
            double corr = (denom_ij > 0.0) ? (G(i, j) / denom_ij) : 0.0;
            double a = std::abs(corr);
            if (a < min_abs) min_abs = a;
        }
    }
    return min_abs;
}

std::vector<CredibleSet> susie_get_cs(const SusieFit& fit, SusieDesign& d,
                                      double coverage, double min_abs_corr) {
    std::vector<CredibleSet> out;
    const int L = (int)fit.alpha.rows();
    const int p = (int)fit.alpha.cols();
    if (L == 0 || p == 0) return out;

    std::vector<std::vector<int32_t> > kept_sets; // for de-duplication

    for (int l = 0; l < L; ++l) {
        // inactive effects (V ~ 0) have uniform alpha and never form a credible
        // set, even in a small high-LD region where purity alone would pass
        // (susieR's susie_get_cs skips V <= prior_tol the same way)
        if (fit.V.size() == L && !(fit.V(l) > SUSIE_PRIOR_TOL)) continue;

        // order variables by descending alpha within this effect
        std::vector<int32_t> order(p);
        std::iota(order.begin(), order.end(), 0);
        const double* arow = fit.alpha.row(l).data();
        std::sort(order.begin(), order.end(),
                  [&](int32_t a, int32_t b) { return fit.alpha(l, a) > fit.alpha(l, b); });

        // accumulate until cumulative alpha reaches the requested coverage
        double cum = 0.0;
        std::vector<int32_t> vars;
        for (int idx = 0; idx < p; ++idx) {
            int32_t j = order[idx];
            vars.push_back(j);
            cum += fit.alpha(l, j);
            if (cum >= coverage) break;
        }
        (void)arow;

        // an effect that explains nothing spreads alpha ~ uniformly; its set will be
        // large with low purity and is filtered out below.
        double purity = set_purity(vars, d);
        if (purity < min_abs_corr) continue;

        // de-duplicate: skip if an identical set was already reported
        std::vector<int32_t> sorted_vars = vars;
        std::sort(sorted_vars.begin(), sorted_vars.end());
        bool dup = false;
        for (const auto& ks : kept_sets) if (ks == sorted_vars) { dup = true; break; }
        if (dup) continue;
        kept_sets.push_back(sorted_vars);

        CredibleSet cs;
        cs.effect_index = l;
        cs.variables = vars;
        cs.pips.reserve(vars.size());
        for (int32_t j : vars) cs.pips.push_back(fit.pip(j));
        cs.coverage = cum;
        cs.requested_coverage = coverage;
        cs.purity = purity;
        cs.lbf = fit.lbf(l);
        cs.V = fit.V(l);
        out.push_back(std::move(cs));
    }
    return out;
}

// ---------------------------------------------------------------------------
// High-level wrapper
// ---------------------------------------------------------------------------

SusieResult fit_susie_trait(SusieDesign& d, const VectorXd& pheno_vec, const VectorXd& Xty,
                            const SusieOptions& opt, double coverage, double min_abs_corr) {
    if (pheno_vec.size() != d.n() || Xty.size() != d.p())
        error("fit_susie_trait(): dimension mismatch (y has %d rows, X'y has %d rows; design is %d x %d)",
              (int)pheno_vec.size(), (int)Xty.size(), d.n(), d.p());
    // center the phenotype (covariates already regressed out by the caller);
    // X'y is unchanged because the design columns are centered
    VectorXd y = pheno_vec.array() - pheno_vec.mean();
    if (!y.allFinite())
        error("fit_susie_trait(): the phenotype vector contains non-finite values");

    SusieResult res;
    if (opt.unmappable_effects == SusieOptions::INF)
        res.fit = fit_susie_inf(d, y, Xty, opt);
    else if (opt.unmappable_effects == SusieOptions::ASH)
        res.fit = fit_susie_ash(d, y, Xty, opt);
    else
        res.fit = fit_susie(d, y, Xty, opt);
    res.cs = susie_get_cs(res.fit, d, coverage, min_abs_corr);
    return res;
}

SusieResult simple_susie_without_missing(const VectorXd& pheno_vec,
                                         const MatrixXd& geno_mat,
                                         const SusieOptions& opt,
                                         double coverage,
                                         double min_abs_corr) {
    SusieDesign d(MatrixXd(geno_mat), opt.standardize);
    VectorXd Xty = d.X().transpose() * pheno_vec;
    return fit_susie_trait(d, pheno_vec, Xty, opt, coverage, min_abs_corr);
}

} // namespace susie
