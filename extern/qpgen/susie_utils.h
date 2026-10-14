#ifndef __SUSIE_UTILS_H
#define __SUSIE_UTILS_H

// SuSiE ("Sum of Single Effects") fine-mapping via the IBSS algorithm, in Eigen.
//
// Core algorithm reimplemented from:
//   G. Wang, A. Sarkar, P. Carbonetto, M. Stephens (2020),
//   "A simple new approach to variable selection in regression, with application
//    to genetic fine mapping", JRSS-B 82(5):1273-1300.
// Output field names mirror susieR's susie() fit object (alpha, mu, mu2,
// lbf_variable, lbf, pip, V, sigma2, elbo, ...) plus credible sets.

#include "Eigen/Dense"
#include <vector>
#include <cstdint>

namespace susie {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

// ---- Fit object, mirroring susieR's susie() ------------------------------
struct SusieFit {
    MatrixXd alpha;        // L x p   posterior inclusion prob within each single effect
    MatrixXd mu;           // L x p   posterior mean | inclusion
    MatrixXd mu2;          // L x p   posterior 2nd moment | inclusion
    MatrixXd lbf_variable; // L x p   log Bayes factor per variable & effect
    VectorXd lbf;          // L       log Bayes factor per single effect
    VectorXd pip;          // p       marginal PIP = 1 - prod_l(1 - alpha_lj)
    VectorXd Xr;           // n       X %*% colSums(alpha*mu)
    VectorXd fitted;       // n       fitted values (= Xr; intercept = 0)
    VectorXd V;            // L       prior variance of each single effect
    double   sigma2;       // residual variance
    double   intercept;    // intercept (0 when data pre-centered)
    std::vector<double> elbo; // objective at each outer iteration
    int      niter;
    bool     converged;
    // SuSiE-inf / SuSiE-ash (unmappable effects) extras.
    //   tau2   : infinitesimal variance (inf); or sum(pi_k*sa2_k)*sigma2 (ash, for reporting)
    //   theta  : p-vector of posterior mean unmappable effects (standardized-X scale)
    //   ash_sa2: K component prior variances on the sa2 grid (ash only; empty otherwise)
    //   ash_pi : K mixture weights (ash only; sum to 1, first entry is null component)
    // For standard SuSiE and inf, ash_sa2/ash_pi are empty and theta may be empty.
    double   tau2 = 0.0;
    VectorXd theta;
    VectorXd ash_sa2;
    VectorXd ash_pi;
};

struct SusieOptions {
    int    L            = 10;
    int    max_iter     = 100;
    double tol          = 1e-3;
    bool   estimate_prior_variance    = true;
    bool   estimate_residual_variance = true;
    double scaled_prior_variance = 0.2;   // prior var = var(y)*this (init/fixed)
    double prior_v_min  = 1e-9;
    double prior_v_max  = 1e3;
    bool   standardize  = true;           // scale X columns to unit variance (SusieDesign)
    enum ConvergenceMethod { ELBO, PIP } convergence_method = ELBO;
    // Unmappable-effects model.
    //   NONE : standard SuSiE.
    //   INF  : SuSiE-inf; theta_j ~ N(0, tau2) with MoM variance components.
    //   ASH  : SuSiE-ash (SuSiE 2.0 manuscript, Algorithm 1); theta_j ~ sum_k pi_k
    //          N(0, sigma2 * sa2_k) on a data-driven grid, fit by Mr.ASH, with the
    //          SER run under Omega = (tau2 XX' + sigma2 I)^{-1}.
    enum UnmappableEffects { NONE, INF, ASH } unmappable_effects = NONE;
    // ASH-only knobs (ignored otherwise): Mr.ASH inner fit, run to convergence in
    // every outer iteration (mr.ash defaults: relative theta change 1e-4).
    int    ash_inner_max    = 1000;      // max coordinate-ascent sweeps per outer iteration
    double ash_inner_tol    = 1e-4;      // relative L2 change of theta
    // Reduction-test hooks.
    //   ash_fix_pi : hold pi at this vector (no EM update of pi); must match
    //                the grid length, so give ash_fix_sa2 as well.
    //   ash_fix_sa2: hold the sa2 grid (first entry should be 0) instead of the
    //                data-driven grid rebuilt every iteration.
    std::vector<double> ash_fix_pi;
    std::vector<double> ash_fix_sa2;
};

// A 95% (or requested coverage) credible set for a single effect.
struct CredibleSet {
    int32_t effect_index;            // 0-based single-effect index (l)
    std::vector<int32_t> variables;  // 0-based variant indices, ordered by descending alpha
    std::vector<double>  pips;        // marginal PIP of each member (parallel to variables)
    double coverage;                 // achieved cumulative alpha within the set
    double requested_coverage;
    double purity;                   // minimum absolute pairwise correlation among members
    double lbf;                      // log Bayes factor of this single effect
    double V;                        // estimated prior variance of this single effect
};

// Thin eigendecomposition X = U D V' of a design (variant-space part only).
struct ThinEigen {
    VectorXd eigval;   // r      d_k^2 (components with ~0 singular value dropped)
    MatrixXd Vmat;     // p x r  variant-space eigenvectors (susieR's svd$v)
};

// Region-level design shared by every trait analyzed on the SAME samples.
//
// All trait-independent work lives here and is done once per region, not once
// per trait: centering/standardizing X, screening unusable columns, X'X, and
// the eigendecomposition used by SuSiE-inf/ash. X'X and the eigendecomposition
// are built on first use and cached. Traits then only need X'y (for a block of
// traits, one GEMM via crossprod()) and their own IBSS iterations.
//
// Traits with different missingness patterns have different sample sets (and,
// after covariate residualization, different X), so they cannot share a design;
// region-assoc therefore allows missing phenotype values only for a single trait.
class SusieDesign {
public:
    // Takes ownership of X (n x p), e.g. std::move(geno_mat). In place, every
    // column is centered and, if standardize, scaled to unit variance (sd with
    // n-1). Columns with non-finite values or zero variance are zeroed and
    // reported once via notice(); they get zero prior weight in every fit.
    SusieDesign(MatrixXd&& X, bool standardize);

    int n() const { return (int)X_.rows(); }
    int p() const { return (int)X_.cols(); }
    const MatrixXd& X() const { return X_; }
    // diag(X'X) of the design (0 for unusable columns)
    const ArrayXd& xtx() const { return xtx_; }
    // input column j (centered) = design column j * col_scale(j): the sd when
    // standardized, 1 otherwise, 0 for unusable columns. Converts design-scale
    // statistics back to the per-allele scale (X'y, x'x, X'X).
    const VectorXd& col_scale() const { return col_scale_; }
    int n_unusable() const { return n_unusable_; }

    // X'Y (p x K) for a block of traits sharing these samples: one GEMM.
    MatrixXd crossprod(const MatrixXd& Y) const;

    // Cached p x p X'X of the design (full symmetric). Built on first call.
    const MatrixXd& XtX();
    bool has_XtX() const { return has_XtX_; }

    // Cached thin eigendecomposition of the design, taken on the smaller of
    // X'X (reusing the cached one) or XX'. Built on first call.
    const ThinEigen& thin_eigen();

private:
    MatrixXd X_;
    ArrayXd  xtx_;
    VectorXd col_scale_;
    int      n_unusable_ = 0;
    bool     has_XtX_ = false;
    MatrixXd XtX_;
    bool     has_eigen_ = false;
    ThinEigen eigen_;
};

// Per-trait fits on a shared design. y is the centered phenotype (no missing
// values, same samples as the design) and Xty = d.X()' y (see crossprod()).
//
// Standard SuSiE: when p <= n it runs on the design's cached X'X (sufficient
// statistics, like susieR's susie_ss); otherwise on X directly.
SusieFit fit_susie(SusieDesign& d, const VectorXd& y, const VectorXd& Xty,
                   const SusieOptions& opt = SusieOptions());

// SuSiE-inf (unmappable_effects = "inf"). Adds an infinitesimal effect via the
// eigenspace Omega = (tau2 XX' + sigma2 I)^{-1} formulation, estimating
// (sigma2, tau2) by method of moments and theta as its BLUP. Uses the design's
// cached eigendecomposition.
SusieFit fit_susie_inf(SusieDesign& d, const VectorXd& y, const VectorXd& Xty,
                       const SusieOptions& opt = SusieOptions());

// SuSiE-ash (simplified port of susieR unmappable_effects = "ash"). theta_j
// has a scale-mixture-of-normals prior sum_k pi_k * N(0, sa2_k * sigma2) fit by
// Mr.ASH coordinate ascent. The grid is fixed and log-spaced (first component
// is the null point mass); pi and sigma2 are estimated by EM. Skips susieR's
// LD-masking and slot-activity heuristics for speed.
SusieFit fit_susie_ash(SusieDesign& d, const VectorXd& y, const VectorXd& Xty,
                       const SusieOptions& opt = SusieOptions());

// Extract credible sets from a fit on design d. Purity (minimum absolute
// correlation) is read from the cached X'X when available, otherwise computed
// from the sampled columns of X. Sets whose purity < min_abs_corr are dropped;
// duplicate sets are removed.
std::vector<CredibleSet> susie_get_cs(const SusieFit& fit, SusieDesign& d,
                                      double coverage = 0.95, double min_abs_corr = 0.5);

// Bundle of everything a caller typically wants.
struct SusieResult {
    SusieFit fit;                     // includes pip, lbf, V, sigma2, elbo, niter, converged
    std::vector<CredibleSet> cs;
};

// Fit one trait on a shared design with the model chosen by
// opt.unmappable_effects, and extract its credible sets. pheno_vec is centered
// internally (Xty is unaffected, as the design columns are centered).
SusieResult fit_susie_trait(SusieDesign& d, const VectorXd& pheno_vec, const VectorXd& Xty,
                            const SusieOptions& opt = SusieOptions(),
                            double coverage = 0.95, double min_abs_corr = 0.5);

// Single-trait convenience entry point: builds a design from a copy of geno_mat
// (n x p, e.g. mean-centered dosages) and fits pheno_vec (covariate-adjusted,
// non-missing). For several traits on the same samples, build one SusieDesign
// and call fit_susie_trait() per trait instead.
SusieResult simple_susie_without_missing(const VectorXd& pheno_vec,
                                         const MatrixXd& geno_mat,
                                         const SusieOptions& opt = SusieOptions(),
                                         double coverage = 0.95,
                                         double min_abs_corr = 0.5);

} // namespace susie

#endif // __SUSIE_UTILS_H
