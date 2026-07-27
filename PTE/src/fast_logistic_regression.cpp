#include "_helper_functions.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

ModelResult fast_logistic_regression_internal(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                              const Eigen::Ref<const Eigen::VectorXd>& y,
                                              const Eigen::Ref<const Eigen::VectorXd>& weights,
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta,
                                              bool smart_cold_start,
                                              int maxit,
                                              double tol,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values,
                                              std::string optimization_alg,
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights,
                                              Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info,
                                              bool estimate_only);

namespace {

inline double plogis_manual(double x) {
    if (x > 20.0) return 1.0;
    if (x < -20.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}

inline double log1pexp_stable(double x) {
    return (x > 0.0) ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

template<typename RDerived, typename WDerived>
inline void score_weighted_crossprod_colwise_assign(const Eigen::MatrixXd& X,
                                                    const Eigen::MatrixBase<RDerived>& residual,
                                                    const Eigen::MatrixBase<WDerived>& w,
                                                    Eigen::VectorXd& score,
                                                    Eigen::MatrixXd& out) {
    const int n = X.rows();
    if (residual.rows() != n || w.rows() != n) {
        Rcpp::stop("score_weighted_crossprod_colwise_assign: vector has incompatible dimensions");
    }
    score.noalias() = X.transpose() * residual;
    out = weighted_crossprod(X, w);
}

class LogisticLbfgsObjective : public Numer::MFuncGrad {
private:
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    const Eigen::Ref<const Eigen::VectorXd> m_y;
    const Eigen::Ref<const Eigen::VectorXd> m_weights;
    const Eigen::Ref<const Eigen::VectorXd> m_eta_fixed;
    bool m_use_weights;
    int m_n;

public:
    LogisticLbfgsObjective(const Eigen::Ref<const Eigen::MatrixXd>& X, const Eigen::Ref<const Eigen::VectorXd>& y, 
                           const Eigen::Ref<const Eigen::VectorXd>& weights, const Eigen::Ref<const Eigen::VectorXd>& eta_fixed, 
                           bool use_weights)
        : m_X(X), m_y(y), m_weights(weights), m_eta_fixed(eta_fixed), 
          m_use_weights(use_weights), m_n(X.rows()) {}

    virtual double f_grad(Numer::Constvec& beta, Numer::Refvec grad) override {
        Eigen::VectorXd eta = m_eta_fixed + m_X * beta;
        Eigen::ArrayXd mu = plogis_array_safe(eta.array());
        Eigen::ArrayXd neg_ll_terms = log1pexp_array_safe(eta.array()) - m_y.array() * eta.array();
        Eigen::VectorXd diff = (mu - m_y.array()).matrix();
        if (m_use_weights) {
            neg_ll_terms *= m_weights.array();
            diff.array() *= m_weights.array();
        }
        const double neg_ll = neg_ll_terms.sum();
        grad.noalias() = m_X.transpose() * diff;
        return neg_ll;
    }
};

inline bool fit_logistic_counted_irls(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                      const Eigen::Ref<const Eigen::VectorXd>& y,
                                      const Eigen::Ref<const Eigen::VectorXd>& counts,
                                      Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta,
                                      Eigen::VectorXd& beta,
                                      int maxit = 100,
                                      double tol = 1e-8) {
    const int n = X.rows();
    const int p = X.cols();
    beta = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        beta = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_beta));
        if (beta.size() != p) stop("warm_start_beta must have length equal to ncol(X)");
    }

    Eigen::VectorXd score(p);
    Eigen::MatrixXd XtWX(p, p);
    Eigen::VectorXd eta(n);
    Eigen::VectorXd mu(n);
    Eigen::VectorXd w(n);
    Eigen::VectorXd diff(n);

    for (int iter = 0; iter < maxit; ++iter) {
        eta.noalias() = X * beta;
        mu = plogis_array_safe(eta.array()).matrix();
        diff.array() = counts.array() * (y.array() - mu.array());
        w.array() = counts.array() * (mu.array() * (1.0 - mu.array())).max(1e-10);
        score.noalias() = X.transpose() * diff;
        XtWX = weighted_crossprod(X, w);

        if (score.norm() < tol) return beta.allFinite();

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) return false;
        Eigen::VectorXd delta = ldlt.solve(score);
        if (!delta.allFinite()) return false;

        beta += delta;
        if (delta.norm() < tol) return beta.allFinite();
    }

    return false;
}

inline Rcpp::NumericVector numeric_vector_from_eigen(const Eigen::VectorXd& x) {
    Rcpp::NumericVector out(x.size());
    for (int i = 0; i < x.size(); ++i) out[i] = x[i];
    return out;
}

inline bool fit_logistic_counted_active_irls(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                             const Eigen::Ref<const Eigen::VectorXd>& y,
                                             const Eigen::Ref<const Eigen::VectorXd>& counts,
                                             const std::vector<int>& active_rows,
                                             const Eigen::VectorXd& start_beta,
                                             Eigen::VectorXd& beta,
                                             int maxit = 100,
                                             double tol = 1e-8) {
    const int p = X.cols();
    beta = start_beta;
    if (beta.size() != p || !beta.allFinite()) beta = Eigen::VectorXd::Zero(p);

    Eigen::VectorXd score(p);
    Eigen::MatrixXd XtWX(p, p);
    bool converged = false;

    for (int iter = 0; iter < maxit; ++iter) {
        score.setZero();
        XtWX.setZero();

        for (const int row : active_rows) {
            const double ci = counts[row];
            if (ci <= 0.0) continue;
            const double eta = X.row(row).dot(beta);
            const double mu = plogis_manual(eta);
            const double diff = ci * (y[row] - mu);
            const double wi = ci * std::max(mu * (1.0 - mu), 1e-10);

            score.noalias() += X.row(row).transpose() * diff;
            XtWX.selfadjointView<Eigen::Upper>().rankUpdate(X.row(row).transpose(), wi);
        }
        XtWX.triangularView<Eigen::Lower>() = XtWX.transpose();

        if (score.norm() < tol) {
            converged = beta.allFinite();
            break;
        }

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) {
            for (int attempt = 0; attempt < 4 && ldlt.info() != Eigen::Success; ++attempt) {
                Eigen::MatrixXd jittered = XtWX;
                jittered.diagonal().array() += std::pow(10.0, -8 + attempt);
                ldlt.compute(jittered);
            }
            if (ldlt.info() != Eigen::Success) break;
        }
        Eigen::VectorXd delta = ldlt.solve(score);
        if (!delta.allFinite()) break;

        beta += delta;
        if (delta.norm() < tol) {
            converged = beta.allFinite();
            break;
        }
    }

    if (converged && beta.allFinite()) return true;

    ModelResult lbfgs_res = fast_logistic_regression_internal(
        X, y, counts, numeric_vector_from_eigen(beta), false, 200, tol,
        R_NilValue, R_NilValue, "lbfgs", R_NilValue, R_NilValue, true);
    if (lbfgs_res.converged && lbfgs_res.b.allFinite()) {
        beta = lbfgs_res.b;
        return true;
    }
    return false;
}

// A (quasi-)separated fold has no real MLE: the log-likelihood keeps improving as some
// coefficient runs off to +-Inf, so IRLS "converges" (score norm < tol) at an essentially
// arbitrary point on that ridge that depends on the starting value. The batch CV/bootstrap
// functions below chain each fold's fit into the next fold's starting point for speed, so an
// undetected separated fold doesn't just report a locally-arbitrary answer, it also pushes every
// subsequent fold's starting point further out, compounding across folds. A fitted linear
// predictor this large (predicted probability within ~1e-7 of 0 or 1) is itself the signature of
// separation for any dataset whose genuinely-identified coefficients are of reasonable magnitude
// -- callers should treat exceeding this the same as non-convergence.
inline double active_rows_max_abs_eta(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                      const std::vector<int>& active_rows,
                                      const Eigen::VectorXd& beta) {
    double max_abs_eta = 0.0;
    for (const int row : active_rows) {
        const double eta = std::abs(X.row(row).dot(beta));
        if (eta > max_abs_eta) max_abs_eta = eta;
    }
    return max_abs_eta;
}

} // namespace

// Internal pure C++ logic
ModelResult fast_logistic_regression_internal(const Eigen::Ref<const Eigen::MatrixXd>& X, 
                                              const Eigen::Ref<const Eigen::VectorXd>& y, 
                                              const Eigen::Ref<const Eigen::VectorXd>& weights,
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                              bool smart_cold_start = false,
                                              int maxit = 100, 
                                              double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls",
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                              bool estimate_only = false) {
    const int n = X.rows();
    const int p = X.cols();
    const bool use_weights = (weights.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    const int p_free = fixed_spec.free_idx.size();
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        beta = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_beta));
    } else if (smart_cold_start) {
        beta = edi_opt::logistic_smart_cold_start(X, y);
    }
    beta = apply_fixed_values(beta, fixed_spec);
    
    Eigen::VectorXd beta_free = subset_vector(beta, fixed_spec.free_idx);
    Eigen::VectorXd eta_fixed = Eigen::VectorXd::Zero(n);
    for(size_t k=0; k<fixed_spec.fixed_idx.size(); k++) {
        eta_fixed.noalias() += X.col(fixed_spec.fixed_idx[k]) * fixed_spec.fixed_values[k];
    }

    if (optimization_alg == "lbfgs") {
        Eigen::MatrixXd X_free(n, p_free);
        for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);

        bool converged = true;
        double fopt = NA_REAL;
        if (p_free > 0) {
            LogisticLbfgsObjective nll(X_free, y, weights, eta_fixed, use_weights);
            int status = Numer::optim_lbfgs(nll, beta_free, fopt, maxit, tol, tol);
            converged = (status >= 0) && beta_free.allFinite();
        }

        ModelResult res;
        res.b = expand_free_params(beta_free, beta, fixed_spec);
        res.neg_ll = fopt;
        if (!estimate_only) {
            Eigen::VectorXd eta = X * res.b;
            res.mu = plogis_array_safe(eta.array()).matrix();
            Eigen::VectorXd w_diag = res.mu.array() * (1.0 - res.mu.array());
            if (use_weights) w_diag.array() *= weights.array();
            w_diag.array() = w_diag.array().max(1e-10);
            res.XtWX = expand_free_covariance(p, fixed_spec, weighted_crossprod(X_free, w_diag), false);
        }
        res.iterations = NA_INTEGER;
        res.converged = converged;
        return res;
    }

    // IRLS Path
    Eigen::MatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);

    Eigen::VectorXd mu(n);
    Eigen::VectorXd w(n);
    Eigen::VectorXd eta(n);
    Eigen::MatrixXd XtWX(p_free, p_free);
    Eigen::VectorXd score_free(p_free);
    Eigen::VectorXd diff(n);
    bool converged = false;
    int iterations = 0;

    for (int iter = 0; iter < maxit; iter++) {
        iterations++;
        eta.noalias() = eta_fixed;
        eta.noalias() += X_free * beta_free;

        mu = plogis_array_safe(eta.array()).matrix();

        if (iter == 0 && warm_start_weights.isNotNull()) {
            Eigen::VectorXd ww = as<Eigen::VectorXd>(warm_start_weights);
            if (ww.size() == n) w = ww;
            else w.array() = mu.array() * (1.0 - mu.array());
        } else {
            w.array() = mu.array() * (1.0 - mu.array());
        }
        if (use_weights) w.array() *= weights.array();
        w.array() = w.array().max(1e-10);
        diff.array() = y.array() - mu.array();
        if (use_weights) diff.array() *= weights.array();

        const bool use_warm_xtwx = (iter == 0) &&
            (warm_start_fisher_info.isNotNull() || (smart_cold_start && warm_start_beta.isNull()));
        if (!use_warm_xtwx) {
            score_weighted_crossprod_colwise_assign(X_free, diff, w, score_free, XtWX);
        } else {
            score_free.noalias() = X_free.transpose() * diff;
            if (warm_start_fisher_info.isNotNull()) {
                XtWX = subset_matrix(as<Eigen::MatrixXd>(warm_start_fisher_info), fixed_spec.free_idx, fixed_spec.free_idx);
            } else {
                XtWX = subset_matrix(edi_opt::logistic_smart_hessian(X, beta), fixed_spec.free_idx, fixed_spec.free_idx);
            }
        }

        if (score_free.norm() < tol) { converged = true; break; }

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) break;
        Eigen::VectorXd delta = ldlt.solve(score_free);
        if (!delta.allFinite()) break;

        beta_free += delta;
        if (delta.norm() < tol) { converged = true; break; }
    }

    ModelResult res;
    res.b = expand_free_params(beta_free, beta, fixed_spec);
    if (!estimate_only) {
        res.mu = mu;
        res.XtWX = expand_free_covariance(p, fixed_spec, XtWX, false);
        Eigen::VectorXd final_eta = eta_fixed + X_free * beta_free;
        Eigen::ArrayXd neg_ll_terms = log1pexp_array_safe(final_eta.array()) - y.array() * final_eta.array();
        if (use_weights) neg_ll_terms *= weights.array();
        res.neg_ll = neg_ll_terms.sum();
        Eigen::VectorXd final_diff = y - mu;
        if (use_weights) final_diff.array() *= weights.array();
        res.score = X.transpose() * final_diff;
    }
    res.iterations = iterations;
    res.converged = converged;
    return res;
}

// [[Rcpp::export]]
Eigen::VectorXd get_logistic_regression_score_cpp(SEXP X_sexp, SEXP y_sexp, SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());

    const int n = X.rows();
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu(n);
    mu = plogis_array_safe(eta.array()).matrix();
    return X.transpose() * (y - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_regression_hessian_cpp(SEXP X_sexp, SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());

    const int n = X.rows();
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd w(n);
    w = plogis_array_safe(eta.array()).matrix();
    w.array() = w.array() * (1.0 - w.array());
    return -weighted_crossprod(X, w);
}

// [[Rcpp::export]]
Eigen::VectorXd get_logistic_regression_weighted_score_cpp(SEXP X_sexp, SEXP y_sexp, SEXP weights_sexp, SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector w_r(weights_sexp);
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> weights(w_r.begin(), w_r.size());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());

    const int n = X.rows();
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu(n);
    mu = plogis_array_safe(eta.array()).matrix();
    Eigen::VectorXd diff = y - mu;
    diff.array() *= weights.array();
    return X.transpose() * diff;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_regression_weighted_hessian_cpp(SEXP X_sexp, SEXP weights_sexp, SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector w_r(weights_sexp);
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> weights(w_r.begin(), w_r.size());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());

    const int n = X.rows();
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd w(n);
    w = plogis_array_safe(eta.array()).matrix();
    w.array() = w.array() * (1.0 - w.array()) * weights.array();
    return -weighted_crossprod(X, w);
}

// [[Rcpp::export]]
List fast_logistic_regression_cpp(SEXP X_sexp, SEXP y_sexp,
                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                  bool smart_cold_start = false,
                                  int maxit = 100, double tol = 1e-8,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "irls",
                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                  bool estimate_only = false) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
    
    if (estimate_only) {
        return List::create(
            Named("b") = res.b,
            Named("converged") = res.converged,
            Named("iterations") = res.iterations
        );
    }
    Eigen::VectorXd weights_vec = res.mu.array() * (1.0 - res.mu.array());
    return List::create(
        Named("b") = res.b,
        Named("w") = weights_vec,
        Named("iterations") = res.iterations,
        Named("fisher_information") = res.XtWX,
        Named("score") = res.score,
        Named("neg_ll") = res.neg_ll,
        Named("converged") = res.converged
    );
}

// [[Rcpp::export]]
List fast_logistic_regression_weighted_cpp(SEXP X_sexp, SEXP y_sexp, SEXP weights_sexp,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                           bool smart_cold_start = false,
                                           int maxit = 100, double tol = 1e-8,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls",
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector w_r(weights_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> weights(w_r.begin(), w_r.size());
    ModelResult res = fast_logistic_regression_internal(X, y, weights, warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    return List::create(
        Named("b") = res.b,
        Named("mu") = res.mu,
        Named("XtWX") = res.XtWX,
        Named("fisher_information") = res.XtWX,
        Named("score") = res.score,
        Named("neg_ll") = res.neg_ll,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations
    );
}

// [[Rcpp::export]]
List fast_logistic_regression_with_var_cpp(SEXP X_sexp, SEXP y_sexp, int j = 2,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                           bool smart_cold_start = false,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls",
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, 100, 1e-8, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    
    Eigen::MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int overall_j) -> int {
        for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
            if (fixed_spec.free_idx[jj] == overall_j) return jj + 1; // 1-based for compute_diagonal_inverse_entry
        return -1;
    };

    int free_j = (j > 0 && j <= X.cols()) ? free_idx_of(j - 1) : -1;
    res.ssq_b_j = (free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;

    int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1;
    res.ssq_b_2 = (free_2 > 0) ? compute_diagonal_inverse_entry(info_free, free_2) : NA_REAL;

    return List::create(
        Named("b") = res.b,
        Named("params") = res.b,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2,
        Named("score") = res.score,
        Named("observed_information") = res.XtWX,
        Named("fisher_information") = res.XtWX,
        Named("information") = res.XtWX,
        Named("information_type") = "fisher",
        Named("hessian") = -res.XtWX,
        Named("neg_loglik") = res.neg_ll,
        Named("neg_ll") = res.neg_ll,
        Named("loglik") = R_finite(res.neg_ll) ? -res.neg_ll : NA_REAL,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations
    );
}

// [[Rcpp::export]]
Rcpp::List fast_default_incidence_cv_run_cpp(SEXP X_obs_sexp,
                                             SEXP y_full_sexp,
                                             SEXP X_tx0_sexp,
                                             SEXP X_tx1_sexp,
                                             SEXP treatment_sexp,
                                             SEXP censored_sexp,
                                             SEXP boot_idx_sexp,
                                             SEXP begin_cutoffs_sexp,
                                             SEXP end_cutoffs_sexp,
                                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                             bool y_higher_is_better = true) {
    NumericMatrix X_obs_r(X_obs_sexp);
    NumericVector y_full_r(y_full_sexp);
    NumericMatrix X_tx0_r(X_tx0_sexp);
    NumericMatrix X_tx1_r(X_tx1_sexp);
    NumericVector treatment_r(treatment_sexp);
    NumericVector censored_r(censored_sexp);
    IntegerVector boot_idx_r(boot_idx_sexp);
    IntegerVector begin_r(begin_cutoffs_sexp);
    IntegerVector end_r(end_cutoffs_sexp);

    Eigen::Map<const Eigen::MatrixXd> X_obs(X_obs_r.begin(), X_obs_r.nrow(), X_obs_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y_full(y_full_r.begin(), y_full_r.size());
    Eigen::Map<const Eigen::MatrixXd> X_tx0(X_tx0_r.begin(), X_tx0_r.nrow(), X_tx0_r.ncol());
    Eigen::Map<const Eigen::MatrixXd> X_tx1(X_tx1_r.begin(), X_tx1_r.nrow(), X_tx1_r.ncol());

    const int n_full = X_obs.rows();
    const int p = X_obs.cols();
    const int n_out = boot_idx_r.size();
    const int n_folds = begin_r.size();
    if (end_r.size() != n_folds) stop("begin_cutoffs and end_cutoffs must have the same length");

    NumericVector est_true(n_out, NA_REAL);
    NumericVector est_counterfactual(n_out, NA_REAL);
    NumericVector given_tx(n_out, NA_REAL);
    NumericVector rec_tx(n_out, NA_REAL);
    NumericVector real_y(n_out, NA_REAL);
    NumericVector censored(n_out, NA_REAL);

    std::vector<int> rows(n_out);
    Eigen::VectorXd full_counts = Eigen::VectorXd::Zero(n_full);
    for (int pos = 0; pos < n_out; ++pos) {
        const int row = boot_idx_r[pos] - 1;
        if (row < 0 || row >= n_full) return List::create(Named("ok") = false);
        rows[pos] = row;
        full_counts[row] += 1.0;
    }

    Eigen::VectorXd fold_start = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        fold_start = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_beta));
        if (fold_start.size() != p) stop("warm_start_beta must have length equal to ncol(X)");
    }

    for (int fold = 0; fold < n_folds; ++fold) {
        const int begin = begin_r[fold] - 1;
        const int end = end_r[fold] - 1;
        if (begin < 0 || end < begin || end >= n_out) return List::create(Named("ok") = false);

        Eigen::VectorXd counts = full_counts;
        for (int pos = begin; pos <= end; ++pos) counts[rows[pos]] -= 1.0;

        std::vector<int> active_rows;
        active_rows.reserve(n_full);
        for (int row = 0; row < n_full; ++row) {
            if (counts[row] > 0.0) active_rows.push_back(row);
        }

        Eigen::VectorXd beta;
        const bool converged = fit_logistic_counted_active_irls(
            X_obs, y_full, counts, active_rows, fold_start, beta, 100, 1e-8);
        if (!converged || !beta.allFinite()) return List::create(Named("ok") = false);
        if (active_rows_max_abs_eta(X_obs, active_rows, beta) > 30.0) return List::create(Named("ok") = false);
        //fold_start is intentionally NOT updated to beta here: every fold starts from the same
        //fixed full-data warm start (see active_rows_max_abs_eta's comment above) so one fold's
        //fit can never drift into the next fold's starting point and compound across folds.

        for (int pos = begin; pos <= end; ++pos) {
            const int row = rows[pos];
            const double yhat0 = X_tx0.row(row).dot(beta);
            const double yhat1 = X_tx1.row(row).dot(beta);
            const double tx = treatment_r[row];
            const double et = yhat0 + tx * (yhat1 - yhat0);
            const double ec = yhat0 + yhat1 - et;
            const bool optimal = y_higher_is_better ? (et > ec) : (et < ec);

            est_true[pos] = et;
            est_counterfactual[pos] = ec;
            given_tx[pos] = tx;
            rec_tx[pos] = tx * optimal + (1.0 - tx) * (1.0 - optimal);
            real_y[pos] = y_full[row];
            censored[pos] = censored_r[row];
        }
    }

    return List::create(
        Named("ok") = true,
        Named("raw_results") = List::create(
            Named("est_true") = est_true,
            Named("est_counterfactual") = est_counterfactual,
            Named("given_tx") = given_tx,
            Named("rec_tx") = rec_tx,
            Named("real_y") = real_y,
            Named("censored") = censored
        )
    );
}

// [[Rcpp::export]]
Rcpp::List fast_default_incidence_bootstrap_q_cpp(SEXP X_obs_sexp,
                                                  SEXP y_full_sexp,
                                                  SEXP X_tx0_sexp,
                                                  SEXP X_tx1_sexp,
                                                  SEXP treatment_sexp,
                                                  int B,
                                                  SEXP begin_cutoffs_sexp,
                                                  SEXP end_cutoffs_sexp,
                                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                                  bool y_higher_is_better = true,
                                                  std::string incidence_metric = "odds_ratio") {
    Rcpp::RNGScope rng_scope;
    NumericMatrix X_obs_r(X_obs_sexp);
    NumericVector y_full_r(y_full_sexp);
    NumericMatrix X_tx0_r(X_tx0_sexp);
    NumericMatrix X_tx1_r(X_tx1_sexp);
    NumericVector treatment_r(treatment_sexp);
    IntegerVector begin_r(begin_cutoffs_sexp);
    IntegerVector end_r(end_cutoffs_sexp);

    Eigen::Map<const Eigen::MatrixXd> X_obs(X_obs_r.begin(), X_obs_r.nrow(), X_obs_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y_full(y_full_r.begin(), y_full_r.size());
    Eigen::Map<const Eigen::MatrixXd> X_tx0(X_tx0_r.begin(), X_tx0_r.nrow(), X_tx0_r.ncol());
    Eigen::Map<const Eigen::MatrixXd> X_tx1(X_tx1_r.begin(), X_tx1_r.nrow(), X_tx1_r.ncol());

    const int n_full = X_obs.rows();
    const int p = X_obs.cols();
    const int n_out = (end_r.size() > 0) ? end_r[end_r.size() - 1] : 0;
    const int n_folds = begin_r.size();
    if (B < 0 || n_out <= 0 || end_r.size() != n_folds) stop("invalid bootstrap q inputs");

    Eigen::VectorXd global_start = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        global_start = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_beta));
        if (global_start.size() != p) stop("warm_start_beta must have length equal to ncol(X)");
    }

    const int metric_code = incidence_metric == "probability_difference" ? 0 :
        (incidence_metric == "risk_ratio" ? 1 : 2);
    NumericVector q_adv(B, NA_REAL), q_avg(B, NA_REAL), q_best(B, NA_REAL);
    int num_bad = 0;
    std::vector<int> rows(n_out);
    std::vector<double> real_y(n_out), given_tx(n_out), rec_tx(n_out);

    for (int b = 0; b < B; ++b) {
        Eigen::VectorXd full_counts = Eigen::VectorXd::Zero(n_full);
        for (int pos = 0; pos < n_out; ++pos) {
            const int row = std::min(n_full - 1, static_cast<int>(std::floor(R::unif_rand() * n_full)));
            rows[pos] = row;
            full_counts[row] += 1.0;
        }

        Eigen::VectorXd fold_start = global_start;
        bool ok = true;
        for (int fold = 0; fold < n_folds && ok; ++fold) {
            const int begin = begin_r[fold] - 1;
            const int end = end_r[fold] - 1;
            if (begin < 0 || end < begin || end >= n_out) { ok = false; break; }

            Eigen::VectorXd counts = full_counts;
            for (int pos = begin; pos <= end; ++pos) counts[rows[pos]] -= 1.0;

            std::vector<int> active_rows;
            active_rows.reserve(n_full);
            for (int row = 0; row < n_full; ++row) {
                if (counts[row] > 0.0) active_rows.push_back(row);
            }

            Eigen::VectorXd beta;
            ok = fit_logistic_counted_active_irls(
                X_obs, y_full, counts, active_rows, fold_start, beta, 100, 1e-8);
            if (!ok || !beta.allFinite()) break;
            if (active_rows_max_abs_eta(X_obs, active_rows, beta) > 30.0) { ok = false; break; }
            //fold_start intentionally not updated -- see the matching comment in
            //fast_default_incidence_cv_run_cpp() above

            for (int pos = begin; pos <= end; ++pos) {
                const int row = rows[pos];
                const double yhat0 = X_tx0.row(row).dot(beta);
                const double yhat1 = X_tx1.row(row).dot(beta);
                const double tx = treatment_r[row];
                const double et = yhat0 + tx * (yhat1 - yhat0);
                const double ec = yhat0 + yhat1 - et;
                const bool optimal = y_higher_is_better ? (et > ec) : (et < ec);
                real_y[pos] = y_full[row];
                given_tx[pos] = tx;
                rec_tx[pos] = tx * optimal + (1.0 - tx) * (1.0 - optimal);
            }
        }

        DefaultQScores q;
        if (ok) {
            if (metric_code == 0) {
                q = default_continuous_or_probability_q_scores_cpp(
                    real_y.data(), given_tx.data(), rec_tx.data(), n_out, y_higher_is_better);
            } else {
                q = default_incidence_ratio_q_scores_cpp(
                    real_y.data(), given_tx.data(), rec_tx.data(), n_out,
                    y_higher_is_better, metric_code);
            }
            q_adv[b] = q.adversarial;
            q_avg[b] = q.average;
            q_best[b] = q.best;
        }
        if (!ok || q.is_bad) ++num_bad;
    }

    return List::create(
        Named("ok") = true,
        Named("q_scores") = List::create(
            Named("adversarial") = q_adv,
            Named("average") = q_avg,
            Named("best") = q_best
        ),
        Named("num_bad") = num_bad
    );
}
