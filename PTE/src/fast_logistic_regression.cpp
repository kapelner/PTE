#include "_helper_functions.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

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
    const int p = X.cols();
    score.setZero();
    out.setZero();
    for (int j = 0; j < p; ++j) {
        const double* xj = X.col(j).data();
        for (int k = j; k < p; ++k) {
            const double* xk = X.col(k).data();
            double acc = 0.0;
            if (k == j) {
                double score_acc = 0.0;
                for (int i = 0; i < n; ++i) {
                    acc += xj[i] * w[i] * xj[i];
                    score_acc += xj[i] * residual[i];
                }
                score[j] = score_acc;
            } else {
                for (int i = 0; i < n; ++i) acc += xj[i] * w[i] * xk[i];
            }
            out(j, k) = acc;
            if (k != j) out(k, j) = acc;
        }
    }
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
        double neg_ll = 0.0;
        Eigen::VectorXd diff(m_n);
        
        for (int i = 0; i < m_n; ++i) {
            double ei = eta[i];
            double prob = plogis_manual(ei);
            double wi = m_use_weights ? m_weights[i] : 1.0;
            neg_ll += wi * (log1pexp_stable(ei) - m_y[i] * ei);
            diff[i] = wi * (prob - m_y[i]);
        }
        grad.noalias() = m_X.transpose() * diff;
        return neg_ll;
    }
};

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

        // Fast vectorized plogis
        mu.array() = 1.0 / (1.0 + (-eta.array()).exp());

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
        double nl = 0.0;
        Eigen::VectorXd final_eta = eta_fixed + X_free * beta_free;
        for (int i = 0; i < n; ++i) {
            double wi = use_weights ? weights[i] : 1.0;
            nl += wi * (log1pexp_stable(final_eta[i]) - y[i] * final_eta[i]);
        }
        res.neg_ll = nl;
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
    mu.array() = 1.0 / (1.0 + (-eta.array()).exp());
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
    w.array() = 1.0 / (1.0 + (-eta.array()).exp()); // mu
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
    mu.array() = 1.0 / (1.0 + (-eta.array()).exp());
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
    w.array() = 1.0 / (1.0 + (-eta.array()).exp()); // mu
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
