#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

class WeibullAFTLikelihood {
private:
    const Eigen::Ref<const Eigen::VectorXd> m_y;
    const Eigen::Ref<const Eigen::VectorXd> m_dead;
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    const int m_n;
    const int m_p;
    const Eigen::VectorXd m_log_y;
    Eigen::VectorXd m_eta;
    Eigen::ArrayXd m_w;
    Eigen::ArrayXd m_exp_w;
    Eigen::VectorXd m_d_eta;
    Eigen::VectorXd m_beta_weights;
    Eigen::VectorXd m_cross_weights;

public:
    WeibullAFTLikelihood(const Eigen::Ref<const Eigen::VectorXd>& y, 
                         const Eigen::Ref<const Eigen::VectorXd>& dead, 
                         const Eigen::Ref<const Eigen::MatrixXd>& X) :
        m_y(y), m_dead(dead), m_X(X), m_n(y.size()), m_p(X.cols()),
        m_log_y(y.array().log().matrix()), m_eta(m_n), m_w(m_n),
        m_exp_w(m_n), m_d_eta(m_n), m_beta_weights(m_n),
        m_cross_weights(m_n) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        // params: [beta (p), log_sigma (1)]
        const auto beta = params.head(m_p);
        double log_sigma = params[m_p];
        double sigma = std::exp(log_sigma);

        m_eta.noalias() = m_X * beta;
        m_w = ((m_log_y - m_eta) / sigma).array().min(700.0);
        m_exp_w = m_w.exp();
        const auto dead = m_dead.array();
        const double loglik = (dead * (m_w - log_sigma - m_log_y.array()) - m_exp_w).sum();
        m_d_eta = ((m_exp_w - dead) / sigma).matrix();
        const double d_ll_d_log_sigma = (m_exp_w * m_w - dead * (m_w + 1.0)).sum();

        grad.setZero();

        grad.head(m_p).noalias() = -m_X.transpose() * m_d_eta;
        grad[m_p] = - d_ll_d_log_sigma;

        return -loglik;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = params.size();
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        const auto beta = params.head(m_p);
        double sigma = std::exp(params[m_p]);
        m_eta.noalias() = m_X * beta;
        m_w = ((m_log_y - m_eta) / sigma).array().min(700.0);
        m_exp_w = m_w.exp();
        m_beta_weights = (m_exp_w / (sigma * sigma)).matrix();
        m_cross_weights =
            ((m_exp_w * (m_w + 1.0) - m_dead.array()) / sigma).matrix();

        H.topLeftCorner(m_p, m_p).noalias() = weighted_crossprod(m_X, m_beta_weights);
        H.topRightCorner(m_p, 1).noalias() = m_X.transpose() * m_cross_weights;
        H(m_p, m_p) = (m_exp_w * (m_w.square() + m_w) - m_dead.array() * m_w).sum();
        H.bottomLeftCorner(1, m_p) = H.topRightCorner(m_p, 1).transpose();
        return H;
    }
};

class WeibullAFTCountedLikelihood {
private:
    const Eigen::Ref<const Eigen::VectorXd> m_y;
    const Eigen::Ref<const Eigen::VectorXd> m_dead;
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    const Eigen::Ref<const Eigen::VectorXd> m_counts;
    const int m_n;
    const int m_p;
    const Eigen::VectorXd m_log_y;
    Eigen::VectorXd m_eta;
    Eigen::ArrayXd m_w;
    Eigen::ArrayXd m_exp_w;
    Eigen::VectorXd m_d_eta;
    Eigen::VectorXd m_beta_weights;
    Eigen::VectorXd m_cross_weights;

public:
    WeibullAFTCountedLikelihood(const Eigen::Ref<const Eigen::VectorXd>& y,
                                const Eigen::Ref<const Eigen::VectorXd>& dead,
                                const Eigen::Ref<const Eigen::MatrixXd>& X,
                                const Eigen::Ref<const Eigen::VectorXd>& counts) :
        m_y(y), m_dead(dead), m_X(X), m_counts(counts), m_n(y.size()), m_p(X.cols()),
        m_log_y(y.array().log().matrix()), m_eta(m_n), m_w(m_n),
        m_exp_w(m_n), m_d_eta(m_n), m_beta_weights(m_n),
        m_cross_weights(m_n) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        const auto beta = params.head(m_p);
        const double log_sigma = params[m_p];
        const double sigma = std::exp(log_sigma);

        m_eta.noalias() = m_X * beta;
        m_w = ((m_log_y - m_eta) / sigma).array().min(700.0);
        m_exp_w = m_w.exp();
        const auto dead = m_dead.array();
        const auto counts = m_counts.array();
        const double loglik = (counts * (dead * (m_w - log_sigma - m_log_y.array()) - m_exp_w)).sum();
        m_d_eta = (counts * ((m_exp_w - dead) / sigma)).matrix();
        const double d_ll_d_log_sigma = (counts * (m_exp_w * m_w - dead * (m_w + 1.0))).sum();

        grad.setZero();
        grad.head(m_p).noalias() = -m_X.transpose() * m_d_eta;
        grad[m_p] = -d_ll_d_log_sigma;

        return -loglik;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        const int total_p = params.size();
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        const auto beta = params.head(m_p);
        const double sigma = std::exp(params[m_p]);

        m_eta.noalias() = m_X * beta;
        m_w = ((m_log_y - m_eta) / sigma).array().min(700.0);
        m_exp_w = m_w.exp();
        const auto counts = m_counts.array();
        m_beta_weights = (counts * (m_exp_w / (sigma * sigma))).matrix();
        m_cross_weights =
            (counts * ((m_exp_w * (m_w + 1.0) - m_dead.array()) / sigma)).matrix();

        H.topLeftCorner(m_p, m_p).noalias() = weighted_crossprod(m_X, m_beta_weights);
        H.topRightCorner(m_p, 1).noalias() = m_X.transpose() * m_cross_weights;
        H(m_p, m_p) = (counts * (m_exp_w * (m_w.square() + m_w) - m_dead.array() * m_w)).sum();
        H.bottomLeftCorner(1, m_p) = H.topRightCorner(m_p, 1).transpose();
        return H;
    }
};

class WeibullAFTActiveCountedLikelihood {
private:
    const Eigen::Ref<const Eigen::VectorXd> m_y;
    const Eigen::Ref<const Eigen::VectorXd> m_dead;
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    const Eigen::Ref<const Eigen::VectorXd> m_counts;
    const std::vector<int>& m_active_rows;
    const int m_p;
    const Eigen::VectorXd m_log_y;

public:
    WeibullAFTActiveCountedLikelihood(const Eigen::Ref<const Eigen::VectorXd>& y,
                                      const Eigen::Ref<const Eigen::VectorXd>& dead,
                                      const Eigen::Ref<const Eigen::MatrixXd>& X,
                                      const Eigen::Ref<const Eigen::VectorXd>& counts,
                                      const std::vector<int>& active_rows) :
        m_y(y), m_dead(dead), m_X(X), m_counts(counts), m_active_rows(active_rows),
        m_p(X.cols()), m_log_y(y.array().log().matrix()) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        const auto beta = params.head(m_p);
        const double log_sigma = params[m_p];
        const double sigma = std::exp(log_sigma);

        double loglik = 0.0;
        grad.setZero();

        for (const int row : m_active_rows) {
            const double count = m_counts[row];
            if (count <= 0.0) continue;

            const double eta = m_X.row(row).dot(beta);
            const double w = std::min((m_log_y[row] - eta) / sigma, 700.0);
            const double exp_w = std::exp(w);
            const double dead = m_dead[row];
            loglik += count * (dead * (w - log_sigma - m_log_y[row]) - exp_w);

            const double d_eta = count * ((exp_w - dead) / sigma);
            for (int j = 0; j < m_p; ++j) {
                grad[j] -= m_X(row, j) * d_eta;
            }
            grad[m_p] -= count * (exp_w * w - dead * (w + 1.0));
        }

        return -loglik;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(params.size(), params.size());
        const auto beta = params.head(m_p);
        const double sigma = std::exp(params[m_p]);

        for (const int row : m_active_rows) {
            const double count = m_counts[row];
            if (count <= 0.0) continue;

            const double eta = m_X.row(row).dot(beta);
            const double w = std::min((m_log_y[row] - eta) / sigma, 700.0);
            const double exp_w = std::exp(w);
            const double dead = m_dead[row];
            const double beta_weight = count * exp_w / (sigma * sigma);
            const double cross_weight = count * ((exp_w * (w + 1.0) - dead) / sigma);

            for (int j = 0; j < m_p; ++j) {
                const double xj = m_X(row, j);
                for (int k = j; k < m_p; ++k) {
                    H(j, k) += beta_weight * xj * m_X(row, k);
                }
                H(j, m_p) += xj * cross_weight;
            }
            H(m_p, m_p) += count * (exp_w * (w * w + w) - dead * w);
        }

        H.triangularView<Eigen::Lower>() = H.transpose();
        return H;
    }
};

// A (quasi-)separated fold has no real MLE: the log-likelihood keeps improving as some
// coefficient runs off to +-Inf, so the optimizer "converges" at an essentially arbitrary point on
// that ridge that depends on the starting value. The batch CV/bootstrap functions below chain each
// fold's fit into the next fold's starting point for speed, so an undetected separated fold
// doesn't just report a locally-arbitrary answer, it also pushes every subsequent fold's starting
// point further out, compounding across folds. A fitted linear predictor this large (implying
// survival times around exp(15) or exp(-15)) is itself the signature of separation for any dataset
// whose genuinely-identified coefficients are of reasonable magnitude -- callers should treat
// exceeding this the same as non-convergence.
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

//' @title Compute Weibull Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a Weibull AFT regression model.
//' @param X_sexp A numeric matrix of predictors.
//' @param y_sexp A numeric vector of survival times.
//' @param dead_sexp A numeric vector of event indicators.
//' @param params_sexp A numeric vector of parameters [beta, log_sigma].
//' @return A numeric vector representing the score.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd get_weibull_regression_score_cpp(SEXP X_sexp,
                                                 SEXP y_sexp,
                                                 SEXP dead_sexp,
                                                 SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector dead_r(dead_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> dead(dead_r.begin(), dead_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    WeibullAFTLikelihood fun(y, dead, X);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return -grad;
}

//' @title Compute Weibull Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a Weibull AFT regression model.
//' @param X_sexp A numeric matrix of predictors.
//' @param y_sexp A numeric vector of survival times.
//' @param dead_sexp A numeric vector of event indicators.
//' @param params_sexp A numeric vector of parameters [beta, log_sigma].
//' @return A numeric matrix representing the Hessian.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd get_weibull_regression_hessian_cpp(SEXP X_sexp,
                                                   SEXP y_sexp,
                                                   SEXP dead_sexp,
                                                   SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector dead_r(dead_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> dead(dead_r.begin(), dead_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    WeibullAFTLikelihood fun(y, dead, X);
    return -fun.hessian(params);
}

// [[Rcpp::export]]
NumericVector compute_weibull_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const NumericMatrix& Xc,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores)
{
    const int n       = i_mat.nrow();
    const int nsim    = i_mat.ncol();
    const int n_full  = y0.size();
    const int p_cov   = Xc.ncol();
    const int p       = 2 + p_cov;   // intercept + treatment + covariates
    const int n_params = p + 1;       // +1 for log_sigma

    const double* y0_ptr   = y0.begin();
    const int*    dead_ptr = dead.begin();
    const double* xc_ptr   = (p_cov > 0) ? Xc.begin() : nullptr;
    const int*    i_ptr    = i_mat.begin();
    const int*    w_ptr    = w_mat.begin();
    const double  mult     = (delta != 0.0) ? std::exp(delta) : 1.0;

    FixedParamSpec fspec = make_fixed_param_spec(n_params);

    std::vector<double> results(nsim, NA_REAL);
    double* res_ptr = results.data();

    const bool has_noise = noise_mat.isNotNull();
    NumericMatrix noise_m;
    const double* noise_ptr = nullptr;
    if (has_noise) {
        noise_m = NumericMatrix(noise_mat);
        noise_ptr = noise_m.begin();
    }

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
    {
        Eigen::VectorXd y_b(n), dead_b(n), params0(n_params);
        Eigen::MatrixXd X_b(n, p);

#pragma omp for schedule(dynamic)
        for (int b = 0; b < nsim; ++b) {
            const int* ic = i_ptr + (size_t)b * n;
            const int* wc = w_ptr + (size_t)b * n;

            int n_t = 0, n_c = 0;
            for (int i = 0; i < n; ++i) {
                const int r  = ic[i] - 1;
                const int wt = (wc[i] == 1) ? 1 : 0;
                X_b(i, 0) = 1.0;
                X_b(i, 1) = static_cast<double>(wt);
                for (int j = 0; j < p_cov; ++j)
                    X_b(i, j + 2) = xc_ptr[(size_t)j * n_full + r];
                double yv = y0_ptr[r];
                if (has_noise) yv += noise_ptr[(size_t)b * n + i];
                y_b(i)    = (wt && mult != 1.0) ? yv * mult : yv;
                dead_b(i) = static_cast<double>(dead_ptr[r]);
                n_t += wt;
                n_c += (1 - wt);
            }
            if (n_t < 2 || n_c < 2) continue;

            params0.setZero();
            WeibullAFTLikelihood fun(y_b, dead_b, X_b);
            LikelihoodFitResult fit = optimize_fixed_likelihood(
                fun, params0, fspec, 100, 1e-8, "lbfgs", "", 0, nullptr);
            if (fit.converged && fit.params.size() > 1 && std::isfinite(fit.params[1]))
                res_ptr[b] = fit.params[1];
        }
    }
    return wrap(results);
}

//' @title Fast Weibull AFT Regression (C++)
//' @description Weibull Accelerated Failure Time model fitting.
//' @param X_sexp A numeric matrix of predictors.
//' @param y_sexp A numeric vector of survival times.
//' @param dead_sexp A numeric vector of event indicators (1=event, 0=censored).
//' @param warm_start_params Optional starting values for coefficients.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess.
//' @param estimate_only Logical. If TRUE, do not compute variance-covariance.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix.
//' @return A list containing coefficients, log_sigma, and convergence status.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List fast_weibull_regression_cpp(SEXP X_sexp,
                                 SEXP y_sexp, 
                                 SEXP dead_sexp,
                                 Nullable<NumericVector> warm_start_params = R_NilValue,
                                 bool smart_cold_start = true,
                                 bool estimate_only = false,
                                 int maxit = 100,
                                 double tol = 1e-8,
                                 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                 Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                 std::string optimization_alg = "lbfgs",
                                 Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector dead_r(dead_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> dead(dead_r.begin(), dead_r.size());

    int p = (int)X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);
    WeibullAFTLikelihood fun(y, dead, X);

    Eigen::VectorXd params = Eigen::VectorXd::Zero(p + 1);
    if (warm_start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(warm_start_params));
        if (params.size() != p + 1) stop("warm_start_params must have length equal to ncol(X) + 1");
    } else {
        WeibullStart legacy_start;
        legacy_start.beta = Eigen::VectorXd::Zero(p);
        legacy_start.log_sigma = 0.0;
        
        WeibullStart start = smart_cold_start ? weibull_aft_start_or_legacy(X, y, dead, legacy_start, fixed_spec) : legacy_start;
        params = weibull_start_to_params(start);
    }
    
    params = apply_fixed_values(params, fixed_spec);

    Eigen::MatrixXd info;
    const Eigen::MatrixXd* info_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_ptr = &info;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, info_ptr);
    
    if (estimate_only) {
        return List::create(
            Named("b") = fit.params.head(p),
            Named("log_sigma") = fit.params[p],
            Named("converged") = fit.converged,
            Named("iterations") = fit.niter
        );
    }

    Eigen::MatrixXd hess = fun.hessian(fit.params);
    return make_uniform_likelihood_fit_result(fit.params, fit.value, fit.converged, -likelihood_score(fun, fit.params), hess, false);
}

// [[Rcpp::export]]
Rcpp::List fast_default_weibull_cv_run_cpp(SEXP X_obs_sexp,
                                           SEXP y_full_sexp,
                                           SEXP dead_full_sexp,
                                           SEXP X_tx0_sexp,
                                           SEXP X_tx1_sexp,
                                           SEXP treatment_sexp,
                                           SEXP boot_idx_sexp,
                                           SEXP begin_cutoffs_sexp,
                                           SEXP end_cutoffs_sexp,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                                           bool y_higher_is_better = true) {
    NumericMatrix X_obs_r(X_obs_sexp);
    NumericVector y_full_r(y_full_sexp);
    NumericVector dead_full_r(dead_full_sexp);
    NumericMatrix X_tx0_r(X_tx0_sexp);
    NumericMatrix X_tx1_r(X_tx1_sexp);
    NumericVector treatment_r(treatment_sexp);
    IntegerVector boot_idx_r(boot_idx_sexp);
    IntegerVector begin_r(begin_cutoffs_sexp);
    IntegerVector end_r(end_cutoffs_sexp);

    Eigen::Map<const Eigen::MatrixXd> X_obs(X_obs_r.begin(), X_obs_r.nrow(), X_obs_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y_full(y_full_r.begin(), y_full_r.size());
    Eigen::Map<const Eigen::VectorXd> dead_full(dead_full_r.begin(), dead_full_r.size());
    Eigen::Map<const Eigen::MatrixXd> X_tx0(X_tx0_r.begin(), X_tx0_r.nrow(), X_tx0_r.ncol());
    Eigen::Map<const Eigen::MatrixXd> X_tx1(X_tx1_r.begin(), X_tx1_r.nrow(), X_tx1_r.ncol());

    const int n_full = X_obs.rows();
    const int p = X_obs.cols();
    const int n_out = boot_idx_r.size();
    const int n_folds = begin_r.size();
    if (end_r.size() != n_folds) stop("begin_cutoffs and end_cutoffs must have the same length");

    Eigen::VectorXd start_params = Eigen::VectorXd::Zero(p + 1);
    if (warm_start_params.isNotNull()) {
        start_params = as<Eigen::VectorXd>(NumericVector(warm_start_params));
        if (start_params.size() != p + 1) stop("warm_start_params must have length equal to ncol(X) + 1");
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1);

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

    Eigen::VectorXd fold_start = start_params;
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

        Eigen::VectorXd params = fold_start;
        LikelihoodFitResult fit;
        if (static_cast<int>(active_rows.size()) * 2 < n_full) {
            WeibullAFTActiveCountedLikelihood fun(y_full, dead_full, X_obs, counts, active_rows);
            fit = optimize_fixed_likelihood(
                fun, params, fixed_spec, 100, 1e-8, "lbfgs", "lbfgs", 0, nullptr);
        } else {
            WeibullAFTCountedLikelihood fun(y_full, dead_full, X_obs, counts);
            fit = optimize_fixed_likelihood(
                fun, params, fixed_spec, 100, 1e-8, "lbfgs", "lbfgs", 0, nullptr);
        }
        if (!fit.converged || !fit.params.allFinite()) return List::create(Named("ok") = false);
        const Eigen::VectorXd beta = fit.params.head(p);
        if (active_rows_max_abs_eta(X_obs, active_rows, beta) > 30.0) return List::create(Named("ok") = false);
        //fold_start is intentionally NOT updated to fit.params here: every fold starts from the
        //same fixed full-data warm start so one fold's fit can never drift into the next fold's
        //starting point and compound across folds (see active_rows_max_abs_eta's comment above).

        for (int pos = begin; pos <= end; ++pos) {
            const int row = rows[pos];
            const double yhat0 = std::exp(X_tx0.row(row).dot(beta));
            const double yhat1 = std::exp(X_tx1.row(row).dot(beta));
            const double tx = treatment_r[row];
            const double et = yhat0 + tx * (yhat1 - yhat0);
            const double ec = yhat0 + yhat1 - et;
            const bool optimal = y_higher_is_better ? (et > ec) : (et < ec);

            est_true[pos] = et;
            est_counterfactual[pos] = ec;
            given_tx[pos] = tx;
            rec_tx[pos] = tx * optimal + (1.0 - tx) * (1.0 - optimal);
            real_y[pos] = y_full[row];
            censored[pos] = dead_full[row];
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
Rcpp::List fast_default_weibull_bootstrap_q_cpp(SEXP X_obs_sexp,
                                                SEXP y_full_sexp,
                                                SEXP dead_full_sexp,
                                                SEXP X_tx0_sexp,
                                                SEXP X_tx1_sexp,
                                                SEXP treatment_sexp,
                                                int B,
                                                SEXP begin_cutoffs_sexp,
                                                SEXP end_cutoffs_sexp,
                                                Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                                                bool y_higher_is_better = true) {
    Rcpp::RNGScope rng_scope;
    NumericMatrix X_obs_r(X_obs_sexp);
    NumericVector y_full_r(y_full_sexp);
    NumericVector dead_full_r(dead_full_sexp);
    NumericMatrix X_tx0_r(X_tx0_sexp);
    NumericMatrix X_tx1_r(X_tx1_sexp);
    NumericVector treatment_r(treatment_sexp);
    IntegerVector begin_r(begin_cutoffs_sexp);
    IntegerVector end_r(end_cutoffs_sexp);

    Eigen::Map<const Eigen::MatrixXd> X_obs(X_obs_r.begin(), X_obs_r.nrow(), X_obs_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y_full(y_full_r.begin(), y_full_r.size());
    Eigen::Map<const Eigen::VectorXd> dead_full(dead_full_r.begin(), dead_full_r.size());
    Eigen::Map<const Eigen::MatrixXd> X_tx0(X_tx0_r.begin(), X_tx0_r.nrow(), X_tx0_r.ncol());
    Eigen::Map<const Eigen::MatrixXd> X_tx1(X_tx1_r.begin(), X_tx1_r.nrow(), X_tx1_r.ncol());

    const int n_full = X_obs.rows();
    const int p = X_obs.cols();
    const int n_out = (end_r.size() > 0) ? end_r[end_r.size() - 1] : 0;
    const int n_folds = begin_r.size();
    if (B < 0 || n_out <= 0 || end_r.size() != n_folds) stop("invalid bootstrap q inputs");

    Eigen::VectorXd global_start = Eigen::VectorXd::Zero(p + 1);
    if (warm_start_params.isNotNull()) {
        global_start = as<Eigen::VectorXd>(NumericVector(warm_start_params));
        if (global_start.size() != p + 1) stop("warm_start_params must have length equal to ncol(X) + 1");
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1);

    NumericVector q_adv(B, NA_REAL), q_avg(B, NA_REAL), q_best(B, NA_REAL);
    int num_bad = 0;
    std::vector<int> rows(n_out);
    std::vector<double> real_y(n_out), dead(n_out), given_tx(n_out), rec_tx(n_out);

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

            LikelihoodFitResult fit;
            if (static_cast<int>(active_rows.size()) * 2 < n_full) {
                WeibullAFTActiveCountedLikelihood fun(y_full, dead_full, X_obs, counts, active_rows);
                fit = optimize_fixed_likelihood(
                    fun, fold_start, fixed_spec, 100, 1e-8, "lbfgs", "lbfgs", 0, nullptr);
            } else {
                WeibullAFTCountedLikelihood fun(y_full, dead_full, X_obs, counts);
                fit = optimize_fixed_likelihood(
                    fun, fold_start, fixed_spec, 100, 1e-8, "lbfgs", "lbfgs", 0, nullptr);
            }
            if (!fit.converged || !fit.params.allFinite()) { ok = false; break; }
            const Eigen::VectorXd beta = fit.params.head(p);
            if (active_rows_max_abs_eta(X_obs, active_rows, beta) > 30.0) { ok = false; break; }
            //fold_start intentionally not updated -- see the matching comment in
            //fast_default_weibull_cv_run_cpp() above

            for (int pos = begin; pos <= end; ++pos) {
                const int row = rows[pos];
                const double yhat0 = std::exp(X_tx0.row(row).dot(beta));
                const double yhat1 = std::exp(X_tx1.row(row).dot(beta));
                const double tx = treatment_r[row];
                const double et = yhat0 + tx * (yhat1 - yhat0);
                const double ec = yhat0 + yhat1 - et;
                const bool optimal = y_higher_is_better ? (et > ec) : (et < ec);
                real_y[pos] = y_full[row];
                dead[pos] = dead_full[row];
                given_tx[pos] = tx;
                rec_tx[pos] = tx * optimal + (1.0 - tx) * (1.0 - optimal);
            }
        }

        DefaultQScores q;
        if (ok) {
            q = default_survival_q_scores_cpp(
                real_y.data(), dead.data(), given_tx.data(), rec_tx.data(), n_out, y_higher_is_better);
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
