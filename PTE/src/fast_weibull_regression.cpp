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

} // namespace

//' @title Compute Weibull Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a Weibull AFT regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of survival times.
//' @param dead A numeric vector of event indicators.
//' @param params A numeric vector of parameters [beta, log_sigma].
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
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of survival times.
//' @param dead A numeric vector of event indicators.
//' @param params A numeric vector of parameters [beta, log_sigma].
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
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of survival times.
//' @param dead A numeric vector of event indicators (1=event, 0=censored).
//' @param warm_start_beta Optional starting values for coefficients.
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
