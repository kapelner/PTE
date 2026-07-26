#ifndef EDI_HELPERS_H
#define EDI_HELPERS_H

// R's CXXFLAGS includes -UNDEBUG which overrides PKG_CXXFLAGS's -DNDEBUG.
// Source-level #define takes effect after command-line flags, so re-defining
// NDEBUG HERE (before ALL includes) disables C assert() for Eigen/LBFGSpp
// headers. Without this, Eigen's BDCSVD and LBFGSpp line search both fire
// abort() on numerical edge cases, crashing R.
#ifndef NDEBUG
#define NDEBUG
#endif

#include "ordinal_fixed_link_helpers.h"
#include "optimization_starts.h"
#include <RcppEigen.h>
#include <R_ext/BLAS.h>
#include <optimization/LBFGS.h>
#include <vector>
#include <set>
#include <limits>
#include <cmath>
#include <string>
#include <type_traits>
#include <Rmath.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Rcpp::as;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::stop;

#include "fast_erfc.h"

// Pure C++ result structure to avoid R List contention
struct ModelResult {
    Eigen::VectorXd b;
    Eigen::VectorXd mu;
    Eigen::MatrixXd XtWX;
    Eigen::VectorXd score;
    double neg_ll;
    double ssq_b_j;
    double ssq_b_2;
    double dispersion;
    double sigma2_hat;
    int iterations;
    bool converged;

    ModelResult() : neg_ll(NA_REAL), ssq_b_j(NA_REAL), ssq_b_2(NA_REAL), dispersion(NA_REAL), sigma2_hat(NA_REAL), iterations(0), converged(false) {}
};

// Pure C++ internal helpers
double compute_diagonal_inverse_entry(const Eigen::Ref<const Eigen::MatrixXd>& M, int j);

struct WeibullStart {
    Eigen::VectorXd beta;
    double log_sigma;

    WeibullStart() : log_sigma(0.0) {}
};

struct OrdinalStart {
    Eigen::VectorXd alpha;
    Eigen::VectorXd beta;
};

// R-facing exports
double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j);
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w);

using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

struct ContiguousGroupLayout {
    std::vector<int> start;
    std::vector<int> size;
    int G;
    int max_size;

    ContiguousGroupLayout() : G(0), max_size(0) {}
};

template <typename GroupIdAccessor>
inline ContiguousGroupLayout build_contiguous_group_layout(int n, GroupIdAccessor group_id_at) {
    ContiguousGroupLayout layout;
    if (n <= 0) return layout;

    int start = 0;
    while (start < n) {
        const auto group_id = group_id_at(start);
        int end = start + 1;
        while (end < n && group_id_at(end) == group_id) ++end;
        const int sz = end - start;
        layout.start.push_back(start);
        layout.size.push_back(sz);
        if (sz > layout.max_size) layout.max_size = sz;
        start = end;
    }

    layout.G = static_cast<int>(layout.start.size());
    return layout;
}

inline bool should_parallelize_replicates(int n_work_items,
                                          int item_size,
                                          int num_cores,
                                          int min_items = 128,
                                          long long min_total_work = 20000) {
    return num_cores > 1 &&
           n_work_items >= min_items &&
           static_cast<long long>(n_work_items) * static_cast<long long>(item_size) >= min_total_work;
}

template<typename Derived, typename WDerived>
inline Eigen::MatrixXd weighted_crossprod(const Eigen::MatrixBase<Derived>& X,
                                          const Eigen::MatrixBase<WDerived>& w) {
    const int n = X.rows();
    const int p = X.cols();
    if (w.rows() != n) {
        Rcpp::stop("weighted_crossprod: weight vector has incompatible dimensions");
    }

    if (Derived::IsRowMajor) {
        Eigen::MatrixXd res = Eigen::MatrixXd::Zero(p, p);
        for (int i = 0; i < n; ++i) {
            double wi = w(i);
            if (wi == 0.0) continue;
            for (int j = 0; j < p; ++j) {
                double xij = X(i, j);
                double w_xij = wi * xij;
                for (int k = j; k < p; ++k) {
                    res(j, k) += w_xij * X(i, k);
                }
            }
        }
        res.triangularView<Eigen::Lower>() = res.transpose();
        return res;
    }

    // col-major: column-pair DSYR — upper triangle only, halves FLOPs vs full GEMM
    Eigen::MatrixXd res(p, p);
    Eigen::VectorXd wXj(n);
    for (int j = 0; j < p; ++j) {
        wXj.noalias() = (w.array() * X.col(j).array()).matrix();
        for (int k = j; k < p; ++k) {
            const double acc = wXj.dot(X.col(k));
            res(j, k) = acc;
            if (k != j) res(k, j) = acc;
        }
    }
    return res;
}

// Overload for Map
template<typename WDerived>
inline Eigen::MatrixXd weighted_crossprod(const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>& X,
                                          const Eigen::MatrixBase<WDerived>& w) {
    const int n = X.rows();
    const int p = X.cols();
    if (w.rows() != n) {
        Rcpp::stop("weighted_crossprod: weight vector has incompatible dimensions");
    }

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(p, p);
    for (int i = 0; i < n; ++i) {
        double wi = w(i);
        if (wi == 0.0) continue;
        for (int j = 0; j < p; ++j) {
            double xij = X(i, j);
            double w_xij = wi * xij;
            for (int k = j; k < p; ++k) {
                res(j, k) += w_xij * X(i, k);
            }
        }
    }
    res.triangularView<Eigen::Lower>() = res.transpose();
    return res;
}

template <typename Derived, typename WDerived, typename YDerived>
inline Eigen::VectorXd weighted_crossprod_rhs(const Eigen::MatrixBase<Derived>& X,
                                              const Eigen::MatrixBase<WDerived>& w,
                                              const Eigen::MatrixBase<YDerived>& y) {
    const int n = X.rows();
    const int p = X.cols();
    if (w.rows() != n || w.cols() != 1 || y.rows() != n || y.cols() != 1) {
        Rcpp::stop("weighted_crossprod_rhs: vectors have incompatible dimensions");
    }

    if (Derived::IsRowMajor) {
        Eigen::VectorXd res = Eigen::VectorXd::Zero(p);
        for (int i = 0; i < n; ++i) {
            double wi_yi = w(i) * y(i);
            if (wi_yi == 0.0) continue;
            for (int j = 0; j < p; ++j) {
                res(j) += X(i, j) * wi_yi;
            }
        }
        return res;
    }

    return X.transpose() * w.cwiseProduct(y);
}

// Overload for Map
template<typename WDerived, typename YDerived>
inline Eigen::VectorXd weighted_crossprod_rhs(const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>& X,
                                              const Eigen::MatrixBase<WDerived>& w,
                                              const Eigen::MatrixBase<YDerived>& y) {
    const int n = X.rows();
    const int p = X.cols();
    if (w.rows() != n || w.cols() != 1 || y.rows() != n || y.cols() != 1) {
        Rcpp::stop("weighted_crossprod_rhs: vectors have incompatible dimensions");
    }

    Eigen::VectorXd res = Eigen::VectorXd::Zero(p);
    for (int i = 0; i < n; ++i) {
        double wi_yi = w(i) * y(i);
        if (wi_yi == 0.0) continue;
        for (int j = 0; j < p; ++j) {
            res(j) += X(i, j) * wi_yi;
        }
    }
    return res;
}

// Compute X^T X using BLAS DSYRK (half FLOPs vs full DGEMM, BLAS-optimized)
template<typename Derived>
inline Eigen::MatrixXd symmetric_crossprod(const Eigen::MatrixBase<Derived>& X) {
    const BLAS_INT n = static_cast<BLAS_INT>(X.rows());
    const BLAS_INT p = static_cast<BLAS_INT>(X.cols());
    Eigen::MatrixXd res(p, p);
    if (p == 0) return res;
    const double alpha = 1.0, beta = 0.0;
    const BLAS_INT lda = static_cast<BLAS_INT>(X.derived().outerStride());
    const BLAS_INT ldc = p;
    F77_CALL(dsyrk)("U", "T", &p, &n, &alpha, X.derived().data(), &lda, &beta, res.data(), &ldc FCONE FCONE);
    res.triangularView<Eigen::Lower>() = res.transpose();
    return res;
}

struct FixedParamSpec {
    Eigen::VectorXi fixed_idx;
    Eigen::VectorXi free_idx;
    Eigen::VectorXd fixed_values;
    bool has_fixed;

    FixedParamSpec() : has_fixed(false) {}
};

inline FixedParamSpec make_fixed_param_spec(
    int n_params,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
    FixedParamSpec spec;

    if (fixed_idx.isNull()) {
        spec.free_idx.resize(n_params);
        for (int i = 0; i < n_params; ++i) spec.free_idx[i] = i;
        return spec;
    }

    Rcpp::IntegerVector fixed_idx_r(fixed_idx);
    if (fixed_idx_r.size() == 0) {
        spec.free_idx.resize(n_params);
        for (int i = 0; i < n_params; ++i) spec.free_idx[i] = i;
        return spec;
    }
    if (fixed_values.isNull()) {
        Rcpp::stop("fixed_values must be supplied when fixed_idx is non-empty");
    }
    Rcpp::NumericVector fixed_values_r(fixed_values);
    if (fixed_values_r.size() != fixed_idx_r.size()) {
        Rcpp::stop("fixed_idx and fixed_values must have the same length");
    }

    std::vector<int> fixed_zero_based(fixed_idx_r.size());
    std::set<int> seen;
    for (int i = 0; i < fixed_idx_r.size(); ++i) {
        int idx = fixed_idx_r[i] - 1;
        if (idx < 0 || idx >= n_params) {
            Rcpp::stop("fixed_idx entries must be one-based indices within the parameter vector");
        }
        if (!seen.insert(idx).second) {
            Rcpp::stop("fixed_idx cannot contain duplicate entries");
        }
        if (!R_finite(fixed_values_r[i])) {
            Rcpp::stop("fixed_values must be finite");
        }
        fixed_zero_based[i] = idx;
    }

    spec.has_fixed = true;
    spec.fixed_idx.resize(fixed_zero_based.size());
    spec.fixed_values.resize(fixed_zero_based.size());
    std::vector<int> is_fixed(n_params, 0);
    for (int i = 0; i < static_cast<int>(fixed_zero_based.size()); ++i) {
        spec.fixed_idx[i] = fixed_zero_based[i];
        spec.fixed_values[i] = fixed_values_r[i];
        is_fixed[fixed_zero_based[i]] = 1;
    }

    int n_free = n_params - fixed_zero_based.size();
    if (n_free <= 0) {
        Rcpp::stop("at least one parameter must remain free");
    }
    spec.free_idx.resize(n_free);
    int k = 0;
    for (int i = 0; i < n_params; ++i) {
        if (!is_fixed[i]) spec.free_idx[k++] = i;
    }

    return spec;
}

inline Eigen::VectorXd subset_vector(const Eigen::VectorXd& x, const Eigen::VectorXi& idx) {
    Eigen::VectorXd out(idx.size());
    for (int i = 0; i < idx.size(); ++i) out[i] = x[idx[i]];
    return out;
}

inline Eigen::MatrixXd subset_matrix(const Eigen::MatrixXd& M, const Eigen::VectorXi& row_idx, const Eigen::VectorXi& col_idx) {
    for (int i = 0; i < row_idx.size(); ++i) {
        if (row_idx[i] < 0 || row_idx[i] >= M.rows()) {
            Rcpp::stop("subset_matrix: row index out of bounds");
        }
    }
    for (int j = 0; j < col_idx.size(); ++j) {
        if (col_idx[j] < 0 || col_idx[j] >= M.cols()) {
            Rcpp::stop("subset_matrix: column index out of bounds");
        }
    }
    Eigen::MatrixXd out(row_idx.size(), col_idx.size());
    for (int i = 0; i < row_idx.size(); ++i) {
        for (int j = 0; j < col_idx.size(); ++j) {
            out(i, j) = M(row_idx[i], col_idx[j]);
        }
    }
    return out;
}

inline Eigen::VectorXd apply_fixed_values(Eigen::VectorXd params, const FixedParamSpec& spec) {
    for (int i = 0; i < spec.fixed_idx.size(); ++i) {
        params[spec.fixed_idx[i]] = spec.fixed_values[i];
    }
    return params;
}

inline Eigen::VectorXd expand_free_params(const Eigen::VectorXd& free_params,
                                          const Eigen::VectorXd& full_template,
                                          const FixedParamSpec& spec) {
    Eigen::VectorXd full = apply_fixed_values(full_template, spec);
    for (int i = 0; i < spec.free_idx.size(); ++i) {
        full[spec.free_idx[i]] = free_params[i];
    }
    return full;
}

inline Eigen::MatrixXd expand_free_covariance(int n_params,
                                              const FixedParamSpec& spec,
                                              const Eigen::MatrixXd& cov_free,
                                              bool fixed_as_na = true) {
    Eigen::MatrixXd cov(n_params, n_params);
    if (fixed_as_na) {
        cov.setConstant(NA_REAL);
    } else {
        cov.setZero();
    }
    for (int i = 0; i < spec.free_idx.size(); ++i) {
        for (int j = 0; j < spec.free_idx.size(); ++j) {
            cov(spec.free_idx[i], spec.free_idx[j]) = cov_free(i, j);
        }
    }
    return cov;
}

inline Eigen::MatrixXd symmetric_pseudo_inverse(const Eigen::MatrixXd& M, double tol = 1e-10) {
    Eigen::MatrixXd Msym = (M + M.transpose()) / 2.0;
    if (!Msym.allFinite()) {
        return Eigen::MatrixXd::Constant(M.rows(), M.cols(), std::numeric_limits<double>::quiet_NaN());
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Msym);
    if (es.info() != Eigen::Success) {
        return Eigen::MatrixXd::Constant(M.rows(), M.cols(), NA_REAL);
    }
    const Eigen::VectorXd evals = es.eigenvalues();
    const double max_abs_eval = evals.cwiseAbs().maxCoeff();
    Eigen::VectorXd inv_evals(evals.size());
    for (int i = 0; i < evals.size(); ++i) {
        inv_evals[i] = (max_abs_eval > 0.0 && std::abs(evals[i]) > tol * max_abs_eval) ? 1.0 / evals[i] : 0.0;
    }
    return es.eigenvectors() * inv_evals.asDiagonal() * es.eigenvectors().transpose();
}

inline double plogis_safe(double x) {
    if (x >= 0.0) { const double z = std::exp(-x); return 1.0 / (1.0 + z); }
    const double z = std::exp(x); return z / (1.0 + z);
}

inline double dplogis_safe(double x) {
    const double p = plogis_safe(x);
    return p * (1.0 - p);
}

inline Eigen::ArrayXd plogis_array_safe(const Eigen::ArrayXd& x) {
    Eigen::ArrayXd res(x.size());
    for (int i = 0; i < x.size(); ++i) {
        if (x[i] >= 0.0) {
            const double z = std::exp(-x[i]);
            res[i] = 1.0 / (1.0 + z);
        } else {
            const double z = std::exp(x[i]);
            res[i] = z / (1.0 + z);
        }
    }
    return res;
}

inline double log1pexp_safe(double x) {
    if (x > 0.0) return x + std::log1p(std::exp(-x));
    return std::log1p(std::exp(x));
}

// Vectorizable accurate log1p for z > -1: built from Eigen's packet .log() so it
// SIMD-vectorizes, unlike Eigen's .log1p() (which falls back to scalar std::log1p).
// Kahan correction keeps full precision even where (1+z) rounds to 1.
inline Eigen::ArrayXd fast_log1p_arr(const Eigen::ArrayXd& z) {
    const Eigen::ArrayXd u = 1.0 + z;
    const Eigen::ArrayXd r = u.log() * (z / (u - 1.0));
    return (u == 1.0).select(z, r);
}

inline Eigen::ArrayXd log1pexp_array_safe(const Eigen::ArrayXd& x) {
    // softplus = max(x,0) + log1p(exp(-|x|)); branchless + fully vectorized (packet exp/log
    // via fast_log1p_arr). Numerically identical to the per-branch std::log1p form.
    const Eigen::ArrayXd z = (-x.abs()).exp();
    return x.max(0.0) + fast_log1p_arr(z);
}

// Scalar log(1+exp(x)) avoiding glibc log1p dispatch entirely.
// atanh identity: log1p(z) = 2*s*(1 + s²/3 + … + s^18/19) where s = z/(2+z).
// For |x| ≤ 37, s = exp(-|x|)/(2+exp(-|x|)) ≤ 1/3 → 10-term Horner error < 5e-12.
inline double fast_log1pexp(double x) {
    if (x > 37.0) return x;
    if (x < -37.0) return std::exp(x);
    const double z  = std::exp(-std::abs(x));
    const double s  = z / (2.0 + z);
    const double s2 = s * s;
    const double p  = 1.0 + s2 * (1.0/3.0 + s2 * (1.0/5.0 + s2 * (1.0/7.0 +
                      s2 * (1.0/9.0 + s2 * (1.0/11.0 + s2 * (1.0/13.0 +
                      s2 * (1.0/15.0 + s2 * (1.0/17.0 + s2/19.0))))))));
    return (x >= 0.0 ? x : 0.0) + 2.0 * s * p;
}

// Vectorized log1pexp — branch-free, fully SIMD-vectorizable (no .log() or .select()).
// Same 10-term atanh Horner as fast_log1pexp; each Eigen operation compiles to one AVX2 instruction.
inline Eigen::ArrayXd log1pexp_array_fast(const Eigen::ArrayXd& x) {
    const Eigen::ArrayXd z  = (-x.abs()).exp();
    const Eigen::ArrayXd s  = z / (2.0 + z);
    const Eigen::ArrayXd s2 = s * s;
    const Eigen::ArrayXd p  = 1.0 + s2 * (1.0/3.0 + s2 * (1.0/5.0 + s2 * (1.0/7.0 +
                               s2 * (1.0/9.0 + s2 * (1.0/11.0 + s2 * (1.0/13.0 +
                               s2 * (1.0/15.0 + s2 * (1.0/17.0 + s2/19.0))))))));
    return x.max(0.0) + 2.0 * s * p;
}

inline bool try_safe_ols_solve(const Eigen::MatrixXd& X,
                               const Eigen::VectorXd& y,
                               Eigen::VectorXd& beta_out) {
    const int p = X.cols();
    if (X.rows() == 0 || p == 0 || y.size() != X.rows() ||
        !X.allFinite() || !y.allFinite()) {
        beta_out = Eigen::VectorXd::Zero(p);
        return false;
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
    beta_out = qr.solve(y);
    if (!beta_out.allFinite()) {
        beta_out = Eigen::VectorXd::Zero(p);
        return false;
    }
    return true;
}

inline Eigen::VectorXd safe_ols_solve(const Eigen::MatrixXd& X,
                                      const Eigen::VectorXd& y) {
    Eigen::VectorXd beta_out;
    try_safe_ols_solve(X, y, beta_out);
    return beta_out;
}

inline bool vector_is_usable_start(const Eigen::VectorXd& x, int expected_size = -1) {
    return x.allFinite() && (expected_size < 0 || x.size() == expected_size);
}

inline Eigen::VectorXd finalize_warm_start_beta(const Eigen::VectorXd& smart_cold_start,
                                           const Eigen::VectorXd& legacy_start,
                                           const FixedParamSpec& fixed_spec,
                                           bool use_smart) {
    Eigen::VectorXd out = use_smart && vector_is_usable_start(smart_cold_start, legacy_start.size()) ?
        smart_cold_start : legacy_start;
    return apply_fixed_values(out, fixed_spec);
}

inline Eigen::VectorXd ols_smart_cold_start_beta(const Eigen::MatrixXd& X,
                                      const Eigen::VectorXd& y) {
    return safe_ols_solve(X, y);
}

inline Eigen::VectorXd ols_smart_cold_start_beta_on_log1p(const Eigen::MatrixXd& X,
                                               const Eigen::VectorXd& y) {
    const int p = X.cols();
    if (y.size() != X.rows() || !y.allFinite()) return Eigen::VectorXd::Zero(p);
    if ((y.array() <= -1.0).any()) return Eigen::VectorXd::Zero(p);
    return safe_ols_solve(X, y.array().log1p().matrix());
}

inline Eigen::VectorXd ols_smart_cold_start_beta_or_legacy(const Eigen::MatrixXd& X,
                                                const Eigen::VectorXd& y,
                                                const Eigen::VectorXd& legacy_start,
                                                const FixedParamSpec& fixed_spec) {
    Eigen::VectorXd beta_out;
    const bool ok = try_safe_ols_solve(X, y, beta_out);
    return finalize_warm_start_beta(beta_out, legacy_start, fixed_spec, ok);
}

inline Eigen::VectorXd ols_smart_cold_start_beta_on_log1p_or_legacy(const Eigen::MatrixXd& X,
                                                         const Eigen::VectorXd& y,
                                                         const Eigen::VectorXd& legacy_start,
                                                         const FixedParamSpec& fixed_spec) {
    Eigen::VectorXd beta_out = Eigen::VectorXd::Zero(X.cols());
    const bool ok = y.size() == X.rows() && y.allFinite() && !(y.array() <= -1.0).any() &&
        try_safe_ols_solve(X, y.array().log1p().matrix(), beta_out);
    return finalize_warm_start_beta(beta_out, legacy_start, fixed_spec, ok);
}

inline Eigen::VectorXd weibull_start_to_params(const WeibullStart& start) {
    Eigen::VectorXd params(start.beta.size() + 1);
    params.head(start.beta.size()) = start.beta;
    params[start.beta.size()] = start.log_sigma;
    return params;
}

inline WeibullStart weibull_start_from_params(const Eigen::VectorXd& params) {
    WeibullStart out;
    if (params.size() == 0) return out;
    out.beta = params.head(params.size() - 1);
    out.log_sigma = params[params.size() - 1];
    return out;
}

inline bool weibull_start_is_usable(const WeibullStart& start, int p) {
    return start.beta.size() == p && start.beta.allFinite() && std::isfinite(start.log_sigma);
}

inline WeibullStart weibull_aft_start(const Eigen::MatrixXd& X,
                                      const Eigen::VectorXd& y,
                                      const Eigen::VectorXd& dead) {
    WeibullStart out;
    const int n = X.rows();
    const int p = X.cols();
    out.beta = Eigen::VectorXd::Zero(p);
    if (y.size() != n || dead.size() != n || n == 0 || p == 0 ||
        !X.allFinite() || !y.allFinite() || !dead.allFinite()) {
        return out;
    }

    std::vector<int> uncensored_rows;
    std::vector<int> positive_rows;
    uncensored_rows.reserve(n);
    positive_rows.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (y[i] > 0.0) {
            positive_rows.push_back(i);
            if (dead[i] > 0.0) uncensored_rows.push_back(i);
        }
    }

    const std::vector<int>& rows_used =
        static_cast<int>(uncensored_rows.size()) > p ? uncensored_rows : positive_rows;
    if (rows_used.empty()) return out;

    Eigen::MatrixXd X_sub(rows_used.size(), p);
    Eigen::VectorXd log_y(rows_used.size());
    for (int i = 0; i < static_cast<int>(rows_used.size()); ++i) {
        const int row = rows_used[i];
        X_sub.row(i) = X.row(row);
        log_y[i] = std::log(y[row]);
    }

    // R survival package (survreg) uses Method of Moments on the log scale:
    // mean(log_y) + 0.572 and var(log_y) / 1.64
    // Here we apply this to the OLS fit on uncensored data.
    Eigen::VectorXd log_y_adj = log_y.array() + 0.5722;
    out.beta = safe_ols_solve(X_sub, log_y_adj);
    if (!out.beta.allFinite()) {
        out.beta = Eigen::VectorXd::Zero(p);
        return out;
    }

    Eigen::VectorXd resid = log_y_adj - X_sub * out.beta;
    const double denom = std::max(1.0, static_cast<double>(rows_used.size() - p));
    const double std_resid_sq = std::max(0.0, resid.squaredNorm() / denom);
    // survreg uses var / 1.645 for the squared scale
    out.log_sigma = 0.5 * std::log(std::max(1e-8, std_resid_sq / 1.6449));
    return out;
}

inline WeibullStart weibull_aft_start_or_legacy(const Eigen::MatrixXd& X,
                                                const Eigen::VectorXd& y,
                                                const Eigen::VectorXd& dead,
                                                const WeibullStart& legacy_start,
                                                const FixedParamSpec& fixed_spec) {
    WeibullStart smart = weibull_aft_start(X, y, dead);
    Eigen::VectorXd params = weibull_start_to_params(
        weibull_start_is_usable(smart, X.cols()) ? smart : legacy_start
    );
    params = apply_fixed_values(params, fixed_spec);
    return weibull_start_from_params(params);
}

inline double ordinal_link_quantile(edi_ordinal::Link link, double p) {
    const double pp = std::min(1.0 - 1e-8, std::max(1e-8, p));
    switch (link) {
    case edi_ordinal::Link::Logit:
        return std::log(pp / (1.0 - pp));
    case edi_ordinal::Link::Probit:
        return R::qnorm5(pp, 0.0, 1.0, 1, 0);
    case edi_ordinal::Link::Cloglog:
        return std::log(-std::log(1.0 - pp));
    case edi_ordinal::Link::Cauchit:
        return std::tan(M_PI * (pp - 0.5));
    }
    return 0.0;
}

inline double ordinal_eta_sign(edi_ordinal::Link link) {
    return link == edi_ordinal::Link::Cloglog ? 1.0 : -1.0;
}

inline Eigen::VectorXd ordinal_start_to_params(const OrdinalStart& start) {
    Eigen::VectorXd params(start.alpha.size() + start.beta.size());
    params.head(start.alpha.size()) = start.alpha;
    params.tail(start.beta.size()) = start.beta;
    return params;
}

inline OrdinalStart ordinal_start_from_params(const Eigen::VectorXd& params, int p) {
    OrdinalStart out;
    const int n_alpha = std::max(0, static_cast<int>(params.size()) - p);
    out.alpha = params.head(n_alpha);
    out.beta = params.tail(p);
    return out;
}

inline bool ordinal_start_is_usable(const OrdinalStart& start, int p, int n_alpha) {
    if (start.beta.size() != p || start.alpha.size() != n_alpha ||
        !start.beta.allFinite() || !start.alpha.allFinite()) {
        return false;
    }
    for (int k = 1; k < n_alpha; ++k) {
        if (!(start.alpha[k] > start.alpha[k - 1])) return false;
    }
    return true;
}

inline OrdinalStart ordinal_smart_cold_start(const Eigen::MatrixXd& X,
                                           const Eigen::VectorXd& y,
                                           edi_ordinal::Link link) {
    OrdinalStart out;
    const int n = X.rows();
    const int p = X.cols();
    out.beta = Eigen::VectorXd::Zero(p);
    if (y.size() != n || n == 0 || !X.allFinite() || !y.allFinite()) return out;

    const std::vector<double> levels = edi_ordinal::init_levels(y);
    const int K = static_cast<int>(levels.size());
    const int n_alpha = std::max(0, K - 1);
    out.alpha = Eigen::VectorXd::Zero(n_alpha);
    if (n_alpha == 0) return out;

    std::vector<int> counts(K, 0);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < K; ++k) {
            if (y[i] <= levels[k]) {
                counts[k]++;
                break;
            }
        }
    }
    int cumulative_count = 0;
    for (int k = 0; k < n_alpha; ++k) {
        cumulative_count += counts[k];
        const double p_k = static_cast<double>(cumulative_count) / static_cast<double>(n);
        out.alpha[k] = ordinal_link_quantile(link, p_k);
    }
    for (int k = 1; k < n_alpha; ++k) {
        if (!(out.alpha[k] > out.alpha[k - 1])) {
            out.alpha[k] = out.alpha[k - 1] + 1e-4;
        }
    }
    return out;
}

inline OrdinalStart ordinal_smart_cold_start_or_legacy(const Eigen::MatrixXd& X,
                                                      const Eigen::VectorXd& y,
                                                      edi_ordinal::Link link,
                                                      const OrdinalStart& legacy_start,
                                                      const FixedParamSpec& fixed_spec) {
    const int n_alpha = legacy_start.alpha.size();
    OrdinalStart smart = ordinal_smart_cold_start(X, y, link);
    Eigen::VectorXd params = ordinal_start_to_params(
        ordinal_start_is_usable(smart, X.cols(), n_alpha) ? smart : legacy_start
    );
    params = apply_fixed_values(params, fixed_spec);
    return ordinal_start_from_params(params, X.cols());
}

template <typename Functor>
inline Eigen::VectorXd numerical_gradient(const Functor& fun, const Eigen::VectorXd& par) {
    const int n = par.size();
    Eigen::VectorXd grad(n);
    for (int i = 0; i < n; ++i) {
        double h = 1e-5 * std::max(1.0, std::abs(par[i]));
        Eigen::VectorXd p_plus = par, p_minus = par;
        p_plus[i] += h; p_minus[i] -= h;
        grad[i] = (fun.value(p_plus) - fun.value(p_minus)) / (2.0 * h);
    }
    return grad;
}

template <typename Functor>
inline Eigen::MatrixXd numerical_hessian(Functor& fun, const Eigen::VectorXd& par) {
    const int n = par.size();
    Eigen::MatrixXd hess(n, n);
    for (int i = 0; i < n; ++i) {
        double h = 1e-4 * std::max(1.0, std::abs(par[i]));
        Eigen::VectorXd p_plus = par, p_minus = par;
        p_plus[i] += h; p_minus[i] -= h;
        
        Eigen::VectorXd g_plus(n), g_minus(n);
        fun(p_plus, g_plus);
        fun(p_minus, g_minus);
        
        hess.col(i) = (g_plus - g_minus) / (2.0 * h);
    }
    return (hess + hess.transpose()) / 2.0;
}

struct LikelihoodFitResult {
    Eigen::VectorXd params;
    double value;
    int niter;
    bool converged;

    LikelihoodFitResult() :
        value(std::numeric_limits<double>::quiet_NaN()),
        niter(0),
        converged(false) {}
};

inline std::string normalize_optimizer_algorithm(const std::string& optimization_alg,
                                                 const std::string& default_optimization_alg,
                                                 bool allow_irls) {
    std::string alg = optimization_alg.empty() ? default_optimization_alg : optimization_alg;
    if (alg == "nr" || alg == "newton" || alg == "newton-raphson") {
        alg = "newton_raphson";
    } else if (alg == "l-bfgs" || alg == "L-BFGS" || alg == "LBFGS") {
        alg = "lbfgs";
    }

    if (alg == "lbfgs" || alg == "newton_raphson" || (allow_irls && alg == "irls")) {
        return alg;
    }
    if (allow_irls) {
        Rcpp::stop("optimization_alg must be one of 'irls', 'lbfgs', or 'newton_raphson'");
    }
    Rcpp::stop("optimization_alg must be one of 'lbfgs' or 'newton_raphson'");
}

template <typename LikelihoodFunctor>
inline double likelihood_value(LikelihoodFunctor& fun,
                               const Eigen::VectorXd& params) {
    Eigen::VectorXd grad(params.size());
    return fun(params, grad);
}

template <typename LikelihoodFunctor>
inline Eigen::VectorXd likelihood_score(LikelihoodFunctor& fun,
                                        const Eigen::VectorXd& params) {
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return -grad;
}

template <typename LikelihoodFunctor>
inline Eigen::MatrixXd likelihood_information(LikelihoodFunctor& fun,
                                              const Eigen::VectorXd& params) {
    return fun.hessian(params);
}

inline Eigen::MatrixXd covariance_from_information(const Eigen::MatrixXd& information) {
    if (!information.allFinite()) {
        return Eigen::MatrixXd::Constant(information.rows(), information.cols(), std::numeric_limits<double>::quiet_NaN());
    }
    Eigen::LDLT<Eigen::MatrixXd> ldlt((information + information.transpose()) / 2.0);
    if (ldlt.info() == Eigen::Success) {
        Eigen::MatrixXd inv = ldlt.solve(Eigen::MatrixXd::Identity(information.rows(), information.cols()));
        if (inv.allFinite()) return (inv + inv.transpose()) / 2.0;
    }
    Eigen::MatrixXd pinv = symmetric_pseudo_inverse(information);
    return (pinv + pinv.transpose()) / 2.0;
}

inline Rcpp::List make_uniform_likelihood_fit_result(const Eigen::VectorXd& params,
                                                     double neg_loglik,
                                                     bool converged,
                                                     const Eigen::VectorXd& score,
                                                     const Eigen::MatrixXd& observed_information,
                                                     bool estimate_only = false,
                                                     const Eigen::MatrixXd* fisher_information = nullptr,
                                                     const std::string& information_type = "observed") {
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("params") = params,
        Rcpp::Named("neg_loglik") = neg_loglik,
        Rcpp::Named("neg_ll") = neg_loglik,
        Rcpp::Named("loglik") = R_finite(neg_loglik) ? -neg_loglik : NA_REAL,
        Rcpp::Named("converged") = converged
    );
    if (!estimate_only) {
        out["score"] = score;
        out["observed_information"] = observed_information;
        out["hessian"] = -observed_information;
        if (fisher_information != nullptr) {
            out["fisher_information"] = *fisher_information;
            out["information"] = *fisher_information;
            out["information_type"] = "fisher";
            out["vcov"] = covariance_from_information(*fisher_information);
        } else {
            out["information"] = observed_information;
            out["information_type"] = information_type;
            out["vcov"] = covariance_from_information(observed_information);
        }
    }
    return out;
}

inline Rcpp::List likelihood_ratio_test_from_negloglik(double unrestricted_neg_loglik,
                                                       double null_neg_loglik,
                                                       int df = 1) {
    double statistic = NA_REAL;
    double p_value = NA_REAL;
    if (R_finite(unrestricted_neg_loglik) && R_finite(null_neg_loglik) && df > 0) {
        statistic = std::max(0.0, 2.0 * (null_neg_loglik - unrestricted_neg_loglik));
        p_value = R::pchisq(statistic, static_cast<double>(df), false, false);
    }
    return Rcpp::List::create(
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("df") = df,
        Rcpp::Named("p_value") = p_value
    );
}

inline Rcpp::List score_test_from_score_information(const Eigen::VectorXd& score,
                                                    const Eigen::MatrixXd& information,
                                                    int tested_idx) {
    const int idx = tested_idx - 1;
    if (idx < 0 || idx >= score.size() || information.rows() != score.size() || information.cols() != score.size()) {
        Rcpp::stop("tested_idx must be a one-based index within the parameter vector");
    }

    std::vector<int> nuisance_idx_v;
    nuisance_idx_v.reserve(score.size() - 1);
    for (int i = 0; i < score.size(); ++i) {
        if (i != idx) nuisance_idx_v.push_back(i);
    }

    double info_eff = information(idx, idx);
    if (!information.allFinite()) {
        return Rcpp::List::create(
            Rcpp::Named("statistic") = NA_REAL,
            Rcpp::Named("df") = 1,
            Rcpp::Named("p_value") = NA_REAL,
            Rcpp::Named("score") = NA_REAL,
            Rcpp::Named("information_effective") = NA_REAL
        );
    }
    if (!nuisance_idx_v.empty()) {
        Eigen::VectorXi nuisance_idx(nuisance_idx_v.size());
        for (int i = 0; i < static_cast<int>(nuisance_idx_v.size()); ++i) nuisance_idx[i] = nuisance_idx_v[i];
        Eigen::MatrixXd I_nn = subset_matrix(information, nuisance_idx, nuisance_idx);
        Eigen::VectorXd I_nT = subset_matrix(information, nuisance_idx, Eigen::VectorXi::Constant(1, idx)).col(0);
        Eigen::VectorXd solved = covariance_from_information(I_nn) * I_nT;
        info_eff -= I_nT.dot(solved);
    }

    double statistic = NA_REAL;
    double p_value = NA_REAL;
    const double score_t = score[idx];
    if (R_finite(score_t) && R_finite(info_eff) && info_eff > 0.0) {
        statistic = score_t * score_t / info_eff;
        p_value = R::pchisq(statistic, 1.0, false, false);
    }

    return Rcpp::List::create(
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("df") = 1,
        Rcpp::Named("p_value") = p_value,
        Rcpp::Named("score") = score_t,
        Rcpp::Named("information_effective") = info_eff
    );
}

inline Rcpp::List gradient_test_from_restricted_score(const Eigen::VectorXd& score,
                                                      double unrestricted_estimate,
                                                      double null_value,
                                                      int tested_idx) {
    const int idx = tested_idx - 1;
    if (idx < 0 || idx >= score.size()) {
        Rcpp::stop("tested_idx must be a one-based index within the parameter vector");
    }

    double statistic = NA_REAL;
    double p_value = NA_REAL;
    const double score_t = score[idx];
    const double estimate_gap = unrestricted_estimate - null_value;

    if (R_finite(score_t) && R_finite(estimate_gap)) {
        statistic = score_t * estimate_gap;
        if (statistic < 0.0) {
            statistic = 0.0;
        }
        if (R_finite(statistic)) {
            p_value = R::pchisq(statistic, 1.0, false, false);
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("df") = 1,
        Rcpp::Named("p_value") = p_value,
        Rcpp::Named("score") = score_t,
        Rcpp::Named("estimate_gap") = estimate_gap
    );
}

template <typename FullFunctor>
class FixedParameterFunctor {
private:
    FullFunctor& m_fun;
    const FixedParamSpec& m_spec;
    const Eigen::VectorXd& m_full_template;

public:
    FixedParameterFunctor(FullFunctor& fun,
                          const FixedParamSpec& spec,
                          const Eigen::VectorXd& full_template) :
        m_fun(fun), m_spec(spec), m_full_template(full_template) {}

    double operator()(const Eigen::VectorXd& free_params, Eigen::VectorXd& grad_free) {
        Eigen::VectorXd full_params = expand_free_params(free_params, m_full_template, m_spec);
        Eigen::VectorXd grad_full(full_params.size());
        double val = m_fun(full_params, grad_full);
        grad_free = subset_vector(grad_full, m_spec.free_idx);
        return val;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& free_params) {
        Eigen::VectorXd full_params = expand_free_params(free_params, m_full_template, m_spec);
        Eigen::MatrixXd H_full = m_fun.hessian(full_params);
        return subset_matrix(H_full, m_spec.free_idx, m_spec.free_idx);
    }

    template <typename F = FullFunctor>
    auto expected_hessian(const Eigen::VectorXd& free_params)
        -> decltype(std::declval<F&>().expected_hessian(std::declval<const Eigen::VectorXd&>())) {
        Eigen::VectorXd full_params = expand_free_params(free_params, m_full_template, m_spec);
        Eigen::MatrixXd H_full = m_fun.expected_hessian(full_params);
        return subset_matrix(H_full, m_spec.free_idx, m_spec.free_idx);
    }

    Eigen::VectorXd expand(const Eigen::VectorXd& free_params) const {
        return expand_free_params(free_params, m_full_template, m_spec);
    }
};

template <typename LikelihoodFunctor>
inline LikelihoodFitResult optimize_likelihood_lbfgs(LikelihoodFunctor& fun,
                                                     Eigen::VectorXd params,
                                                     int maxit,
                                                     double tol,
                                                     int max_linesearch = 0) {
    LBFGSpp::LBFGSParam<double> lbfgs_params;
    lbfgs_params.epsilon = tol;
    lbfgs_params.epsilon_rel = tol;
    lbfgs_params.past = 1;
    lbfgs_params.delta = tol;
    lbfgs_params.max_iterations = maxit;
    lbfgs_params.max_linesearch = (max_linesearch > 0) ? max_linesearch : 100;
    lbfgs_params.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;

    LBFGSpp::LBFGSSolver<double> solver(lbfgs_params);
    LikelihoodFitResult fit;
    fit.params = params;
    try {
        fit.niter = solver.minimize(fun, fit.params, fit.value);
        fit.converged = (fit.niter < maxit);
    } catch (...) {
        fit.value = std::numeric_limits<double>::quiet_NaN();
        fit.converged = false;
        fit.niter = maxit;
    }
    return fit;
}

template <class T, class = void>
struct has_expected_hessian : std::false_type {};

template <class T>
struct has_expected_hessian<T, decltype(void(
    std::declval<T&>().expected_hessian(std::declval<const Eigen::VectorXd&>())
))> : std::true_type {};

template <class F>
inline Eigen::MatrixXd hessian_for_opt(F& fun, const Eigen::VectorXd& params) {
    // Newton-Raphson uses local observed curvature. Expected information belongs
    // to Fisher scoring/IRLS and inference, not generic Newton updates.
    return fun.hessian(params);
}

inline bool is_valid_warm_start_information(const Eigen::MatrixXd& H, int expected_dim) {
    if (H.rows() != expected_dim || H.cols() != expected_dim || !H.allFinite()) {
        return false;
    }
    Eigen::MatrixXd sym = (H + H.transpose()) / 2.0;
    Eigen::LLT<Eigen::MatrixXd> llt(sym);
    return llt.info() == Eigen::Success;
}

template <typename LikelihoodFunctor>
inline LikelihoodFitResult optimize_likelihood_newton(LikelihoodFunctor& fun,
                                                      Eigen::VectorXd params,
                                                      int maxit,
                                                      double tol,
                                                      const Eigen::MatrixXd* warm_start_hessian = nullptr) {
    (void)warm_start_hessian;
    LikelihoodFitResult fit;
    fit.params = params;

    for (int iter = 0; iter < maxit; ++iter) {
        Eigen::VectorXd grad(params.size());
        double current_value = fun(params, grad);
        if (!std::isfinite(current_value) || !grad.allFinite()) break;
        if (grad.norm() < tol) {
            fit.value = current_value;
            fit.niter = iter;
            fit.converged = true;
            fit.params = params;
            return fit;
        }

        Eigen::MatrixXd H = hessian_for_opt(fun, params);
        
        if (!H.allFinite()) break;
        Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
        if (!lu.isInvertible()) break;
        Eigen::VectorXd step = lu.solve(grad);
        if (!step.allFinite()) break;

        double step_scale = 1.0;
        bool accepted = false;
        while (step_scale > 1e-4) {
            Eigen::VectorXd candidate = params - step_scale * step;
            double candidate_value = likelihood_value(fun, candidate);
            if (std::isfinite(candidate_value) && candidate_value < current_value) {
                params = candidate;
                fit.value = candidate_value;
                accepted = true;
                break;
            }
            step_scale *= 0.5;
        }
        if (!accepted) break;

        fit.niter = iter + 1;
        if ((step_scale * step).norm() < tol) {
            fit.converged = true;
            fit.params = params;
            return fit;
        }
    }

    fit.params = params;
    fit.value = likelihood_value(fun, params);
    fit.converged = false;
    return fit;
}

template <typename LikelihoodFunctor>
inline LikelihoodFitResult optimize_likelihood_newton_then_lbfgs(LikelihoodFunctor& fun,
                                                                 Eigen::VectorXd params,
                                                                 int maxit,
                                                                 double tol,
                                                                 int max_linesearch = 0,
                                                                 const Eigen::MatrixXd* warm_start_hessian = nullptr) {
    LikelihoodFitResult newton_fit = optimize_likelihood_newton(fun, params, maxit, tol, warm_start_hessian);
    if (newton_fit.converged) {
        return newton_fit;
    }

    try {
        LikelihoodFitResult lbfgs_fit = optimize_likelihood_lbfgs(fun, newton_fit.params, maxit, tol, max_linesearch);
        if (lbfgs_fit.converged) {
            return lbfgs_fit;
        }
        if (std::isfinite(lbfgs_fit.value) &&
            (!std::isfinite(newton_fit.value) || lbfgs_fit.value < newton_fit.value)) {
            return lbfgs_fit;
        }
    } catch (...) {
        // Keep the best damped-Newton result when the fallback optimizer also fails.
    }

    return newton_fit;
}

template <typename FullFunctor>
inline LikelihoodFitResult optimize_fixed_likelihood_lbfgs(FullFunctor& fun,
                                                           Eigen::VectorXd params,
                                                           const FixedParamSpec& fixed_spec,
                                                           int maxit,
                                                           double tol,
                                                           int max_linesearch = 0) {
    params = apply_fixed_values(params, fixed_spec);
    Eigen::VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
    FixedParameterFunctor<FullFunctor> fixed_fun(fun, fixed_spec, params);
    LikelihoodFitResult fit = optimize_likelihood_lbfgs(fixed_fun, params_free, maxit, tol, (max_linesearch > 0 ? max_linesearch : 100));
    fit.params = fixed_fun.expand(fit.params);
    return fit;
}

template <typename FullFunctor>
inline LikelihoodFitResult optimize_fixed_likelihood_newton(FullFunctor& fun,
                                                            Eigen::VectorXd params,
                                                            const FixedParamSpec& fixed_spec,
                                                            int maxit,
                                                            double tol,
                                                            const Eigen::MatrixXd* warm_start_hessian = nullptr) {
    params = apply_fixed_values(params, fixed_spec);
    Eigen::VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
    FixedParameterFunctor<FullFunctor> fixed_fun(fun, fixed_spec, params);
    
    Eigen::MatrixXd H_free;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_hessian != nullptr) {
        H_free = subset_matrix(*warm_start_hessian, fixed_spec.free_idx, fixed_spec.free_idx);
        if (is_valid_warm_start_information(H_free, params_free.size())) {
            h_ptr = &H_free;
        }
    }
    
    LikelihoodFitResult fit = optimize_likelihood_newton(fixed_fun, params_free, maxit, tol, h_ptr);
    fit.params = fixed_fun.expand(fit.params);
    return fit;
}

template <typename FullFunctor>
inline LikelihoodFitResult optimize_fixed_likelihood(FullFunctor& fun,
                                                     Eigen::VectorXd params,
                                                     const FixedParamSpec& fixed_spec,
                                                     int maxit,
                                                     double tol,
                                                     const std::string& optimization_alg,
                                                     const std::string& default_optimization_alg,
                                                     int max_linesearch = 0,
                                                     const Eigen::MatrixXd* warm_start_hessian = nullptr) {
    std::string alg = normalize_optimizer_algorithm(optimization_alg, default_optimization_alg, false);
    if (alg == "lbfgs") {
        return optimize_fixed_likelihood_lbfgs(fun, params, fixed_spec, maxit, tol, max_linesearch);
    }
    params = apply_fixed_values(params, fixed_spec);
    Eigen::VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
    FixedParameterFunctor<FullFunctor> fixed_fun(fun, fixed_spec, params);
    
    Eigen::MatrixXd H_free;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_hessian != nullptr) {
        H_free = subset_matrix(*warm_start_hessian, fixed_spec.free_idx, fixed_spec.free_idx);
        if (is_valid_warm_start_information(H_free, params_free.size())) {
            h_ptr = &H_free;
        }
    }

    LikelihoodFitResult fit = optimize_likelihood_newton_then_lbfgs(fixed_fun, params_free, maxit, tol, max_linesearch, h_ptr);
    fit.params = fixed_fun.expand(fit.params);
    return fit;
}

template <typename LikelihoodFunctor>
inline LikelihoodFitResult optimize_likelihood(LikelihoodFunctor& fun,
                                               Eigen::VectorXd params,
                                               int maxit,
                                               double tol,
                                               const std::string& optimization_alg,
                                               const std::string& default_optimization_alg,
                                               int max_linesearch = 0,
                                               const Eigen::MatrixXd* warm_start_hessian = nullptr) {
    std::string alg = normalize_optimizer_algorithm(optimization_alg, default_optimization_alg, false);
    if (alg == "lbfgs") {
        return optimize_likelihood_lbfgs(fun, params, maxit, tol, max_linesearch);
    }
    return optimize_likelihood_newton_then_lbfgs(fun, params, maxit, tol, max_linesearch, warm_start_hessian);
}

// Fast digamma via A&S 6.3.18 asymptotic expansion + recurrence shift.
// Accurate to ≤ 4e-12 relative error for x > 0; falls back to R::digamma for x <= 0.
inline double fast_digamma(double x) {
    if (x <= 0.0) return R::digamma(x);
    double r = 0.0;
    while (x < 8.0) { r -= 1.0 / x; x += 1.0; }
    const double ix = 1.0 / x, ix2 = ix * ix;
    r += std::log(x) - 0.5 * ix
       - ix2 * (1.0/12.0 - ix2 * (1.0/120.0 - ix2 * (1.0/252.0 - ix2 * (1.0/240.0 - ix2/132.0))));
    return r;
}

inline double fast_lgamma_stirling(double x) {
    const double inv = 1.0 / x;
    const double inv2 = inv * inv;
    const double series = inv * (
        1.0 / 12.0 + inv2 * (
        -1.0 / 360.0 + inv2 * (
         1.0 / 1260.0 + inv2 * (
        -1.0 / 1680.0 + inv2 * (
         1.0 / 1188.0 + inv2 * (
        -691.0 / 360360.0 + inv2 * (1.0 / 156.0)))))));

    return (x - 0.5) * std::log(x) - x
        + 0.91893853320467274178032973640562 + series;
}

inline double fast_lgamma_lanczos(double x) {
    const double z = x - 1.0;
    double a = 0.99999999999980993;
    a += 676.5203681218851 / (z + 1.0);
    a += -1259.1392167224028 / (z + 2.0);
    a += 771.32342877765313 / (z + 3.0);
    a += -176.61502916214059 / (z + 4.0);
    a += 12.507343278686905 / (z + 5.0);
    a += -0.13857109526572012 / (z + 6.0);
    a += 9.9843695780195716e-6 / (z + 7.0);
    a += 1.5056327351493116e-7 / (z + 8.0);
    const double t = z + 7.5;
    return 0.91893853320467274178032973640562 + (z + 0.5) * std::log(t) - t + std::log(a);
}

// Fast log-gamma for positive arguments. Uses a Lanczos rational approximation
// for moderate x and Stirling for large x; falls back for nonpositive/nonfinite.
inline double fast_lgamma(double x) {
    if (x <= 0.0 || !std::isfinite(x)) return std::lgamma(x);
    if (x < 0.5) return fast_lgamma_lanczos(x + 1.0) - std::log(x);
    if (x < 8.0) return fast_lgamma_lanczos(x);
    return fast_lgamma_stirling(x);
}

#endif
