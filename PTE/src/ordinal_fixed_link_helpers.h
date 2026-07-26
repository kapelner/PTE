#ifndef EDI_ORDINAL_FIXED_LINK_HELPERS_H
#define EDI_ORDINAL_FIXED_LINK_HELPERS_H

#include "fast_erfc.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace edi_ordinal {

enum class Link {
    Logit,
    Probit,
    Cloglog,
    Cauchit
};

inline std::vector<double> init_levels(const Eigen::VectorXd& y) {
    std::vector<double> levels(y.data(), y.data() + y.size());
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
    return levels;
}

// Cephes-style range reduction followed by a degree-5/5 rational minimax
// approximation on |x| <= tan(pi/8). Dense validation spanning the finite
// double range found a maximum absolute error of one ulp against libm atan.
inline double fast_atan(double x) {
    constexpr double tan_pi_over_8 = 0.41421356237309504880;
    constexpr double tan_3pi_over_8 = 2.41421356237309504880;
    constexpr double pi_over_4 = 0.78539816339744830962;
    constexpr double pi_over_2 = 1.57079632679489661923;

    double reduced = std::abs(x);
    double offset = 0.0;
    if (reduced > tan_3pi_over_8) {
        offset = pi_over_2;
        reduced = -1.0 / reduced;
    } else if (reduced > tan_pi_over_8) {
        offset = pi_over_4;
        reduced = (reduced - 1.0) / (reduced + 1.0);
    }

    const double z = reduced * reduced;
    const double p = ((((-0.87506086000319041228 * z
                         - 16.157537187333650637) * z
                         - 75.008557923146046673) * z
                         - 122.88666844901361724) * z
                         - 64.850219049420253718);
    const double q = (((((z + 24.858464901423062980) * z
                          + 165.02700983169885444) * z
                          + 432.88106049129026689) * z
                          + 485.39039963591369627) * z
                          + 194.55065714826139644);
    const double result = offset + reduced + reduced * z * p / q;
    return std::copysign(result, x);
}

inline double cdf(Link link, double z) {
    switch (link) {
    case Link::Logit:
        if (z >= 0.0) {
            const double e = std::exp(-z);
            return 1.0 / (1.0 + e);
        } else {
            const double e = std::exp(z);
            return e / (1.0 + e);
        }
    case Link::Probit:
        return pnorm_fast(z);
    case Link::Cloglog:
        if (z > 5.0) return 1.0;
        if (z < -37.0) return 0.0;
        return 1.0 - std::exp(-std::exp(z));
    case Link::Cauchit:
        return 0.5 + fast_atan(z) / M_PI;
    }
    return NA_REAL;
}

inline double pdf(Link link, double z) {
    switch (link) {
    case Link::Logit: {
        const double F = cdf(link, z);
        return F * (1.0 - F);
    }
    case Link::Probit:
        return dnorm_fast(z);
    case Link::Cloglog:
        if (z > 5.0 || z < -37.0) return 0.0;
        return std::exp(z - std::exp(z));
    case Link::Cauchit:
        return 1.0 / (M_PI * (1.0 + z * z));
    }
    return NA_REAL;
}

inline double pdf_derivative(Link link, double z) {
    switch (link) {
    case Link::Logit: {
        const double F = cdf(link, z);
        const double f = F * (1.0 - F);
        return f * (1.0 - 2.0 * F);
    }
    case Link::Probit: {
        const double f = dnorm_fast(z);
        return -z * f;
    }
    case Link::Cloglog: {
        if (z > 5.0 || z < -37.0) return 0.0;
        const double ez = std::exp(z);
        const double f = std::exp(z - ez);
        return f * (1.0 - ez);
    }
    case Link::Cauchit: {
        const double denom = 1.0 + z * z;
        return -2.0 * z / (M_PI * denom * denom);
    }
    }
    return NA_REAL;
}

// Fused cloglog CDF+PDF: computes exp(z) and exp(-exp(z)) once — 2 exp total vs 4.
inline void cdf_and_pdf(Link link, double z, double& F, double& f) {
    if (link == Link::Cloglog) {
        if (z > 5.0)  { F = 1.0; f = 0.0; return; }
        if (z < -37.0) { F = 0.0; f = 0.0; return; }
        const double e_z  = std::exp(z);
        const double emez = std::exp(-e_z);
        F = 1.0 - emez;
        f = e_z * emez;
    } else {
        F = cdf(link, z);
        f = pdf(link, z);
    }
}

// Fused cloglog CDF+PDF+PDF': 2 exp total vs 6 for hessian path.
inline void cdf_pdf_fpdf(Link link, double z, double& F, double& f, double& fp) {
    if (link == Link::Cloglog) {
        if (z > 5.0)  { F = 1.0; f = 0.0; fp = 0.0; return; }
        if (z < -37.0) { F = 0.0; f = 0.0; fp = 0.0; return; }
        const double e_z  = std::exp(z);
        const double emez = std::exp(-e_z);
        F  = 1.0 - emez;
        f  = e_z * emez;
        fp = f * (1.0 - e_z);
    } else {
        F  = cdf(link, z);
        f  = pdf(link, z);
        fp = pdf_derivative(link, z);
    }
}

inline int level_index(const std::vector<double>& levels, double y) {
    for (int k = 0; k < static_cast<int>(levels.size()); ++k) {
        if (y == levels[k]) return k;
    }
    return -1;
}

class FixedOrdinalRegression {
private:
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    const Eigen::Ref<const Eigen::VectorXd> m_y;
    const Eigen::Ref<const Eigen::VectorXd> m_weights;
    const std::vector<double> m_levels;
    const int m_n;
    const int m_p;
    const int m_K;
    const Link m_link;
    const double m_eta_sign;
    const bool m_use_weights;
    mutable Eigen::VectorXd m_scratch_dq;
    mutable Eigen::MatrixXd m_scratch_d2q;
    mutable Eigen::VectorXd m_scratch_v;

    inline double obs_weight(int i) const {
        return m_use_weights ? std::max(m_weights[i], 0.0) : 1.0;
    }

    bool validate_params(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        const int n_alpha = m_K - 1;
        if (params.size() != n_alpha + m_p) return false;
        for (int k = 1; k < n_alpha; ++k) {
            if (params[k] <= params[k - 1]) return false;
        }
        return true;
    }

    void add_endpoint_derivatives(int alpha_idx,
                                  double z,
                                  double endpoint_sign,
                                  const Eigen::Ref<const Eigen::RowVectorXd>& x,
                                  Eigen::Ref<Eigen::VectorXd> dq,
                                  Eigen::Ref<Eigen::MatrixXd> d2q) const {
        m_scratch_v.setZero();
        m_scratch_v[alpha_idx] = 1.0;
        m_scratch_v.tail(m_p) = m_eta_sign * x.transpose();

        const double f = pdf(m_link, z);
        const double fp = pdf_derivative(m_link, z);
        dq.noalias() += endpoint_sign * f * m_scratch_v;
        d2q.noalias() += endpoint_sign * fp * (m_scratch_v * m_scratch_v.transpose());
    }

    void add_endpoint_gradient(int alpha_idx,
                               double z,
                               double endpoint_sign,
                               const Eigen::Ref<const Eigen::RowVectorXd>& x,
                               Eigen::Ref<Eigen::VectorXd> dq) const {
        const double f = pdf(m_link, z);
        dq[alpha_idx] += endpoint_sign * f;
        dq.tail(m_p).noalias() += endpoint_sign * f * m_eta_sign * x.transpose();
    }

public:
    FixedOrdinalRegression(const Eigen::Ref<const Eigen::MatrixXd>& X,
                           const Eigen::Ref<const Eigen::VectorXd>& y,
                           Link link,
                           double eta_sign,
                           const Eigen::Ref<const Eigen::VectorXd>& weights = Eigen::VectorXd()) :
        m_X(X), m_y(y), m_weights(weights), m_levels(init_levels(y)), m_n(X.rows()), m_p(X.cols()),
        m_K(m_levels.size()), m_link(link), m_eta_sign(eta_sign), m_use_weights(weights.size() == X.rows()),
        m_scratch_dq((m_K - 1) + m_p),
        m_scratch_d2q((m_K - 1) + m_p, (m_K - 1) + m_p),
        m_scratch_v((m_K - 1) + m_p) {}

    static std::vector<double> init_levels_static(const Eigen::Ref<const Eigen::VectorXd>& y) {
        return init_levels(y);
    }

    int n_levels() const {
        return m_K;
    }

    double neg_log_likelihood(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        if (!validate_params(params)) return 1e10;
        const int n_alpha = m_K - 1;
        const Eigen::VectorXd alpha = params.head(n_alpha);
        const Eigen::VectorXd beta = params.tail(m_p);
        const Eigen::VectorXd eta = m_X * beta;
        double nll = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const int yi_idx = level_index(m_levels, m_y[i]);
            if (yi_idx < 0) return 1e10;
            const double p_upper = (yi_idx == m_K - 1) ? 1.0 : cdf(m_link, alpha[yi_idx] + m_eta_sign * eta[i]);
            const double p_lower = (yi_idx == 0) ? 0.0 : cdf(m_link, alpha[yi_idx - 1] + m_eta_sign * eta[i]);
            const double prob = std::max(1e-12, p_upper - p_lower);
            nll -= obs_weight(i) * std::log(prob);
        }
        return nll;
    }

    double operator()(const Eigen::Ref<const Eigen::VectorXd>& params, Eigen::Ref<Eigen::VectorXd> grad) const {
        const int n_params = params.size();
        grad.setZero();
        if (!validate_params(params)) return 1e10;

        const int n_alpha = m_K - 1;
        const Eigen::VectorXd alpha = params.head(n_alpha);
        const Eigen::VectorXd beta = params.tail(m_p);
        const Eigen::VectorXd eta = m_X * beta;
        double nll = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const int yi_idx = level_index(m_levels, m_y[i]);
            if (yi_idx < 0) return 1e10;

            double p_upper = 1.0, f_upper = 0.0;
            double p_lower = 0.0, f_lower = 0.0;
            if (yi_idx < m_K - 1)
                cdf_and_pdf(m_link, alpha[yi_idx] + m_eta_sign * eta[i], p_upper, f_upper);
            if (yi_idx > 0)
                cdf_and_pdf(m_link, alpha[yi_idx - 1] + m_eta_sign * eta[i], p_lower, f_lower);
            const double prob = std::max(1e-12, p_upper - p_lower);
            m_scratch_dq.setZero();
            const double wi = obs_weight(i);
            if (yi_idx < m_K - 1) {
                m_scratch_dq[yi_idx] += f_upper;
                m_scratch_dq.tail(m_p).noalias() += f_upper * m_eta_sign * m_X.row(i).transpose();
            }
            if (yi_idx > 0) {
                m_scratch_dq[yi_idx - 1] -= f_lower;
                m_scratch_dq.tail(m_p).noalias() -= f_lower * m_eta_sign * m_X.row(i).transpose();
            }
            nll -= wi * std::log(prob);
            grad.noalias() -= wi * m_scratch_dq / prob;
        }
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        const int n_params = params.size();
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n_params, n_params);
        if (!validate_params(params)) {
            H.diagonal().array() = 1e10;
            return H;
        }

        const int n_alpha = m_K - 1;
        const Eigen::VectorXd alpha = params.head(n_alpha);
        const Eigen::VectorXd beta = params.tail(m_p);
        const Eigen::VectorXd eta = m_X * beta;
        double* H_data = H.data();

        // Local arrays are faster than Eigen scratch vectors for these small
        // parameter dimensions and are allocated only on Hessian calls.
        std::vector<double> dq(n_params, 0.0);
        std::vector<double> v_upper(n_params, 0.0);
        std::vector<double> v_lower(n_params, 0.0);

        for (int i = 0; i < m_n; ++i) {
            const int yi_idx = level_index(m_levels, m_y[i]);
            if (yi_idx < 0) continue;

            const double eta_i = eta[i];
            const double wi = obs_weight(i);

            // Compute endpoint quantities — fused helper saves 4 exp per threshold vs separate calls
            double f_upper = 0.0, fp_upper = 0.0, p_upper = 1.0;
            double f_lower = 0.0, fp_lower = 0.0, p_lower = 0.0;
            const bool has_upper = (yi_idx < n_alpha);
            const bool has_lower = (yi_idx > 0);

            if (has_upper) {
                const double z_upper = alpha[yi_idx] + m_eta_sign * eta_i;
                cdf_pdf_fpdf(m_link, z_upper, p_upper, f_upper, fp_upper);
            }
            if (has_lower) {
                const double z_lower = alpha[yi_idx - 1] + m_eta_sign * eta_i;
                cdf_pdf_fpdf(m_link, z_lower, p_lower, f_lower, fp_lower);
            }
            const double prob = std::max(1e-12, p_upper - p_lower);

            // Build v_upper and v_lower as plain arrays
            const double* xi = m_X.data() + i;  // xi[j * m_n] == X(i,j)
            if (has_upper) {
                std::fill(v_upper.begin(), v_upper.end(), 0.0);
                v_upper[yi_idx] = 1.0;
                for (int j = 0; j < m_p; ++j)
                    v_upper[n_alpha + j] = m_eta_sign * xi[j * m_n];
            }
            if (has_lower) {
                std::fill(v_lower.begin(), v_lower.end(), 0.0);
                v_lower[yi_idx - 1] = 1.0;
                for (int j = 0; j < m_p; ++j)
                    v_lower[n_alpha + j] = m_eta_sign * xi[j * m_n];
            }

            // Build dq = f_upper * v_upper - f_lower * v_lower
            std::fill(dq.begin(), dq.end(), 0.0);
            if (has_upper) {
                for (int j = 0; j < n_params; ++j)
                    dq[j] += f_upper * v_upper[j];
            }
            if (has_lower) {
                for (int j = 0; j < n_params; ++j)
                    dq[j] -= f_lower * v_lower[j];
            }

            // Accumulate into H (upper triangle only)
            const double inv_prob = 1.0 / prob;
            const double inv_prob2 = inv_prob * inv_prob;
            for (int c = 0; c < n_params; ++c) {
                for (int r = 0; r <= c; ++r) {
                    double val = wi * inv_prob2 * dq[r] * dq[c];
                    if (has_upper)
                        val -= wi * inv_prob * fp_upper * v_upper[r] * v_upper[c];
                    if (has_lower)
                        val += wi * inv_prob * fp_lower * v_lower[r] * v_lower[c];
                    H_data[r + c * n_params] += val;
                }
            }
        }
        // Reflect upper triangle to lower
        for (int c = 0; c < n_params; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * n_params] = H_data[r + c * n_params];
        return H;
    }

    Eigen::MatrixXd expected_hessian(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        const int n_params = params.size();
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n_params, n_params);
        if (!validate_params(params)) return H;

        const int n_alpha = m_K - 1;
        const Eigen::VectorXd alpha = params.head(n_alpha);
        const Eigen::VectorXd beta = params.tail(m_p);
        const Eigen::VectorXd eta = m_X * beta;
        double* H_data = H.data();

        // Cached per-category gradients for this observation.
        std::vector<std::vector<double>> grad_pi(m_K, std::vector<double>(n_params, 0.0));
        std::vector<double> f(n_alpha, 0.0);
        std::vector<double> p(n_alpha, 0.0);

        for (int i = 0; i < m_n; ++i) {
            const double eta_i = eta[i];
            const double wi = obs_weight(i);
            const double* xi = m_X.data() + i;

            // 1. Compute f_k and p_k (cdf) for all thresholds — fused saves 2 exp per threshold
            for (int k = 0; k < n_alpha; ++k) {
                const double z = alpha[k] + m_eta_sign * eta_i;
                cdf_and_pdf(m_link, z, p[k], f[k]);
            }

            // 2. Compute grad_pi_k for each category k
            // grad_pi_k = \nabla (p_k - p_{k-1}) = f_k * v_k - f_{k-1} * v_{k-1}
            // where v_k = [0...1...0, sign * x]
            for (int k = 0; k < m_K; ++k) {
                std::fill(grad_pi[k].begin(), grad_pi[k].end(), 0.0);
                
                // Contribution from upper threshold p_k
                if (k < n_alpha) {
                    grad_pi[k][k] += f[k];
                    for (int j = 0; j < m_p; ++j)
                        grad_pi[k][n_alpha + j] += f[k] * m_eta_sign * xi[j * m_n];
                }
                
                // Contribution from lower threshold p_{k-1}
                if (k > 0) {
                    grad_pi[k][k - 1] -= f[k - 1];
                    for (int j = 0; j < m_p; ++j)
                        grad_pi[k][n_alpha + j] -= f[k - 1] * m_eta_sign * xi[j * m_n];
                }
            }

            // 3. Accumulate: I += \sum_k (1/pi_k) * grad_pi_k * grad_pi_k^T
            for (int k = 0; k < m_K; ++k) {
                double pi_k = (k == 0) ? p[0] : ((k == m_K - 1) ? (1.0 - p[k - 1]) : (p[k] - p[k - 1]));
                pi_k = std::max(pi_k, 1e-12);
                const double inv_pi = 1.0 / pi_k;
                const double* gk = grad_pi[k].data();
                
                for (int c = 0; c < n_params; ++c) {
                    if (gk[c] == 0.0) continue;
                    for (int r = 0; r <= c; ++r) {
                        if (gk[r] == 0.0) continue;
                        H_data[r + c * n_params] += wi * inv_pi * gk[r] * gk[c];
                    }
                }
            }
        }

        // Reflect
        for (int c = 0; c < n_params; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * n_params] = H_data[r + c * n_params];
        
        return H;
    }
};

} // namespace edi_ordinal

#endif
