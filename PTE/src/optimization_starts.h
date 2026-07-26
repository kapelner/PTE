#ifndef EDI_OPTIMIZATION_STARTS_H
#define EDI_OPTIMIZATION_STARTS_H

#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <algorithm>

namespace edi_opt {

// Robust OLS solver using Complete Orthogonal Decomposition to handle rank-deficiency.
inline bool robust_ols_solve(const Eigen::MatrixXd& X,
                             const Eigen::VectorXd& y,
                             Eigen::VectorXd& beta_out) {
    if (X.rows() == 0 || X.cols() == 0 || y.size() != X.rows() ||
        !X.allFinite() || !y.allFinite()) {
        beta_out = Eigen::VectorXd::Zero(X.cols());
        return false;
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
    beta_out = qr.solve(y);
    return beta_out.allFinite();
}

// Logistic smart start: OLS on logit-transformed means.
inline Eigen::VectorXd logistic_smart_cold_start(const Eigen::MatrixXd& X,
                                           const Eigen::VectorXd& y) {
    const int n = X.rows();
    const int p = X.cols();
    Eigen::VectorXd y_adj(n);
    for (int i = 0; i < n; ++i) {
        // mustart = (y + 0.5) / 2
        // logit(mustart) = log((y + 0.5) / (1.5 - y))
        y_adj[i] = std::log((y[i] + 0.5) / (1.5 - y[i]));
    }
    Eigen::VectorXd beta;
    if (robust_ols_solve(X, y_adj, beta)) return beta;
    return Eigen::VectorXd::Zero(p);
}

// Poisson smart start: OLS on log-transformed counts.
inline Eigen::VectorXd poisson_smart_cold_start(const Eigen::MatrixXd& X,
                                          const Eigen::VectorXd& y) {
    const int n = X.rows();
    const int p = X.cols();
    // mustart = y + 0.1
    Eigen::VectorXd y_adj = (y.array() + 0.1).log().matrix();
    Eigen::VectorXd beta;
    if (robust_ols_solve(X, y_adj, beta)) return beta;
    return Eigen::VectorXd::Zero(p);
}

// Smart Hessian for Logistic: X^T * diag(mu*(1-mu)) * X
inline Eigen::MatrixXd logistic_smart_hessian(const Eigen::MatrixXd& X,
                                              const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    // plogis stable
    Eigen::VectorXd mu = (1.0 / (1.0 + (-eta.array().max(-20.0).min(20.0)).exp())).matrix();
    Eigen::VectorXd w = (mu.array() * (1.0 - mu.array())).matrix();
    return X.transpose() * w.asDiagonal() * X;
}

// Smart Hessian for Poisson: X^T * diag(mu) * X
inline Eigen::MatrixXd poisson_smart_hessian(const Eigen::MatrixXd& X,
                                             const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().max(-20.0).min(700.0).exp().matrix();
    return X.transpose() * mu.asDiagonal() * X;
}

// Smart Hessian for Negative Binomial
inline Eigen::MatrixXd negbin_smart_hessian(const Eigen::MatrixXd& X,
                                            const Eigen::VectorXd& params,
                                            const Eigen::VectorXi& y) {
    const int n = X.rows();
    const int p = X.cols();
    Eigen::VectorXd beta = params.head(p);
    double theta = std::exp(params[p]);
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp();
    
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(p + 1, p + 1);
    for (int i = 0; i < n; ++i) {
        double mu_i = mu[i];
        double yi = static_cast<double>(y[i]);
        double denom = theta + mu_i;
        Eigen::VectorXd x = X.row(i).transpose();

        double w_beta = mu_i * theta * (yi + theta) / (denom * denom);
        H.topLeftCorner(p, p).noalias() += w_beta * (x * x.transpose());

        double d_score_beta_d_log_theta = theta * mu_i * (yi - mu_i) / (denom * denom);
        H.topRightCorner(p, 1).noalias() -= d_score_beta_d_log_theta * x;

        double A = R::digamma(yi + theta) - R::digamma(theta) +
            std::log(theta) - std::log(denom) + 1.0 - (yi + theta) / denom;
        double dA_dtheta = R::trigamma(yi + theta) - R::trigamma(theta) +
            1.0 / theta - 1.0 / denom + (yi - mu_i) / (denom * denom);
        H(p, p) -= theta * A + theta * theta * dA_dtheta;
    }
    H.bottomLeftCorner(1, p) = H.topRightCorner(p, 1).transpose();
    return H;
}

// ZINB smart start: Combine component-wise smart starts.
inline Eigen::VectorXd zinb_smart_cold_start(const Eigen::MatrixXd& Xc,
                                        const Eigen::MatrixXd& Xz,
                                        const Eigen::VectorXd& y) {
    const int pc = Xc.cols();
    const int pz = Xz.cols();
    Eigen::VectorXd par(pc + pz + 1);
    
    // Cond component: Poisson start
    par.head(pc) = poisson_smart_cold_start(Xc, y);
    
    // ZI component: Logistic start on (y == 0)
    Eigen::VectorXd y_is_zero = (y.array() == 0.0).cast<double>();
    par.segment(pc, pz) = logistic_smart_cold_start(Xz, y_is_zero);
    
    // log_theta = 0
    par[pc + pz] = 0.0;
    return par;
}

// ZAP smart start: Combine component-wise smart starts.
inline Eigen::VectorXd zap_smart_cold_start(const Eigen::MatrixXd& Xc,
                                       const Eigen::MatrixXd& Xz,
                                       const Eigen::VectorXd& y) {
    const int pc = Xc.cols();
    const int pz = Xz.cols();
    Eigen::VectorXd par(pc + pz);
    
    // Cond component: Poisson start on non-zeros
    std::vector<int> nz_idx;
    for (int i = 0; i < y.size(); ++i) if (y[i] > 0) nz_idx.push_back(i);
    
    if (nz_idx.size() > (size_t)pc) {
        Eigen::MatrixXd Xc_nz(nz_idx.size(), pc);
        Eigen::VectorXd y_nz(nz_idx.size());
        for (size_t i = 0; i < nz_idx.size(); ++i) {
            Xc_nz.row(i) = Xc.row(nz_idx[i]);
            y_nz[i] = y[nz_idx[i]];
        }
        par.head(pc) = poisson_smart_cold_start(Xc_nz, y_nz);
    } else {
        par.head(pc) = poisson_smart_cold_start(Xc, y);
    }
    
    // Hurdle component: Logistic start on (y == 0)
    Eigen::VectorXd y_is_zero = (y.array() == 0.0).cast<double>();
    par.tail(pz) = logistic_smart_cold_start(Xz, y_is_zero);
    
    return par;
}

// Weibull AFT smart start: OLS on uncensored rows.
inline Eigen::VectorXd weibull_smart_cold_start(const Eigen::MatrixXd& X,
                                          const Eigen::VectorXd& y,
                                          const Eigen::VectorXd& dead) {
    const int n = X.rows();
    const int p = X.cols();
    std::vector<int> uc_idx;
    for (int i = 0; i < n; ++i) if (dead[i] > 0) uc_idx.push_back(i);
    
    Eigen::VectorXd par = Eigen::VectorXd::Zero(p + 1);
    if (uc_idx.size() > (size_t)p) {
        Eigen::MatrixXd X_uc(uc_idx.size(), p);
        Eigen::VectorXd y_uc(uc_idx.size());
        for (size_t i = 0; i < uc_idx.size(); ++i) {
            X_uc.row(i) = X.row(uc_idx[i]);
            y_uc[i] = std::log(std::max(y[uc_idx[i]], 1e-8));
        }
        Eigen::VectorXd beta;
        if (robust_ols_solve(X_uc, y_uc, beta)) {
            par.head(p) = beta;
            // log_sigma start
            Eigen::VectorXd res = y_uc - X_uc * beta;
            double s2 = res.squaredNorm() / static_cast<double>(uc_idx.size() - p);
            par[p] = std::log(std::sqrt(std::max(1e-4, s2)));
        }
    } else {
        par[0] = std::log(y.mean());
    }
    return par;
}

} // namespace edi_opt

#endif
