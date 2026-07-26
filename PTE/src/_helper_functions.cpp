#include "_helper_functions.h"
using namespace Rcpp;


//typedef Eigen::Map<Eigen::MatrixXd> MapMat;
//typedef Eigen::Map<Eigen::VectorXd> MapVec;

double compute_diagonal_inverse_entry(const Eigen::Ref<const Eigen::MatrixXd>& M, int j) {
	Eigen::VectorXd b = Eigen::VectorXd::Unit(M.rows(), j - 1);

	// Prefer a direct LDLT decomposition for well-behaved symmetric systems.
	Eigen::LDLT<Eigen::MatrixXd> ldlt(M);

	if (ldlt.info() == Eigen::Success) {
		Eigen::VectorXd x = ldlt.solve(b);
		if (x.allFinite()) {
			return x(j - 1);
		}
	}

	// Fall back to a rank-revealing solve for near-singular or indefinite systems.
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(M);
	if (qr.rank() == 0) {
		return NA_REAL;
	}

	Eigen::VectorXd x = qr.solve(b);
	if (!x.allFinite()) {
		return NA_REAL;
	}

	return x(j - 1);
}

// [[Rcpp::export]]
Eigen::MatrixXd eigen_Xt_times_X_cpp(SEXP X_sexp) {
	NumericMatrix X_r(X_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	return X.transpose() * X;
}

// [[Rcpp::export]]
double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(SEXP M_sexp, int j) {
	NumericMatrix M_r(M_sexp);
	Eigen::Map<const Eigen::MatrixXd> M(M_r.begin(), M_r.nrow(), M_r.ncol());
    return compute_diagonal_inverse_entry(M, j);
}

// [[Rcpp::export]]
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(SEXP X_sexp, SEXP w_sexp) {
	NumericMatrix X_r(X_sexp);
	NumericVector w_r(w_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> w(w_r.begin(), w_r.size());
	return weighted_crossprod(X, w);
}

// [[Rcpp::export]]
Rcpp::List likelihood_ratio_test_from_negloglik_cpp(double unrestricted_neg_loglik,
                                                    double null_neg_loglik,
                                                    int df = 1) {
	return likelihood_ratio_test_from_negloglik(unrestricted_neg_loglik, null_neg_loglik, df);
}

// [[Rcpp::export]]
Rcpp::List score_test_from_score_information_cpp(SEXP score_sexp,
                                                 SEXP information_sexp,
                                                 int tested_idx) {
	NumericVector score_r(score_sexp);
	NumericMatrix information_r(information_sexp);
	Eigen::Map<const Eigen::VectorXd> score(score_r.begin(), score_r.size());
	Eigen::Map<const Eigen::MatrixXd> information(information_r.begin(), information_r.nrow(), information_r.ncol());
	return score_test_from_score_information(score, information, tested_idx);
}

// [[Rcpp::export]]
Rcpp::List gradient_test_from_restricted_score_cpp(SEXP score_sexp,
                                                   double unrestricted_estimate,
                                                   double null_value,
                                                   int tested_idx) {
	NumericVector score_r(score_sexp);
	Eigen::Map<const Eigen::VectorXd> score(score_r.begin(), score_r.size());
	return gradient_test_from_restricted_score(score, unrestricted_estimate, null_value, tested_idx);
}

// [[Rcpp::export]]
double mean_cpp(SEXP x_sexp) {
	NumericVector x_r(x_sexp);
	Eigen::Map<const Eigen::VectorXd> x(x_r.begin(), x_r.size());
	if (x.size() == 0) {
	return NA_REAL;
	}
	return x.mean();
}

// [[Rcpp::export]]
double var_cpp(SEXP x_sexp) {
	NumericVector x_r(x_sexp);
	Eigen::Map<const Eigen::VectorXd> x(x_r.begin(), x_r.size());
	if (x.size() <= 1) {
	return NA_REAL;
	}
	return (x.array() - x.mean()).square().sum() / (x.size() - 1);
}
