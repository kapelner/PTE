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

// [[Rcpp::export]]
Rcpp::List fast_default_continuous_cv_run_cpp(SEXP X_obs_sexp,
                                              SEXP y_full_sexp,
                                              SEXP X_tx0_sexp,
                                              SEXP X_tx1_sexp,
                                              SEXP treatment_sexp,
                                              SEXP censored_sexp,
                                              SEXP boot_idx_sexp,
                                              SEXP begin_cutoffs_sexp,
                                              SEXP end_cutoffs_sexp,
                                              bool y_higher_is_better) {
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
	Eigen::MatrixXd full_xtx = Eigen::MatrixXd::Zero(p, p);
	Eigen::VectorXd full_xty = Eigen::VectorXd::Zero(p);
	for (int pos = 0; pos < n_out; ++pos) {
		const int row = boot_idx_r[pos] - 1;
		if (row < 0 || row >= n_full) return List::create(Named("ok") = false);
		rows[pos] = row;
		const double yi = y_full[row];
		full_xty.noalias() += X_obs.row(row).transpose() * yi;
		full_xtx.selfadjointView<Eigen::Upper>().rankUpdate(X_obs.row(row).transpose());
	}
	full_xtx.triangularView<Eigen::Lower>() = full_xtx.transpose();

	for (int fold = 0; fold < n_folds; ++fold) {
		const int begin = begin_r[fold] - 1;
		const int end = end_r[fold] - 1;
		if (begin < 0 || end < begin || end >= n_out) return List::create(Named("ok") = false);

		Eigen::MatrixXd xtx = full_xtx;
		Eigen::VectorXd xty = full_xty;
		for (int pos = begin; pos <= end; ++pos) {
			const int row = rows[pos];
			const double yi = y_full[row];
			xty.noalias() -= X_obs.row(row).transpose() * yi;
			xtx.selfadjointView<Eigen::Upper>().rankUpdate(X_obs.row(row).transpose(), -1.0);
		}
		xtx.triangularView<Eigen::Lower>() = xtx.transpose();

		Eigen::LLT<Eigen::MatrixXd> llt(xtx);
		Eigen::VectorXd beta;
		if (llt.info() == Eigen::Success) {
			beta = llt.solve(xty);
		} else {
			Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(xtx);
			beta = cod.solve(xty);
		}
		if (!beta.allFinite()) return List::create(Named("ok") = false);

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
Rcpp::List default_continuous_q_scores_cpp(SEXP y_sexp,
                                           SEXP given_tx_sexp,
                                           SEXP rec_tx_sexp,
                                           bool y_higher_is_better) {
	NumericVector y_r(y_sexp);
	NumericVector given_tx_r(given_tx_sexp);
	NumericVector rec_tx_r(rec_tx_sexp);
	if (given_tx_r.size() != y_r.size() || rec_tx_r.size() != y_r.size()) {
		stop("q-score vectors must have the same length");
	}
	DefaultQScores q = default_continuous_or_probability_q_scores_cpp(
		y_r.begin(), given_tx_r.begin(), rec_tx_r.begin(), y_r.size(), y_higher_is_better);
	return List::create(
		Named("q_adversarial") = q.adversarial,
		Named("q_average") = q.average,
		Named("q_best") = q.best,
		Named("is_bad") = q.is_bad
	);
}

// [[Rcpp::export]]
Rcpp::List default_incidence_q_scores_cpp(SEXP y_sexp,
                                          SEXP given_tx_sexp,
                                          SEXP rec_tx_sexp,
                                          bool y_higher_is_better,
                                          std::string incidence_metric) {
	NumericVector y_r(y_sexp);
	NumericVector given_tx_r(given_tx_sexp);
	NumericVector rec_tx_r(rec_tx_sexp);
	if (given_tx_r.size() != y_r.size() || rec_tx_r.size() != y_r.size()) {
		stop("q-score vectors must have the same length");
	}

	DefaultQScores q;
	if (incidence_metric == "probability_difference") {
		q = default_continuous_or_probability_q_scores_cpp(
			y_r.begin(), given_tx_r.begin(), rec_tx_r.begin(), y_r.size(), y_higher_is_better);
	} else {
		q = default_incidence_ratio_q_scores_cpp(
			y_r.begin(), given_tx_r.begin(), rec_tx_r.begin(), y_r.size(),
			y_higher_is_better, incidence_metric == "risk_ratio" ? 1 : 2);
	}
	return List::create(
		Named("q_adversarial") = q.adversarial,
		Named("q_average") = q.average,
		Named("q_best") = q.best,
		Named("is_bad") = q.is_bad
	);
}

// [[Rcpp::export]]
Rcpp::List fast_default_continuous_bootstrap_q_cpp(SEXP X_obs_sexp,
                                                   SEXP y_full_sexp,
                                                   SEXP X_tx0_sexp,
                                                   SEXP X_tx1_sexp,
                                                   SEXP treatment_sexp,
                                                   SEXP censored_sexp,
                                                   int B,
                                                   SEXP begin_cutoffs_sexp,
                                                   SEXP end_cutoffs_sexp,
                                                   bool y_higher_is_better) {
	Rcpp::RNGScope rng_scope;
	NumericMatrix X_obs_r(X_obs_sexp);
	NumericVector y_full_r(y_full_sexp);
	NumericMatrix X_tx0_r(X_tx0_sexp);
	NumericMatrix X_tx1_r(X_tx1_sexp);
	NumericVector treatment_r(treatment_sexp);
	NumericVector censored_r(censored_sexp);
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

	NumericVector q_adv(B, NA_REAL), q_avg(B, NA_REAL), q_best(B, NA_REAL);
	int num_bad = 0;
	std::vector<int> rows(n_out);
	std::vector<double> real_y(n_out), given_tx(n_out), rec_tx(n_out), censored(n_out);

	for (int b = 0; b < B; ++b) {
		Eigen::MatrixXd full_xtx = Eigen::MatrixXd::Zero(p, p);
		Eigen::VectorXd full_xty = Eigen::VectorXd::Zero(p);
		for (int pos = 0; pos < n_out; ++pos) {
			const int row = std::min(n_full - 1, static_cast<int>(std::floor(R::unif_rand() * n_full)));
			rows[pos] = row;
			const double yi = y_full[row];
			full_xty.noalias() += X_obs.row(row).transpose() * yi;
			full_xtx.selfadjointView<Eigen::Upper>().rankUpdate(X_obs.row(row).transpose());
		}
		full_xtx.triangularView<Eigen::Lower>() = full_xtx.transpose();

		bool ok = true;
		for (int fold = 0; fold < n_folds && ok; ++fold) {
			const int begin = begin_r[fold] - 1;
			const int end = end_r[fold] - 1;
			if (begin < 0 || end < begin || end >= n_out) { ok = false; break; }

			Eigen::MatrixXd xtx = full_xtx;
			Eigen::VectorXd xty = full_xty;
			for (int pos = begin; pos <= end; ++pos) {
				const int row = rows[pos];
				const double yi = y_full[row];
				xty.noalias() -= X_obs.row(row).transpose() * yi;
				xtx.selfadjointView<Eigen::Upper>().rankUpdate(X_obs.row(row).transpose(), -1.0);
			}
			xtx.triangularView<Eigen::Lower>() = xtx.transpose();

			Eigen::LLT<Eigen::MatrixXd> llt(xtx);
			Eigen::VectorXd beta;
			if (llt.info() == Eigen::Success) {
				beta = llt.solve(xty);
			} else {
				Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(xtx);
				beta = cod.solve(xty);
			}
			if (!beta.allFinite()) { ok = false; break; }

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
				censored[pos] = censored_r[row];
			}
		}

		DefaultQScores q;
		if (ok) {
			q = default_continuous_or_probability_q_scores_cpp(
				real_y.data(), given_tx.data(), rec_tx.data(), n_out, y_higher_is_better);
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
