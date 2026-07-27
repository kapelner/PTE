#include "_helper_functions.h"
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

//' Kaplan-Meier median survival time for a single group
//'
//' A minimal, compiled replacement for \code{summary(survival::survfit(Surv(y, dead) ~ 1))$table['median']}
//' -- computes only the median (no confidence intervals, no full survival table), which is all
//' \code{create_PTE_results_object} ever needs. See \code{get_survival_stat_diff} for the
//' treatment-vs-control difference used directly by the package.
//'
//' @param y_sexp Numeric vector of survival/censoring times.
//' @param dead_sexp Integer vector of event indicators (1 = event/death, 0 = censored).
//' @return The KM median survival time, or \code{Inf} if the curve never drops below 0.5.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double get_survival_stat_for_group_cpp(SEXP y_sexp, SEXP dead_sexp) {
	NumericVector y_r(y_sexp);
	IntegerVector dead_r(dead_sexp);
	std::vector<double> y(y_r.begin(), y_r.end());
	std::vector<int> dead(dead_r.begin(), dead_r.end());
	return km_median_for_group_cpp_internal(y, dead);
}

//' Difference in Kaplan-Meier median survival time between two groups
//'
//' A minimal, compiled replacement for the package's previous
//' \code{summary(survival::survfit(Surv(y, dead) ~ w))$table[, 'median']} pattern, which paid for
//' the full \code{survfit}/\code{summary.survfit} machinery (confidence intervals, the full
//' survival table, formula/factor bookkeeping) purely to extract two numbers. This computes just
//' the two group medians directly and returns their difference.
//'
//' @param y_sexp Numeric vector of survival/censoring times (both groups concatenated).
//' @param dead_sexp Integer vector of event indicators (1 = event/death, 0 = censored).
//' @param w_sexp Integer vector of group membership (1 = first/"treatment" group, 0 = second/"control" group).
//' @return \code{median(group w == 1) - median(group w == 0)}. \code{NA} if either group's median
//' is not reached (the KM curve never drops below 0.5) or either group is empty, matching the
//' \code{NA} that \code{summary.survfit} would have produced in the same situation.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double get_survival_stat_diff_cpp(SEXP y_sexp, SEXP dead_sexp, SEXP w_sexp) {
	NumericVector y_r(y_sexp);
	IntegerVector dead_r(dead_sexp);
	IntegerVector w_r(w_sexp);
	int n = y_r.size();

	std::vector<double> y0, y1;
	std::vector<int> dead0, dead1;
	for (int i = 0; i < n; ++i) {
		if (w_r[i] == 0) {
			y0.push_back(y_r[i]);
			dead0.push_back(dead_r[i]);
		} else {
			y1.push_back(y_r[i]);
			dead1.push_back(dead_r[i]);
		}
	}

	double stat0 = km_median_for_group_cpp_internal(y0, dead0);
	double stat1 = km_median_for_group_cpp_internal(y1, dead1);

	if (!R_FINITE(stat0) || !R_FINITE(stat1)) return NA_REAL; //either median unreached (Inf) or group empty (NA) -> NA, matching summary.survfit
	return stat1 - stat0;
}

// [[Rcpp::export]]
Rcpp::List default_survival_q_scores_cpp(SEXP y_sexp,
                                         SEXP dead_sexp,
                                         SEXP given_tx_sexp,
                                         SEXP rec_tx_sexp,
                                         bool y_higher_is_better) {
	NumericVector y_r(y_sexp);
	NumericVector dead_r(dead_sexp);
	NumericVector given_tx_r(given_tx_sexp);
	NumericVector rec_tx_r(rec_tx_sexp);
	if (dead_r.size() != y_r.size() || given_tx_r.size() != y_r.size() || rec_tx_r.size() != y_r.size()) {
		stop("q-score vectors must have the same length");
	}
	DefaultQScores q = default_survival_q_scores_cpp(
		y_r.begin(), dead_r.begin(), given_tx_r.begin(), rec_tx_r.begin(), y_r.size(), y_higher_is_better);
	return List::create(
		Named("q_adversarial") = q.adversarial,
		Named("q_average") = q.average,
		Named("q_best") = q.best,
		Named("is_bad") = q.is_bad
	);
}
