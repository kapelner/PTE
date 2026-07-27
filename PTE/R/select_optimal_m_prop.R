#' Select an optimal m-out-of-n bootstrap resample size via the minimum volatility method
#'
#' \code{\link{PTE_bootstrap_inference}}'s \code{m_pow_of_n} argument sets each bootstrap
#' replicate's resample size to \eqn{m = \lfloor n^{\gamma} \rfloor}{m = round(n^gamma)}, which is
#' required to be less than \eqn{n} (i.e. \eqn{\gamma < 1}{gamma < 1}) for the bootstrap to be
#' consistent for our non-smooth policy-value estimands (see \code{\link{PTE_bootstrap_inference}}).
#' There is, however, no single population-optimal \eqn{\gamma}{gamma} -- asymptotic theory only
#' requires \eqn{m \to \infty}{m -> infinity} and \eqn{m / n \to 0}{m / n -> 0}. Politis, Romano &
#' Wolf (\emph{Subsampling}, 1999, Ch. 9) and Bickel & Sakov (2008) address this by proposing the
#' "minimum volatility" (MV) method: compute the bootstrap answer (here, the width of the
#' rate-corrected percentile confidence interval for one of the three PTE estimands) over a grid of
#' candidate \eqn{\gamma}{gamma} (equivalently \eqn{m}{m}) values, then, for each grid point, measure
#' how much that answer varies over a small centered neighborhood of adjacent grid points ("volatility").
#' The chosen \eqn{\gamma}{gamma} is the one whose neighborhood is most stable, i.e. the region of the
#' grid where the bootstrap's answer has stopped changing in response to further tuning of \eqn{m}{m} --
#' this is the practical, data-driven proxy for "\eqn{m}{m} large enough for the CLT-type approximation
#' underlying the bootstrap to have kicked in, but still small enough relative to \eqn{n}{n} for
#' consistency."
#'
#' This function is a calibration step: it runs \code{\link{PTE_bootstrap_inference}} once per grid
#' point (so its total cost is roughly \code{length(m_pow_of_n_grid)} times that of a single call) and
#' does not itself return final inference -- rerun \code{\link{PTE_bootstrap_inference}} at the
#' selected \code{m_pow_of_n_optimal} (with a larger \code{B}) to get the confidence intervals and
#' p-values you report.
#'
#' @references
#' Politis, D.N., Romano, J.P. and Wolf, M. (1999) \emph{Subsampling}. Springer Series in Statistics.
#'
#' Bickel, P.J. and Sakov, A. (2008) On the choice of m in the m out of n bootstrap and confidence
#' bounds for extrema. \emph{Statistica Sinica}, 18(3), 967-985.
#'
#' @param X 						A \eqn{n \times p}{n x p} dataframe of covariates where one column is labeled "treatment" and it
#' 									is a binary vector of treatment allocations in the study. See \code{\link{PTE_bootstrap_inference}}.
#' @param y 						An \eqn{n}-length numeric vector which is the response. See \code{\link{PTE_bootstrap_inference}}.
#' @param regression_type			See \code{\link{PTE_bootstrap_inference}}. Default \code{"continuous"}.
#' @param incidence_metric			See \code{\link{PTE_bootstrap_inference}}. Default \code{"odds_ratio"}.
#' @param personalized_model_build_function 	See \code{\link{PTE_bootstrap_inference}}. Default \code{NULL} (use the built-in default for \code{regression_type}).
#' @param censored					See \code{\link{PTE_bootstrap_inference}}. Required if \code{regression_type} is \code{"survival"}.
#' @param predict_function 			See \code{\link{PTE_bootstrap_inference}}. Default \code{NULL} (use the built-in default).
#' @param difference_function		See \code{\link{PTE_bootstrap_inference}}. Default \code{NULL} (use the built-in default for \code{regression_type}).
#' @param cleanup_mod_function 		See \code{\link{PTE_bootstrap_inference}}. Default \code{NULL} (no cleanup).
#' @param y_higher_is_better 		See \code{\link{PTE_bootstrap_inference}}. Default \code{TRUE}.
#' @param pct_leave_out 			See \code{\link{PTE_bootstrap_inference}}. Default \code{0.10}.
#' @param m_pow_of_n_grid			The grid of candidate \eqn{\gamma}{gamma} exponents to evaluate (each entry is a candidate value
#' 									for \code{\link{PTE_bootstrap_inference}}'s \code{m_pow_of_n} argument). Must have at least
#' 									\code{volatility_window} entries so every interior grid point has a complete neighborhood. The
#' 									default, \code{seq(0.5, 1, by = 0.05)}, brackets the Politis-Romano-Wolf-recommended
#' 									\eqn{\gamma = 0.75}{gamma = 0.75} on both sides, up to (but not below) \eqn{\gamma = 0.5}{gamma = 0.5}
#' 									(\eqn{m = \sqrt{n}}{m = sqrt(n)}, a common lower-rate benchmark) and \eqn{\gamma = 1}{gamma = 1}
#' 									(the classical, here potentially inconsistent, \eqn{n}-out-of-\eqn{n} bootstrap, included only as a
#' 									reference point -- see \code{volatility_window}).
#' @param estimand					Which of the three PTE estimands' confidence interval width to stabilize: \code{"adversarial"},
#' 									\code{"average"} (the default) or \code{"best"}; see \code{\link{PTE_bootstrap_inference}}.
#' @param B 						The number of bootstrap samples to take \emph{at each grid point}. Since the total cost is
#' 									\code{length(m_pow_of_n_grid) * B} bootstrap replicates, this defaults to a much smaller
#' 									\code{250} than \code{\link{PTE_bootstrap_inference}}'s own default of \code{3000} -- this
#' 									function is for calibrating \code{m_pow_of_n}, not for the final reported inference.
#' @param alpha 					Confidence interval size (1 - alpha) used to measure CI width at each grid point. Defaults to 0.05.
#' @param volatility_window			The number of adjacent grid points (must be odd and \eqn{\geq 3}{>= 3}) averaged into the
#' 									volatility measure at each interior grid point. Only grid points with a complete, centered
#' 									window are eligible to be selected (this excludes the \code{(volatility_window - 1) / 2} grid
#' 									points at each end of \code{m_pow_of_n_grid} from selection, since their neighborhoods would
#' 									otherwise be one-sided and their apparent volatility artificially low). Default \code{3}.
#' @param num_cores					See \code{\link{PTE_bootstrap_inference}}. Passed through unchanged to each grid point's call
#' 									(grid points are evaluated sequentially, not in additional parallel, to avoid nested
#' 									parallelism on top of \code{\link{PTE_bootstrap_inference}}'s own).
#'
#' @return 							A list with class \code{"PTE_optimal_m_selection"} containing:
#' 									\describe{
#' 										\item{\code{m_pow_of_n_optimal}}{The selected \eqn{\gamma}{gamma}, ready to pass as
#' 											\code{m_pow_of_n} to \code{\link{PTE_bootstrap_inference}}.}
#' 										\item{\code{m_optimal}}{The corresponding resample size, \code{round(n^m_pow_of_n_optimal)}.}
#' 										\item{\code{grid_table}}{A \code{data.frame} with one row per grid point: \code{m_pow_of_n},
#' 											\code{m}, the point \code{estimate} and \code{ci_width} for \code{estimand}, and its
#' 											\code{volatility} (\code{NA} for the non-eligible boundary points described above).}
#' 										\item{\code{estimand}, \code{B}, \code{alpha}, \code{volatility_window}}{Echoed back for reference.}
#' 									}
#'
#' @author Adam Kapelner
#'
#' @examples
#' \dontrun{
#' 	library(PTE)
#' 	data(continuous_example)
#' 	X = continuous_example$X
#' 	y = continuous_example$y
#'
#' 	m_selection = select_optimal_m_prop(
#' 		X, y,
#' 		regression_type = "continuous",
#' 		m_pow_of_n_grid = seq(0.5, 1, by = 0.05),
#' 		B = 250,
#' 		num_cores = 1
#' 	)
#' 	m_selection$grid_table
#' 	m_selection$m_pow_of_n_optimal
#'
#' 	#now rerun inference at the selected m_pow_of_n with a production-sized B
#' 	pte_results = PTE_bootstrap_inference(
#' 		X, y,
#' 		regression_type = "continuous",
#' 		m_pow_of_n = m_selection$m_pow_of_n_optimal,
#' 		B = 3000
#' 	)
#' 	pte_results
#' }
#' @export
select_optimal_m_prop = function(X, y,
		regression_type = "continuous",
		incidence_metric = "odds_ratio",
		personalized_model_build_function = NULL,
		censored = NULL,
		predict_function = NULL,
		difference_function = NULL,
		cleanup_mod_function = NULL,
		y_higher_is_better = TRUE,
		pct_leave_out = 0.10,
		m_pow_of_n_grid = seq(0.5, 1, by = 0.05),
		estimand = "average",
		B = 250,
		alpha = 0.05,
		volatility_window = 3,
		num_cores = NULL
	){

	#check validity of all values the user input
	checkmate::assertDataFrame(X)
	n = nrow(X)
	checkmate::assertNumeric(y, len = n, .var.name = "y (response vector)")
	checkmate::assertChoice(regression_type, c("continuous", "incidence", "survival"))
	#incidence_metric is only meaningful (and thus only validated) when it will actually be used --
	#it's documented as ignored once a custom difference_function is supplied (mirrors PTE_bootstrap_inference).
	if (regression_type == "incidence" && is.null(difference_function)){
		checkmate::assertChoice(incidence_metric, c("probability_difference", "risk_ratio", "odds_ratio"))
	}
	checkmate::assertFunction(personalized_model_build_function, null.ok = TRUE)
	checkmate::assertFunction(predict_function, null.ok = TRUE)
	checkmate::assertFunction(difference_function, null.ok = TRUE)
	checkmate::assertFunction(cleanup_mod_function, null.ok = TRUE)
	checkmate::assertFlag(y_higher_is_better)
	checkmate::assertNumber(pct_leave_out, lower = 0, upper = 1)
	checkmate::assertNumeric(m_pow_of_n_grid, lower = 0, upper = 1, min.len = 3, any.missing = FALSE, unique = TRUE,
		.var.name = "m_pow_of_n_grid (grid of candidate m-out-of-n bootstrap exponents)")
	checkmate::assertChoice(estimand, c("adversarial", "average", "best"))
	checkmate::assertCount(B, positive = TRUE)
	checkmate::assertNumber(alpha, lower = 0, upper = 1)
	checkmate::assertCount(volatility_window, positive = TRUE)
	if (volatility_window < 3 || volatility_window %% 2 == 0){
		stop("\"volatility_window\" must be an odd integer >= 3 so every interior grid point has a well-defined, centered neighborhood.")
	}
	checkmate::assertNumber(volatility_window, upper = length(m_pow_of_n_grid),
		.var.name = "volatility_window (must be <= length(m_pow_of_n_grid), so at least one grid point has a complete neighborhood)")
	checkmate::assertCount(num_cores, positive = TRUE, null.ok = TRUE)

	m_pow_of_n_grid = sort(unique(m_pow_of_n_grid))
	K = length(m_pow_of_n_grid)

	ci_field = paste0("ci_q_", estimand)
	est_field = paste0("est_q_", estimand)

	#evaluate the bootstrap CI width at each candidate m -- sequential over the grid (not parallelized
	#across grid points) since PTE_bootstrap_inference() already parallelizes across its B replicates
	#via num_cores; parallelizing both loops would just create nested parallelism for no benefit.
	m_vals = round(n ^ m_pow_of_n_grid)
	ci_widths = rep(NA_real_, K)
	estimates = rep(NA_real_, K)
	for (k in 1 : K){
		res_k = PTE_bootstrap_inference(
			X, y,
			regression_type = regression_type,
			incidence_metric = incidence_metric,
			personalized_model_build_function = personalized_model_build_function,
			censored = censored,
			predict_function = predict_function,
			difference_function = difference_function,
			cleanup_mod_function = cleanup_mod_function,
			y_higher_is_better = y_higher_is_better,
			pct_leave_out = pct_leave_out,
			m_pow_of_n = m_pow_of_n_grid[k],
			B = B,
			alpha = alpha,
			num_cores = num_cores
		)
		ci_widths[k] = diff(res_k[[ci_field]])
		estimates[k] = res_k[[est_field]]
	}

	#minimum volatility (MV) method (Politis, Romano & Wolf, 1999, Ch. 9; Bickel & Sakov, 2008): at each
	#interior grid point, measure the spread (sd) of the CI width over a complete, centered window of
	#neighboring grid points, then pick the grid point whose neighborhood is most stable. Boundary points
	#that cannot have a complete centered window are left ineligible (volatility = NA) since a one-sided
	#window's apparent volatility is not comparable to a two-sided one and would bias selection toward the
	#edges of the grid.
	half_window = (volatility_window - 1) %/% 2
	volatility = rep(NA_real_, K)
	eligible = (half_window + 1) : (K - half_window)
	for (k in eligible){
		volatility[k] = sd(ci_widths[(k - half_window) : (k + half_window)])
	}
	if (all(is.na(volatility[eligible]))){
		#every eligible grid point had at least one NA CI width in its window (e.g. bootstrap
		#replicates were mostly/all bad at every candidate m -- a very aggressive m_pow_of_n_grid
		#combined with rare categorical covariate levels can do this). which.min() on an all-NA
		#vector would otherwise fail with an opaque "argument is of length 0" error.
		stop("Every eligible grid point's bootstrap CI width was NA, so no m_pow_of_n could be selected. ",
			"This usually means m = round(n^m_pow_of_n) was too small (relative to n and/or the number ",
			"of levels of a categorical covariate) for most bootstrap replicates to produce a valid fit ",
			"at one or more grid points. Try a m_pow_of_n_grid biased toward larger values, or increase B.")
	}
	k_star = eligible[which.min(volatility[eligible])]

	grid_table = data.frame(
		m_pow_of_n = m_pow_of_n_grid,
		m = m_vals,
		estimate = estimates,
		ci_width = ci_widths,
		volatility = volatility
	)

	return_obj = list(
		m_pow_of_n_optimal = m_pow_of_n_grid[k_star],
		m_optimal = m_vals[k_star],
		grid_table = grid_table,
		estimand = estimand,
		B = B,
		alpha = alpha,
		volatility_window = volatility_window
	)
	class(return_obj) = "PTE_optimal_m_selection"
	return_obj
}
