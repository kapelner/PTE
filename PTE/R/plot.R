
#' Plots a summary of the bootstrap samples.
#' 
#' @param x		 					A \code{PTE_bootstrap_results} model object built via
#' 									running the \code{PTE_bootstrap_inference} function.
#' @param ... 						Other methods passed to plot
#' @method plot PTE_bootstrap_results
#' 
#' @author Adam Kapelner
#' @export
plot.PTE_bootstrap_results = function(x, ...){
	if (x$regression_type == "continuous"){
		xlab = "I (average response difference)"
	} else if (x$regression_type == "survival"){
		xlab = "I (average median survival difference)"
	} else if (x$incidence_metric == "probability_difference"){
		xlab = "I (average probability difference)"
	} else if (x$incidence_metric == "risk_ratio"){
		xlab = "I (average risk ratio)"
	} else if (x$incidence_metric == "odds_ratio"){
		xlab = "I (average odds ratio)"
	}
	#display params
	min_q = min(x$q_scores$average, x$q_scores$best)
	max_q = max(x$q_scores$average, x$q_scores$best)
	par(mfrow = c(2, 1))
	#first plot
	hist(x$q_scores$average, br = x$B / 3, xlab = xlab, xlim = c(min_q, max_q), main = "Average I's")
	abline(v = x$est_q_average, col = "forestgreen", lwd = 3)
	abline(v = x$ci_q_average[1], col = "firebrick3", lwd = 1)
	abline(v = x$ci_q_average[2], col = "firebrick3", lwd = 1)
	if (x$run_bca_bootstrap){
		abline(v = x$bca_ci_q_average[1], col = "dodgerblue3", lwd = 1)
		abline(v = x$bca_ci_q_average[2], col = "dodgerblue3", lwd = 1)
	}
	abline(v = x$H_0_mu_equals, col = "gray")
	#second plot
	hist(x$q_scores$best, br = x$B / 3, xlab = xlab, xlim = c(min_q, max_q), main = "Best I's")
	abline(v = x$est_q_best, col = "forestgreen", lwd = 3)
	abline(v = x$ci_q_best[1], col = "firebrick3", lwd = 1)
	abline(v = x$ci_q_best[2], col = "firebrick3", lwd = 1)
	if (x$run_bca_bootstrap){
		abline(v = x$bca_ci_q_best[1], col = "dodgerblue3", lwd = 1)
		abline(v = x$bca_ci_q_best[2], col = "dodgerblue3", lwd = 1)
	}
	abline(v = x$H_0_mu_equals, col = "gray")
}