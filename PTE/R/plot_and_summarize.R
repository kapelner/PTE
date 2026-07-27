
#' Plots a summary of the bootstrap samples
#'
#' @param x		 					A \code{PTE_bootstrap_results} model object built via
#' 									running the \code{PTE_bootstrap_inference} function.
#' @param ... 						Other methods passed to plot
#' @method plot PTE_bootstrap_results
#'
#' @examples
#' library(PTE)
#' data(continuous_example)
#' pte_results = PTE_bootstrap_inference(continuous_example$X, continuous_example$y,
#'     regression_type = "continuous", B = 5, num_cores = 1)
#' plot(pte_results)
#'
#' @author Adam Kapelner
#' @export
plot.PTE_bootstrap_results = function(x, ...){
	checkmate::assertClass(x, "PTE_bootstrap_results")
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

	q_df = data.frame(
		value = c(x$q_scores$average, x$q_scores$best),
		panel = factor(
			rep(c("Average I's", "Best I's"), c(length(x$q_scores$average), length(x$q_scores$best))),
			levels = c("Average I's", "Best I's")
		)
	)

	vline_df = data.frame(
		panel = factor(rep(c("Average I's", "Best I's"), each = 6 + 2 * x$run_bca_bootstrap), levels = c("Average I's", "Best I's")),
		xintercept = c(
			x$est_q_average, x$basic_ci_q_average[1], x$basic_ci_q_average[2], x$ci_q_average[1], x$ci_q_average[2],
			if (x$run_bca_bootstrap) c(x$bca_ci_q_average[1], x$bca_ci_q_average[2]),
			x$H_0_mu_equals,
			x$est_q_best, x$basic_ci_q_best[1], x$basic_ci_q_best[2], x$ci_q_best[1], x$ci_q_best[2],
			if (x$run_bca_bootstrap) c(x$bca_ci_q_best[1], x$bca_ci_q_best[2]),
			x$H_0_mu_equals
		),
		type = factor(
			rep(c(
				"Point estimate", "Basic CI", "Basic CI", "Percentile CI", "Percentile CI",
				if (x$run_bca_bootstrap) c("BCa CI", "BCa CI"),
				"Null hypothesis"
			), 2),
			levels = c("Point estimate", "Basic CI", "Percentile CI", "BCa CI", "Null hypothesis")
		)
	)

	line_colors = c(
		"Point estimate" = "forestgreen",
		"Basic CI" = "purple",
		"Percentile CI" = "firebrick3",
		"BCa CI" = "dodgerblue3",
		"Null hypothesis" = "gray"
	)
	line_widths = c(
		"Point estimate" = 1,
		"Basic CI" = 0.5,
		"Percentile CI" = 0.5,
		"BCa CI" = 0.5,
		"Null hypothesis" = 0.5
	)

	p = ggplot2::ggplot(q_df, ggplot2::aes(x = .data$value)) +
		ggplot2::geom_histogram(bins = round(x$B / 3)) +
		ggplot2::geom_vline(
			data = vline_df,
			ggplot2::aes(xintercept = .data$xintercept, color = .data$type, linewidth = .data$type)
		) +
		ggplot2::scale_color_manual(name = NULL, values = line_colors, drop = TRUE) +
		ggplot2::scale_linewidth_manual(name = NULL, values = line_widths, drop = TRUE, guide = "none") +
		ggplot2::facet_wrap(~ panel, ncol = 1) +
		ggplot2::labs(x = xlab, y = "Count") +
		ggplot2::theme_minimal()

	print(p)
	invisible(NULL)
}

#' Prints a summary of the model to the console
#'
#' @param x 								A \code{PTE_bootstrap_results} model object built via
#' 											running the \code{PTE_bootstrap_inference} function.
#' @param ...								Other methods passed to print
#' @method print PTE_bootstrap_results
#'
#' @examples
#' library(PTE)
#' data(continuous_example)
#' pte_results = PTE_bootstrap_inference(continuous_example$X, continuous_example$y,
#'     regression_type = "continuous", B = 5, num_cores = 1)
#' print(pte_results)
#'
#' @author Adam Kapelner
#' @export
print.PTE_bootstrap_results = function(x, ...){
	checkmate::assertClass(x, "PTE_bootstrap_results")
	if (x$display_adversarial_score){
		cat("    I_adversarial observed est = ",
				round(x$observed_q_scores$adversarial, 3),
				",  pctile p-val = ",
				round(x$p_val_adversarial, 3),
				", \n      ",
				round(100 * (1 - x$alpha), 1),
				"% CI's: basic = [",
				round(x$basic_ci_q_adversarial[1], 3),
				", ",
				round(x$basic_ci_q_adversarial[2], 3),
				"], ",
				"pctile = [",
				round(x$ci_q_adversarial[1], 3),
				", ",
				round(x$ci_q_adversarial[2], 3),
				"], ",
				ifelse(x$run_bca_bootstrap,
						paste(
								"BCa = [",
								round(x$bca_ci_q_adversarial[1], 3),
								", ",
								round(x$bca_ci_q_adversarial[2], 3),
								"],", sep = ""),
						""),
				"\n", sep = "")
	}
	cat("    I_random observed_est = ",
			round(x$observed_q_scores$average, 3),
			",  pctile p-val = ",
			round(x$p_val_average, 3),
			", \n      ",
			round(100 * (1 - x$alpha), 1),
			"% CI's: basic = [",
			round(x$basic_ci_q_average[1], 3),
			", ",
			round(x$basic_ci_q_average[2], 3),
			"], ",
			"pctile = [",
			round(x$ci_q_average[1], 3),
			", ",
			round(x$ci_q_average[2], 3),
			"], ",
			ifelse(x$run_bca_bootstrap,
				paste(
					"BCa = [",
					round(x$bca_ci_q_average[1], 3),
					", ",
					round(x$bca_ci_q_average[2], 3),
					"],", sep = ""),
				""),
			sep = "")
	cat("\n    I_best observed_est = ",
			round(x$observed_q_scores$best, 3),
			",  p val = ",
			round(x$p_val_best, 3),
			", \n      ",
			round(100 * (1 - x$alpha), 1),
			"% CI's: basic = [",
			round(x$basic_ci_q_best[1], 3),
			", ",
			round(x$basic_ci_q_best[2], 3),
			"], ",
			"pctile = [",
			round(x$ci_q_best[1], 3),
			", ",
			round(x$ci_q_best[2], 3),
			"], ",
			ifelse(x$run_bca_bootstrap,
				paste(
					"BCa = [",
					round(x$bca_ci_q_best[1], 3),
					", ",
					round(x$bca_ci_q_best[2], 3),
					"]", sep = ""),
				""),
			sep = "")
	cat("\n")
}

#' Prints a summary of the model to the console
#'
#' @param object 					A \code{PTE_bootstrap_results} model object built via
#' 									running the \code{PTE_bootstrap_inference} function.
#' @param ... 						Other methods passed to summary
#' @method summary PTE_bootstrap_results
#'
#' @examples
#' library(PTE)
#' data(continuous_example)
#' pte_results = PTE_bootstrap_inference(continuous_example$X, continuous_example$y,
#'     regression_type = "continuous", B = 5, num_cores = 1)
#' summary(pte_results)
#'
#' @author Adam Kapelner
#' @export
summary.PTE_bootstrap_results = function(object, ...){
	checkmate::assertClass(object, "PTE_bootstrap_results")
	print(object)
}
