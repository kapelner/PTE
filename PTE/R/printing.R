#' Prints a summary of the model to the console
#' 
#' @param x 								A \code{PTE_bootstrap_results} model object built via
#' 											running the \code{PTE_bootstrap_inference} function.
#' @param ...								Other methods passed to print
#' @method print PTE_bootstrap_results
#' 
#' @author Adam Kapelner
#' @export
print.PTE_bootstrap_results = function(x, ...){
	if (x$display_adversarial_score){
		cat("    I_adversarial observed est = ", 
				round(x$observed_q_scores$adversarial, 3), 
				",  p val = ",
				round(x$p_val_adversarial, 3), 
				", \n      ",			
				round(100 * (1 - x$alpha), 1),
				"% CI's: pctile = [",
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
			",  p val = ", 
			round(x$p_val_average, 3),
			", \n      ",
			round(100 * (1 - x$alpha), 1),
			"% CI's: pctile = [",
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
			"% CI's: pctile = [",
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
#' @author Adam Kapelner
#' @export
summary.PTE_bootstrap_results = function(object, ...){
	print(object)
}