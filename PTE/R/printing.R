#' Prints a summary of the model to the console
#' 
#' @param x 								A \code{PTE_bootstrap_results} model object built via
#' 											running the \code{PTE_bootstrap_inference} function.
#' @param ...								Other methods passed to print
#' @examples beta0 = 1
#'	beta1 = -1
#'	gamma0 = 0
#'	gamma1 = sqrt(2 * pi)
#'	mu_x = 0
#'	sigsq_x = 1
#'	sigsq_e = 1
#'	num_boot = 20 #for speed only
#'	n = 50 #for speed only
#'	
#'	x = sort(rnorm(n, mu_x, sigsq_x))
#'	noise = rnorm(n, 0, sigsq_e)
#'	
#'	treatment = sample(c(rep(1, n / 2), rep(0, n / 2)))
#'	y = beta0 + beta1 * x + treatment * (gamma0 + gamma1 * x) + noise
#'	
#'	X = data.frame(treatment, x)
#'	
#'	res = bootstrap_inference(X, y,
#'			"lm(y ~ . + treatment * ., data = Xyleft)",
#'			num_cores = 1,
#'			B = num_boot, 
#'			plot = FALSE)
#'	print(res)
#' 
#' @method print PTE_bootstrap_results
#' 
#' @author Adam Kapelner
#' @export
print.PTE_bootstrap_results = function(x, ...){
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
					"],", 
					""), sep = ""),  
			sep = "")
	cat("\n    I_random observed_est = ", 
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
#' @examples beta0 = 1
#'	beta1 = -1
#'	gamma0 = 0
#'	gamma1 = sqrt(2 * pi)
#'	mu_x = 0
#'	sigsq_x = 1
#'	sigsq_e = 1
#'	num_boot = 20 #for speed only
#'	n = 50 #for speed only
#'	
#'	x = sort(rnorm(n, mu_x, sigsq_x))
#'	noise = rnorm(n, 0, sigsq_e)
#'	
#'	treatment = sample(c(rep(1, n / 2), rep(0, n / 2)))
#'	y = beta0 + beta1 * x + treatment * (gamma0 + gamma1 * x) + noise
#'	
#'	X = data.frame(treatment, x)
#'	
#'	res = bootstrap_inference(X, y,
#'			"lm(y ~ . + treatment * ., data = Xyleft)",
#'			num_cores = 1,
#'			B = num_boot, 
#'			plot = FALSE)
#'	print(res)
#' 
#' @method summary PTE_bootstrap_results
#' 
#' @author Adam Kapelner
#' @export
summary.PTE_bootstrap_results = function(object, ...){
	print(object)
}