create_raw_results_matrix = function(n){
	raw_results = as.data.frame(matrix(NA, nrow = n, ncol = 6))
	colnames(raw_results) = c("est_true", "est_counterfactual", "given_tx", "rec_tx", "real_y", "censored")
	raw_results
}