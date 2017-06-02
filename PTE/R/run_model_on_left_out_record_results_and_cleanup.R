run_model_on_left_out_record_results_and_cleanup = function(
		Xy, 
		regression_type,
		y_higher_is_better,
		leave_outs_to_be_predicted, 
		train_on_all_except_these,
		personalized_model_build_function,
		predict_function,
		cleanup_mod_function,
		full_verbose,
		verbose){
	
	#the left one out matrix has n-1 rows and will be considered the "training data"
	Xytrain = Xy[-train_on_all_except_these, ]
	
	#pull out the record of the left-one-out subject
	Xyleftout = Xy[leave_outs_to_be_predicted, 1 : (ncol(Xy) - 1)] #leave out y
	
	if (regression_type != "survival"){
		Xytrain$censored = NULL
	}
	#build the model via the user-specified string
	mod = personalized_model_build_function(Xytrain) #this function makes use of the "Xyleft" object
	if (full_verbose){
		print(summary(mod))
	}
	
	#also take note of what actually happened to this subject in the experiment
	real_ys = Xy[leave_outs_to_be_predicted, ncol(Xy)]
	orig_trts = Xyleftout$treatment
	
	#now evaluate the left-one-out subject on the model for both his true treatment and his counterfactual
	Xyleftout$treatment = 0
	yhatTx0s = predict_function(mod, Xyleftout)
	Xyleftout$treatment = 1
	yhatTx1s = predict_function(mod, Xyleftout)
#	cat("yhatTx0s", yhatTx0s, "\nyhatTx1s", yhatTx1s, "\n")
	
	#give the user some indication of progress if they want to see it
	if (full_verbose){
		cat("model #", leave_outs_to_be_predicted, "/", nrow(Xy), " yhatTx0:1 = ", round(yhatTx0s, 2), " : ", round(yhatTx1s, 2), "\n", sep = "")
	} else if (verbose){
		cat(".")
	}
	
	#if the models need to be cleaned up in some way, do it now before the next iteration of the leave-one-out
	if (!is.null(cleanup_mod_function)){
		cleanup_mod_function()
	}
	
	#tabulate the result for the prediction on this left one out model
	res = matrix(NA, nrow = length(orig_trts), ncol = 6)
	
	for (i_left_out in 1 : length(leave_outs_to_be_predicted)){
		orig_trt_i = orig_trts[i_left_out]
		est_true = ifelse(orig_trt_i == 0, yhatTx0s[i_left_out], yhatTx1s[i_left_out])
		est_counterfactual = ifelse(orig_trt_i == 0, yhatTx1s[i_left_out], yhatTx0s[i_left_out])
		if (y_higher_is_better){
			optimal = est_true > est_counterfactual
		} else {
			optimal = est_true < est_counterfactual
		}
#		cat("i_left_out", i_left_out, "orig_trt_i", orig_trt_i, "est_true", est_true, "est_counterfactual", est_counterfactual, "y0", yhatTx0s[i_left_out], "y1", yhatTx1s[i_left_out], "\n")
		res[i_left_out, ] = c(
			est_true, 
			est_counterfactual, 
			orig_trt_i, 
			ifelse(optimal, orig_trt_i, 1 - orig_trt_i),
			real_ys[i_left_out],
			Xyleftout$censored[i_left_out]
		)		
	}
#	colnames(res) = c("est_true", "est_counterfactual", "given_tx", "rec_tx", "real_y", "censored") #for debugging only
	res
}