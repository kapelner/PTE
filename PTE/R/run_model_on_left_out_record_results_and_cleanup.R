run_model_on_left_out_record_results_and_cleanup = function( 
		leave_outs_to_be_predicted, 
		train_on_all_except_these){
	
	#the left one out matrix has n-1 rows and will be considered the "training data"
	Xytrain = Xy[-train_on_all_except_these, ]
	
	#pull out the record of the left-one-out subject
	Xyleftout = Xy[leave_outs_to_be_predicted, 1 : (ncol(Xy) - 1)]
	
	if (regression_type != "survival"){
		Xytrain$censored = NULL
	}
	#build the model via the user-specified string
	mod = personalized_model_build_function(Xytrain) #this function makes use of the "Xyleft" object
	print(summary(mod))
	
	
	#also take note of what actually happened to this subject in the experiment
	real_ys = Xy[leave_outs_to_be_predicted, ncol(Xy)]
	orig_trts = Xyleftout$treatment
	
	#now evaluate the left-one-out subject on the model for both his true treatment and his counterfactual
	Xyleftout$treatment = 0
	yhatTx0s = predict_function(mod, Xyleftout)
	Xyleftout$treatment = 1
	yhatTx1s = predict_function(mod, Xyleftout)
	
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
	tabulate_results_for_left_one_out_subject(orig_trts, yhatTx0s, yhatTx1s, real_ys, Xyleftout$censored)
}