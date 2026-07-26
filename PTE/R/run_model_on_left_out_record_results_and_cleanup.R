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
		verbose,
		fast_lm_objects = NULL,
		boot_idx = NULL){

	if (!is.null(fast_lm_objects)){
		#fast path for the default continuous (OLS) model: Xy here is the ORIGINAL (unresampled)
		#master data frame and boot_idx maps fold positions back to its rows, so we can reuse the
		#design matrices precomputed once in create_fast_lm_objects() instead of paying for
		#model.frame/model.matrix/lm()/predict.lm() on every single fold of every bootstrap sample
		train_rows = boot_idx[-train_on_all_except_these]
		test_rows = boot_idx[leave_outs_to_be_predicted]

		X_train = fast_lm_objects$X_obs[train_rows, , drop = FALSE]
		fit = .lm.fit(X_train, fast_lm_objects$y_full[train_rows])

		if (fit$rank == ncol(X_train)){
			#the common case: full rank, so we can skip lm()/predict.lm() entirely
			beta = fit$coefficients
			yhatTx0s = as.numeric(fast_lm_objects$X_tx0[test_rows, , drop = FALSE] %*% beta)
			yhatTx1s = as.numeric(fast_lm_objects$X_tx1[test_rows, , drop = FALSE] %*% beta)
		} else {
			#rare rank-deficient fold (e.g. a bootstrap resample drops a factor level combination):
			#.lm.fit() does NOT reconstruct which coefficient belongs to which column the way
			#lm.fit()/lm() do (see ?.lm.fit, "the coefficients ... may be lost"), so a naive
			#X %*% coefficients here can silently diverge from predict.lm()'s answer. Fall back to
			#the exact original computation for this fold only, to guarantee identical behavior.
			Xytrain = Xy[train_rows, ]
			Xytrain$censored = NULL
			Xyleftout = Xy[test_rows, 1 : (ncol(Xy) - 1)]
			mod = personalized_model_build_function(Xytrain)
			Xyleftout$treatment = 0
			yhatTx0s = predict_function(mod, Xyleftout)
			Xyleftout$treatment = 1
			yhatTx1s = predict_function(mod, Xyleftout)
		}

		real_ys = Xy$y[test_rows]
		orig_trts = Xy$treatment[test_rows]
		censored_col = Xy$censored[test_rows]

		if (verbose && !full_verbose){
			cat(".")
		}
	} else {
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

		censored_col = Xyleftout$censored
	}

	#tabulate the result for the prediction on this left one out model (vectorized over all left-out subjects)
	est_true = ifelse(orig_trts == 0, yhatTx0s, yhatTx1s)
	est_counterfactual = ifelse(orig_trts == 0, yhatTx1s, yhatTx0s)
	optimal = if (y_higher_is_better){
		est_true > est_counterfactual
	} else {
		est_true < est_counterfactual
	}
	rec_tx = ifelse(optimal, orig_trts, 1 - orig_trts)

	matrix(
		c(est_true, est_counterfactual, orig_trts, rec_tx, real_ys, censored_col),
		ncol = 6,
		dimnames = NULL
	)
}
