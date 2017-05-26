THRESHOLD_FOR_BOOTSTRAP_WARNING_MESSAGE = 0.01

#' Bootstrap inference for a prespecified personalization / recommendation model
#' 
#' Runs B bootstrap samples using a prespecified model then computes the two I estimates based on cross validation. 
#' p values of the two I estimates are computed for a given \eqn{H_0: \mu_{I_0} = \mu_0}{H_0: mu_I_0 = mu_0} and 
#' confidence intervals are provided.
#' 
#' @param X 						A \eqn{n \times p}{n x p} dataframe of covariates.
#' @param y 						An \eqn{n}-length numeric vector which is the response
#' @param model_string 				A string of R code that will be evaluated to construct the leave one out model. Make sure the covariate data is
#' 									referred to as \code{Xyleft}.
#' @param regression_type			A string indicating the regression problem. Legal values are "continous" (the response \code{y} is
#' 									a real number with no missing data, the default), "incidence" (the reponse \code{y} is
#' 									either 0 or 1) and "survival" (the response is a time value with NA's for all uncensored
#' 									responses). If the type is "survival", the user must also supply additional data via the 
#' 									parameter \code{censored}.
#' @param censored					Only required if the \code{regression_type} is "survival". In this case, this vector is of length
#' 									\eqn{n} and is binary where 0 indicates censorship (e.g. the patient died).   
#' @param predict_string 			A string of R code that will be evaluated on left out data after the model is built with the training data. Make sure
#' 									the forecast data (the left one out data) is referred to as \code{obs_left_out} and the model is referred to as \code{mod}.
#' @param cleanup_mod_function 		A function that is called at the end of a cross validation iteration to cleanup the model 
#' 									in some way. This is used for instance if you would like to release the memory your model is using but generall does not apply.
#' 									The default is \code{NA} for "no function."
#' @param y_higher_is_better 		True if a response value being higher is clinically "better" than one that is lower (e.g. cognitive ability in a drug trial for the 
#' 									mentally ill). False if the response value being lower is clinically "better" than one that is higher (e.g. amount of weight lost 
#' 									in a weight-loss trial). Default is \code{TRUE}.
#' @param verbose 					Prints out a dot for each bootstrap sample. This only works on some platforms.
#' @param full_verbose 				Prints out full information for each cross validation model for each bootstrap sample. This only works on some platforms.
#' @param H_0_mu_equals 			The \eqn{\mu_{I_0}}{mu_I_0} value in \eqn{H_0}{H_0}. Default is 0 which answers the question: does my allocation procedure do better than a naive
#' 									allocation procedure.
#' @param pct_leave_out 			In the cross-validation, the proportion of the original dataset left out to estimate out-of-sample metrics. The default is 0.1
#' 									which corresponds to 10-fold cross validation.
#' @param B 						The number of bootstrap samples to take. We recommend making this as high as you can tolerate given speed considerations.
#' 									The default is 3000.
#' @param m_prop 					Within each bootstrap sample, the proportion of the total number of rows of \code{X} to sample without replacement. \code{m_prop < 1} ensures
#' 									the number of rows sampled is less than \code{n} which fixes the consistency of the bootstrap estimator of a non-smooth functional. The default 
#' 									is 1 since non-smoothness may not be a common issue.
#' @param alpha 					Defines the confidence interval size (1 - alpha). Defaults to 0.05.
#' @param plot 						Illustrates the estimate, the bootstrap samples and the confidence intervals on a histogram plot. Default to TRUE.
#' @param num_cores					The number of cores to use in parallel to run the bootstrap samples more rapidly. Defaults to serial by using 1 core.   
#' @param ... 						Additional parameters to be sent to the model constructor. Note that if you wish to pass these parameters, 
#' 									"..." must be specified in model_string.  
#' 
#' @return 
#' 
#' @author Adam Kapelner
#' @export
PTE_bootstrap_inference = function(X, y,  
		regression_type = "response",
		personalized_model_build_function = NULL,
		censored = NULL,
		predict_function = function(){predict(mod, obs_left_out);},
		cleanup_mod_function = NA,
		y_higher_is_better = TRUE,		
		verbose = TRUE,
		full_verbose = FALSE,
		H_0_mu_equals = 0,
		pct_leave_out = 0.10,
		m_prop = 1,
		B = 3000,
		alpha = 0.05,
		plot = TRUE,
        num_cores = 1, 
		...){
	
	#check validity of all values that user input
	if (!(regression_type %in% c("continuous", "incidence", "survival"))){
		stop("The \"regression_type\" argument must be one of the following three:\n  continuous, incidence, survival.\n")
	}
	if (regression_type == "survival" && is.null(censored)){
		stop("If you are doing a survival regression, you must pass in a binary \"censored\" vector.")
	}
	
	#data shared throughout all bootstrap simulations
	Xy = cbind(X, y)
	n = nrow(Xy)
		
	if (!is.null(censored) && length(censored) != n){
		stop("The binary \"censored\" vector must be the same length as the number of observations.")
	}
		
	#create default for model building function - always first order model with interactions
	if (is.null(personalized_model_build_function)){
		switch(regression_type,
				continuous = function(){
					lm(y ~ . + treatment * ., 
						data = Xyleft)
				},
				incidence = function(){
					glm(y ~ . + treatment * ., 
						data = Xyleft, 
						family = "binomial")
				},
				survival = function(){
					survreg(Surv(y, censored) ~ . + treatment, 
						data = Xyleft, 
						dist = "weibull")	
				}
		)
	}

	#take care of cutoffs for leave out windows
	cutoff_obj = create_cutoffs_for_K_fold_cv(pct_leave_out, n)	
	
    ##run actual model to get observed score
    observed_run_results = list()
    observed_q_scores = list()
  
	#run oos results
	observed_raw_results = create_raw_results_matrix(n)
	for (l_test in 1 : cutoff_obj$num_windows){		
		left_out_window_test = cutoff_obj$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj$end_cutoffs_for_leave_outs[l_test]
	  	observed_raw_results[left_out_window_test, ] = run_model_on_left_out_record_results_and_cleanup(Xy, 
				left_out_window_test,
				left_out_window_test,
		     	model_string, 
		     	predict_string, 
		    	cleanup_mod_function, 
		    	y_higher_is_better,
		     	full_verbose = full_verbose,
		     	...) #all other arguments head on into the model building which is done in this function
	}
	observed_run_results = create_PTE_results_object(observed_raw_results, y_higher_is_better)
	observed_q_scores$adversarial = observed_run_results$q_adversarial
	observed_q_scores$average = observed_run_results$q_average
	observed_q_scores$best = observed_run_results$q_best
  
    ##Now move on to bootstrap sampling...
	
	#place to store results from all bootstrap samples
	#raw_results = list()
	run_results = list()
	q_scores = list()
	q_scores[["adversarial"]] = array(NA, B)
	q_scores[["average"]] = array(NA, B)
	q_scores[["best"]] = array(NA, B)
	
    #will work on windows -- not sure about unix/mac
	cluster = makeCluster(num_cores)
	registerDoParallel(cluster)
  
	boot_list = foreach(b = 1 : B) %dopar% {
    	iter_list = list()
    	iter_list$q_scores = list()
		raw_results = create_raw_results_matrix(n)
		
		#pull a bootstrap sample
		Xyb = Xy[sample(1 : n, round(m_prop * n), replace = TRUE), ]
		
		for (l_test in 1 : cutoff_obj$num_windows){
			left_out_window_test = cutoff_obj$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj$end_cutoffs_for_leave_outs[l_test]
			raw_results[left_out_window_test, ] = run_model_on_left_out_record_results_and_cleanup(Xyb, 
					left_out_window_test,
					left_out_window_test,
					model_string, 
					predict_string, 
					cleanup_mod_function, 
					y_higher_is_better,
					full_verbose = full_verbose,
					...) #all other arguments head on into the model building which is done in this function
		}
		#iter_list$raw_results = raw_results
		iter_list$run_results = create_PTE_results_object(raw_results, y_higher_is_better)
		iter_list$q_scores$adversarial = ifelse(length(iter_list$run_results$q_adversarial) == 1, iter_list$run_results$q_adversarial, NA)
		iter_list$q_scores$average = ifelse(length(iter_list$run_results$q_average) == 1, iter_list$run_results$q_average, NA)
		iter_list$q_scores$best = ifelse(length(iter_list$run_results$q_best) == 1, iter_list$run_results$q_best, NA)
		if (verbose){
			cat(".")
		}
    	iter_list ##doParallel makes a list of these iter_lists by returning this object to the function
	}
	if (verbose){
		cat("\n")
	}	
  
##now populate existing vecs to proceed
  num_bad = 0
  for (b in 1 : B){
    run_results[[b]] = boot_list[[b]]$run_results
	#raw_results[[b]] = boot_list[[b]]$raw_results
    q_scores$adversarial[b] = boot_list[[b]]$q_scores$adversarial
    q_scores$average[b] = boot_list[[b]]$q_scores$average
    q_scores$best[b] = boot_list[[b]]$q_scores$best
	num_bad = num_bad + ifelse(boot_list[[b]]$run_results$is_bad, 1, 0)
  }
  
  

	#now see what happened with the hypothesis tests
	if (y_higher_is_better){
		#we're looking for a q score GREATER than H_0_mu_equals so the p value will be the proportion below
		p_val_adversarial = sum(q_scores$adversarial < H_0_mu_equals) / (B + 1)
		p_val_average = sum(q_scores$average < H_0_mu_equals) / (B + 1)
		p_val_best = sum(q_scores$best < H_0_mu_equals) / (B + 1)
	} else {
		#we're looking for a q score LESS than H_0_mu_equals so the pval will be the proportion above
		p_val_adversarial = sum(q_scores[["adversarial"]] > H_0_mu_equals) / (B + 1)
		p_val_average = sum(q_scores[["average"]] > H_0_mu_equals) / (B + 1)
		p_val_best = sum(q_scores[["best"]] > H_0_mu_equals) / (B + 1)
	}
	
	est_q_adversarial = mean(q_scores$adversarial)
	est_q_average = mean(q_scores$average)
	est_q_best = mean(q_scores$best)
  
    ##percentile method
	ci_q_adversarial = c(quantile(q_scores$adversarial, alpha / 2), quantile(q_scores$adversarial, 1 - alpha / 2))
	ci_q_average = c(quantile(q_scores$average, alpha / 2), quantile(q_scores$average, 1 - alpha / 2))
	ci_q_best = c(quantile(q_scores$best, alpha / 2), quantile(q_scores$best, 1 - alpha / 2))
  
  
    ##BCA CIS
    ##compute acceleration
	bca_run_results = list()
	bca_q_scores = list()
	bca_q_scores[["adversarial"]] = array(NA, n)
	bca_q_scores[["average"]] = array(NA, n)
	bca_q_scores[["best"]] = array(NA, n)
	
	#due to the jacknife, we now need new beginning and endpoints in the sliding window
	cutoff_obj = create_cutoffs_for_K_fold_cv(pct_leave_out, n - 1)

    ###leave out data point out and run procedure
    full_list = foreach(i = 1 : n) %dopar%{
	    iter_list = list() 
	    iter_list$bca_q_scores = list()
	    
	    bca_raw_results = create_raw_results_matrix(n - 1)
	    Xy_minus_i = Xy[-i, ]
			    
		for (l_test in 1 : cutoff_obj$num_windows){
			left_out_window_test = cutoff_obj$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj$end_cutoffs_for_leave_outs[l_test]
	    	bca_raw_results[left_out_window_test, ] = run_model_on_left_out_record_results_and_cleanup(Xy_minus_i, 
					left_out_window_test,
					left_out_window_test,
			        model_string, 
			        predict_string, 
			        cleanup_mod_function, 
			        y_higher_is_better,
			        full_verbose = full_verbose,
			        ...) #all other arguments head on into the model building which is done in this function
	    }
	    iter_list$bca_run_results = create_PTE_results_object(bca_raw_results, y_higher_is_better)
	    iter_list$bca_q_scores$adversarial = iter_list$bca_run_results$q_adversarial
	    iter_list$bca_q_scores$average = iter_list$bca_run_results$q_average
	    iter_list$bca_q_scores$best = iter_list$bca_run_results$q_best
	
	    if (verbose){
	      cat(".")
	    }
	    iter_list
  	}
	


    #fill in vecs
	for (i in 1 : n){
	  bca_run_results[[i]] = full_list[[i]]$run_results ##do we use this?
	  bca_q_scores$adversarial[i] = full_list[[i]]$bca_q_scores$adversarial
	  bca_q_scores$average[i] = full_list[[i]]$bca_q_scores$average
  	  bca_q_scores$best[i] = full_list[[i]]$bca_q_scores$best
  	}
  
    stopCluster(cluster) ## remove cluster

    ##need to deal with the reversal of signs -- just flip sign if y_higher_is_better is FALSE
    if (!y_higher_is_better){
      bca_q_scores$adversarial = -bca_q_scores$adversarial
      bca_q_scores$average = -bca_q_scores$average
      bca_q_scores$best = -bca_q_scores$best
      observed_q_scores$adversarial =  -observed_q_scores$adversarial
      observed_q_scores$average = -observed_q_scores$average
      observed_q_scores$best = -observed_q_scores$best
      q_scores$adversarial = -q_scores$adversarial
      q_scores$average = -q_scores$average
      q_scores$best = -q_scores$best
    }

    ##now compute the accelerations
    diff_adversarial = mean(bca_q_scores$adversarial) - bca_q_scores$adversarial
    a_adversarial = sum(diff_adversarial^3) / (6*sum(diff_adversarial^2)^1.5)
  
	diff_average = mean(bca_q_scores$average) - bca_q_scores$average
	a_average = sum(diff_average^3) / (6*sum(diff_average^2)^1.5)

    diff_best = mean(bca_q_scores$best) - bca_q_scores$best
    a_best = sum(diff_best^3) / (6*sum(diff_best^2)^1.5)
    
    ##z0 values
    z0_adversarial =  qnorm(sum(q_scores$adversarial <= observed_q_scores$adversarial)/length(q_scores$adversarial)) ## proportion less than estimate
	z0_average =  qnorm(sum(q_scores$average <= observed_q_scores$average)/length(q_scores$average)) ## proportion less than estimate
	z0_best =  qnorm(sum(q_scores$best <= observed_q_scores$best)/length(q_scores$best)) ## proportion less than estimate

    ##Now compute bca CIs
    left_adversarial = z0_adversarial + qnorm(alpha/2)
    right_adversarial = z0_adversarial + qnorm(1 - alpha/2)
	bca_ci_q_adversarial_quantiles = c(pnorm(z0_adversarial + (left_adversarial)/(1 - a_adversarial * left_adversarial)),
	                         pnorm(z0_adversarial + (right_adversarial)/(1 - a_adversarial * right_adversarial)))

    left_average = z0_average + qnorm(alpha/2)
    right_average = z0_average + qnorm(1 - alpha/2)
    bca_ci_q_average_quantiles = c(pnorm(z0_average + (left_average)/(1 - a_average * left_average)),
                         pnorm(z0_average + (right_average)/(1 - a_average * right_average)))

    left_best = z0_best + qnorm(alpha/2)
    right_best = z0_best + qnorm(1 - alpha/2)
    bca_ci_q_best_quantiles = c(pnorm(z0_best + (left_best)/(1 - a_best * left_best)),
                         pnorm(z0_best + (right_best)/(1 - a_best * right_best)))

    bca_ci_q_adversarial = quantile(q_scores$adversarial, probs = bca_ci_q_adversarial_quantiles)
	bca_ci_q_average = quantile(q_scores$average, probs = bca_ci_q_average_quantiles)
	bca_ci_q_best = quantile(q_scores$best, probs = bca_ci_q_best_quantiles)  


    ##convert back to correctly signed units if y_higher_is_better is false
    if (!y_higher_is_better){
      bca_q_scores$adversarial = -bca_q_scores$adversarial
      bca_q_scores$average = -bca_q_scores$average
      bca_q_scores$best = -bca_q_scores$best
      observed_q_scores$adversarial = -observed_q_scores$adversarial
      observed_q_scores$average = -observed_q_scores$average
      observed_q_scores$best = -observed_q_scores$best   
      q_scores$adversarial = -q_scores$adversarial
      q_scores$average = -q_scores$average
      q_scores$best = -q_scores$best    
    
      bca_ci_q_adversarial = -bca_ci_q_adversarial[2:1]
      bca_ci_q_average = -bca_ci_q_average[2:1]
      bca_ci_q_best = -bca_ci_q_best[2:1]
    }
	
	if (plot){
		min_q = min(q_scores$average, q_scores$best)
		max_q = max(q_scores$average, q_scores$best)
		par(mfrow = c(2, 1))
		hist(q_scores$average, br = B / 3, xlab = "I", xlim = c(min_q, max_q), main = "Average I's")
		abline(v = est_q_average, col = "forestgreen", lwd = 3)
		abline(v = ci_q_average[1], col = "firebrick3", lwd = 1)
		abline(v = ci_q_average[2], col = "firebrick3", lwd = 1)
		abline(v = bca_ci_q_average[1], col = "dodgerblue3", lwd = 1)
		abline(v = bca_ci_q_average[2], col = "dodgerblue3", lwd = 1)
		hist(q_scores$best, br = B / 3, xlab = "I", xlim = c(min_q, max_q), main = "Best I's")
		abline(v = est_q_best, col = "forestgreen", lwd = 3)
		abline(v = ci_q_best[1], col = "firebrick3", lwd = 1)
		abline(v = ci_q_best[2], col = "firebrick3", lwd = 1)	
		abline(v = bca_ci_q_best[1], col = "dodgerblue3", lwd = 1)
		abline(v = bca_ci_q_best[2], col = "dodgerblue3", lwd = 1)
	}
	
	#print a warning message if need be
	if (num_bad  / B > THRESHOLD_FOR_BOOTSTRAP_WARNING_MESSAGE){
		warning("This inference may be suspect since ", num_bad, " bootstrap samples were invalid (", round(num_bad  / B * 100, 2), "%).", sep = "")
	}
	
	return_obj = list()
	return_obj$Xy = Xy
	return_obj$model_string = model_string
	return_obj$predict_string = predict_string
	return_obj$y_higher_is_better = y_higher_is_better
	#return_obj$raw_results = raw_results
	return_obj$run_results = run_results
	return_obj$num_bad = num_bad
    return_obj$observed_q_adversarial = observed_q_scores$adversarial
    return_obj$observed_q_average = observed_q_scores$average
    return_obj$observed_q_best = observed_q_scores$best
	return_obj$q_scores = q_scores
	return_obj$p_val_adversarial = p_val_adversarial
	return_obj$p_val_average = p_val_average
	return_obj$p_val_best = p_val_best
	return_obj$alpha = alpha
	return_obj$est_q_adversarial = est_q_adversarial
	return_obj$est_q_average = est_q_average
	return_obj$est_q_best = est_q_best
	return_obj$ci_q_adversarial = ci_q_adversarial
	return_obj$ci_q_average = ci_q_average
	return_obj$ci_q_best = ci_q_best
    return_obj$bca_ci_q_adversarial = bca_ci_q_adversarial
    return_obj$bca_ci_q_average = bca_ci_q_average
    return_obj$bca_ci_q_best = bca_ci_q_best
	class(return_obj) = "PTE_bootstrap_results"
	return_obj
}