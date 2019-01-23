THRESHOLD_FOR_BOOTSTRAP_WARNING_MESSAGE = 0.01

#' Bootstrap inference for a prespecified personalization / recommendation model
#' 
#' Runs B bootstrap samples using a prespecified model then computes the two I estimates based on cross validation. 
#' p values of the two I estimates are computed for a given \eqn{H_0: \mu_{I_0} = \mu_0}{H_0: mu_I_0 = mu_0} and 
#' confidence intervals are provided.
#' 
#' @param X 						A \eqn{n \times p}{n x p} dataframe of covariates where one column is labeled "treatment" and it
#' 									is a binary vector of treatment allocations in the study.
#' @param y 						An \eqn{n}-length numeric vector which is the response
#' @param personalized_model_build_function 	An R function that will be evaluated to construct the personalized medicine / recommendation 
#' 												model. In the formula for the model, the response is "y", the treatment vector is 
#' 												"treatment" and the data is "Xytrain". This function must return some type of object
#' 												that can be used for prediction later via \code{predict_function}. Here are the defaults
#' 												for each \code{regression_type}. They are linear models with first order interactions:
#' 
#' 											 		personalized_model_build_function = switch(regression_type,
#' 														continuous = function(Xytrain){ #defalt is OLS regression
#' 															lm(y ~ . * treatment, 
#' 																data = Xytrain)
#' 														},
#' 														incidence = function(Xytrain){ #default is logistic regression
#' 															glm(y ~ . * treatment, 
#' 																data = Xytrain, 
#' 																family = "binomial")
#' 														},
#' 														survival = function(Xytrain){ #default is Weibull regression
#' 															survreg(Surv(Xytrain$y, Xytrain$censored) ~ (. - censored) * treatment, 
#' 																data = Xytrain, 
#' 																dist = "weibull")
#' 														}
#' 													)
#' 
#' @param regression_type			A string indicating the regression problem. Legal values are "continous" (the response \code{y} is
#' 									a real number with no missing data, the default), "incidence" (the reponse \code{y} is
#' 									either 0 or 1) and "survival". If the type is "survival", the user must also supply additional data via the 
#' 									parameter \code{censored}.
#' @param incidence_metric			Ignored unless the \code{regression_type} is "incidence" and \code{difference_function} is set to \code{NULL} (in 
#' 									the latter case, you have specified a more custom metric). Then, this parameter allows the user to select which 
#' 									of the three standard metrics to use for comparison: "probability_difference", "risk_ratio", "odds_ratio" where 
#' 									the default is "odds_ratio". 
#' @param censored					Only required if the \code{regression_type} is "survival". In this case, this vector is of length \eqn{n} and is binary 
#' 									where 0 indicates censored and 1 indicates uncensored. In a clinical trial, someone who is still alive 
#' 									at the end of the study or was lost to follow up will receive a censor value of 0, while someone who died during the study 
#'                  				will receive a censor value of 1. 
#' 									\eqn{n} and is binary where 0 indicates censorship (e.g. the patient died).   
#' @param predict_function 			An R function that will be evaluated on left out data after the model is built with the training data. This function
#' 									uses the object "mod" that is the result of the \code{personalized_model_build_function} and it must make use of
#' 									"Xyleftout", a subset of observations from \code{X}. This function must return a 
#' 									scalar numeric quantity for comparison. The default function is \code{predict(mod, obs_left_out)} e.g. the default looks like:
#' 
#' 										function(mod, Xyleftout){
#' 											predict(mod, Xyleftout)
#' 										}
#' 
#' @param difference_function		A function which takes the result of one out of sample experiment (boostrap or not) of all n samples and converts it into a difference that
#' 									will be used as a score in a score distribution to determine if the personalization model is statistically significantly able to distinguish
#' 									subjects. The function looks as follows:
#' 									
#' 										function(results, indices_1_1, indices_0_0, indices_0_1, indices_1_0){
#' 											...
#' 											c(rec_vs_non_rec_diff_score, rec_vs_all_score, rec_vs_best_score)
#' 										} 
#' 
#' 
#' 									where \code{results} is a matrix consisting of columns of the estimated response of the treatment administered,
#' 									the estimated response of the counterfactual treatment, the administered treatment, the recommended treatment based on the personalization
#' 									model, the real response, and if this subject was censored (0 if so). Here are a couple of example entries:
#' 
#' 									      	est_true 	est_counterfactual 	given_tx 	rec_tx 	real_y 	censored
#' 											166.8       152.2    			1       	1    	324     1
#' 											1679.1     	2072.0    			1       	0    	160     0
#' 
#' 									
#' 									The arguments \code{indices_1_1, indices_0_0, indices_0_1, indices_1_0} give the indices of the subjects whose treatment was administered 
#' 									as 1 and whose optimal was 1, whose treatment was administered 0 and whose optimal was 0, etc.
#' 
#' 									This function should return three numeric scores: the recommend vs. the non-recommended (adversarial), the recommended 
#' 									vs. all (all) and the recommended vs. the best average treatment (best) as a 3-dimensional vector as illustrated above.
#' 
#' 									By default, this parameter is \code{NULL} which means for continuous and incidence the average difference is used and
#' 									for survival, the median Kaplan-Meier survival is used.
#' 
#' @param cleanup_mod_function 		A function that is called at the end of a cross validation iteration to cleanup the model 
#' 									in some way. This is used for instance if you would like to release the memory your model is using but generally does not apply.
#' 									The default is \code{NA} for "no function."
#' @param y_higher_is_better 		True if a response value being higher is clinically "better" than one that is lower (e.g. cognitive ability in a drug trial for the 
#' 									mentally ill). False if the response value being lower is clinically "better" than one that is higher (e.g. amount of weight lost 
#' 									in a weight-loss trial). Default is \code{TRUE}.
#' @param verbose 					Prints out a dot for each bootstrap sample. This only works on some platforms.
#' @param full_verbose 				Prints out full information for each cross validation model for each bootstrap sample. This only works on some platforms.
#' @param H_0_mu_equals 			The \eqn{\mu_{I_0}}{mu_I_0} value in \eqn{H_0}{H_0}. Default is \code{NULL} which specifies 0 for regression types continuous,
#' 									survival and incidence (with incidence metric "probability_difference") or 1 if the regression type is incidence and the incidence
#' 									metric is "risk_ratio" or "odds_ratio". These defaults essentially answer the question: does my allocation procedure do better 
#' 									than the business-as-usual / naive allocation procedure?
#' @param pct_leave_out 			In the cross-validation, the proportion of the original dataset left out to estimate out-of-sample metrics. The default is 0.1
#' 									which corresponds to 10-fold cross validation.
#' @param B 						The number of bootstrap samples to take. We recommend making this as high as you can tolerate given speed considerations.
#' 									The default is 3000.
#' @param m_prop 					Within each bootstrap sample, the proportion of the total number of rows of \code{X} to sample without replacement. \code{m_prop < 1} ensures
#' 									the number of rows sampled is less than \code{n} which fixes the consistency of the bootstrap estimator of a non-smooth functional. The default 
#' 									is 1 since non-smoothness may not be a common issue.
#' @param alpha 					Defines the confidence interval size (1 - alpha). Defaults to 0.05.
#' @param run_bca_bootstrap			Do the BCA bootstrap as well. This takes double the time. It defaults to \code{FALSE}.
#' @param display_adversarial_score	The adversarial score records the personalization metric versus the deliberate opposite of the personalization. This does not correspond
#' 									to any practical situation but it is useful for debugging. Default is \code{FALSE}.
#' @param num_cores					The number of cores to use in parallel to run the bootstrap samples more rapidly. 
#' 									Defaults to \code{NULL} which automatically sets it to one if there is one available processor or
#' 									if there are multiple available processors, the number of available processors save one.   
#' 
#' @return 							A results object of type "PTE_bootstrap_results" that contains much information about the observed results
#' 									and the bootstrap runs, including hypothesis testing and confidence intervals.
#' 
#' @author Adam Kapelner
#' 
#' @examples
#' \dontrun{
#' 	library(PTE)
#' 	B = 1000 #lower this for quicker demos
#' 
#' 	##response: continuous
#' 	data(continuous_example)
#' 	X = continuous_example$X
#'  y = continuous_example$y
#' 	pte_results = PTE_bootstrap_inference(X, y, regression_type = "continuous", B = B)
#' 	pte_results
#' 
#' 	##response: incidence
#' 	data(continuous_example)
#' 	X = continuous_example$X
#'  y = continuous_example$y
#' 	y = ifelse(y > quantile(y, 0.75), 1, 0) #force incidence and pretend y came to you this way
#'	#there are three ways to assess incidence effects below: 
#' 	#	odds ratio, risk ratio and probability difference 
#' 	pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", B = B)
#' 	pte_results
#' 	pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", B = B, 
#'                                       incidence_metric = "risk_ratio")
#' 	pte_results
#' 	pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", B = B, 
#' 	                                      incidence_metric = "probability_difference")
#' 	pte_results
#' 
#' 	##response: survival
#' 	data(survival_example)
#' 	X = survival_example$X
#' 	y = survival_example$y
#'  censored = survival_example$censored
#' 	pte_results = PTE_bootstrap_inference(X, y, censored = censored, 
#'     	regression_type = "survival", 
#'         B = 1000)
#' 	pte_results
#' }
#' @export
PTE_bootstrap_inference = function(X, y,  
		regression_type = "continuous",
		incidence_metric = "odds_ratio",
		personalized_model_build_function = NULL,
		censored = NULL,
		predict_function = function(mod, Xyleftout){predict(mod, Xyleftout)},
		difference_function = NULL,
		cleanup_mod_function = NULL,
		y_higher_is_better = TRUE,		
		verbose = FALSE,
		full_verbose = FALSE,
		H_0_mu_equals = NULL,
		pct_leave_out = 0.10,
		m_prop = 1,
		B = 3000,
		alpha = 0.05,
		run_bca_bootstrap = FALSE,
		display_adversarial_score = FALSE,
        num_cores = NULL
	){
	
	#check validity of all values that user input
	if (!(regression_type %in% c("continuous", "incidence", "survival"))){
		stop("The \"regression_type\" argument must be one of the following three:\n  continuous, incidence, survival.\n")
	}
	if (regression_type == "survival" && is.null(censored)){
		stop("If you are doing a survival comparison, you must pass in a binary \"censored\" vector.")
	}
#	if (regression_type == "survival" && !y_higher_is_better){
#		warning("You have a survival regression where y_higher_is_better is set to FALSE indicating lower survival times are better. Is this in error?")
#	}
	if (regression_type == "incidence" && is.null(difference_function) && !(incidence_metric %in% c("probability_difference", "risk_ratio", "odds_ratio"))){
		stop("If you are doing an incidence comparison, the \"incidence_metric\" parameter must be one of the following: \"probability_difference\", \"risk_ratio\" or \"odds_ratio\".")
	}
	
	if (is.null(H_0_mu_equals)){
		if (regression_type != "incidence"){
			H_0_mu_equals = 0
		} else if (incidence_metric == "risk_ratio" || incidence_metric == "odds_ratio"){
			H_0_mu_equals = 1
		} else {
			H_0_mu_equals = 0
		}
	}
	
	#ensure we have a treatment column in X
	if (!("treatment" %in% colnames(X))){
		stop("Your data frame must have a column \"\treatment\" which is an indicator vector of the allocation in the RCT.")
	}
	#ensure treatment is a factor variable with levels zero and one
	if (!(class(X$treatment) %in% c("numeric", "integer")) && identical(names(table(X$treatment)), c("0", "1"))){
		stop("Your data frame must have a column \"\treatment\" which is a numeric variable with only two values: \"0\" and \"1\".")
	}
	
	n = nrow(X)
	if (length(y) != n){
		stop("The response vector must have the same length as the number of observations in X.")
	}
	
	if (!is.null(censored) && length(censored) != n){
		stop("The binary \"censored\" vector must be the same length as the number of observations.")
	}
	
	if (is.null(censored)){
		censored = rep(NA, n)
	}
		
	#create default for model building function - always first order model with interactions
	if (is.null(personalized_model_build_function)){
		personalized_model_build_function = switch(regression_type,
			continuous = function(Xytrain){ #defalt is OLS regression
				lm(y ~ . * treatment, 
					data = Xytrain)
			},
			incidence = function(Xytrain){ #default is logistic regression
				glm(y ~ . * treatment, 
					data = Xytrain, 
					family = "binomial")
			},
			survival = function(Xytrain){ #default is Weibull regression
				survreg(Surv(Xytrain$y, Xytrain$censored) ~ (. - censored) * treatment, 
					data = Xytrain, 
					dist = "weibull")
			}
		)
	}

	#create master dataframe for convenience
	Xy = cbind(X, censored, y)

	#take care of cutoffs for leave out windows
	cutoff_obj = create_cutoffs_for_K_fold_cv(pct_leave_out, n)	
	
    ##run actual model to get observed score
    observed_run_results = list()
    observed_q_scores = list()

	#run oos results
	observed_raw_results = create_raw_results_matrix(n)
	for (l_test in 1 : cutoff_obj$num_windows){		
		left_out_window_test = cutoff_obj$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj$end_cutoffs_for_leave_outs[l_test]
#		print(left_out_window_test)
	  	observed_raw_results[left_out_window_test, ] = 
			run_model_on_left_out_record_results_and_cleanup(
					Xy,
					regression_type,
					y_higher_is_better,
					left_out_window_test, 
					left_out_window_test,
					personalized_model_build_function,
					predict_function,
					cleanup_mod_function,
					full_verbose,
					verbose)
	}
	observed_raw_results
	observed_run_results = create_PTE_results_object(observed_raw_results, regression_type, y_higher_is_better, difference_function, incidence_metric)
	observed_run_results
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
	if (is.null(num_cores)){
		num_cores = as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
		num_cores = max(num_cores - 1, 1)
	}
	cluster = makeCluster(num_cores)
	registerDoParallel(cluster)
  
	boot_list = foreach(b = 1 : B, .packages = (.packages())) %dopar% {
		
    	iter_list = list()
    	iter_list$q_scores = list()
		raw_results = create_raw_results_matrix(n)
		
		#pull a bootstrap sample
		Xyb = Xy[sample(1 : n, round(m_prop * n), replace = TRUE), ]
		
		for (l_test in 1 : cutoff_obj$num_windows){
			left_out_window_test = cutoff_obj$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj$end_cutoffs_for_leave_outs[l_test]
			raw_results[left_out_window_test, ] = run_model_on_left_out_record_results_and_cleanup(
					Xyb,
					regression_type,
					y_higher_is_better,
					left_out_window_test,
					left_out_window_test,
					personalized_model_build_function,
					predict_function,
					cleanup_mod_function,
					full_verbose,
					verbose)
		}
		#iter_list$raw_results = raw_results
		iter_list$run_results = create_PTE_results_object(raw_results, regression_type, y_higher_is_better, difference_function, incidence_metric)
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
		p_val_adversarial = sum(q_scores$adversarial > H_0_mu_equals) / (B + 1)
		p_val_average = sum(q_scores$average > H_0_mu_equals) / (B + 1)
		p_val_best = sum(q_scores$best > H_0_mu_equals) / (B + 1)
	}
	
	est_q_adversarial = mean(q_scores$adversarial)
	est_q_average = mean(q_scores$average)
	est_q_best = mean(q_scores$best)
  
    ##percentile method
	ci_q_adversarial = c(quantile(q_scores$adversarial, alpha / 2, na.rm = TRUE), quantile(q_scores$adversarial, 1 - alpha / 2, na.rm = TRUE))
	ci_q_average = c(quantile(q_scores$average, alpha / 2, na.rm = TRUE), quantile(q_scores$average, 1 - alpha / 2, na.rm = TRUE))
	ci_q_best = c(quantile(q_scores$best, alpha / 2, na.rm = TRUE), quantile(q_scores$best, 1 - alpha / 2, na.rm = TRUE))
  
  

	
	
	if (run_bca_bootstrap){
		
		##need to deal with the reversal of signs -- just flip sign if y_higher_is_better is FALSE
		if (!y_higher_is_better){
			observed_q_scores$adversarial =  -observed_q_scores$adversarial
			observed_q_scores$average = -observed_q_scores$average
			observed_q_scores$best = -observed_q_scores$best
			q_scores$adversarial = -q_scores$adversarial
			q_scores$average = -q_scores$average
			q_scores$best = -q_scores$best
		}
		
		
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
		full_list = foreach(i = 1 : n, .packages = (.packages())) %dopar%{
			iter_list = list() 
			iter_list$bca_q_scores = list()
			
			bca_raw_results = create_raw_results_matrix(n - 1)
			Xy_minus_i = Xy[-i, ]
			
			for (l_test in 1 : cutoff_obj$num_windows){
				left_out_window_test = cutoff_obj$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj$end_cutoffs_for_leave_outs[l_test]
				bca_raw_results[left_out_window_test, ] = run_model_on_left_out_record_results_and_cleanup(
						Xy_minus_i,
						regression_type,
						y_higher_is_better,
						left_out_window_test,
						left_out_window_test,
						personalized_model_build_function,
						predict_function,
						cleanup_mod_function,
						full_verbose,
						verbose)
			}
			iter_list$bca_run_results = create_PTE_results_object(bca_raw_results, regression_type, y_higher_is_better, difference_function, incidence_metric)
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
		
		bca_ci_q_adversarial = quantile(q_scores$adversarial, probs = bca_ci_q_adversarial_quantiles, na.rm = TRUE)
		bca_ci_q_average = quantile(q_scores$average, probs = bca_ci_q_average_quantiles, na.rm = TRUE)
		bca_ci_q_best = quantile(q_scores$best, probs = bca_ci_q_best_quantiles, na.rm = TRUE)  
		
		
		##convert back to correctly signed units if y_higher_is_better is false
		if (!y_higher_is_better){
			
			observed_q_scores$adversarial = -observed_q_scores$adversarial
			observed_q_scores$average = -observed_q_scores$average
			observed_q_scores$best = -observed_q_scores$best 
			
			q_scores$adversarial = -q_scores$adversarial
			q_scores$average = -q_scores$average
			q_scores$best = -q_scores$best
			
			
			bca_q_scores$adversarial = -bca_q_scores$adversarial
			bca_q_scores$average = -bca_q_scores$average
			bca_q_scores$best = -bca_q_scores$best
			
			bca_ci_q_adversarial = -bca_ci_q_adversarial[2:1]
			bca_ci_q_average = -bca_ci_q_average[2:1]
			bca_ci_q_best = -bca_ci_q_best[2:1]
		}
	}
	
	#print a warning message if need be
	if (num_bad  / B > THRESHOLD_FOR_BOOTSTRAP_WARNING_MESSAGE){
		warning("This inference may be suspect since ", num_bad, " bootstrap samples were invalid (", round(num_bad  / B * 100, 2), "%).", sep = "")
	}
	
	
	#return the final model too
	if (regression_type != "survival"){
		Xy$censored = NULL
	}
	personalization_model = personalized_model_build_function(Xy)
	
	return_obj = list()
	return_obj$Xy = Xy
	return_obj$regression_type = regression_type
	return_obj$incidence_metric = incidence_metric
	return_obj$personalized_model_build_function = personalized_model_build_function
	return_obj$predict_function = predict_function
	return_obj$y_higher_is_better = y_higher_is_better
	return_obj$difference_function = difference_function
	return_obj$cleanup_mod_function = cleanup_mod_function
	return_obj$alpha = alpha
	return_obj$H_0_mu_equals = H_0_mu_equals
	return_obj$personalization_model = personalization_model
	return_obj$run_bca_bootstrap = run_bca_bootstrap
	return_obj$run_results = run_results
	return_obj$num_bad = num_bad
    return_obj$observed_q_scores = observed_q_scores
	return_obj$q_scores = q_scores
	return_obj$p_val_adversarial = p_val_adversarial
	return_obj$p_val_average = p_val_average
	return_obj$p_val_best = p_val_best
	return_obj$est_q_adversarial = est_q_adversarial
	return_obj$est_q_average = est_q_average
	return_obj$est_q_best = est_q_best
	return_obj$ci_q_adversarial = ci_q_adversarial
	return_obj$ci_q_average = ci_q_average
	return_obj$ci_q_best = ci_q_best
	return_obj$display_adversarial_score = display_adversarial_score
	return_obj$B = B
	if (run_bca_bootstrap){
		return_obj$bca_q_scores = bca_q_scores
	    return_obj$bca_ci_q_adversarial = bca_ci_q_adversarial
	    return_obj$bca_ci_q_average = bca_ci_q_average
	    return_obj$bca_ci_q_best = bca_ci_q_best
	}
	class(return_obj) = "PTE_bootstrap_results"
	return_obj
}