THRESHOLD_FOR_BOOTSTRAP_WARNING_MESSAGE = 0.01

utils::globalVariables(c("boot_chunk", "bca_chunk"))

initialize_pte_bootstrap_backend = function(num_cores,
		os_type = .Platform$OS.type,
		make_cluster = parallel::makeCluster,
		make_mirai_cluster = NULL,
		register_parallel = doParallel::registerDoParallel,
		register_seq = foreach::registerDoSEQ,
		require_namespace = requireNamespace){
	if (num_cores == 1){
		register_seq()
		return(list(cluster = NULL, backend = "sequential"))
	}
	if (os_type == "windows"){
		if (!require_namespace("mirai", quietly = TRUE)){
			stop("The \"mirai\" package is required to parallelize on Windows. Please install it via install.packages(\"mirai\").")
		}
		if (is.null(make_mirai_cluster)){
			make_mirai_cluster = getExportedValue("mirai", "make_cluster")
		}
		cluster = make_mirai_cluster(num_cores)
		register_parallel(cluster)
		return(list(cluster = cluster, backend = "mirai"))
	}
	cluster = make_cluster(num_cores, type = "FORK")
	register_parallel(cluster)
	list(cluster = cluster, backend = "fork")
}

#' Bootstrap inference for a prespecified personalization / recommendation model
#' 
#' Runs B bootstrap samples using a prespecified model then computes the two I estimates based on cross validation. 
#' p values of the two I estimates are computed for a given \eqn{H_0: \mu_{I_0} = \mu_0}{H_0: mu_I_0 = mu_0} and 
#' confidence intervals are provided using the basic, percentile methods by default and the BCa method as well 
#' if desired.
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
#' 									\preformatted{
#' 									personalized_model_build_function = switch(regression_type,
#' 										continuous = function(Xytrain) \{ #defalt is OLS regression
#' 											lm(y ~ . * treatment,
#' 												data = Xytrain)
#' 										\},
#' 										incidence = function(Xytrain) \{ #default is logistic regression
#' 											glm(y ~ . * treatment,
#' 												data = Xytrain,
#' 												family = "binomial")
#' 										\},
#' 										survival = function(Xytrain) \{ #default is Weibull regression
#' 											survreg(Surv(Xytrain$y, Xytrain$censored) ~ (. - censored) * treatment,
#' 												data = Xytrain,
#' 												dist = "weibull")
#' 										\}
#' 									)
#' 									}
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
#' 									\preformatted{
#' 									function(mod, Xyleftout) \{
#' 										predict(mod, Xyleftout)
#' 									\}
#' 									}
#'
#' @param difference_function		A function which takes the result of one out of sample experiment (boostrap or not) of all n samples and converts it into a difference that
#' 									will be used as a score in a score distribution to determine if the personalization model is statistically significantly able to distinguish
#' 									subjects. The function looks as follows:
#' 									
#' 									\preformatted{
#' 									function(results, indices_1_1, indices_0_0, indices_0_1, indices_1_0) \{
#' 										...
#' 										c(rec_vs_non_rec_diff_score, rec_vs_all_score, rec_vs_best_score)
#' 									\}
#' 									}
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
#' @param m_pow_of_n 				Within each bootstrap sample, the exponent \eqn{\gamma}{gamma} used to set the resample size (sampled with replacement)
#' 									to \eqn{m = \lfloor n^{\gamma} \rfloor}{m = round(n^gamma)}, implementing the "m-out-of-n bootstrap." Politis, Romano
#' 									& Wolf (\emph{Subsampling}, 1999) show that bootstrap consistency for non-smooth (non-regular) functionals -- as arises
#' 									here because the treatment recommendation is a hard threshold on the estimated treatment effect -- requires
#' 									\eqn{m \to \infty}{m -> infinity} and \eqn{m / n \to 0}{m / n -> 0}, i.e. \eqn{\gamma < 1}{gamma < 1}; \eqn{\gamma = 1}{gamma = 1}
#' 									recovers the classical (and here potentially inconsistent) \eqn{n}-out-of-\eqn{n} bootstrap. The default
#' 									\eqn{\gamma = 0.75}{gamma = 0.75} follows the Politis-Romano-Wolf recommended rate as a practical compromise between
#' 									bias (\eqn{\gamma}{gamma} too close to 1, reintroducing the non-regularity) and excess bootstrap Monte Carlo variance
#' 									(\eqn{\gamma}{gamma} too small); see their book for the data-driven "minimum volatility" method for calibrating
#' 									\eqn{\gamma}{gamma} to a specific dataset. When \eqn{\gamma < 1}{gamma < 1}, the reported confidence intervals are
#' 									rescaled by \eqn{\sqrt{m / n}}{sqrt(m / n)} as required by m-out-of-n bootstrap theory; this rescaling is a no-op when
#' 									\eqn{\gamma = 1}{gamma = 1} (i.e. \eqn{m = n}{m = n}).
#' @param alpha 					Defines the confidence interval size (1 - alpha). Defaults to 0.05.
#' @param run_bca_bootstrap			Do the BCA bootstrap as well. This takes double the time. It defaults to \code{FALSE}.
#' @param display_adversarial_score	The adversarial score records the personalization metric versus the deliberate opposite of the personalization. This does not correspond
#' 									to any practical situation but it is useful for debugging. Default is \code{FALSE}.
#' @param num_cores					The number of cores to use in parallel to run the bootstrap samples more rapidly.
#' 									Defaults to \code{NULL} which uses \code{max(cores - 1, 1)} where \code{cores} is detected
#' 									once when the PTE package is loaded and cached in \code{options(mc.cores = ...)} (see
#' 									\code{?options}). Pass an explicit value here to override the cached default for a single call.
#' 
#' @return 							A results object of type "PTE_bootstrap_results" that contains much information about the observed results
#' 									and the bootstrap runs, including hypothesis testing and confidence intervals.
#' 
#' @author Adam Kapelner
#' 
#' @examples
#' library(PTE)
#'
#' ##B and num_cores are kept small here so the examples run quickly; in practice make B as
#' ##large as you can tolerate and omit num_cores to use all-but-one core.
#'
#' ##response: continuous
#' data(continuous_example)
#' X = continuous_example$X
#' y = continuous_example$y
#' pte_results = PTE_bootstrap_inference(X, y, regression_type = "continuous", B = 5, num_cores = 1)
#' pte_results
#'
#' ##response: incidence. incidence_metric also accepts "risk_ratio" and "probability_difference"
#' ##(see the incidence_metric parameter above); odds_ratio is illustrated here for speed.
#' ##force incidence and pretend y came to you this way
#' y_incidence = ifelse(y > quantile(y, 0.75), 1, 0)
#' pte_results = PTE_bootstrap_inference(X, y_incidence,
#'     regression_type = "incidence",
#'     B = 5, num_cores = 1)
#' pte_results
#'
#' ##response: survival
#' data(survival_example)
#' pte_results = PTE_bootstrap_inference(survival_example$X, survival_example$y,
#'     censored = survival_example$censored,
#'     regression_type = "survival",
#'     B = 5, num_cores = 1)
#' pte_results
#' @export
PTE_bootstrap_inference = function(X, y,  
		regression_type = "continuous",
		incidence_metric = "odds_ratio",
		personalized_model_build_function = NULL,
		censored = NULL,
		predict_function = NULL,
		difference_function = NULL,
		cleanup_mod_function = NULL,
		y_higher_is_better = TRUE,		
		verbose = FALSE,
		full_verbose = FALSE,
		H_0_mu_equals = NULL,
		pct_leave_out = 0.10,
		m_pow_of_n = 0.75,
		B = 3000,
		alpha = 0.05,
		run_bca_bootstrap = FALSE,
		display_adversarial_score = FALSE,
        num_cores = NULL
	){

	#check validity of all values that user input
	checkmate::assertDataFrame(X)
	n = nrow(X)
	checkmate::assertNumeric(y, len = n, .var.name = "y (response vector)")
	checkmate::assertChoice(regression_type, c("continuous", "incidence", "survival"))
	if (regression_type == "survival" && is.null(censored)){
		stop("If you are doing a survival comparison, you must pass in a binary \"censored\" vector.")
	}
#	if (regression_type == "survival" && !y_higher_is_better){
#		warning("You have a survival regression where y_higher_is_better is set to FALSE indicating lower survival times are better. Is this in error?")
#	}
	#incidence_metric is only meaningful (and thus only validated) when it will actually be used --
	#it's documented as ignored once a custom difference_function is supplied.
	if (regression_type == "incidence" && is.null(difference_function)){
		checkmate::assertChoice(incidence_metric, c("probability_difference", "risk_ratio", "odds_ratio"))
	}
	checkmate::assertNumber(H_0_mu_equals, null.ok = TRUE)

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
	checkmate::assertNames(colnames(X), must.include = "treatment",
		.var.name = "X (must contain a \"treatment\" column, the indicator vector of the allocation in the RCT)")
	#ensure treatment is numeric with levels zero and one
	checkmate::assertNumeric(X$treatment, .var.name = "X$treatment (must be a numeric variable)")
	checkmate::assertSetEqual(sort(unique(X$treatment)), c(0, 1),
		.var.name = "X$treatment (must be a numeric variable with only two values: 0 and 1)")

	checkmate::assertNumber(pct_leave_out, lower = 0, upper = 1)
	checkmate::assertNumber(m_pow_of_n, lower = 0, upper = 1)
	#m-out-of-n bootstrap resample size (Politis, Romano & Wolf, "Subsampling", 1999): m_pow_of_n < 1
	#ensures m = round(n^m_pow_of_n) < n which restores bootstrap consistency for the non-smooth
	#treatment-recommendation functional; m_pow_of_n = 1 recovers the classical m = n bootstrap.
	m_boot = round(n^m_pow_of_n)

	checkmate::assertIntegerish(censored, len = n, null.ok = TRUE,
		.var.name = "censored (binary censoring indicator vector)")
	if (is.null(censored)){
		censored = rep(NA, n)
	}

	checkmate::assertFunction(personalized_model_build_function, null.ok = TRUE)
	checkmate::assertFunction(predict_function, null.ok = TRUE)
	checkmate::assertFunction(difference_function, null.ok = TRUE)
	checkmate::assertFunction(cleanup_mod_function, null.ok = TRUE)
	checkmate::assertFlag(y_higher_is_better)
	checkmate::assertFlag(verbose)
	checkmate::assertFlag(full_verbose)
	checkmate::assertCount(B, positive = TRUE)
	checkmate::assertNumber(alpha, lower = 0, upper = 1)
	checkmate::assertFlag(run_bca_bootstrap)
	checkmate::assertFlag(display_adversarial_score)
	checkmate::assertCount(num_cores, positive = TRUE, null.ok = TRUE)
		
	#create default for model building function - always first order model with interactions
	if (is.null(personalized_model_build_function)){
		personalized_model_build_function_default = TRUE
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
	} else {
		personalized_model_build_function_default = FALSE
	}

	predict_function_default = is.null(predict_function)
	if (predict_function_default){
		predict_function = function(mod, Xyleftout){predict(mod, Xyleftout)}
	}

	#create master dataframe for convenience
	Xy = cbind(X, censored, y)

	#the default continuous (OLS), incidence (logistic) and survival (Weibull AFT) model/predict
	#pairs can each bypass their R-level fit/predict function's per-fold model.frame/model.matrix/
	#factor-level bookkeeping by precomputing the design matrices once up front and fitting via a
	#compiled Rcpp routine + matrix multiplication per fold instead (see fast_lm_default.R,
	#fast_glm_default.R, fast_weibull_default.R). This only applies to the *default* functions since
	#custom user functions are opaque to us, and only when there's no missing data in the columns
	#actually used by the default model since NA-dropping would otherwise need to happen per-fold.
	#For continuous/incidence, censored is an ignored placeholder column and may be all NA.
	fast_path_has_no_missing = !anyNA(X) && !anyNA(y) && (regression_type != "survival" || !anyNA(censored))
	fast_path_eligible = personalized_model_build_function_default && predict_function_default && fast_path_has_no_missing
	use_fast_lm = fast_path_eligible && regression_type == "continuous"
	use_fast_glm = fast_path_eligible && regression_type == "incidence"
	use_fast_weibull = fast_path_eligible && regression_type == "survival"
	use_fast_path = use_fast_lm || use_fast_glm || use_fast_weibull
	fast_lm_objects = if (use_fast_lm) create_fast_lm_objects(Xy) else NULL
	fast_glm_objects = if (use_fast_glm) create_fast_glm_objects(Xy) else NULL
	fast_weibull_objects = if (use_fast_weibull) create_fast_weibull_objects(Xy) else NULL

	#take care of cutoffs for leave out windows
	cutoff_obj = create_cutoffs_for_K_fold_cv(pct_leave_out, n)
	#separate cutoffs sized to the m-out-of-n bootstrap resample: each bootstrap replicate's index
	#vector boot_idx_b has length m_boot (not n) when m_pow_of_n < 1, and the fold-building code
	#below treats cutoff windows as POSITIONS WITHIN boot_idx_b, not absolute subject ids -- so the
	#windows must partition 1:m_boot, not 1:n, or folds silently run off the end of boot_idx_b.
	cutoff_obj_boot = if (m_boot == n) cutoff_obj else create_cutoffs_for_K_fold_cv(pct_leave_out, m_boot)

    ##run actual model to get observed score
    observed_run_results = list()
    observed_q_scores = list()

	#run oos results
	observed_boot_idx = 1 : n #no resampling for the observed (non-bootstrap) run
	observed_raw_results = run_default_fast_cv_results(
		Xy, y_higher_is_better, cutoff_obj, observed_boot_idx,
		fast_lm_objects, fast_glm_objects, fast_weibull_objects,
		verbose, full_verbose
	)
	if (is.null(observed_raw_results)){
		observed_raw_results = create_raw_results_list(n)
		for (l_test in 1 : cutoff_obj$num_windows){
			left_out_window_test = cutoff_obj$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj$end_cutoffs_for_leave_outs[l_test]
	#		print(left_out_window_test)
			observed_raw_results = assign_fold_results(observed_raw_results, left_out_window_test,
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
						verbose,
						fast_lm_objects,
						observed_boot_idx,
						fast_glm_objects,
						fast_weibull_objects))
		}
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
	
	if (is.null(num_cores)){
		num_cores = getOption("mc.cores")
	}
	backend = initialize_pte_bootstrap_backend(num_cores)
	cluster = backend$cluster
	on.exit(if (!is.null(cluster)) stopCluster(cluster), add = TRUE)

	#push the objects that are invariant across all B bootstrap iterations (the master data frame,
	#the precomputed fast-path design matrices, the model/predict/cleanup functions, etc.) to each
	#worker's global environment ONCE here, rather than letting foreach's automatic per-statement
	#export re-serialize them; .noexport below tells foreach not to bundle them again itself.
	pte_cluster_export_vars = c(
		"n", "m_boot", "use_fast_path", "Xy", "cutoff_obj_boot",
		"regression_type", "y_higher_is_better",
		"personalized_model_build_function", "predict_function", "cleanup_mod_function",
		"full_verbose", "verbose",
		"fast_lm_objects", "fast_glm_objects", "fast_weibull_objects",
		"difference_function", "incidence_metric"
	)
	if (!is.null(cluster)){
		clusterExport(cluster, pte_cluster_export_vars, envir = environment())
	}

	use_fast_bootstrap_q = use_fast_path && is.null(difference_function) && m_boot == n && !verbose && !full_verbose

	#chunk the B bootstrap iterations into one task per worker instead of one task per iteration --
	#foreach/doParallel dispatches one blocking IPC round trip per task, so with B often in the
	#thousands and num_cores in the single digits, per-task round-trip overhead (not computation)
	#can dominate wall-clock time. splitIndices() partitions 1:B into num_cores contiguous, balanced
	#chunks; each task runs its whole chunk locally via lapply() and returns a list of iter_lists,
	#which .combine = "c" concatenates back into a single flat, B-length, order-preserving list --
	#byte-identical in shape to the un-chunked version, so nothing downstream needs to change.
	boot_chunks = splitIndices(B, num_cores)
	if (use_fast_bootstrap_q){
		boot_chunk_seeds = sample.int(.Machine$integer.max, length(boot_chunks))
		boot_list = foreach(chunk_i = seq_along(boot_chunks), .packages = (.packages()), .noexport = pte_cluster_export_vars, .combine = "c") %dopar% {
			boot_chunk = boot_chunks[[chunk_i]]
			set.seed(boot_chunk_seeds[chunk_i])
			list(run_default_fast_bootstrap_q_scores(
				Xy, y_higher_is_better, cutoff_obj_boot, length(boot_chunk), incidence_metric,
				fast_lm_objects, fast_glm_objects, fast_weibull_objects
			))
		}
	} else {
		boot_list = foreach(boot_chunk = boot_chunks, .packages = (.packages()), .noexport = pte_cluster_export_vars, .combine = "c") %dopar% {
			lapply(boot_chunk, function(b){
				iter_list = list()
				iter_list$q_scores = list()

				#pull a bootstrap sample. When a fast path is active, avoid materializing the
				#resampled data frame entirely -- just keep the index vector and index into the
				#precomputed design matrices per-fold instead (see run_model_on_left_out_record_results_and_cleanup.R)
				boot_idx_b = sample(1 : n, m_boot, replace = TRUE)
				raw_results = run_default_fast_cv_results(
					Xy, y_higher_is_better, cutoff_obj_boot, boot_idx_b,
					fast_lm_objects, fast_glm_objects, fast_weibull_objects,
					verbose, full_verbose
				)
				if (is.null(raw_results)){
					raw_results = create_raw_results_list(m_boot)
					Xyb = if (use_fast_path) Xy else fast_row_subset(Xy, boot_idx_b)

					for (l_test in 1 : cutoff_obj_boot$num_windows){
						left_out_window_test = cutoff_obj_boot$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj_boot$end_cutoffs_for_leave_outs[l_test]
						raw_results = assign_fold_results(raw_results, left_out_window_test,
							run_model_on_left_out_record_results_and_cleanup(
								Xyb,
								regression_type,
								y_higher_is_better,
								left_out_window_test,
								left_out_window_test,
								personalized_model_build_function,
								predict_function,
								cleanup_mod_function,
								full_verbose,
								verbose,
								fast_lm_objects,
								boot_idx_b,
								fast_glm_objects,
								fast_weibull_objects))
					}
				}
				#iter_list$raw_results = raw_results
				iter_list$run_results = create_PTE_q_score_object(raw_results, regression_type, y_higher_is_better, difference_function, incidence_metric)
				iter_list$q_scores$adversarial = ifelse(length(iter_list$run_results$q_adversarial) == 1, iter_list$run_results$q_adversarial, NA)
				iter_list$q_scores$average = ifelse(length(iter_list$run_results$q_average) == 1, iter_list$run_results$q_average, NA)
				iter_list$q_scores$best = ifelse(length(iter_list$run_results$q_best) == 1, iter_list$run_results$q_best, NA)
				if (verbose){
					cat(".")
				}
				iter_list ##doParallel makes a list of these iter_lists by returning this object to the function
			})
		}
	}
	if (verbose){
		cat("\n")
	}	
  
##now populate existing vecs to proceed
  num_bad = 0
  if (use_fast_bootstrap_q){
	offset = 0
	for (chunk_i in seq_along(boot_list)){
		chunk_res = boot_list[[chunk_i]]
		chunk_n = length(chunk_res$q_scores$adversarial)
		chunk_range = (offset + 1) : (offset + chunk_n)
		q_scores$adversarial[chunk_range] = chunk_res$q_scores$adversarial
		q_scores$average[chunk_range] = chunk_res$q_scores$average
		q_scores$best[chunk_range] = chunk_res$q_scores$best
		num_bad = num_bad + chunk_res$num_bad
		offset = offset + chunk_n
	}
  } else {
	  store_run_results = !is.null(difference_function)
	  for (b in 1 : B){
	    if (store_run_results){
			run_results[[b]] = boot_list[[b]]$run_results
	    }
	    q_scores$adversarial[b] = boot_list[[b]]$q_scores$adversarial
	    q_scores$average[b] = boot_list[[b]]$q_scores$average
	    q_scores$best[b] = boot_list[[b]]$q_scores$best
		num_bad = num_bad + ifelse(boot_list[[b]]$run_results$is_bad, 1, 0)
	  }
  }
  
  

	#now see what happened with the hypothesis tests
	if (y_higher_is_better){
		#we're looking for a q score GREATER than H_0_mu_equals so the p value will be the proportion below.
		#na.rm = TRUE on both sum() here and mean() below: a "bad" bootstrap replicate (rank-deficient
		#fold, non-convergent fit, or a (quasi-)separated fit rejected by the fast-path guards) leaves
		#an NA in q_scores, and NA propagates through sum()/mean() by default -- without na.rm, a
		#single bad replicate out of B would silently turn the point estimate and p-value NA even
		#though the CI computation below (which already uses na.rm = TRUE) still succeeds on the
		#remaining valid replicates. The (B + 1) denominator is intentionally NOT reduced to the valid
		#count: bad replicates still count against the denominator, matching how num_bad's "may be
		#suspect" warning already treats them as diluting the evidence rather than being ignored.
		p_val_adversarial = sum(q_scores$adversarial < H_0_mu_equals, na.rm = TRUE) / (B + 1)
		p_val_average = sum(q_scores$average < H_0_mu_equals, na.rm = TRUE) / (B + 1)
		p_val_best = sum(q_scores$best < H_0_mu_equals, na.rm = TRUE) / (B + 1)
	} else {
		#we're looking for a q score LESS than H_0_mu_equals so the pval will be the proportion above
		p_val_adversarial = sum(q_scores$adversarial > H_0_mu_equals, na.rm = TRUE) / (B + 1)
		p_val_average = sum(q_scores$average > H_0_mu_equals, na.rm = TRUE) / (B + 1)
		p_val_best = sum(q_scores$best > H_0_mu_equals, na.rm = TRUE) / (B + 1)
	}

	est_q_adversarial = mean(q_scores$adversarial, na.rm = TRUE)
	est_q_average = mean(q_scores$average, na.rm = TRUE)
	est_q_best = mean(q_scores$best, na.rm = TRUE)
  
	##rate correction for the m-out-of-n bootstrap (Politis, Romano & Wolf, "Subsampling", 1999 /
	##Bickel-Gotze-van Zwet): the replicates in q_scores were computed on resamples of size
	##m = round(n^m_pow_of_n) rather than n, so sqrt(m) * (q_scores - observed) approximates the sampling
	##distribution of sqrt(n) * (observed - truth), not sqrt(n) * (q_scores - observed) as the classical
	##(m = n) bootstrap assumes. Both CI methods below must therefore rescale the bootstrap deviations by
	##sqrt(m / n) to remain consistent when m_pow_of_n < 1. This is a no-op when m_boot = n (m_pow_of_n = 1),
	##so it exactly reproduces prior behavior in that case.
	rate_adj = sqrt(m_boot / n)

    ##percentile method
	alpha_over_two_quantile_adversarial = quantile(q_scores$adversarial, alpha / 2, na.rm = TRUE)
	one_minus_alpha_over_two_quantile_adversarial = quantile(q_scores$adversarial, 1 - alpha / 2, na.rm = TRUE)
	ci_q_adversarial = c(
		observed_q_scores$adversarial + rate_adj * (alpha_over_two_quantile_adversarial - observed_q_scores$adversarial),
		observed_q_scores$adversarial + rate_adj * (one_minus_alpha_over_two_quantile_adversarial - observed_q_scores$adversarial)
	)

	alpha_over_two_quantile_average = quantile(q_scores$average, alpha / 2, na.rm = TRUE)
	one_minus_alpha_over_two_quantile_average = quantile(q_scores$average, 1 - alpha / 2, na.rm = TRUE)
	ci_q_average = c(
		observed_q_scores$average + rate_adj * (alpha_over_two_quantile_average - observed_q_scores$average),
		observed_q_scores$average + rate_adj * (one_minus_alpha_over_two_quantile_average - observed_q_scores$average)
	)

	alpha_over_two_quantile_best = quantile(q_scores$best, alpha / 2, na.rm = TRUE)
	one_minus_alpha_over_two_quantile_best = quantile(q_scores$best, 1 - alpha / 2, na.rm = TRUE)
	ci_q_best = c(
		observed_q_scores$best + rate_adj * (alpha_over_two_quantile_best - observed_q_scores$best),
		observed_q_scores$best + rate_adj * (one_minus_alpha_over_two_quantile_best - observed_q_scores$best)
	)

	##basic method
	basic_ci_q_adversarial = c(
		observed_q_scores$adversarial - rate_adj * (one_minus_alpha_over_two_quantile_adversarial - observed_q_scores$adversarial),
		observed_q_scores$adversarial - rate_adj * (alpha_over_two_quantile_adversarial - observed_q_scores$adversarial)
	)
	basic_ci_q_average = c(
		observed_q_scores$average - rate_adj * (one_minus_alpha_over_two_quantile_average - observed_q_scores$average),
		observed_q_scores$average - rate_adj * (alpha_over_two_quantile_average - observed_q_scores$average)
	)
	basic_ci_q_best = c(
		observed_q_scores$best - rate_adj * (one_minus_alpha_over_two_quantile_best - observed_q_scores$best),
		observed_q_scores$best - rate_adj * (alpha_over_two_quantile_best - observed_q_scores$best)
	)
	
	
	if (run_bca_bootstrap){

		#the BCa acceleration/bias-correction (z0, a) below is Efron's classical construction, derived
		#assuming q_scores comes from an n-out-of-n bootstrap -- it has no established m-out-of-n
		#generalization analogous to the sqrt(m / n) rate correction applied to the percentile/basic
		#CIs above. At m_pow_of_n < 1, q_scores is drawn from m-sized (not n-sized) resamples, so the
		#BCa interval below is not on the same theoretical footing as the percentile/basic ones.
		if (m_boot != n){
			warning("run_bca_bootstrap = TRUE was requested with m_pow_of_n = ", m_pow_of_n, " (m = ", m_boot, " < n = ", n, "). ",
				"The BCa interval's bias-correction/acceleration is only theoretically justified for the classical n-out-of-n bootstrap ",
				"(m_pow_of_n = 1); it has no established m-out-of-n analogue here and is not rate-corrected like the percentile/basic ",
				"intervals are. Treat the reported BCa interval with caution, or rerun with m_pow_of_n = 1 for a theoretically grounded BCa CI.")
		}

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
		#cutoff_obj was just reassigned for the n - 1 jackknife folds -- refresh the worker copy
		#(everything else in pte_cluster_export_vars is unchanged since the export above)
		if (!is.null(cluster)){
			clusterExport(cluster, "cutoff_obj", envir = environment())
		}

		###leave out data point out and run procedure -- chunked one task per worker, see the
		#bootstrap loop above for why (splitIndices() + lapply() + .combine = "c")
		bca_chunks = splitIndices(n, num_cores)
		full_list = foreach(bca_chunk = bca_chunks, .packages = (.packages()), .noexport = pte_cluster_export_vars, .combine = "c") %dopar%{
			lapply(bca_chunk, function(i){
				iter_list = list()
				iter_list$bca_q_scores = list()

				bca_boot_idx = (1 : n)[-i] #maps positions in the "leave subject i out entirely" set back to Xy's rows
				bca_raw_results = run_default_fast_cv_results(
					Xy, y_higher_is_better, cutoff_obj, bca_boot_idx,
					fast_lm_objects, fast_glm_objects, fast_weibull_objects,
					verbose, full_verbose
				)
				if (is.null(bca_raw_results)){
					bca_raw_results = create_raw_results_list(n - 1)
					Xy_minus_i = if (use_fast_path) Xy else fast_row_subset(Xy, -i)

					for (l_test in 1 : cutoff_obj$num_windows){
						left_out_window_test = cutoff_obj$begin_cutoffs_for_leave_outs[l_test] : cutoff_obj$end_cutoffs_for_leave_outs[l_test]
						bca_raw_results = assign_fold_results(bca_raw_results, left_out_window_test,
							run_model_on_left_out_record_results_and_cleanup(
								Xy_minus_i,
								regression_type,
								y_higher_is_better,
								left_out_window_test,
								left_out_window_test,
								personalized_model_build_function,
								predict_function,
								cleanup_mod_function,
								full_verbose,
								verbose,
								fast_lm_objects,
								bca_boot_idx,
								fast_glm_objects,
								fast_weibull_objects))
					}
				}
				iter_list$bca_run_results = create_PTE_q_score_object(bca_raw_results, regression_type, y_higher_is_better, difference_function, incidence_metric)
				iter_list$bca_q_scores$adversarial = iter_list$bca_run_results$q_adversarial
				iter_list$bca_q_scores$average = iter_list$bca_run_results$q_average
				iter_list$bca_q_scores$best = iter_list$bca_run_results$q_best

				if (verbose){
					cat(".")
				}
				iter_list
			})
		}
		
		
		
		#fill in vecs
		store_bca_run_results = !is.null(difference_function)
		for (i in 1 : n){
			if (store_bca_run_results){
				bca_run_results[[i]] = full_list[[i]]$bca_run_results ##do we use this?
			}
			bca_q_scores$adversarial[i] = full_list[[i]]$bca_q_scores$adversarial
			bca_q_scores$average[i] = full_list[[i]]$bca_q_scores$average
			bca_q_scores$best[i] = full_list[[i]]$bca_q_scores$best
		}
		
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
	return_obj$basic_ci_q_adversarial = basic_ci_q_adversarial
	return_obj$basic_ci_q_average = basic_ci_q_average
	return_obj$basic_ci_q_best = basic_ci_q_best
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
