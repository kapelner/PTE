## Internal helper functions used by PTE_bootstrap_inference() and its bootstrap/cross-validation
## loop. Consolidated here (rather than one file per function) since each is small and none is
## part of the public API.

create_cutoffs_for_K_fold_cv = function(pct_leave_out, n){
	if (n * pct_leave_out < 1){
		#fewer than one observation per fold at this pct_leave_out -- the seq()-based partition below
		#assumes a step size >= 1 and silently produces degenerate/overlapping windows otherwise (e.g.
		#when the m-out-of-n bootstrap's resample size m = round(n^m_pow_of_n) is small). Fall back to
		#leave-one-out cross-validation, which is well-defined for any n >= 1.
		return(list(begin_cutoffs_for_leave_outs = 1 : n, end_cutoffs_for_leave_outs = 1 : n, num_windows = n))
	}
	begin_cutoffs_for_leave_outs = seq(1, n)[seq(1, n, n * pct_leave_out)] - 1
	begin_cutoffs_for_leave_outs[1] = 1
	end_cutoffs_for_leave_outs = seq(1, n)[seq(1, n, n * pct_leave_out)] - 2
	end_cutoffs_for_leave_outs = c(end_cutoffs_for_leave_outs[2 : length(end_cutoffs_for_leave_outs)], n)
	num_windows = length(begin_cutoffs_for_leave_outs)
	list(begin_cutoffs_for_leave_outs = begin_cutoffs_for_leave_outs, end_cutoffs_for_leave_outs = end_cutoffs_for_leave_outs, num_windows = num_windows)
}

#' Preallocate the per-subject raw results accumulator
#'
#' A plain named list of six equal-length numeric vectors, not a data.frame or matrix: this gets
#' filled in fold-by-fold (potentially thousands of times per bootstrap sample) via
#' \code{assign_fold_results()}, and list-element/vector-slice assignment is a direct memory write,
#' unlike \code{[<-.data.frame} which re-checks/re-coerces column types and rebuilds the underlying
#' list on every call -- a large, avoidable cost at this scale. Downstream code
#' (\code{create_PTE_results_object()}) reads columns back out via \code{results$colname}.
#'
#' @param n Number of subjects.
#' @keywords internal
#' @noRd
create_raw_results_list = function(n){
	list(
		est_true = rep(NA_real_, n),
		est_counterfactual = rep(NA_real_, n),
		given_tx = rep(NA_real_, n),
		rec_tx = rep(NA_real_, n),
		real_y = rep(NA_real_, n),
		censored = rep(NA_real_, n)
	)
}

#' Write one cross-validation fold's results into the accumulator at the given subject indices
#'
#' @param raw_results 	An accumulator list from \code{create_raw_results_list()} (or a prior call
#' to this function).
#' @param window 		The subject indices this fold covers.
#' @param fold_result 	The named list returned by \code{run_model_on_left_out_record_results_and_cleanup()}
#' for those same subjects -- same names as \code{raw_results}.
#' @keywords internal
#' @noRd
assign_fold_results = function(raw_results, window, fold_result){
	for (nm in names(fold_result)){
		raw_results[[nm]][window] = fold_result[[nm]]
	}
	raw_results
}

#' Subset rows (and optionally columns) of a data frame without \code{[.data.frame}'s overhead
#'
#' \code{Xy[idx, ]} (and especially \code{Xy[sample(1:n, n, replace = TRUE), ]}) is expensive
#' when called tens of thousands of times across the bootstrap x cross-validation loop: it
#' re-derives factor levels, re-dispatches through \code{[[.data.frame}/\code{[.data.frame} for
#' every column, and -- worst of all when there are repeated indices from bootstrap resampling --
#' calls \code{make.unique()} to disambiguate duplicated row names via \code{paste}/\code{unique}.
#' None of that bookkeeping is needed here since nothing downstream looks up results by row name.
#' This does the same per-column extraction (still dispatching to \code{[.factor} etc. so factor
#' levels are preserved correctly) but skips the data-frame-level machinery, and is safe to use
#' regardless of which \code{personalized_model_build_function}/\code{predict_function} the caller
#' plugs in since the result is a normal, valid data frame (just with plain integer row names).
#'
#' @param df 	A data frame.
#' @param idx 	A row index vector (positive or negative, exactly as you'd pass to \code{df[idx, ]}).
#' @param cols 	A column index vector (defaults to all columns).
#' @keywords internal
#' @noRd
fast_row_subset = function(df, idx, cols = seq_along(df)){
	n_out = length(attr(df, "row.names")[idx])
	out = vector("list", length(cols))
	for (k in seq_along(cols)){
		col = .subset2(df, cols[k])
		#a data frame column can itself be an embedded matrix (e.g. from scale() or poly());
		#plain vector-style col[idx] would silently flatten/corrupt that, so subset rows only
		out[[k]] = if (is.matrix(col)) col[idx, , drop = FALSE] else col[idx]
	}
	names(out) = names(df)[cols]
	attr(out, "row.names") = seq_len(n_out)
	class(out) = "data.frame"
	out
}

#' Precompute design matrices for the default continuous (OLS) personalization model
#'
#' Builds the model matrix for \code{y ~ . * treatment} once against the FULL master data frame,
#' along with the two counterfactual design matrices (treatment forced to 0 and to 1 for every row).
#' This lets \code{run_model_on_left_out_record_results_and_cleanup} fit each fold via \code{.lm.fit}
#' and predict via a plain matrix multiply, instead of calling \code{lm()}/\code{predict.lm()} (which
#' each re-run \code{model.frame}/\code{model.matrix} formula and factor-level bookkeeping) tens of
#' thousands of times across the bootstrap x cross-validation loop.
#'
#' Only used when the user has NOT supplied a custom \code{personalized_model_build_function} or
#' \code{predict_function} and \code{regression_type == "continuous"}. Assumes no missing data
#' (callers must check \code{anyNA(Xy)} beforehand), since row-deletion for NAs would otherwise
#' need to happen per-fold, defeating the point of precomputing the matrices once.
#'
#' @param Xy A data frame with columns "treatment", the covariates and "y" (and possibly "censored",
#' which is dropped since it is never a predictor for the continuous regression type).
#' @keywords internal
#' @noRd
create_fast_lm_objects = function(Xy){
	Xy$censored = NULL #never a predictor for the continuous regression type

	mf_full = model.frame(y ~ . * treatment, data = Xy)
	lm_terms = attr(mf_full, "terms")
	X_obs = model.matrix(lm_terms, mf_full)
	y_full = model.response(mf_full)

	Xy_tx0 = Xy
	Xy_tx0$treatment = 0
	Xy_tx1 = Xy
	Xy_tx1$treatment = 1

	list(
		X_obs = X_obs,
		y_full = y_full,
		X_tx0 = model.matrix(lm_terms, data = Xy_tx0),
		X_tx1 = model.matrix(lm_terms, data = Xy_tx1)
	)
}

#' Precompute design matrices for the default incidence (logistic) personalization model
#'
#' Mirrors \code{create_fast_lm_objects()} but for the default incidence model: builds the
#' model matrix for \code{y ~ . * treatment} once against the FULL master data frame, along with
#' the two counterfactual design matrices (treatment forced to 0 and to 1 for every row). This
#' lets \code{run_model_on_left_out_record_results_and_cleanup} fit each fold via
#' \code{fast_logistic_regression_cpp} and predict via a plain matrix multiply on the linear
#' predictor scale (matching \code{predict.glm}'s default \code{type = "link"}), instead of
#' calling \code{glm()}/\code{predict.glm()} (which each re-run \code{model.frame}/\code{model.matrix}
#' formula and factor-level bookkeeping, plus IRLS via the general-purpose \code{glm.fit}) tens of
#' thousands of times across the bootstrap x cross-validation loop.
#'
#' Only used when the user has NOT supplied a custom \code{personalized_model_build_function} or
#' \code{predict_function} and \code{regression_type == "incidence"}. Assumes no missing data
#' (callers must check \code{anyNA(Xy)} beforehand), since row-deletion for NAs would otherwise
#' need to happen per-fold, defeating the point of precomputing the matrices once.
#'
#' @param Xy A data frame with columns "treatment", the covariates and "y" (and possibly "censored",
#' which is dropped since it is never a predictor for the incidence regression type).
#' @keywords internal
#' @noRd
create_fast_glm_objects = function(Xy){
	Xy$censored = NULL #never a predictor for the incidence regression type

	mf_full = model.frame(y ~ . * treatment, data = Xy)
	glm_terms = attr(mf_full, "terms")
	X_obs = model.matrix(glm_terms, mf_full)
	y_full = model.response(mf_full)

	Xy_tx0 = Xy
	Xy_tx0$treatment = 0
	Xy_tx1 = Xy
	Xy_tx1$treatment = 1

	#a fit on the FULL (unresampled) data is a good, cheap-to-compute warm start for every
	#fold's IRLS: bootstrap/CV training sets are large overlapping subsets of the same n rows, so
	#their MLE is almost always close to the full-data MLE. Using this fixed full-data fit (rather
	#than chaining each fold's estimate into the next) keeps every fold's fit independent of
	#iteration order, so the bootstrap %dopar% loop stays embarrassingly parallel.
	full_fit = fast_logistic_regression_cpp(X_obs, y_full, estimate_only = TRUE)
	warm_start_beta = if (full_fit$converged) full_fit$b else NULL

	list(
		X_obs = X_obs,
		y_full = y_full,
		X_tx0 = model.matrix(glm_terms, data = Xy_tx0),
		X_tx1 = model.matrix(glm_terms, data = Xy_tx1),
		warm_start_beta = warm_start_beta
	)
}

#' Precompute design matrices for the default survival (Weibull AFT) personalization model
#'
#' Mirrors \code{create_fast_lm_objects()} but for the default survival model: builds the model
#' matrix for \code{(. - censored) * treatment} once against the FULL master data frame (matching
#' the RHS of \code{survreg(Surv(Xytrain$y, Xytrain$censored) ~ (. - censored) * treatment, ...)}),
#' along with the two counterfactual design matrices. This lets
#' \code{run_model_on_left_out_record_results_and_cleanup} fit each fold via
#' \code{fast_weibull_regression_cpp} and predict via \code{exp(linear predictor)}, which is
#' exactly what \code{predict.survreg(type = "response")} computes for \code{dist = "weibull"}
#' (its \code{itrans} is \code{exp()}), instead of calling \code{survreg()}/\code{predict.survreg()}
#' tens of thousands of times across the bootstrap x cross-validation loop.
#'
#' Note the event/censoring indicator: the default model builder fits
#' \code{survreg(Surv(Xytrain$y, Xytrain$censored) ~ ...)}, i.e. it feeds the \code{censored}
#' column directly as \code{Surv()}'s event-status argument (1 = event/death, 0 = censored per the
#' \code{survival} package convention), so despite its name the \code{censored} column IS the
#' "dead" indicator expected by \code{fast_weibull_regression_cpp} -- it must NOT be flipped here.
#'
#' Only used when the user has NOT supplied a custom \code{personalized_model_build_function} or
#' \code{predict_function} and \code{regression_type == "survival"}. Assumes no missing data
#' (callers must check \code{anyNA(Xy)} beforehand), since row-deletion for NAs would otherwise
#' need to happen per-fold, defeating the point of precomputing the matrices once.
#'
#' @param Xy A data frame with columns "treatment", "censored", the covariates and "y".
#' @keywords internal
#' @noRd
create_fast_weibull_objects = function(Xy){
	mf_full = model.frame(y ~ (. - censored) * treatment, data = Xy)
	weibull_terms = attr(mf_full, "terms")
	X_obs = model.matrix(weibull_terms, mf_full)
	y_full = model.response(mf_full)
	dead_full = Xy$censored #see note above: "censored" is fed directly as Surv()'s event/dead argument

	Xy_tx0 = Xy
	Xy_tx0$treatment = 0
	Xy_tx1 = Xy
	Xy_tx1$treatment = 1

	#see the matching note in create_fast_glm_objects(): a fixed full-data fit is a cheap, good
	#warm start for every fold, and keeps folds independent of each other / of iteration order.
	full_fit = fast_weibull_regression_cpp(X_obs, y_full, dead_full, estimate_only = TRUE)
	warm_start_params = if (full_fit$converged) c(full_fit$b, full_fit$log_sigma) else NULL

	list(
		X_obs = X_obs,
		y_full = y_full,
		dead_full = dead_full,
		X_tx0 = model.matrix(weibull_terms, data = Xy_tx0),
		X_tx1 = model.matrix(weibull_terms, data = Xy_tx1),
		warm_start_params = warm_start_params
	)
}

run_default_fast_cv_results = function(
		Xy,
		y_higher_is_better,
		cutoff_obj,
		boot_idx,
		fast_lm_objects = NULL,
		fast_glm_objects = NULL,
		fast_weibull_objects = NULL,
		verbose = FALSE,
		full_verbose = FALSE){

	if (verbose || full_verbose){
		return(NULL)
	}

	begin_cutoffs = as.integer(cutoff_obj$begin_cutoffs_for_leave_outs)
	end_cutoffs = as.integer(cutoff_obj$end_cutoffs_for_leave_outs)
	if (length(boot_idx) != max(end_cutoffs)){
		return(NULL)
	}
	boot_idx = as.integer(boot_idx)

	if (!is.null(fast_lm_objects)){
		res = fast_default_continuous_cv_run_cpp(
			fast_lm_objects$X_obs,
			fast_lm_objects$y_full,
			fast_lm_objects$X_tx0,
			fast_lm_objects$X_tx1,
			Xy$treatment,
			Xy$censored,
			boot_idx,
			begin_cutoffs,
			end_cutoffs,
			y_higher_is_better
		)
	} else if (!is.null(fast_glm_objects)){
		res = fast_default_incidence_cv_run_cpp(
			fast_glm_objects$X_obs,
			fast_glm_objects$y_full,
			fast_glm_objects$X_tx0,
			fast_glm_objects$X_tx1,
			Xy$treatment,
			Xy$censored,
			boot_idx,
			begin_cutoffs,
			end_cutoffs,
			fast_glm_objects$warm_start_beta,
			y_higher_is_better
		)
	} else if (!is.null(fast_weibull_objects)){
		res = fast_default_weibull_cv_run_cpp(
			fast_weibull_objects$X_obs,
			fast_weibull_objects$y_full,
			fast_weibull_objects$dead_full,
			fast_weibull_objects$X_tx0,
			fast_weibull_objects$X_tx1,
			Xy$treatment,
			boot_idx,
			begin_cutoffs,
			end_cutoffs,
			fast_weibull_objects$warm_start_params,
			y_higher_is_better
		)
	} else {
		return(NULL)
	}

	if (isTRUE(res$ok)){
		res$raw_results
	} else {
		NULL
	}
}

run_default_fast_bootstrap_q_scores = function(
		Xy,
		y_higher_is_better,
		cutoff_obj,
		B,
		incidence_metric,
		fast_lm_objects = NULL,
		fast_glm_objects = NULL,
		fast_weibull_objects = NULL){

	begin_cutoffs = as.integer(cutoff_obj$begin_cutoffs_for_leave_outs)
	end_cutoffs = as.integer(cutoff_obj$end_cutoffs_for_leave_outs)

	if (!is.null(fast_lm_objects)){
		fast_default_continuous_bootstrap_q_cpp(
			fast_lm_objects$X_obs,
			fast_lm_objects$y_full,
			fast_lm_objects$X_tx0,
			fast_lm_objects$X_tx1,
			Xy$treatment,
			Xy$censored,
			as.integer(B),
			begin_cutoffs,
			end_cutoffs,
			y_higher_is_better
		)
	} else if (!is.null(fast_glm_objects)){
		fast_default_incidence_bootstrap_q_cpp(
			fast_glm_objects$X_obs,
			fast_glm_objects$y_full,
			fast_glm_objects$X_tx0,
			fast_glm_objects$X_tx1,
			Xy$treatment,
			as.integer(B),
			begin_cutoffs,
			end_cutoffs,
			fast_glm_objects$warm_start_beta,
			y_higher_is_better,
			incidence_metric
		)
	} else if (!is.null(fast_weibull_objects)){
		fast_default_weibull_bootstrap_q_cpp(
			fast_weibull_objects$X_obs,
			fast_weibull_objects$y_full,
			fast_weibull_objects$dead_full,
			fast_weibull_objects$X_tx0,
			fast_weibull_objects$X_tx1,
			Xy$treatment,
			as.integer(B),
			begin_cutoffs,
			end_cutoffs,
			fast_weibull_objects$warm_start_params,
			y_higher_is_better
		)
	} else {
		NULL
	}
}

mean_from_sum_count = function(sum_x, n_x){
	sum_x / n_x
}

#' Fit + predict both counterfactual treatment arms, treating a hard fit/predict error as a bad fold
#'
#' A small bootstrap fold (e.g. a low m-out-of-n resample size) can leave a factor covariate constant
#' within that fold's training data even though the covariate has 2+ declared levels overall (a rare
#' category simply wasn't drawn) -- \code{personalized_model_build_function()}/\code{predict_function()},
#' which may be raw \code{lm()}/\code{glm()}/\code{survreg()}, then throw a hard error (e.g. "contrasts
#' can be applied only to factors with 2 or more levels") rather than returning a degenerate fit. The
#' fast paths' own rank-deficiency/non-convergence/separation guards (see
#' \code{run_model_on_left_out_record_results_and_cleanup()}) only catch problems detectable from the
#' precomputed design matrices, not this. Catching the error here and returning \code{NA} predictions
#' instead lets the existing \code{NA}/\code{is_bad} machinery downstream flag the replicate as bad
#' (matching how rank-deficient/non-convergent folds are already handled) instead of aborting the
#' entire bootstrap run over one degenerate fold.
#'
#' @keywords internal
#' @noRd
safe_build_and_predict_counterfactuals = function(Xytrain, Xyleftout, personalized_model_build_function, predict_function){
	tryCatch({
		mod = personalized_model_build_function(Xytrain)
		Xyleftout$treatment = 0
		yhatTx0s = predict_function(mod, Xyleftout)
		Xyleftout$treatment = 1
		yhatTx1s = predict_function(mod, Xyleftout)
		list(yhatTx0s = yhatTx0s, yhatTx1s = yhatTx1s)
	}, error = function(e){
		n_test = nrow(Xyleftout)
		list(yhatTx0s = rep(NA_real_, n_test), yhatTx1s = rep(NA_real_, n_test))
	})
}

#' Build only the q-score fields needed inside bootstrap resampling loops
#'
#' @param results A named list from \code{create_raw_results_list()}/\code{assign_fold_results()}.
#' @keywords internal
#' @noRd
create_PTE_q_score_object = function(results, regression_type, y_higher_is_better, difference_function, incidence_metric){
	if (!is.null(difference_function)){
		return(create_PTE_results_object(results, regression_type, y_higher_is_better, difference_function, incidence_metric))
	}

	if (regression_type == "continuous"){
		default_continuous_q_scores_cpp(results$real_y, results$given_tx, results$rec_tx, y_higher_is_better)
	} else if (regression_type == "incidence"){
		default_incidence_q_scores_cpp(results$real_y, results$given_tx, results$rec_tx, y_higher_is_better, incidence_metric)
	} else {
		default_survival_q_scores_cpp(results$real_y, results$censored, results$given_tx, results$rec_tx, y_higher_is_better)
	}
}

#' Fit and evaluate one cross-validation fold's left-out subjects
#'
#' Trains on all subjects except \code{train_on_all_except_these} and evaluates the subjects in
#' \code{leave_outs_to_be_predicted} under both counterfactual treatments, using whichever of the
#' three precomputed \code{fast_*_objects} is non-\code{NULL} (see \code{create_fast_lm_objects()},
#' \code{create_fast_glm_objects()}, \code{create_fast_weibull_objects()} above) to skip
#' \code{personalized_model_build_function}/\code{predict_function} entirely for the default
#' continuous/incidence/survival models, falling back to those exact user-supplied functions
#' whenever no fast path applies or a fold's fast fit doesn't converge/is rank-deficient.
#'
#' @keywords internal
#' @noRd
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
		boot_idx = NULL,
		fast_glm_objects = NULL,
		fast_weibull_objects = NULL){

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
			Xytrain = fast_row_subset(Xy, train_rows)
			Xytrain$censored = NULL
			Xyleftout = fast_row_subset(Xy, test_rows, 1 : (ncol(Xy) - 1))
			pred = safe_build_and_predict_counterfactuals(Xytrain, Xyleftout, personalized_model_build_function, predict_function)
			yhatTx0s = pred$yhatTx0s
			yhatTx1s = pred$yhatTx1s
		}

		real_ys = Xy$y[test_rows]
		orig_trts = Xy$treatment[test_rows]
		censored_col = Xy$censored[test_rows]

		if (verbose && !full_verbose){
			cat(".")
		}
	} else if (!is.null(fast_glm_objects)){
		#fast path for the default incidence (logistic) model: mirrors the fast_lm_objects branch
		#above but fits via fast_logistic_regression_cpp() (IRLS in C++/Eigen) instead of glm()/
		#glm.fit(), and predicts on the link (log-odds) scale to match predict.glm()'s default
		#type = "link" -- see create_fast_glm_objects() above
		train_rows = boot_idx[-train_on_all_except_these]
		test_rows = boot_idx[leave_outs_to_be_predicted]

		X_train = fast_glm_objects$X_obs[train_rows, , drop = FALSE]
		#warm-starting from the full-data fit cuts IRLS iterations substantially since every fold
		#trains on a large overlapping subset of the same n rows (see create_fast_glm_objects())
		fit = fast_logistic_regression_cpp(X_train, fast_glm_objects$y_full[train_rows],
			warm_start_beta = fast_glm_objects$warm_start_beta, estimate_only = TRUE)

		#a (quasi-)separated fold has no real MLE -- the log-likelihood keeps improving as some
		#coefficient runs off to +-Inf, so IRLS "converges" (score norm < tol) at an essentially
		#arbitrary point on that ridge, one that depends on the starting value. Since we warm-start
		#from the full-data fit rather than glm()'s own beta=0 start, that arbitrary point can land
		#meaningfully far from what glm()/predict.glm() would have reported for this exact fold.
		#None of this dataset's genuinely-identified coefficients exceed a magnitude of a few units,
		#so a fitted linear predictor this large (predicted probability within ~1e-7 of 0 or 1) is
		#itself the signature of separation -- treat it the same as non-convergence and fall back to
		#the exact computation rather than report a fit-path-dependent (and therefore not
		#reproducible via the documented glm() default) coefficient.
		separated = fit$converged && max(abs(X_train %*% fit$b)) > 30

		if (fit$converged && !separated){
			beta = fit$b
			yhatTx0s = as.numeric(fast_glm_objects$X_tx0[test_rows, , drop = FALSE] %*% beta)
			yhatTx1s = as.numeric(fast_glm_objects$X_tx1[test_rows, , drop = FALSE] %*% beta)
		} else {
			#rare non-convergent or (quasi-)separated fold (e.g. a bootstrap resample induces
			#near-perfect separation): fall back to the exact original computation for this fold
			#only, to guarantee identical behavior to glm()/predict.glm() rather than risk silently
			#divergent coefficients.
			Xytrain = fast_row_subset(Xy, train_rows)
			Xytrain$censored = NULL
			Xyleftout = fast_row_subset(Xy, test_rows, 1 : (ncol(Xy) - 1))
			pred = safe_build_and_predict_counterfactuals(Xytrain, Xyleftout, personalized_model_build_function, predict_function)
			yhatTx0s = pred$yhatTx0s
			yhatTx1s = pred$yhatTx1s
		}

		real_ys = Xy$y[test_rows]
		orig_trts = Xy$treatment[test_rows]
		censored_col = Xy$censored[test_rows]

		if (verbose && !full_verbose){
			cat(".")
		}
	} else if (!is.null(fast_weibull_objects)){
		#fast path for the default survival (Weibull AFT) model: mirrors the fast_lm_objects branch
		#above but fits via fast_weibull_regression_cpp() (Newton/L-BFGS in C++/Eigen) instead of
		#survreg(), and predicts on the response scale via exp(linear predictor), which is exactly
		#what predict.survreg(type = "response") computes for dist = "weibull" (its itrans is exp())
		#-- see create_fast_weibull_objects() above
		train_rows = boot_idx[-train_on_all_except_these]
		test_rows = boot_idx[leave_outs_to_be_predicted]

		X_train = fast_weibull_objects$X_obs[train_rows, , drop = FALSE]
		#warm-starting from the full-data fit (see create_fast_weibull_objects()) for the same
		#reason as the logistic path above
		fit = fast_weibull_regression_cpp(X_train, fast_weibull_objects$y_full[train_rows],
			fast_weibull_objects$dead_full[train_rows],
			warm_start_params = fast_weibull_objects$warm_start_params, estimate_only = TRUE)

		#same (quasi-)separation guard as the logistic path above: a linear predictor this extreme
		#(implying survival times on the order of exp(15) or exp(-15)) only happens when the warm
		#start has pushed a fold's fit to an arbitrary, fit-path-dependent point on a ridge rather
		#than a real optimum -- fall back rather than report it.
		separated = fit$converged && max(abs(X_train %*% fit$b)) > 30

		if (fit$converged && !separated){
			beta = fit$b
			yhatTx0s = exp(as.numeric(fast_weibull_objects$X_tx0[test_rows, , drop = FALSE] %*% beta))
			yhatTx1s = exp(as.numeric(fast_weibull_objects$X_tx1[test_rows, , drop = FALSE] %*% beta))
		} else {
			#rare non-convergent or (quasi-)separated fold: fall back to the exact original
			#computation for this fold only, to guarantee identical behavior to
			#survreg()/predict.survreg().
			Xytrain = fast_row_subset(Xy, train_rows)
			Xyleftout = fast_row_subset(Xy, test_rows, 1 : (ncol(Xy) - 1))
			pred = safe_build_and_predict_counterfactuals(Xytrain, Xyleftout, personalized_model_build_function, predict_function)
			yhatTx0s = pred$yhatTx0s
			yhatTx1s = pred$yhatTx1s
		}

		real_ys = Xy$y[test_rows]
		orig_trts = Xy$treatment[test_rows]
		censored_col = Xy$censored[test_rows]

		if (verbose && !full_verbose){
			cat(".")
		}
	} else {
		#the left one out matrix has n-1 rows and will be considered the "training data"
		Xytrain = fast_row_subset(Xy, -train_on_all_except_these)

		#pull out the record of the left-one-out subject
		Xyleftout = fast_row_subset(Xy, leave_outs_to_be_predicted, 1 : (ncol(Xy) - 1)) #leave out y

		if (regression_type != "survival"){
			Xytrain$censored = NULL
		}

		#also take note of what actually happened to this subject in the experiment (independent of
		#whether the model build/predict below succeeds)
		real_ys = .subset2(Xy, ncol(Xy))[leave_outs_to_be_predicted]
		orig_trts = Xyleftout$treatment

		#build the model via the user-specified function (makes use of the "Xytrain" object) and
		#evaluate the left-one-out subject on it for both his true treatment and his counterfactual.
		#A fold that leaves a factor covariate constant (e.g. a rare category absent from a small
		#m-out-of-n resample) can make personalized_model_build_function()/predict_function() throw a
		#hard error (e.g. raw lm()'s "contrasts can be applied only to factors with 2 or more levels")
		#rather than return a degenerate fit -- see safe_build_and_predict_counterfactuals().
		pred = tryCatch({
			mod = personalized_model_build_function(Xytrain)
			if (full_verbose){
				print(summary(mod))
			}
			Xyleftout$treatment = 0
			yhatTx0s = predict_function(mod, Xyleftout)
			Xyleftout$treatment = 1
			yhatTx1s = predict_function(mod, Xyleftout)
			list(yhatTx0s = yhatTx0s, yhatTx1s = yhatTx1s)
		}, error = function(e){
			n_test = nrow(Xyleftout)
			list(yhatTx0s = rep(NA_real_, n_test), yhatTx1s = rep(NA_real_, n_test))
		})
		yhatTx0s = pred$yhatTx0s
		yhatTx1s = pred$yhatTx1s
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

	#tabulate the result for the prediction on this left one out model (vectorized over all left-out subjects).
	#orig_trts/optimal are 0/1 (resp. logical 0/1), so plain arithmetic reproduces ifelse()'s
	#per-element selection without its dispatch overhead
	est_true = yhatTx0s + orig_trts * (yhatTx1s - yhatTx0s)
	est_counterfactual = yhatTx0s + yhatTx1s - est_true
	optimal = if (y_higher_is_better){
		est_true > est_counterfactual
	} else {
		est_true < est_counterfactual
	}
	rec_tx = orig_trts * optimal + (1 - orig_trts) * (1 - optimal)

	list(
		est_true = est_true,
		est_counterfactual = est_counterfactual,
		given_tx = orig_trts,
		rec_tx = rec_tx,
		real_y = real_ys,
		censored = censored_col
	)
}

#' Build the \code{PTE_bootstrap_results} q-score/summary object for one run's raw results
#'
#' @param results A named list from \code{create_raw_results_list()}/\code{assign_fold_results()}
#' (or, for a custom \code{difference_function}, whatever shape that function expects).
#' @keywords internal
#' @noRd
create_PTE_results_object = function(results, regression_type, y_higher_is_better, difference_function, incidence_metric){
	n = length(results$real_y)

	#build return object for the user
	return_obj = list()
	return_obj$results = results

	#get real y's (results is a plain named list of equal-length vectors -- see create_raw_results_list())
	y_all = results$real_y
	indices_1_1 = results$given_tx == 1 & results$rec_tx == 1 #optimal
	indices_0_0 = results$given_tx == 0 & results$rec_tx == 0 #optimal
	indices_0_1 = results$given_tx == 0 & results$rec_tx == 1 #non-optimal
	indices_1_0 = results$given_tx == 1 & results$rec_tx == 0 #non-optimal

	if (is.null(difference_function)){
		Y00 = y_all[indices_0_0]
		Y01 = y_all[indices_0_1]
		Y10 = y_all[indices_1_0]
		Y11 = y_all[indices_1_1]
		return_obj$Y00 = Y00
		return_obj$Y01 = Y01
		return_obj$Y10 = Y10
		return_obj$Y11 = Y11

		#four groups summary square - less useful for survival, but still nice to see
		summary_square = matrix(NA, 3, 3)
		rownames(summary_square) = c("received Trt 0", "received Trt 1", "Avg Opt")
		colnames(summary_square) = c("optimal Trt 0", "optimal Trt 1", "Avg Trt")
		summary_square[1, 1] = mean(Y00, na.rm = TRUE)
		summary_square[1, 2] = mean(Y01, na.rm = TRUE)
		summary_square[2, 1] = mean(Y10, na.rm = TRUE)
		summary_square[2, 2] = mean(Y11, na.rm = TRUE)
		summary_square[1, 3] = mean(c(Y00, Y01), na.rm = TRUE)
		summary_square[2, 3] = mean(c(Y10, Y11), na.rm = TRUE)
		summary_square[3, 1] = mean(c(Y00, Y10), na.rm = TRUE)
		summary_square[3, 2] = mean(c(Y01, Y11), na.rm = TRUE)
		summary_square[3, 3] = mean(y_all, na.rm = TRUE)
		return_obj$summary_square = summary_square
		return_obj$is_bad = is.nan(sum(summary_square)) #if any of them are NaN this will pick it up

		#four groups
		ns = matrix(NA, 3, 3)
		rownames(ns) = c("received Trt 0", "received Trt 1", "Totals")
		colnames(ns) = c("optimal Trt 0", "optimal Trt 1", "Totals")
		ns[1, 1] = length(Y00[!is.na(Y00)])
		ns[1, 2] = length(Y01[!is.na(Y01)])
		ns[2, 1] = length(Y10[!is.na(Y10)])
		ns[2, 2] = length(Y11[!is.na(Y11)])
		#now add em up to get tots
		ns[1, 3] = ns[1, 1] + ns[1, 2]
		ns[2, 3] = ns[2, 1] + ns[2, 2]
		ns[3, 1] = ns[1, 1] + ns[2, 1]
		ns[3, 2] = ns[1, 2] + ns[2, 2]
		ns[3, 3] = ns[1, 1] + ns[1, 2] + ns[2, 1] + ns[2, 2]
		return_obj$ns = ns
		return_obj$pct_data_used = round(ns[3, 3] / n * 100, 3)




		if (regression_type == "continuous" || (regression_type == "incidence" && incidence_metric == "probability_difference")){
			return_obj$pred_differences_avg = mean(abs(results$est_true - results$est_counterfactual), na.rm = TRUE)
			return_obj$pred_differences_sd = sd(abs(results$est_true - results$est_counterfactual), na.rm = TRUE)
			return_obj$avg_rec = mean(c(Y00, Y11), na.rm = TRUE)
			return_obj$avg_non_rec = mean(c(Y01, Y10), na.rm = TRUE)
			return_obj$avg_all = mean(y_all, na.rm = TRUE)
			avg_ys_tx_1 = mean(c(Y00, Y01), na.rm = TRUE)
			avg_ys_tx_2 = mean(c(Y10, Y11), na.rm = TRUE)
			if (avg_ys_tx_1 >= avg_ys_tx_2 && y_higher_is_better){
				return_obj$avg_best = avg_ys_tx_1
			} else if (avg_ys_tx_1 >= avg_ys_tx_2 && !y_higher_is_better){
				return_obj$avg_best = avg_ys_tx_2
			} else if (avg_ys_tx_1 < avg_ys_tx_2 && y_higher_is_better){
				return_obj$avg_best = avg_ys_tx_2
			} else if (avg_ys_tx_1 < avg_ys_tx_2 && !y_higher_is_better){
				return_obj$avg_best = avg_ys_tx_1
			}
			return_obj$q_adversarial = return_obj$avg_rec - return_obj$avg_non_rec
			return_obj$q_average = return_obj$avg_rec - return_obj$avg_all
			return_obj$q_best = return_obj$avg_rec - return_obj$avg_best
			#some more metrics that may be of use for the continuous case
			if (regression_type == "continuous" ){
				sse = sum((results$real_y - results$est_true)^2, na.rm = TRUE)
				return_obj$oos_rmse = sqrt(sse / n)
				sst = sum((results$real_y - mean(results$real_y, na.rm = TRUE))^2, na.rm = TRUE)
				return_obj$out_of_sample_Rsq = 1 - sse / sst
			}
		} else if (regression_type == "incidence"){
			#set up all the data
			p_all = mean(y_all, na.rm = TRUE)
			p_rec = mean(c(Y00, Y11), na.rm = TRUE)
			p_non_rec = mean(c(Y01, Y10), na.rm = TRUE)
			p_1 = mean(c(Y10, Y11), na.rm = TRUE)
			p_0 = mean(c(Y00, Y01), na.rm = TRUE)

			if (p_1 >= p_0 && y_higher_is_better){
				p_best = p_1
			} else if (p_1 < p_0 && y_higher_is_better){
				p_best = p_0
			} else if (p_1 <= p_0 && !y_higher_is_better){
				p_best = p_0
			} else {
				p_best = p_1
			}

			if (incidence_metric == "risk_ratio"){
				return_obj$q_adversarial = p_rec / p_non_rec
				return_obj$q_average = p_rec / p_all
				return_obj$q_best = p_rec / p_best
			} else if (incidence_metric == "odds_ratio"){
				p_rec_odds = (p_rec / (1 - p_rec))
				return_obj$q_adversarial = p_rec_odds / (p_non_rec / (1 - p_non_rec))
				return_obj$q_average = p_rec_odds / (p_all / (1 - p_all))
				return_obj$q_best = p_rec_odds / (p_best / (1 - p_best))
			}
		} else if (regression_type == "survival"){ #this set should be collectively exhaustive
			c_all = results$censored

			#get_survival_stat_diff_cpp(y, dead, w) computes the KM median directly for the two
			#groups indicated by w (1 vs 0) and returns median(w==1) - median(w==0), NA if either
			#group's median is unreached or empty -- exactly what
			#summary(survfit(Surv(y, dead) ~ w))$table[, 'median'][2] - [...][1] computed before,
			#but without paying for survfit()/summary.survfit()'s full table/CI/formula machinery
			#(which now dominates this path since the Weibull fit itself got fast -- see fast_survival_stats.cpp)

			#first do adversarial
			y_rec = c(Y00, Y11)
			y_non_rec = c(Y01, Y10)
			c_rec = c(c_all[indices_0_0], c_all[indices_1_1])
			c_non_rec = c(c_all[indices_0_1], c_all[indices_1_0])

			return_obj$q_adversarial = get_survival_stat_diff_cpp(
				c(y_rec, y_non_rec), c(c_rec, c_non_rec), c(rep(1L, length(y_rec)), rep(0L, length(y_non_rec)))
			)

			#now do average
			return_obj$q_average = get_survival_stat_diff_cpp(
				c(y_rec, y_all), c(c_rec, c_all), c(rep(1L, length(y_rec)), rep(0L, length(y_all)))
			)

			#now do best... this requires running two models... the first to see who's better on "average" (or
			#at the median) and the second to run the best tx vs the recommended tx.
			y_1 = c(Y10, Y11)
			c_1 = c(c_all[indices_1_0], c_all[indices_1_1])
			y_0 = c(Y01, Y00)
			c_0 = c(c_all[indices_0_1], c_all[indices_0_0])
			median_diff = get_survival_stat_diff_cpp(c(y_1, y_0), c(c_1, c_0), c(rep(1L, length(y_1)), rep(0L, length(y_0))))
			#if median_diff is greater than 0, tx1 is the better treatment
			if (median_diff == 0){ #measure 0 event in theory... but in practice - who knows?
				random_bernoulli = runif(1) < 0.5
				y_best = if (random_bernoulli){y_0} else {y_1}
				c_best = if (random_bernoulli){c_0} else {c_1}
			} else if (median_diff > 0 && y_higher_is_better){
				y_best = y_1
				c_best = c_1
			} else if (median_diff < 0 && y_higher_is_better){
				y_best = y_0
				c_best = c_0
			} else if (median_diff > 0 && !y_higher_is_better){
				y_best = y_0
				c_best = c_0
			} else {
				y_best = y_1
				c_best = c_1
			}
			return_obj$q_best = get_survival_stat_diff_cpp(
				c(y_rec, y_best), c(c_rec, c_best), c(rep(1L, length(y_rec)), rep(0L, length(y_best)))
			)
		}
	} else { ###user custom function output
		all_diffs = difference_function(results, indices_1_1, indices_0_0, indices_0_1, indices_1_0)
		return_obj$q_adversarial = all_diffs[1]
		return_obj$q_average = all_diffs[2]
		return_obj$q_best = all_diffs[3]
		return_obj$is_bad = is.nan(sum(all_diffs)) || any(is.na(all_diffs))
	}
	return_obj
}
