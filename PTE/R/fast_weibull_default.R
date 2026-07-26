#' Precompute design matrices for the default survival (Weibull AFT) personalization model
#'
#' Mirrors \code{create_fast_lm_objects()} but for the default survival model: builds the model
#' matrix for \code{(. - censored) * treatment} once against the FULL master data frame (matching
#' the RHS of \code{survreg(Surv(Xytrain$y, Xytrain$censored) ~ (. - censored) * treatment, ...)}),
#' along with the two counterfactual design matrices. This lets
#' \code{run_model_on_left_out_record_results_and_cleanup} fit each fold via
#' \code{fast_weibull_regression_cpp} and predict via a plain matrix multiply on the linear
#' predictor scale, instead of calling \code{survreg()}/\code{predict.survreg()} tens of thousands
#' of times across the bootstrap x cross-validation loop.
#'
#' The linear predictor (rather than \code{predict.survreg}'s default \code{type = "response"},
#' which additionally scales by \code{exp()} and a \code{Gamma(1 + scale)} factor common to both
#' counterfactuals) is sufficient here: \code{create_PTE_results_object} determines the survival
#' q-scores from the real (not predicted) response and censoring status, so only the *sign* of the
#' difference between the two counterfactual predictions -- which decides the recommended treatment
#' -- matters, and that sign is unchanged by a monotonic increasing transform shared across both
#' counterfactuals.
#'
#' Only used when the user has NOT supplied a custom \code{personalized_model_build_function} or
#' \code{predict_function} and \code{regression_type == "survival"}. Assumes no missing data
#' (callers must check \code{anyNA(Xy)} beforehand), since row-deletion for NAs would otherwise
#' need to happen per-fold, defeating the point of precomputing the matrices once.
#'
#' @param Xy A data frame with columns "treatment", "censored", the covariates and "y".
create_fast_weibull_objects = function(Xy){
	mf_full = model.frame(y ~ (. - censored) * treatment, data = Xy)
	weibull_terms = attr(mf_full, "terms")
	X_obs = model.matrix(weibull_terms, mf_full)
	y_full = model.response(mf_full)
	dead_full = 1 - Xy$censored

	Xy_tx0 = Xy
	Xy_tx0$treatment = 0
	Xy_tx1 = Xy
	Xy_tx1$treatment = 1

	list(
		X_obs = X_obs,
		y_full = y_full,
		dead_full = dead_full,
		X_tx0 = model.matrix(weibull_terms, data = Xy_tx0),
		X_tx1 = model.matrix(weibull_terms, data = Xy_tx1)
	)
}
