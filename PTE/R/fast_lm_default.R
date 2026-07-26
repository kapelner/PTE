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
