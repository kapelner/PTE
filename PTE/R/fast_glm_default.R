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

	list(
		X_obs = X_obs,
		y_full = y_full,
		X_tx0 = model.matrix(glm_terms, data = Xy_tx0),
		X_tx1 = model.matrix(glm_terms, data = Xy_tx1)
	)
}
