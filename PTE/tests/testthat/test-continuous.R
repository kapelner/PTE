# regression_type = "continuous": y is a real-valued response, default model is
# lm(y ~ . * treatment). Fixture has a strong, deterministic treatment*x interaction so that
# even with a small B the sign/rough magnitude of the estimated effect is stable across the
# runs (the bootstrap resampling itself is not seedable -- see helper-fixtures.R).

dat = make_continuous_data()
B = 30

res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "continuous", B = B, num_cores = 2)

test_that("basic structure and contract hold", {
	expect_valid_pte_result(res, B)
})

test_that("H_0_mu_equals defaults to 0 for continuous", {
	expect_identical(res$H_0_mu_equals, 0)
})

test_that("no bootstrap samples are flagged invalid for a well-behaved dataset", {
	expect_identical(res$num_bad, 0)
})

test_that("a strong true effect is detected (low p-value, CI excludes 0)", {
	expect_lt(res$p_val_average, 0.2)
	expect_lt(res$p_val_best, 0.2)
	expect_gt(res$ci_q_average[1], 0)
	expect_gt(res$est_q_average, 0)
})

test_that("the fitted personalization_model is an lm fit on the full data", {
	expect_s3_class(res$personalization_model, "lm")
	expect_identical(nobs(res$personalization_model), nrow(dat$X))
})

test_that("returned Xy no longer carries a censored column for non-survival types", {
	expect_false("censored" %in% colnames(res$Xy))
})

test_that("custom personalized_model_build_function and predict_function are honored", {
	custom_build = function(Xytrain) lm(y ~ . * treatment, data = Xytrain)
	custom_predict = function(mod, Xyleftout) predict(mod, Xyleftout)

	res_custom = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = B, num_cores = 2,
		personalized_model_build_function = custom_build,
		predict_function = custom_predict
	)

	expect_valid_pte_result(res_custom, B)
	expect_identical(res_custom$num_bad, 0)
	# should reproduce essentially the same statistical conclusion as the built-in default
	expect_lt(res_custom$p_val_average, 0.2)
	expect_gt(res_custom$est_q_average, 0)
})

test_that("a custom difference_function is honored and produces a valid result", {
	custom_diff = function(results, indices_1_1, indices_0_0, indices_0_1, indices_1_0){
		y_all = results$real_y
		rec = c(y_all[indices_1_1], y_all[indices_0_0])
		non_rec = c(y_all[indices_0_1], y_all[indices_1_0])
		best = max(mean(y_all[results$given_tx == 0]), mean(y_all[results$given_tx == 1]))
		c(mean(rec) - mean(non_rec), mean(rec) - mean(y_all), mean(rec) - best)
	}

	res_diff = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = B, num_cores = 2,
		difference_function = custom_diff
	)

	expect_valid_pte_result(res_diff, B)
	expect_gt(res_diff$est_q_average, 0)
})

test_that("H_0_mu_equals can be overridden", {
	res_h0 = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = 5, num_cores = 1,
		H_0_mu_equals = 0.5
	)
	expect_identical(res_h0$H_0_mu_equals, 0.5)
})

test_that("y_higher_is_better = FALSE flips the direction of the estimated effect", {
	res_flip = PTE_bootstrap_inference(
		dat$X, -dat$y, regression_type = "continuous", B = B, num_cores = 2,
		y_higher_is_better = FALSE
	)
	expect_valid_pte_result(res_flip, B)
	# negating y and flipping the "better" direction should recover the same sign of effect
	expect_gt(res_flip$est_q_average, 0)
	expect_lt(res_flip$p_val_average, 0.2)
})

test_that("run_bca_bootstrap = TRUE additionally produces valid BCa confidence intervals", {
	dat_small = make_continuous_data(n = 40)
	res_bca = PTE_bootstrap_inference(
		dat_small$X, dat_small$y, regression_type = "continuous", B = 15, num_cores = 2,
		run_bca_bootstrap = TRUE
	)
	expect_valid_pte_result(res_bca, 15, expect_bca = TRUE)
	expect_true(res_bca$run_bca_bootstrap)
})

test_that("num_bad is reported and a warning is issued when most bootstrap samples are degenerate", {
	# a single bootstrap replicate of size 4 with one covariate is prone to producing a
	# degenerate (constant) recommendation, which create_PTE_results_object flags as "bad"
	tiny_X = data.frame(treatment = rep(0:1, each = 2), x = c(-1, 1, -1, 1))
	tiny_y = c(-1, 1, 1, -1)
	expect_warning(
		res_tiny <- PTE_bootstrap_inference(
			tiny_X, tiny_y, regression_type = "continuous", B = 10, num_cores = 1,
			pct_leave_out = 0.5
		),
		"invalid"
	)
	expect_gt(res_tiny$num_bad, 0)
})
