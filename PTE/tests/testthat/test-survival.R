# regression_type = "survival": y is a time-to-event response, censored is required, default
# model is survreg(Surv(y, censored) ~ (. - censored) * treatment, dist = "weibull").

dat = make_survival_data()
B = 15

test_that("basic structure and contract hold, H_0_mu_equals defaults to 0", {
	res = suppressWarnings(
		PTE_bootstrap_inference(
			dat$X, dat$y, censored = dat$censored, regression_type = "survival", B = B, num_cores = 1
		)
	)
	expect_valid_pte_result(res, B)
	expect_identical(res$H_0_mu_equals, 0)
})

test_that("the fitted personalization_model is a survreg fit with a censored column retained", {
	res = suppressWarnings(
		PTE_bootstrap_inference(
			dat$X, dat$y, censored = dat$censored, regression_type = "survival", B = 5, num_cores = 1
		)
	)
	expect_s3_class(res$personalization_model, "survreg")
	expect_true("censored" %in% colnames(res$Xy))
})

test_that("H_0_mu_equals can be overridden for survival", {
	res_h0 = suppressWarnings(
		PTE_bootstrap_inference(
			dat$X, dat$y, censored = dat$censored, regression_type = "survival", B = 5, num_cores = 1,
			H_0_mu_equals = 1
		)
	)
	expect_identical(res_h0$H_0_mu_equals, 1)
})

test_that("y_higher_is_better = FALSE runs and produces a valid result for survival data", {
	# lower time-to-event being "better" is a documented, realistic use case (e.g. time to
	# symptom resolution) -- see the y_higher_is_better docs
	res_flip = suppressWarnings(
		PTE_bootstrap_inference(
			dat$X, dat$y, censored = dat$censored, regression_type = "survival", B = B, num_cores = 1,
			y_higher_is_better = FALSE
		)
	)
	expect_valid_pte_result(res_flip, B)
})

test_that("a custom difference_function is honored for survival data", {
	custom_diff = function(results, indices_1_1, indices_0_0, indices_0_1, indices_1_0){
		y_all = results$real_y
		rec = c(y_all[indices_1_1], y_all[indices_0_0])
		non_rec = c(y_all[indices_0_1], y_all[indices_1_0])
		c(mean(rec) - mean(non_rec), mean(rec) - mean(y_all), mean(rec) - max(mean(y_all[results$given_tx == 0]), mean(y_all[results$given_tx == 1])))
	}

	res_diff = suppressWarnings(
		PTE_bootstrap_inference(
			dat$X, dat$y, censored = dat$censored, regression_type = "survival", B = B, num_cores = 1,
			difference_function = custom_diff
		)
	)
	expect_valid_pte_result(res_diff, B)
})
