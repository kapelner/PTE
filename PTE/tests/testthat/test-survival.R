# regression_type = "survival": y is a time-to-event response, censored is required, default
# model is survreg(Surv(y, censored) ~ (. - censored) * treatment, dist = "weibull").

dat = make_survival_data()
B = 15

test_that("basic structure and contract hold, H_0_mu_equals defaults to 0", {
	res = suppressWarnings(
		PTE_bootstrap_inference(
			dat$X, dat$y, censored = dat$censored, regression_type = "survival", B = B, num_cores = 2
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
