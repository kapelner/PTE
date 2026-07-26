# regression_type = "incidence": y is binary, default model is glm(y ~ . * treatment, family
# = "binomial"). Three incidence_metric options change both H_0_mu_equals's default and how
# the q-scores are computed (probability_difference, risk_ratio, odds_ratio).

dat = make_incidence_data()
B = 30

test_that("default incidence_metric is odds_ratio with H_0_mu_equals defaulting to 1", {
	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "incidence", B = B, num_cores = 2)
	expect_valid_pte_result(res, B)
	expect_identical(res$incidence_metric, "odds_ratio")
	expect_identical(res$H_0_mu_equals, 1)
})

test_that("incidence_metric = 'risk_ratio' defaults H_0_mu_equals to 1", {
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", incidence_metric = "risk_ratio", B = B, num_cores = 2
	)
	expect_valid_pte_result(res, B)
	expect_identical(res$H_0_mu_equals, 1)
})

test_that("incidence_metric = 'probability_difference' defaults H_0_mu_equals to 0", {
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", incidence_metric = "probability_difference", B = B, num_cores = 2
	)
	expect_valid_pte_result(res, B)
	expect_identical(res$H_0_mu_equals, 0)
})

test_that("a strong true effect is detected for the default odds_ratio metric", {
	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "incidence", B = B, num_cores = 2)
	# odds ratios are positive by construction and a strong recommendation effect should push
	# the estimate comfortably above the null of 1
	expect_gt(res$est_q_average, 1)
	expect_lt(res$p_val_average, 0.3)
})

test_that("the fitted personalization_model is a glm with a binomial family", {
	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "incidence", B = 5, num_cores = 1)
	expect_s3_class(res$personalization_model, "glm")
	expect_identical(family(res$personalization_model)$family, "binomial")
})
