# regression_type = "incidence": y is binary, default model is glm(y ~ . * treatment, family
# = "binomial"). Three incidence_metric options change both H_0_mu_equals's default and how
# the q-scores are computed (probability_difference, risk_ratio, odds_ratio).

dat = make_incidence_data()
B = 30

test_that("default incidence_metric is odds_ratio with H_0_mu_equals defaulting to 1", {
	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "incidence", B = B, num_cores = 1)
	expect_valid_pte_result(res, B)
	expect_identical(res$incidence_metric, "odds_ratio")
	expect_identical(res$H_0_mu_equals, 1)
})

test_that("incidence_metric = 'risk_ratio' defaults H_0_mu_equals to 1", {
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", incidence_metric = "risk_ratio", B = B, num_cores = 1
	)
	expect_valid_pte_result(res, B)
	expect_identical(res$H_0_mu_equals, 1)
})

test_that("incidence_metric = 'probability_difference' defaults H_0_mu_equals to 0", {
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", incidence_metric = "probability_difference", B = B, num_cores = 1
	)
	expect_valid_pte_result(res, B)
	expect_identical(res$H_0_mu_equals, 0)
})

test_that("a strong true effect is detected for the default odds_ratio metric", {
	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "incidence", B = B, num_cores = 1)
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

test_that("probability_difference metric detects a strong true effect with a positive difference", {
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", incidence_metric = "probability_difference", B = B, num_cores = 1
	)
	expect_gt(res$est_q_average, 0)
	expect_lt(res$p_val_average, 0.3)
})

test_that("a custom difference_function is honored for incidence data", {
	custom_diff = function(results, indices_1_1, indices_0_0, indices_0_1, indices_1_0){
		y_all = results$real_y
		rec = c(y_all[indices_1_1], y_all[indices_0_0])
		non_rec = c(y_all[indices_0_1], y_all[indices_1_0])
		c(mean(rec) - mean(non_rec), mean(rec) - mean(y_all), mean(rec) - max(mean(y_all[results$given_tx == 0]), mean(y_all[results$given_tx == 1])))
	}

	res_diff = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", B = B, num_cores = 1,
		difference_function = custom_diff
	)
	expect_valid_pte_result(res_diff, B)
	expect_gt(res_diff$est_q_average, 0)
})

test_that("y_higher_is_better = FALSE produces a valid, differently-directed result", {
	res_default = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "incidence", B = B, num_cores = 1)
	res_flip = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", B = B, num_cores = 1, y_higher_is_better = FALSE
	)
	expect_valid_pte_result(res_flip, B)
	# flipping which outcome is "better" reverses which treatment is recommended, which should
	# pull the odds ratio estimate to the opposite side of 1 relative to the default direction
	expect_true((res_default$est_q_average - 1) * (res_flip$est_q_average - 1) <= 0)
})
