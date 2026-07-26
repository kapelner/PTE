# All of these should error out during input validation, before the parallel bootstrap cluster
# is ever started, so these tests are fast.

base_X = data.frame(treatment = rep(0:1, each = 5), x = rnorm(10))
base_y = rnorm(10)

test_that("invalid regression_type errors", {
	expect_error(
		PTE_bootstrap_inference(base_X, base_y, regression_type = "bogus"),
		"regression_type"
	)
})

test_that("survival regression_type without censored errors", {
	expect_error(
		PTE_bootstrap_inference(base_X, base_y, regression_type = "survival"),
		"censored"
	)
})

test_that("invalid incidence_metric errors when no custom difference_function is given", {
	expect_error(
		PTE_bootstrap_inference(base_X, base_y, regression_type = "incidence", incidence_metric = "bogus"),
		"incidence_metric"
	)
})

test_that("invalid incidence_metric is tolerated when a custom difference_function is given", {
	# incidence_metric is documented as ignored once difference_function is supplied
	custom_diff = function(results, indices_1_1, indices_0_0, indices_0_1, indices_1_0) c(0, 0, 0)
	binary_y = rep(0:1, each = 5)
	expect_no_error(
		res <- PTE_bootstrap_inference(
			base_X, binary_y,
			regression_type = "incidence", incidence_metric = "bogus",
			difference_function = custom_diff,
			B = 5, num_cores = 1
		)
	)
	expect_s3_class(res, "PTE_bootstrap_results")
})

test_that("missing treatment column errors", {
	X_no_treatment = data.frame(x = rnorm(10))
	expect_error(
		PTE_bootstrap_inference(X_no_treatment, base_y),
		"indicator vector of the allocation"
	)
})

test_that("non-numeric treatment column errors", {
	X_factor_treatment = data.frame(treatment = factor(rep(c("A", "B"), each = 5)), x = rnorm(10))
	expect_error(
		PTE_bootstrap_inference(X_factor_treatment, base_y),
		"treatment"
	)
})

test_that("treatment column with values other than 0/1 errors", {
	X_bad_treatment = data.frame(treatment = rep(1:2, each = 5), x = rnorm(10))
	expect_error(
		PTE_bootstrap_inference(X_bad_treatment, base_y),
		"treatment"
	)
})

test_that("mismatched y length errors", {
	expect_error(
		PTE_bootstrap_inference(base_X, rnorm(9)),
		"response vector"
	)
})

test_that("mismatched censored length errors", {
	expect_error(
		PTE_bootstrap_inference(
			base_X, base_y,
			regression_type = "survival",
			censored = rep(0:1, each = 4) # length 8 != 10
		),
		"censored"
	)
})
