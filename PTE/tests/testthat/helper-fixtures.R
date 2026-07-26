# Shared, deterministic-data fixtures for the PTE test suite.
#
# The bootstrap resampling and cross-validation folds inside PTE_bootstrap_inference are NOT
# reproducible via set.seed() -- each parallel worker gets its own fresh RNG stream, so exact
# numeric values differ run to run even with the same input data. To get tests that are both
# fast and non-flaky, these fixtures use a single strongly-signalled covariate (mirroring the
# simulation design in paper_duplication/section_4.1.R) so that the *direction* and rough
# *magnitude* of the estimated effect are stable across runs, while B is kept small since nearly
# all wall-clock time is parallel cluster startup, not iteration count.

make_continuous_data = function(n = 80, seed = 1984, gamma1 = sqrt(2 * pi)){
	set.seed(seed)
	x = sort(rnorm(n))
	noise = rnorm(n)
	treatment = sample(rep(0:1, each = n / 2))
	y = 1 - x + treatment * (gamma1 * x) + noise
	list(X = data.frame(treatment, x), y = y)
}

make_incidence_data = function(n = 100, seed = 1984, gamma1 = 3){
	set.seed(seed)
	x = sort(rnorm(n))
	treatment = sample(rep(0:1, each = n / 2))
	lin = -x + treatment * (gamma1 * x)
	p = 1 / (1 + exp(-lin))
	y = rbinom(n, 1, p)
	list(X = data.frame(treatment, x), y = y)
}

make_survival_data = function(n = 150, seed = 1984, gamma1 = 1.5){
	set.seed(seed)
	x = sort(rnorm(n))
	treatment = sample(rep(0:1, each = n / 2))
	lin = 1 - 0.5 * x + treatment * (gamma1 * x)
	y = rexp(n, rate = exp(-lin))
	censored = rbinom(n, 1, 0.8)
	list(X = data.frame(treatment, x), y = y, censored = censored)
}

# Structural/contract assertions common to every PTE_bootstrap_results object, regardless of
# regression_type or which optional features (BCA, custom functions) were used.
expect_valid_pte_result = function(res, B, expect_bca = FALSE){
	expect_s3_class(res, "PTE_bootstrap_results")

	for (nm in c("adversarial", "average", "best")){
		expect_length(res$q_scores[[nm]], B)
	}

	expect_true(is.numeric(res$p_val_average) && res$p_val_average >= 0 && res$p_val_average <= 1)
	expect_true(is.numeric(res$p_val_best) && res$p_val_best >= 0 && res$p_val_best <= 1)

	expect_true(res$ci_q_average[1] <= res$ci_q_average[2])
	expect_true(res$ci_q_best[1] <= res$ci_q_best[2])
	expect_true(res$basic_ci_q_average[1] <= res$basic_ci_q_average[2])
	expect_true(res$basic_ci_q_best[1] <= res$basic_ci_q_best[2])

	if (expect_bca){
		expect_true(res$bca_ci_q_average[1] <= res$bca_ci_q_average[2])
		expect_true(res$bca_ci_q_best[1] <= res$bca_ci_q_best[2])
	}

	expect_identical(res$B, B)
	expect_true(is.numeric(res$num_bad) && res$num_bad >= 0)
}
