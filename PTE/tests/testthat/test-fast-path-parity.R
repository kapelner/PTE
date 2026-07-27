# Regression tests for the compiled fast paths (create_fast_lm/glm/weibull_objects and the
# fast_default_*_cv_run_cpp / fast_default_*_bootstrap_q_cpp batch functions) agreeing with the
# exact R-level reference computation (lm()/glm()/survreg() + predict()). Unlike the bootstrap
# loop, the *observed* pass (B bootstrap replicates aside) is fully deterministic -- it always
# uses boot_idx = 1:n with no sample() call -- so it's the only part of PTE_bootstrap_inference's
# output that can be compared for exact numeric equality across two separate calls. See
# helper-fixtures.R's header comment for why the bootstrap-sampled q_scores themselves can't be
# compared this way.
#
# These tests exist because the fast paths have twice silently diverged from the documented R
# defaults during development: once when IRLS/optimizer warm-starting chained each
# cross-validation fold's fit into the next fold's starting point (only caught by comparing the
# observed pass to glm()/survreg() on a dataset with a quasi-separated coefficient), and once when
# a hand-written Kaplan-Meier median replaced survival::survfit() but used continuous
# interpolation between step heights instead of survfit()'s actual (mostly step-function, with a
# specific tie-breaking rule) definition.

expect_observed_q_scores_match = function(res_fast, res_slow, tolerance = 1e-6){
	expect_equal(
		unlist(res_fast$observed_q_scores[c("adversarial", "average", "best")]),
		unlist(res_slow$observed_q_scores[c("adversarial", "average", "best")]),
		tolerance = tolerance
	)
}

test_that("fast and slow (explicit default-equivalent) paths agree exactly for continuous", {
	dat = make_continuous_data()
	custom_build = function(Xytrain) lm(y ~ . * treatment, data = Xytrain)
	custom_predict = function(mod, Xyleftout) predict(mod, Xyleftout)

	res_fast = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "continuous", B = 1, num_cores = 1)
	res_slow = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = 1, num_cores = 1,
		personalized_model_build_function = custom_build, predict_function = custom_predict
	)
	expect_observed_q_scores_match(res_fast, res_slow)
})

test_that("fast and slow (explicit default-equivalent) paths agree exactly for incidence", {
	dat = make_incidence_data()
	custom_build = function(Xytrain) glm(y ~ . * treatment, data = Xytrain, family = "binomial")
	custom_predict = function(mod, Xyleftout) predict(mod, Xyleftout)

	res_fast = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "incidence", B = 1, num_cores = 1)
	res_slow = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", B = 1, num_cores = 1,
		personalized_model_build_function = custom_build, predict_function = custom_predict
	)
	expect_observed_q_scores_match(res_fast, res_slow)
})

test_that("fast and slow (explicit default-equivalent) paths agree exactly for survival", {
	dat = make_survival_data()
	custom_build = function(Xytrain) survival::survreg(survival::Surv(Xytrain$y, Xytrain$censored) ~ (. - censored) * treatment, data = Xytrain, dist = "weibull")
	custom_predict = function(mod, Xyleftout) predict(mod, Xyleftout)

	res_fast = suppressWarnings(PTE_bootstrap_inference(dat$X, dat$y, censored = dat$censored, regression_type = "survival", B = 1, num_cores = 1))
	res_slow = suppressWarnings(PTE_bootstrap_inference(
		dat$X, dat$y, censored = dat$censored, regression_type = "survival", B = 1, num_cores = 1,
		personalized_model_build_function = custom_build, predict_function = custom_predict
	))
	expect_observed_q_scores_match(res_fast, res_slow)
})

make_separated_incidence_data = function(n = 100, seed = 7){
	set.seed(seed)
	treatment = sample(rep(0:1, each = n / 2))
	x1 = rnorm(n)
	# x2 almost perfectly predicts y, inducing (quasi-)separation on its coefficient -- this is
	# the failure mode that silently broke the fast incidence path during development (see the
	# file header comment)
	x2 = factor(rbinom(n, 1, 0.5))
	y = ifelse(x2 == 1, rbinom(n, 1, 0.98), rbinom(n, 1, 0.02))
	list(X = data.frame(treatment, x1, x2), y = y)
}

test_that("the fast incidence path does not diverge from glm() when a covariate is (quasi-)separated", {
	dat = make_separated_incidence_data()
	custom_build = function(Xytrain) glm(y ~ . * treatment, data = Xytrain, family = "binomial")
	custom_predict = function(mod, Xyleftout) predict(mod, Xyleftout)

	res_fast = suppressWarnings(PTE_bootstrap_inference(dat$X, dat$y, regression_type = "incidence", B = 1, num_cores = 1))
	res_slow = suppressWarnings(PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "incidence", B = 1, num_cores = 1,
		personalized_model_build_function = custom_build, predict_function = custom_predict
	))
	# a looser tolerance than the well-behaved-data tests above: near-separated MLEs are
	# inherently sensitive to solver/starting-point details, so the bar here is "doesn't diverge
	# wildly" (the historical bug produced differences of ~15-20 on this scale), not bit-exactness
	expect_observed_q_scores_match(res_fast, res_slow, tolerance = 1)
})

test_that("get_survival_stat_diff_cpp matches survival::survfit()'s Kaplan-Meier median difference", {
	km_diff_r = function(y, dead, w){
		mod = survival::survfit(survival::Surv(y, dead) ~ w)
		tab = summary(mod)$table
		if (is.null(dim(tab))) return(NA_real_)
		med = tab[, "median"]
		unname(med[2] - med[1])
	}

	set.seed(11)
	n_checked = 0
	for (i in 1:40){
		n = sample(10:60, 1)
		y = round(rexp(n, 0.3), 2)
		dead = rbinom(n, 1, 0.7)
		w = rbinom(n, 1, 0.5)
		if (sum(w == 0) == 0 || sum(w == 1) == 0) next
		n_checked = n_checked + 1

		expected = suppressWarnings(km_diff_r(y, dead, w))
		actual = PTE:::get_survival_stat_diff_cpp(y, dead, as.integer(w))
		if (is.na(expected)){
			expect_true(is.na(actual), info = paste("seed draw", i))
		} else {
			expect_equal(actual, expected, tolerance = 1e-6, info = paste("seed draw", i))
		}
	}
	expect_gt(n_checked, 20)
})
