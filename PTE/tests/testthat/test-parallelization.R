# Covers PTE/R/zzz.R's package-load-time mc.cores default and the OS-conditional
# cluster handling (FORK on Unix, mirai on Windows) in PTE_bootstrap_inference --
# see docs/superpowers/specs/2026-07-26-mirai-fork-parallelization-design.md.

test_that(".onLoad sets mc.cores to max(detectCores() - 1, 1) when unset", {
	old = getOption("mc.cores")
	on.exit(options(mc.cores = old), add = TRUE)

	options(mc.cores = NULL)
	PTE:::.onLoad(NULL, "PTE")

	expect_identical(getOption("mc.cores"), max(parallel::detectCores() - 1, 1))
})

test_that(".onLoad does not clobber a pre-existing mc.cores option", {
	old = getOption("mc.cores")
	on.exit(options(mc.cores = old), add = TRUE)

	options(mc.cores = 999L)
	PTE:::.onLoad(NULL, "PTE")

	expect_identical(getOption("mc.cores"), 999L)
})

test_that("num_cores defaults to getOption('mc.cores') when left NULL", {
	dat = make_continuous_data()
	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "continuous", B = 5)
	expect_valid_pte_result(res, 5)
})

test_that("cluster connections are cleaned up after a call without BCA", {
	baseline = nrow(showConnections(all = TRUE))

	dat = make_continuous_data()
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = 5, num_cores = 2,
		run_bca_bootstrap = FALSE
	)

	expect_valid_pte_result(res, 5)
	expect_equal(nrow(showConnections(all = TRUE)), baseline)
})

test_that("cluster connections are cleaned up after a call with BCA", {
	baseline = nrow(showConnections(all = TRUE))

	dat = make_continuous_data(n = 40)
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = 5, num_cores = 2,
		run_bca_bootstrap = TRUE
	)

	expect_valid_pte_result(res, 5, expect_bca = TRUE)
	expect_equal(nrow(showConnections(all = TRUE)), baseline)
})

test_that("num_cores = 1 uses a sequential backend with no cluster overhead", {
	baseline = nrow(showConnections(all = TRUE))

	dat = make_continuous_data()
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = 5, num_cores = 1
	)

	expect_valid_pte_result(res, 5)
	expect_identical(foreach::getDoParName(), "doSEQ")
	expect_equal(nrow(showConnections(all = TRUE)), baseline)
})
