dat = make_continuous_data()

test_that("returns a valid PTE_optimal_m_selection object with a sensible grid_table", {
	grid = seq(0.6, 1, by = 0.1)
	m_selection = select_optimal_m_prop(
		dat$X, dat$y, regression_type = "continuous",
		m_pow_of_n_grid = grid, B = 15, num_cores = 1
	)

	expect_s3_class(m_selection, "PTE_optimal_m_selection")
	expect_true(m_selection$m_pow_of_n_optimal %in% grid)
	expect_identical(m_selection$m_optimal, round(nrow(dat$X) ^ m_selection$m_pow_of_n_optimal))
	expect_identical(m_selection$estimand, "average")
	expect_identical(m_selection$B, 15)

	gt = m_selection$grid_table
	expect_s3_class(gt, "data.frame")
	expect_identical(nrow(gt), length(grid))
	expect_identical(gt$m_pow_of_n, sort(grid))
	expect_identical(gt$m, round(nrow(dat$X) ^ sort(grid)))

	#volatility_window = 3 (the default) means the first and last grid points can never have a
	#complete centered neighborhood, so they must be ineligible (NA volatility) and can never be
	#selected as optimal
	expect_true(is.na(gt$volatility[1]))
	expect_true(is.na(gt$volatility[nrow(gt)]))
	expect_false(m_selection$m_pow_of_n_optimal %in% c(min(grid), max(grid)))
})

test_that("an unsorted, duplicate-free grid is sorted internally and still selects an interior point", {
	grid = c(1, 0.6, 0.8, 0.7, 0.9)
	m_selection = select_optimal_m_prop(
		dat$X, dat$y, regression_type = "continuous",
		m_pow_of_n_grid = grid, B = 15, num_cores = 1
	)
	expect_identical(m_selection$grid_table$m_pow_of_n, sort(grid))
	expect_false(m_selection$m_pow_of_n_optimal %in% c(min(grid), max(grid)))
})

test_that("estimand can be switched to adversarial or best", {
	m_selection = select_optimal_m_prop(
		dat$X, dat$y, regression_type = "continuous",
		m_pow_of_n_grid = seq(0.6, 1, by = 0.1), B = 15, num_cores = 1,
		estimand = "best"
	)
	expect_identical(m_selection$estimand, "best")
})

test_that("invalid regression_type errors", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, regression_type = "bogus"),
		"regression_type"
	)
})

test_that("m_pow_of_n_grid with fewer than 3 points errors", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, m_pow_of_n_grid = c(0.7, 0.9)),
		"m_pow_of_n_grid"
	)
})

test_that("m_pow_of_n_grid values outside [0, 1] error", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, m_pow_of_n_grid = c(0.5, 0.8, 1.2)),
		"m_pow_of_n_grid"
	)
})

test_that("duplicate m_pow_of_n_grid values error", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, m_pow_of_n_grid = c(0.5, 0.8, 0.8, 1)),
		"m_pow_of_n_grid"
	)
})

test_that("invalid estimand errors", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, estimand = "bogus"),
		"estimand"
	)
})

test_that("an even volatility_window errors", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, volatility_window = 4),
		"volatility_window"
	)
})

test_that("a volatility_window below 3 errors", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, volatility_window = 1),
		"volatility_window"
	)
})

test_that("a volatility_window larger than the grid errors", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, m_pow_of_n_grid = seq(0.6, 1, by = 0.1), volatility_window = 7),
		"volatility_window"
	)
})

test_that("a non-positive B errors", {
	expect_error(
		select_optimal_m_prop(dat$X, dat$y, B = 0),
		"B"
	)
})

test_that("an all-bad grid (every eligible point NA) errors with a clear message instead of which.min's opaque error", {
	testthat::local_mocked_bindings(
		PTE_bootstrap_inference = function(...) list(ci_q_average = c(NA_real_, NA_real_), est_q_average = NA_real_),
		.package = "PTE"
	)
	expect_error(
		select_optimal_m_prop(
			dat$X, dat$y, regression_type = "continuous",
			m_pow_of_n_grid = seq(0.6, 1, by = 0.1), B = 5, num_cores = 1
		),
		"no m_pow_of_n could be selected"
	)
})
