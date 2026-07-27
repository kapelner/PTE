# Covers PTE/R/zzz.R's package-load-time mc.cores default and the OS-conditional
# cluster handling (FORK on Unix, mirai on Windows) in PTE_bootstrap_inference --
# see docs/superpowers/specs/2026-07-26-mirai-fork-parallelization-design.md.

can_create_fork_cluster = function(){
	if (.Platform$OS.type == "windows"){
		return(FALSE)
	}
	cluster = NULL
	ok = tryCatch({
		cluster <<- parallel::makeCluster(1, type = "FORK")
		TRUE
	}, error = function(e) FALSE)
	if (!is.null(cluster)){
		parallel::stopCluster(cluster)
	}
	ok
}

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
	old = getOption("mc.cores")
	on.exit(options(mc.cores = old), add = TRUE)

	options(mc.cores = 1)
	dat = make_continuous_data()
	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "continuous", B = 5)
	expect_valid_pte_result(res, 5)
})

test_that("cluster connections are cleaned up after a call without BCA", {
	skip_if_not(can_create_fork_cluster(), "FORK clusters are unavailable in this check environment")

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
	skip_if_not(can_create_fork_cluster(), "FORK clusters are unavailable in this check environment")

	baseline = nrow(showConnections(all = TRUE))

	dat = make_continuous_data(n = 40)
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = 5, num_cores = 2,
		run_bca_bootstrap = TRUE
	)

	expect_valid_pte_result(res, 5, expect_bca = TRUE)
	expect_equal(nrow(showConnections(all = TRUE)), baseline)
})

test_that("serial operation uses a sequential backend with no cluster overhead", {
	baseline = nrow(showConnections(all = TRUE))

	dat = make_continuous_data()
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = 5, num_cores = 1
	)

	expect_valid_pte_result(res, 5)
	expect_identical(foreach::getDoParName(), "doSEQ")
	expect_equal(nrow(showConnections(all = TRUE)), baseline)
})

test_that("backend initialization uses sequential operation for one core", {
	register_seq_called = FALSE
	register_parallel_called = FALSE

	backend = PTE:::initialize_pte_bootstrap_backend(
		num_cores = 1,
		register_seq = function(){
			register_seq_called <<- TRUE
		},
		register_parallel = function(cluster){
			register_parallel_called <<- TRUE
		}
	)

	expect_null(backend$cluster)
	expect_identical(backend$backend, "sequential")
	expect_true(register_seq_called)
	expect_false(register_parallel_called)
})

test_that("backend initialization uses FORK clusters on non-Windows platforms", {
	make_cluster_args = NULL
	registered_cluster = NULL
	fake_cluster = structure(list(), class = "cluster")

	backend = PTE:::initialize_pte_bootstrap_backend(
		num_cores = 2,
		os_type = "unix",
		make_cluster = function(...){
			make_cluster_args <<- list(...)
			fake_cluster
		},
		register_parallel = function(cluster){
			registered_cluster <<- cluster
		}
	)

	expect_identical(backend$cluster, fake_cluster)
	expect_identical(backend$backend, "fork")
	expect_identical(make_cluster_args[[1]], 2)
	expect_identical(make_cluster_args$type, "FORK")
	expect_identical(registered_cluster, fake_cluster)
})

test_that("forking operation completes on Unix-like platforms", {
	skip_on_os("windows")
	skip_if_not(can_create_fork_cluster(), "FORK clusters are unavailable in this check environment")

	dat = make_continuous_data()
	res = PTE_bootstrap_inference(
		dat$X, dat$y, regression_type = "continuous", B = 5, num_cores = 2
	)

	expect_valid_pte_result(res, 5)
})

test_that("backend initialization uses mirai clusters on Windows", {
	require_namespace_args = NULL
	make_mirai_cluster_arg = NULL
	registered_cluster = NULL
	fake_cluster = structure(list(), class = "miraiCluster")

	backend = PTE:::initialize_pte_bootstrap_backend(
		num_cores = 2,
		os_type = "windows",
		make_mirai_cluster = function(num_cores){
			make_mirai_cluster_arg <<- num_cores
			fake_cluster
		},
		register_parallel = function(cluster){
			registered_cluster <<- cluster
		},
		require_namespace = function(package, quietly){
			require_namespace_args <<- list(package = package, quietly = quietly)
			TRUE
		}
	)

	expect_identical(backend$cluster, fake_cluster)
	expect_identical(backend$backend, "mirai")
	expect_identical(require_namespace_args, list(package = "mirai", quietly = TRUE))
	expect_identical(make_mirai_cluster_arg, 2)
	expect_identical(registered_cluster, fake_cluster)
})

test_that("mirai operation reports the missing optional dependency clearly", {
	expect_error(
		PTE:::initialize_pte_bootstrap_backend(
			num_cores = 2,
			os_type = "windows",
			require_namespace = function(package, quietly) FALSE
		),
		"mirai.*required"
	)
})
