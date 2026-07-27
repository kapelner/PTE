test_that("v2.0 default results are statistically similar to v1.7", {
	skip_if_not(
		identical(Sys.getenv("PTE_RUN_V17_COMPARISON"), "true"),
		"Set PTE_RUN_V17_COMPARISON=true to run the v1.7 statistical comparison integration test."
	)
	skip_if_not(nzchar(Sys.which("git")), "git is required to materialize the v1.7 reference package.")
	repo_root = system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE)
	skip_if_not(length(repo_root) == 1 && dir.exists(file.path(repo_root, ".git")), "the v1.7 comparison must be run from inside the git repository.")

	v17_lib = tempfile("pte-v17-lib-")
	v17_src = tempfile("pte-v17-src-")
	dir.create(v17_lib)
	dir.create(v17_src)

	tarfile = tempfile(fileext = ".tar")
	archive_status = system2("git", c("-C", repo_root, "archive", "--format=tar", paste0("--output=", tarfile), "503f350"), stdout = TRUE, stderr = TRUE)
	skip_if_not(identical(archive_status, character(0)), paste("could not archive v1.7 reference:", paste(archive_status, collapse = "\n")))
	untar(tarfile, exdir = v17_src)
	description_files = list.files(v17_src, pattern = "^DESCRIPTION$", recursive = TRUE, full.names = TRUE)
	skip_if_not(length(description_files) >= 1, "could not locate v1.7 package DESCRIPTION after extracting git archive.")
	v17_pkg_dir = dirname(description_files[[1]])

	install_output = suppressWarnings(system2(
		file.path(R.home("bin"), "R"),
		c("CMD", "INSTALL", "-l", v17_lib, v17_pkg_dir),
		stdout = TRUE,
		stderr = TRUE
	))
	install_status = attr(install_output, "status")
	skip_if(is.numeric(install_status) && install_status != 0, paste("could not install v1.7 reference:", paste(install_output, collapse = "\n")))

	run_reference_case = function(lib_path, scenario, out_file, B, patch_parallel = FALSE){
		script_file = tempfile(fileext = ".R")
		parallel_patch_lines = if (patch_parallel) {
			c(
				"pte_ns = asNamespace('PTE')",
				"patch_pte_import = function(name, value){",
				"	for (env in list(pte_ns, parent.env(pte_ns))){",
				"		if (exists(name, envir = env, inherits = FALSE)){",
				"			if (bindingIsLocked(name, env)) unlockBinding(name, env)",
				"			assign(name, value, envir = env)",
				"			lockBinding(name, env)",
				"			return(invisible(NULL))",
				"		}",
				"	}",
				"	stop('could not patch PTE import: ', name)",
				"}",
				"patch_pte_import('makeCluster', function(...) NULL)",
				"patch_pte_import('registerDoParallel', function(...) foreach::registerDoSEQ())",
				"patch_pte_import('stopCluster', function(...) invisible(NULL))"
			)
		} else {
			character(0)
		}
		writeLines(c(
			sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib_path)),
			"suppressPackageStartupMessages(library(PTE))",
			parallel_patch_lines,
			"make_continuous_data = function(n = 80, seed = 1984, gamma1 = sqrt(2 * pi)){",
			"	set.seed(seed)",
			"	x = sort(rnorm(n))",
			"	noise = rnorm(n)",
			"	treatment = sample(rep(0:1, each = n / 2))",
			"	y = 1 - x + treatment * (gamma1 * x) + noise",
			"	list(X = data.frame(treatment, x), y = y)",
			"}",
			"make_incidence_data = function(n = 100, seed = 1984, gamma1 = 3){",
			"	set.seed(seed)",
			"	x = sort(rnorm(n))",
			"	treatment = sample(rep(0:1, each = n / 2))",
			"	lin = -x + treatment * (gamma1 * x)",
			"	p = 1 / (1 + exp(-lin))",
			"	y = rbinom(n, 1, p)",
			"	list(X = data.frame(treatment, x), y = y)",
			"}",
			"make_survival_data = function(n = 150, seed = 1984, gamma1 = 1.5){",
			"	set.seed(seed)",
			"	x = sort(rnorm(n))",
			"	treatment = sample(rep(0:1, each = n / 2))",
			"	lin = 1 - 0.5 * x + treatment * (gamma1 * x)",
			"	y = rexp(n, rate = exp(-lin))",
			"	censored = rbinom(n, 1, 0.8)",
			"	list(X = data.frame(treatment, x), y = y, censored = censored)",
			"}",
			sprintf("scenario = %s", deparse(scenario)),
			sprintf("B = %d", B),
			"set.seed(20260726)",
			"if (scenario == 'continuous'){",
			"	dat = make_continuous_data()",
			"	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = 'continuous', B = B, num_cores = 1)",
			"} else if (scenario == 'incidence'){",
			"	dat = make_incidence_data()",
			"	res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = 'incidence', B = B, num_cores = 1)",
			"} else {",
			"	dat = make_survival_data()",
			"	res = PTE_bootstrap_inference(dat$X, dat$y, censored = dat$censored, regression_type = 'survival', B = B, num_cores = 1)",
			"}",
			"payload = list(",
			"	version = as.character(utils::packageVersion('PTE')),",
			"	observed_q_scores = unlist(res$observed_q_scores[c('adversarial', 'average', 'best')]),",
			"	q_scores = res$q_scores[c('adversarial', 'average', 'best')],",
			"	est = c(adversarial = res$est_q_adversarial, average = res$est_q_average, best = res$est_q_best),",
			"	ci = list(adversarial = res$ci_q_adversarial, average = res$ci_q_average, best = res$ci_q_best),",
			"	num_bad = res$num_bad",
			")",
			sprintf("saveRDS(payload, %s)", deparse(out_file))
		), script_file)
		output = suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script_file, stdout = TRUE, stderr = TRUE))
		status = attr(output, "status")
		if (is.numeric(status) && status != 0){
			stop(paste(output, collapse = "\n"), call. = FALSE)
		}
		readRDS(out_file)
	}

	assert_statistically_similar = function(current, reference, scenario){
		expect_named(current$q_scores, c("adversarial", "average", "best"))
		expect_named(reference$q_scores, c("adversarial", "average", "best"))

		observed_tol = switch(scenario, continuous = 1e-7, incidence = 1e-6, survival = 1e-7)
		expect_equal(unname(current$observed_q_scores), unname(reference$observed_q_scores), tolerance = observed_tol)

		for (score_name in c("adversarial", "average", "best")){
			x = current$q_scores[[score_name]]
			y = reference$q_scores[[score_name]]
			x = x[is.finite(x)]
			y = y[is.finite(y)]
			expect_true(length(x) > 0, info = paste(scenario, score_name, "current scores are all non-finite"))
			expect_true(length(y) > 0, info = paste(scenario, score_name, "reference scores are all non-finite"))

			pooled_se = sqrt(stats::var(x) / length(x) + stats::var(y) / length(y))
			mean_tol = max(0.05 * max(1, abs(mean(y))), 4 * pooled_se)
			expect_true(abs(mean(x) - mean(y)) <= mean_tol, info = paste(scenario, score_name, "bootstrap means differ"))

			current_ci = current$ci[[score_name]]
			reference_ci = reference$ci[[score_name]]
			expect_true(
				current_ci[1] <= reference_ci[2] && reference_ci[1] <= current_ci[2],
				info = paste(scenario, score_name, "percentile CIs do not overlap")
			)
		}
	}

	current_lib = dirname(system.file(package = "PTE"))
	scenario_B = c(continuous = 250, incidence = 250, survival = 150)
	for (scenario in names(scenario_B)){
		current_out = tempfile(fileext = ".rds")
		reference_out = tempfile(fileext = ".rds")
		current = run_reference_case(current_lib, scenario, current_out, scenario_B[[scenario]])
		reference = run_reference_case(v17_lib, scenario, reference_out, scenario_B[[scenario]], patch_parallel = TRUE)

		expect_equal(current$version, as.character(utils::packageVersion("PTE")))
		expect_equal(reference$version, "1.7")
		assert_statistically_similar(current, reference, scenario)
	}
})
