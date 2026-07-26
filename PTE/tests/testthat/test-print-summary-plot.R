# S3 methods for PTE_bootstrap_results: print(), summary(), plot().

dat = make_continuous_data(n = 60)
res = PTE_bootstrap_inference(dat$X, dat$y, regression_type = "continuous", B = 15, num_cores = 2)

test_that("print.PTE_bootstrap_results reports I_random and I_best but not I_adversarial by default", {
	out = capture.output(print(res))
	txt = paste(out, collapse = "\n")
	expect_match(txt, "I_random")
	expect_match(txt, "I_best")
	expect_false(grepl("I_adversarial", txt))
})

test_that("print.PTE_bootstrap_results reports I_adversarial when display_adversarial_score is TRUE", {
	res_adv = res
	res_adv$display_adversarial_score = TRUE
	out = capture.output(print(res_adv))
	txt = paste(out, collapse = "\n")
	expect_match(txt, "I_adversarial")
})

test_that("summary.PTE_bootstrap_results produces the same output as print", {
	print_out = capture.output(print(res))
	summary_out = capture.output(summary(res))
	expect_identical(print_out, summary_out)
})

test_that("plot.PTE_bootstrap_results runs without error and returns invisibly", {
	grDevices::pdf(grDevices::nullfile())
	on.exit(grDevices::dev.off(), add = TRUE)

	expect_no_error(return_value <- withVisible(plot(res)))
	expect_false(return_value$visible)
})
