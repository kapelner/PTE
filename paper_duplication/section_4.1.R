library(PTE)

set.seed(1984)


ns = c(100, 200, 500, 1000)

beta0 = 1
beta1 = -1
gamma0 = 0
gamma1 = sqrt(2 * pi)
mu_x = 0
sigsq_x = 1
sigsq_e = 1
num_boot = 3000


graphics.off()
par(mfrow = c(length(ns), 2), mar = c(3, 0.5, 0.5, 0.5))

#Figure 2
for (i_n in 1 : length(ns)){
	n = ns[i_n]
	
	x = sort(rnorm(n, mu_x, sigsq_x))
	noise = rnorm(n, 0, sigsq_e)
	
	treatment = sample(c(rep(1, n / 2), rep(0, n / 2)))
	y = beta0 + beta1 * x + treatment * (gamma0 + gamma1 * x) + noise
	
	X = data.frame(treatment, x)
	
	res = bootstrap_inference(X, y,
			"lm(y ~ . + treatment * ., data = Xyleft)",
			num_cores = 4,
			B = num_boot, 
			plot = FALSE)
#	print(n)
#	print(res)
	
	min_q = 0.25
	max_q = 1.5

	hist(res$q_scores$average, br = num_boot / 50, xlab = "", xlim = c(min_q, max_q), main = "", yaxt = "n", xaxt = "n", ylab = "", col = rgb(0.9, 0.9, 0.9), border = rgb(0.8, 0.8, 0.8))
	axis(1, at = c(0.5, 1, 1.5), labels = c(0.5, 1, 1.5), cex.axis = 1.3)
	abline(v = res$observed_q_average, col = "black", lwd = 4)
	abline(v = res$ci_q_average[1], col = "black", lwd = 1)
	abline(v = res$ci_q_average[2], col = "black", lwd = 1)
	hist(res$q_scores$best, br = num_boot / 50, xlab = "", xlim = c(min_q, max_q), main = "", yaxt = "n", xaxt = "n", ylab = "", col = rgb(0.9, 0.9, 0.9), border = rgb(0.8, 0.8, 0.8))
	axis(1, at = c(0.5, 1, 1.5), labels = c(0.5, 1, 1.5), cex.axis = 1.3)
	abline(v = res$observed_q_best, col = "black", lwd = 4)
	abline(v = res$ci_q_best[1], col = "black", lwd = 1)
	abline(v = res$ci_q_best[2], col = "black", lwd = 1)
}

