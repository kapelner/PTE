library(PTE)
options(mc.cores = 4)
set.seed(1984)
graphics.off()


sigsq_e = 1
beta_0 = 1
beta_1 = -1
beta_2 = 0.5
gamma_0 = 0
gamma_1 = 1
gamma_2 = -3

#we first need to find the true muIo*

ntest = 1e5
xtest = rnorm(ntest)
ztest = rnorm(ntest)
epsilons_test = rnorm(ntest, mean = 0, sd = sqrt(sigsq_e))

y_d_star = array(NA, ntest)
for (i in 1 : ntest){
	#evaluate d(x,z) - the real d, not the estimated d
	xtrain = xtest[i]
	z = ztest[i]
	eps = epsilons_test[i]
	if (gamma_0 + gamma_1 * xtrain^3 + gamma_2 * z > 0){
		#administer treatment
		y_d_star[i] = beta_0 + beta_1 * xtrain + beta_2 * z + gamma_0 + gamma_1 * xtrain^3 + gamma_2 * z + eps
	} else {
		#do not administer treatment
		y_d_star[i] = beta_0 + beta_1 * xtrain + beta_2 * z + eps
	}
	
}

par(mfrow = c(3,1))
hist(y_d_star, br = 1000, xlim = c(-5, 10))
abline(v = mean(y_d_star), col = "blue", lwd = 4)

xtest = rnorm(ntest)
ztest = rnorm(ntest)
epsilons_test = rnorm(ntest, mean = 0, sd = sqrt(sigsq_e))

y_d_0 = array(NA, ntest)
A_rand = rbinom(ntest, 1, 0.5)
sum(A_rand) / ntest
for (i in 1 : ntest){
	xtrain = xtest[i]
	z = ztest[i]
	eps = epsilons_test[i]
	#flip a coin
	if (A_rand[i] == 1){
		#administer treatment
		y_d_0[i] = beta_0 + beta_1 * xtrain + beta_2 * z + gamma_0 + gamma_1 * xtrain^3 + gamma_2 * z + eps
	} else {
		#do not administer treatment
		y_d_0[i] = beta_0 + beta_1 * xtrain + beta_2 * z + eps
	}
	
}
hist(y_d_0, br = 1000, xlim = c(-5, 10))
abline(v = mean(y_d_0), col = "blue", lwd = 4)

boxplot(y_d_star, y_d_0)
t.test(y_d_star, y_d_0)

muI0star = mean(y_d_star) - mean(y_d_0)



#we first need to find the true muIo*

ntest = 1e5
xtest = rnorm(ntest)
ztest = rnorm(ntest)
treatment_test = rbinom(ntest, 1, 0.5)
epsilons_test = rnorm(ntest, mean = 0, sd = sqrt(sigsq_e))

#We need to fish up the best guesses now as to all of the coefficients
#so generate the responses under the true DGP
ytest = beta_0 +
		beta_1 * xtest +
		beta_2 * ztest + 
		treatment_test * (
			gamma_0 + 
			gamma_1 * xtest^3 +
			gamma_2 * ztest
			) +
		epsilons_test

#now regress on what our f is

mod = lm(ytest ~ xtest * treatment_test)
gamma_0_hat = coef(mod)[3]
gamma_1_hat = coef(mod)[4]

#now simulate our E[Y] under d over X

xtest = rnorm(ntest)
ztest = rnorm(ntest)
epsilons_test = rnorm(ntest, mean = 0, sd = sqrt(sigsq_e))

y_d = array(NA, ntest)
for (i in 1 : ntest){
	#evaluate d(x,z) - the real d, not the estimated d
	xtrain = xtest[i]
	z = ztest[i]
	eps = epsilons_test[i]
	if (gamma_0_hat + gamma_1_hat * xtrain > 0){
		#administer treatment
		y_d[i] = beta_0 + beta_1 * xtrain + beta_2 * z + gamma_0 + gamma_1 * xtrain^3 + gamma_2 * z + eps
	} else {
		#do not administer treatment
		y_d[i] = beta_0 + beta_1 * xtrain + beta_2 * z + eps
	}
	
}

par(mfrow = c(3,1))
hist(y_d, br = 1000, xlim = c(-5, 10))
abline(v = mean(y_d), col = "blue", lwd = 4)

xtest = rnorm(ntest)
ztest = rnorm(ntest)
epsilons_test = rnorm(ntest, mean = 0, sd = sqrt(sigsq_e))

y_d_0 = array(NA, ntest)
A_rand = rbinom(ntest, 1, 0.5)
sum(A_rand) / ntest
for (i in 1 : ntest){
	xtrain = xtest[i]
	z = ztest[i]
	eps = epsilons_test[i]
	#flip a coin
	if (A_rand[i] == 1){
		#administer treatment
		y_d_0[i] = beta_0 + beta_1 * xtrain + beta_2 * z + gamma_0 + gamma_1 * xtrain^3 + gamma_2 * z + eps
	} else {
		#do not administer treatment
		y_d_0[i] = beta_0 + beta_1 * xtrain + beta_2 * z + eps
	}
	
}
hist(y_d_0, br = 1000, xlim = c(-5, 10))
abline(v = mean(y_d_0), col = "blue", lwd = 4)

boxplot(y_d, y_d_0)
t.test(y_d, y_d_0)

muI0 = mean(y_d) - mean(y_d_0)



muI0star #This is the 1.65 on p19 of the paper
muI0 #This is the 0.79 on p19 of the paper
#To find the 85% and 15% breakdown, you can alter the code above and recalculate muI0

#now we can go on simulating

ntrains = c(100, 200, 500, 1000)
num_boot = 3000


graphics.off()
par(mfrow = c(length(ntrains), 2), mar = c(3, 0.5, 0.5, 0.5))


set.seed(1984)

#Figure 3
for (i_n in 1 : length(ntrains)){
	ntrain = ntrains[i_n]
	
	xtrain = rnorm(ntrain)
	ztrain = rnorm(ntrain)
	treatment = sample(c(rep(1, ntrain / 2), rep(0, ntrain / 2)))
	epsilon_trains = rnorm(ntrain, mean = 0, sd = sqrt(sigsq_e))
	
	
	ytrain = beta_0 +
			beta_1 * xtrain +
			beta_2 * ztrain + 
			treatment * (
				gamma_0 + 
				gamma_1 * xtrain^3 +
				gamma_2 * ztrain
				) +
			epsilon_trains
	
	Xtrain = data.frame(treatment, xtrain)
	
	res = PTE_bootstrap_inference(Xtrain, ytrain, B = num_boot)
	print(i_n)
	print(res)
	
	min_q = -0.2#min(res$q_scores$average)
	max_q = 1.67#max(res$q_scores$average)
	
	hist(res$q_scores$average, br = num_boot / 3, xlab = "", xlim = c(min_q, max_q), main = "", yaxt = "n", xaxt = "n", ylab = "", col = rgb(0.9, 0.9, 0.9), border = rgb(0.8, 0.8, 0.8))  
	axis(1, at = c(0, 0.5, muI0, 1, 1.5), labels = c(0, 0.5, 0.79, 1, 1.5), cex.axis = 1.3)
	abline(v = res$observed_q_average, col = "black", lwd = 4)
	abline(v = res$ci_q_average[1], col = "black", lwd = 1)
	abline(v = res$ci_q_average[2], col = "black", lwd = 1)
	abline(v = muI0star, col = "black", lwd = 1, lty = 2)
	hist(res$q_scores$best, br = num_boot / 3, xlab = "", xlim = c(min_q, max_q), main = "", yaxt = "n", xaxt = "n", ylab = "", col = rgb(0.9, 0.9, 0.9), border = rgb(0.8, 0.8, 0.8)) # 
	axis(1, at = c(0, 0.5, muI0, 1, 1.5), labels = c(0, 0.5, 0.79, 1, 1.5), cex.axis = 1.3)
	abline(v = res$observed_q_best, col = "black", lwd = 4)
	abline(v = res$ci_q_best[1], col = "black", lwd = 1)
	abline(v = res$ci_q_best[2], col = "black", lwd = 1)
	abline(v = muI0star, col = "black", lwd = 1, lty = 2)
}
