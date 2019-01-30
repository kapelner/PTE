set.seed(1984)

library(PTE)
options(mc.cores = 4)
Xy = read.csv("cpt_original.csv") #this dataset is not released due to potentially identifying information
y = Xy$wk15_HRSD
X = Xy[, c("treatment", "age", "chronic", "life_stressors", "pdstatus", "unemployed", "married")]
#recode a couple of variables
X$pdstatus = ifelse(X$pdstatus == 0.5, 1, 0)
X$treatment = ifelse(X$treatment == 1, 1, 0)
head(X)

# s = Sys.time()
res = PTE_bootstrap_inference(X, y, y_higher_is_better = FALSE, B = 3000, run_bca_bootstrap = TRUE)
# print(Sys.time() - s)
res


# par(mar = c(2, 0, 0.5, 0))
# B = 3000
# min_q = min(res$q_scores$average, res$q_scores$best)
# max_q = max(res$q_scores$average, res$q_scores$best)
# 
# #Figure 6a
# hist(res$q_scores$average, br = B / 25, xlab = "", xlim = c(min_q, max_q), main = "", ylab = "", yaxt = "n", col = rgb(0.9, 0.9, 0.9), border = rgb(0.8, 0.8, 0.8)) #, main = "I_Avg's"
# abline(v = res$est_q_average, col = "black", lwd = 4)
# abline(v = res$ci_q_average[1], col = "black", lwd = 1)
# abline(v = res$ci_q_average[2], col = "black", lwd = 1)
# 
# #Figure 6b
# hist(res$q_scores$best, br = B / 25, xlab = "", xlim = c(min_q, max_q), main = "", ylab = "", yaxt = "n", col = rgb(0.9, 0.9, 0.9), border = rgb(0.8, 0.8, 0.8)) #, main = "I_Best's"
# abline(v = res$est_q_best, col = "black", lwd = 4)
# abline(v = res$ci_q_best[1], col = "black", lwd = 1)
# abline(v = res$ci_q_best[2], col = "black", lwd = 1)
