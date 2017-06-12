######## SURVIVAL
rm(list = ls())


library(PTE)
data(survival)
X$site = NULL
pte_results = PTE_bootstrap_inference(X, y, censored = censored, regression_type = "survival", B = 200, y_higher_is_better = FALSE)
pte_results

######## CONTINUOUS
rm(list = ls())
data(cpt1)
X$site = NULL

library(PTE)
pte_results = PTE_bootstrap_inference(X, y, regression_type = "continuous", B = 200, y_higher_is_better = FALSE)
pte_results


######## INCIDENCE
rm(list = ls())
data(cpt1)
X$site = NULL
y = ifelse(y > 15, 0, 1) #force incidence here and y=1 is better (not depressed)

library(PTE)
pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", B = 400)
pte_results
pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", B = 200, 
                                      incidence_metric = "risk_ratio")
pte_results
pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", B = 200, 
                                      incidence_metric = "probability_difference")
pte_results











