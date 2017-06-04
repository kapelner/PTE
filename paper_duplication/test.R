######## SURVIVAL

#Alina's code to get data loaded
rm(censoredtimes)
rm(noncensoredtimes)
rm(censor_indicator)
head(sampledata)
sampledata = na.omit(sampledata)
sampledata$AD = NULL

#recode
sampledata = sampledata[!is.na(sampledata$TIME), ]
y = sampledata$TIME
sampledata$TIME = NULL
censored = sampledata$censor_indicator
sampledata$censor_indicator = NULL
colnames(sampledata)[3] = "treatment"
X = sampledata
X$treatment = X$treatment - 1
rm(sampledata)
head(X)

library(PTE)
pte_results = PTE_bootstrap_inference(X, y, censored = censored, regression_type = "survival", B = 200)


######## CONTINUOUS
rm(list = ls())
load("../PTE/data/cpt1.RData")
X$site = NULL

library(PTE)
pte_results = PTE_bootstrap_inference(X, y, regression_type = "continuous", B = 200, y_higher_is_better = FALSE)
pte_results


######## INCIDENCE

rm(list = ls())
load("../PTE/data/cpt1.RData")
X$site = NULL
y = ifelse(y > 15, 0, 1) #force incidence here and y=1 is better (not depressed)

library(PTE)
pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", B = 200)
pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", incidence_metric = "risk_ratio", B = 200)
pte_results = PTE_bootstrap_inference(X, y, regression_type = "incidence", incidence_metric = "probability_difference", B = 200)
pte_results











