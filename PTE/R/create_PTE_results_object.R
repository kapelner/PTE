create_PTE_results_object = function(results, regression_type, y_higher_is_better, difference_function, incidence_metric){
	n = nrow(results)
	
	#build return object for the user
	return_obj = list()
	return_obj$results = results
	
	#get real y's
	y_all = results$real_y
	indices_1_1 = results$given_tx == 1 & results$rec_tx == 1 #optimal
	indices_0_0 = results$given_tx == 0 & results$rec_tx == 0 #optimal
	indices_0_1 = results$given_tx == 0 & results$rec_tx == 1 #non-optimal
	indices_1_0 = results$given_tx == 1 & results$rec_tx == 0 #non-optimal
	
	if (is.null(difference_function)){	
		Y00 = y_all[indices_0_0] 
		Y01 = y_all[indices_0_1]
		Y10 = y_all[indices_1_0]
		Y11 = y_all[indices_1_1]
		return_obj$Y00 = Y00
		return_obj$Y01 = Y01
		return_obj$Y10 = Y10
		return_obj$Y11 = Y11
		
		#four groups summary square - less useful for survival, but still nice to see
		summary_square = matrix(NA, 3, 3)
		rownames(summary_square) = c("received Trt 0", "received Trt 1", "Avg Opt")
		colnames(summary_square) = c("optimal Trt 0", "optimal Trt 1", "Avg Trt")
		summary_square[1, 1] = mean(Y00, na.rm = TRUE)
		summary_square[1, 2] = mean(Y01, na.rm = TRUE)
		summary_square[2, 1] = mean(Y10, na.rm = TRUE)
		summary_square[2, 2] = mean(Y11, na.rm = TRUE)
		summary_square[1, 3] = mean(c(Y00, Y01), na.rm = TRUE)
		summary_square[2, 3] = mean(c(Y10, Y11), na.rm = TRUE)
		summary_square[3, 1] = mean(c(Y00, Y10), na.rm = TRUE)
		summary_square[3, 2] = mean(c(Y01, Y11), na.rm = TRUE)
		summary_square[3, 3] = mean(y_all, na.rm = TRUE)
		return_obj$summary_square = summary_square
		return_obj$is_bad = is.nan(sum(summary_square)) #if any of them are NaN this will pick it up
		
		#four groups
		ns = matrix(NA, 3, 3)
		rownames(ns) = c("received Trt 0", "received Trt 1", "Totals")
		colnames(ns) = c("optimal Trt 0", "optimal Trt 1", "Totals")
		ns[1, 1] = length(Y00[!is.na(Y00)])
		ns[1, 2] = length(Y01[!is.na(Y01)])
		ns[2, 1] = length(Y10[!is.na(Y10)])
		ns[2, 2] = length(Y11[!is.na(Y11)])
		#now add em up to get tots
		ns[1, 3] = ns[1, 1] + ns[1, 2]
		ns[2, 3] = ns[2, 1] + ns[2, 2]
		ns[3, 1] = ns[1, 1] + ns[2, 1]
		ns[3, 2] = ns[1, 2] + ns[2, 2]
		ns[3, 3] = ns[1, 1] + ns[1, 2] + ns[2, 1] + ns[2, 2]
		return_obj$ns = ns
		return_obj$pct_data_used = round(ns[3, 3] / n * 100, 3)
	
		

		
		if (regression_type == "continuous" || (regression_type == "incidence" && incidence_metric == "probability_difference")){
			return_obj$pred_differences_avg = mean(abs(results[, 1] - results[, 2]), na.rm = TRUE)
			return_obj$pred_differences_sd = sd(abs(results[, 1] - results[, 2]), na.rm = TRUE)
			return_obj$avg_rec = mean(c(Y00, Y11), na.rm = TRUE)
			return_obj$avg_non_rec = mean(c(Y01, Y10), na.rm = TRUE)
			return_obj$avg_all = mean(y_all, na.rm = TRUE)
			avg_ys_tx_1 = mean(c(Y00, Y01), na.rm = TRUE)
			avg_ys_tx_2 = mean(c(Y10, Y11), na.rm = TRUE)
			if (avg_ys_tx_1 >= avg_ys_tx_2 && y_higher_is_better){
				return_obj$avg_best = avg_ys_tx_1		
			} else if (avg_ys_tx_1 >= avg_ys_tx_2 && !y_higher_is_better){
				return_obj$avg_best = avg_ys_tx_2
			} else if (avg_ys_tx_1 < avg_ys_tx_2 && y_higher_is_better){
				return_obj$avg_best = avg_ys_tx_2			
			} else if (avg_ys_tx_1 < avg_ys_tx_2 && !y_higher_is_better){
				return_obj$avg_best = avg_ys_tx_1
			}
			return_obj$q_adversarial = return_obj$avg_rec - return_obj$avg_non_rec
			return_obj$q_average = return_obj$avg_rec - return_obj$avg_all
			return_obj$q_best = return_obj$avg_rec - return_obj$avg_best
			#some more metrics that may be of use for the continuous case
			if (regression_type == "continuous" ){
				sse = sum((results$real_y - results$est_true)^2, na.rm = TRUE)
				return_obj$oos_rmse = sqrt(sse / n)
				sst = sum((results$real_y - mean(results$real_y, na.rm = TRUE))^2, na.rm = TRUE)
				return_obj$out_of_sample_Rsq = 1 - sse / sst
			}
		} else if (regression_type == "incidence"){
			#set up all the data
			p_all = mean(y_all, na.rm = TRUE)
			p_rec = mean(c(Y00, Y11), na.rm = TRUE)
			p_non_rec = mean(c(Y01, Y10), na.rm = TRUE)
			p_1 = mean(c(Y10, Y11), na.rm = TRUE)
			p_0 = mean(c(Y00, Y01), na.rm = TRUE)
			
			if (p_1 >= p_0 && y_higher_is_better){
				p_best = p_1
			} else if (p_1 < p_0 && y_higher_is_better){
				p_best = p_0			
			} else if (p_1 <= p_0 && !y_higher_is_better){
				p_best = p_0
			} else {
				p_best = p_1
			}
			
			if (incidence_metric == "risk_ratio"){
				return_obj$q_adversarial = p_rec / p_non_rec
				return_obj$q_average = p_rec / p_all
				return_obj$q_best = p_rec / p_best
			} else if (incidence_metric == "odds_ratio"){
				p_rec_odds = (p_rec / (1 - p_rec))
				return_obj$q_adversarial = p_rec_odds / (p_non_rec / (1 - p_non_rec))
				return_obj$q_average = p_rec_odds / (p_all / (1 - p_all))
				return_obj$q_best = p_rec_odds / (p_best / (1 - p_best))
			}
		} else if (regression_type == "survival"){ #this set should be collectively exhaustive
			c_all = results$censored
			
			#first do adversarial
			y_rec = c(Y00, Y11)
			y_non_rec = c(Y01, Y10)			
			c_rec = c(c_all[indices_0_0], c_all[indices_1_1])
			c_non_rec = c(c_all[indices_0_1], c_all[indices_1_0])
			rec_indicator = c(rep(1, length(y_rec)), rep(0, length(y_non_rec)))
						
			mod = survfit(Surv(c(y_rec, y_non_rec), c(c_rec, c_non_rec)) ~ rec_indicator)
			median_estimates = summary(mod)$table[, 'median']
			median_diff = median_estimates[2] - median_estimates[1]
#			plot(mod, col = c("blue", "red"))					
			return_obj$q_adversarial = median_diff
			
			#now do average
			rec_indicator = c(rep(1, length(y_rec)), rep(0, length(y_all)))
			
			mod = survfit(Surv(c(y_rec, y_all), c(c_rec, c_all)) ~ rec_indicator)
			median_estimates = summary(mod)$table[, 'median']
			median_diff = median_estimates[2] - median_estimates[1]
#			plot(mod, col = c("blue", "red"))			
			return_obj$q_average = median_diff
			
			#now do best... this requires running two models... the first to see who's better on "average" (or
			#at the median) and the second to run the best tx vs the recommended tx.
			y_1 = c(Y10, Y11)
			c_1 = c(c_all[indices_1_0], c_all[indices_1_1])
			y_0 = c(Y01, Y00)
			c_0 = c(c_all[indices_0_1], c_all[indices_0_0])
			rec_indicator = c(rep(1, length(y_1)), rep(0, length(y_0))) #call tx 1 the 1 level
			mod = survfit(Surv(c(y_1, y_0), c(c_1, c_0)) ~ rec_indicator)
#			plot(mod, col = c("blue", "red"))
			median_estimates = summary(mod)$table[, 'median']
			median_diff = median_estimates[2] - median_estimates[1]
			#if median_diff is greater than 0, tx1 is the better treatment
			if (median_diff == 0){ #measure 0 event in theory... but in practice - who knows?
				random_bernoulli = runif(1) < 0.5
				y_best = if (random_bernoulli){y_0} else {y_1}
				c_best = if (random_bernoulli){c_0} else {c_1}
			} else if (median_diff > 0){
				y_best = y_1
				c_best = c_1
			} else {
				y_best = y_0
				c_best = c_0	
			}
			
			rec_indicator = c(rep(1, length(y_rec)), rep(0, length(y_best)))
			
			mod = survfit(Surv(c(y_rec, y_best), c(c_rec, c_best)) ~ rec_indicator)
			median_estimates = summary(mod)$table[, 'median']
			median_diff = median_estimates[2] - median_estimates[1]
#			plot(mod, col = c("blue", "red"))			
			return_obj$q_best = median_diff
			
		}
	} else { ###user custom function output
		all_diffs = difference_function(results, indices_1_1, indices_0_0, indices_0_1, indices_1_0)
		return_obj$q_adversarial = all_diffs[1]
		return_obj$q_average = all_diffs[2]
		return_obj$q_best = all_diffs[3]	
	}
	return_obj	
}
