create_PTE_results_object = function(raw_results, y_higher_is_better){
	n = nrow(raw_results)
	
	#build return object for the user
	return_obj = list()
	return_obj$results = raw_results
	
	#get real y's
	indices_1_1 = raw_results$trt0 == 1 & raw_results$optimal == 1 #optimal
	indices_0_0 = raw_results$trt0 == 0 & raw_results$optimal == 0 #optimal
	indices_0_1 = raw_results$trt0 == 0 & raw_results$optimal == 1 #non-optimal
	indices_1_0 = raw_results$trt0 == 1 & raw_results$optimal == 0 #non-optimal
	
	Ya = raw_results[indices_0_0, "real_y"] 
	Yb = raw_results[indices_0_1, "real_y"]
	Yc = raw_results[indices_1_0, "real_y"]
	Yd = raw_results[indices_1_1, "real_y"]
	return_obj$Ya = Ya
	return_obj$Yb = Yb
	return_obj$Yc = Yc
	return_obj$Yd = Yd
	
	#four groups
	summary_square = matrix(NA, 3, 3)
	rownames(summary_square) = c("received Trt 0", "received Trt 1", "Avg Opt")
	colnames(summary_square) = c("optimal Trt 0", "optimal Trt 1", "Avg Trt")
	summary_square[1, 1] = mean(Ya, na.rm = TRUE)
	summary_square[1, 2] = mean(Yb, na.rm = TRUE)
	summary_square[2, 1] = mean(Yc, na.rm = TRUE)
	summary_square[2, 2] = mean(Yd, na.rm = TRUE)
	summary_square[1, 3] = mean(c(Ya, Yb), na.rm = TRUE)
	summary_square[2, 3] = mean(c(Yc, Yd), na.rm = TRUE)
	summary_square[3, 1] = mean(c(Ya, Yc), na.rm = TRUE)
	summary_square[3, 2] = mean(c(Yb, Yd), na.rm = TRUE)
	summary_square[3, 3] = mean(c(Ya, Yb, Yc, Yd), na.rm = TRUE)
	return_obj$summary_square = summary_square
	return_obj$is_bad = is.nan(sum(summary_square)) #if any of them are NaN this will pick it up
	
	#four groups
	ns = matrix(NA, 3, 3)
	rownames(ns) = c("received Trt 0", "received Trt 1", "Totals")
	colnames(ns) = c("optimal Trt 0", "optimal Trt 1", "Totals")
	ns[1, 1] = length(Ya[!is.na(Ya)])
	ns[1, 2] = length(Yb[!is.na(Yb)])
	ns[2, 1] = length(Yc[!is.na(Yc)])
	ns[2, 2] = length(Yd[!is.na(Yd)])
	#now add em up to get tots
	ns[1, 3] = ns[1, 1] + ns[1, 2]
	ns[2, 3] = ns[2, 1] + ns[2, 2]
	ns[3, 1] = ns[1, 1] + ns[2, 1]
	ns[3, 2] = ns[1, 2] + ns[2, 2]
	ns[3, 3] = ns[1, 1] + ns[1, 2] + ns[2, 1] + ns[2, 2]
	return_obj$ns = ns
	return_obj$pct_data_used = round(ns[3, 3] / n, 3)

	
	if (regression_type == "continuous"){
		sse = sum((raw_results$real_y - raw_results$est_true)^2, na.rm = TRUE)
		return_obj$oos_rmse = sqrt(sse / n)
		sst = sum((raw_results$real_y - mean(raw_results$real_y, na.rm = TRUE))^2, na.rm = TRUE)
		return_obj$out_of_sample_Rsq = 1 - sse / sst
	}

	if (regression_type == "survival"){
		#first do adversarial
		y_opt = c(Ya, Yd)
		censor_opt = c()
		
		mod = survfit(Surv() ~ sampledata$CONDAC)
		median_diff = summary(mod)$table[,'median'][1] -
						summary(mod)$table[,'median'][2]
		
		
		return_obj$q_adversarial = 
		return_obj$q_average = 
		return_obj$q_best = 
		
	} else {
		return_obj$pred_differences_avg = mean(abs(raw_results[, 1] - raw_results[, 2]), na.rm = TRUE)
		return_obj$pred_differences_sd = sd(abs(raw_results[, 1] - raw_results[, 2]), na.rm = TRUE)
		return_obj$avg_optimals = mean(c(Ya, Yd), na.rm = TRUE)
		return_obj$avg_non_optimals = mean(c(Yb, Yc), na.rm = TRUE)
		return_obj$avg_all = mean(c(Ya, Yb, Yc, Yd), na.rm = TRUE)
		avg_ys_tx_1 = mean(c(Ya, Yb), na.rm = TRUE)
		avg_ys_tx_2 = mean(c(Yc, Yd), na.rm = TRUE)
		if (avg_ys_tx_1 >= avg_ys_tx_2 && y_higher_is_better){ #sometimes continuous data aint continuous and you can have a "measure 0" event of equality here - at equality should pick group 1 or 2 with equal prob (not done)
			return_obj$avg_best = avg_ys_tx_1		
		} else if (avg_ys_tx_1 >= avg_ys_tx_2 && !y_higher_is_better){
			return_obj$avg_best = avg_ys_tx_2
		} else if (avg_ys_tx_1 < avg_ys_tx_2 && y_higher_is_better){
			return_obj$avg_best = avg_ys_tx_2			
		} else if (avg_ys_tx_1 < avg_ys_tx_2 && !y_higher_is_better){
			return_obj$avg_best = avg_ys_tx_1
		}
		return_obj$q_adversarial = return_obj$avg_optimals - return_obj$avg_non_optimals
		return_obj$q_average = return_obj$avg_optimals - return_obj$avg_all
		return_obj$q_best = return_obj$avg_optimals - return_obj$avg_best		
	}

	return_obj	
}
