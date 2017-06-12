#' Personalized Medicine...
#'
#' @name 		PTE
#' @docType 	package
#' @title 		Personalized Medicine Inference
#' @author 		Adam Kapelner \email{kapelner@@qc.cuny.edu}, Alina Levine and Justin Bleich
#' @references 	Kapelner, A, Bleich, J, Cohen, ZD, DeRubeis, RJ and Berk, R (2014) Inference for Treatment Regime Models in Personalized Medicine, arXiv
#' @keywords 	Personalized medicine, bootstrap
#' @import      foreach parallel doParallel survival stats graphics
##### Run "library(roxygen2); roxygenise("PTE", clean = TRUE)" to regenerate all Rd files and NAMESPACE and DESCRIPTION file
##### but make sure you are in the root directory of the project
NULL

#' Mock RCT data with a continuous endpoint.
#' 
#' A list with two objects (a) \code{X}, a dataframe with n rows representing clinical subjects and columns: 
#' treatment, x1, x2, x3, x4 and x5 where treatment is binary indicating the two arms of the clinical trial
#' and x1, ..., x5 are covariates that were collected about each subject and (b) \code{y}, a length n vector storing 
#' the continuous response values where, in this mock dataset, larger values indicate "better" outcomes for the 
#' subjects.
#'
#' @name continuous_example
#' @docType data
#' @author My Name \email{kapelner@@qc.cuny.edu}
NULL

#' Mock RCT data with a survival endpoint.
#' 
#' A list with three objects (a) \code{X}, a dataframe with n rows representing clinical subjects and columns: 
#' treatment, x1, x2, x3 and x4 where treatment is binary indicating the two arms of the clinical trial
#' and x1, ..., x4 are covariates that were collected about each subject (b) \code{y}, a length n vector storing 
#' the survival response values (a time measurement) where, in this mock dataset, smaller values indicate "better" 
#' survival outcomes for the subjects and (c) \code{censored}, a length n vector storing the censor dummies where 
#' c_16 = 1 means the response y_16 was censored and thus the truth value of y_16 is unknown and y_16 only represents
#' the moment it was censored (and c_16 = 0 means it was uncensored and y_16 is the true response value). 
#'
#' @name survival_example
#' @docType data
#' @author My Name \email{kapelner@@qc.cuny.edu}
NULL