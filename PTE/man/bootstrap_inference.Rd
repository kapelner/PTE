\name{bootstrap_inference}
\alias{bootstrap_inference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap inference for prespecified models
}
\description{
Runs B bootstrap samples using a prespecified model then computes the two I estimates based on cross validation. 
p values of the two I estimates are computed for a given \eqn{H_0: \mu_{I_0} = \mu_0}{H_0: mu_I_0 = mu_0} and 
confidence intervals are provided.
}
\usage{
bootstrap_inference(X, y, 
		model_string,
		predict_string = "predict(mod, obs_left_out)",
		cleanup_mod_function = NA,
		y_higher_is_better = TRUE,
		verbose = TRUE,
		full_verbose = FALSE,
		H_0_mu_equals = 0,
		pct_leave_out = 0.10,
		B = 3000,
		m_prop = 1,
		alpha = 0.05,
		plot = TRUE,
        num_cores = 1,
        ...)
}
		
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A \eqn{n \times p}{n x p} dataframe of covariates.
}
	\item{y}{
An \eqn{n}-length numeric vector which is the response
}
  \item{model_string}{
A string of R code that will be evaluated to construct the leave one out model. Make sure the covariate data is
referred to as \code{Xyleft}.
}
  \item{predict_string}{
A string of R code that will be evaluated on left out data after the model is built with the training data. Make sure
the forecast data (the left one out data) is referred to as \code{obs_left_out} and the model is referred to as \code{mod}.
}
  \item{cleanup_mod_function}{
A function that is called at the end of a cross validation iteration to cleanup the model 
in some way. This is used for instance if you would like to release the memory your model is using but generall does not apply.
The default is \code{NA} for "no function."
}
  \item{y_higher_is_better}{
True if a response value being higher is clinically "better" than one that is lower (e.g. cognitive ability in a drug trial for the 
mentally ill). False if the response value being lower is clinically "better" than one that is higher (e.g. amount of weight lost 
in a weight-loss trial). Default is \code{TRUE}.
}
  \item{verbose}{
Prints out a dot for each bootstrap sample. This only works on some platforms.
}
  \item{full_verbose}{
Prints out full information for each cross validation model for each bootstrap sample. This only works on some platforms.
}
  \item{H_0_mu_equals}{
The \eqn{\mu_{I_0}}{mu_I_0} value in \eqn{H_0}{H_0}. Default is 0 which answers the question: does my allocation procedure do better than a naive
allocation procedure.
}
  \item{pct_leave_out}{
In the cross-validation, the proportion of the original dataset left out to estimate out-of-sample metrics. The default is 0.1
which corresponds to 10-fold cross validation.
}
  \item{B}{
The number of bootstrap samples to take. We recommend making this as high as you can tolerate given speed considerations.
The default is 3000.
}
  \item{m_prop}{
Within each bootstrap sample, the proportion of the total number of rows of \code{X} to sample without replacement. \code{m_prop < 1} ensures
the number of rows sampled is less than \code{n} which fixes the consistency of the bootstrap estimator of a non-smooth functional. The default 
is 1 since non-smoothness may not be a common issue. 
}
  \item{alpha}{
Defines the confidence interval size (1 - alpha). Defaults to 0.05.
}
  \item{plot}{
Illustrates the estimate, the bootstrap samples and the confidence intervals on a histogram plot. Default to TRUE.
}
  \item{num_cores}{
The number of cores to use in parallel to run the bootstrap samples more rapidly. Defaults to serial by using 1 core.  
}
  \item{...}{
Additional parameters to be sent to the model constructor. Note that if you wish to pass these parameters, 
"..." must be specified in model_string.  
}
}
\value{
Returns a list object containing results of the procedure.
}

\author{
Adam Kapelner and Justin Bleich
}
\references{
Kapelner, A, Bleich, J, Cohen, ZD, DeRubeis, RJ and Berk, R (2014) Inference for Treatment Regime Models in Personalized Medicine, arXiv
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
	beta0 = 1
	beta1 = -1
	gamma0 = 0
	gamma1 = sqrt(2 * pi)
	mu_x = 0
	sigsq_x = 1
	sigsq_e = 1
	num_boot = 20 #for speed only
	n = 50 #for speed only
	
	x = sort(rnorm(n, mu_x, sigsq_x))
	noise = rnorm(n, 0, sigsq_e)
	
	treatment = sample(c(rep(1, n / 2), rep(0, n / 2)))
	y = beta0 + beta1 * x + treatment * (gamma0 + gamma1 * x) + noise
	
	X = data.frame(treatment, x)
	
	res = bootstrap_inference(X, y,
			"lm(y ~ . + treatment * ., data = Xyleft)",
			num_cores = 1,
			B = num_boot)
}
