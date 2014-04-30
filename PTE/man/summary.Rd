\name{summary.PTE_bootstrap_results}
\alias{summary.PTE_bootstrap_results}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print the results of a bootstrap model run to screen
}
\description{
Prints p-values and confidence intervals to the screen for both the random and best business-as-usual allocation procedures.
}
\usage{
\method{summary}{PTE_bootstrap_results}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	An object of class ``PTE_bootstrap_results''.
}
  \item{...}{
	Parameters that are ignored.
}
}

\author{
Adam Kapelner and Justin Bleich
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
bootstrap_inference
}
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
			B = num_boot, 
			plot = FALSE)
	summary(res)
}
\references{
Kapelner, A, Bleich, J, Cohen, ZD, DeRubeis, RJ and Berk, R (2014) Inference for Treatment Regime Models in Personalized Medicine, arXiv
}
