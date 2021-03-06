% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIMpoisson_glmm_regr.R
\name{SIMpoisson_glmm_regr}
\alias{SIMpoisson_glmm_regr}
\alias{SIMpoisson_glmm_beta}
\alias{SIMpoisson_glmm_sigma}
\title{Simulation function for a Poisson model with one random component}
\usage{
SIMpoisson_glmm_regr(nk, theta, x)

SIMpoisson_glmm_beta(nk, theta, sigma, x)

SIMpoisson_glmm_sigma(nk, theta, beta, x)
}
\arguments{
\item{nk}{integer, number of datasets to be simulated}

\item{theta}{numeric, parameter vector. theta = c(beta, sigma)}

\item{x}{numeric matrix, explanatory variables}
}
\description{
The function simulates data under a Poisson model with given explanatory variables \code{x} and parameter vector \code{theta}, consisting of the coefficients of the fixed effects, beta, and sigma, the standard deviation of a single random effect with distribution \code{N(0, sigma^2)}. For each of the \code{nk} simulated datasets, a vector of summary statistics is computed, consisting of the regression coefficients estimated in a Poisson model without random effects, and an estimate of the dispersion based on the Pearson statistic (sum of the squared Pearson residuals) of this model. The function is designed to estimate theta by the approximate maximum likelihood algorithm in \code{KDKW.FD} or \code{KDKW.SP}.
}
\details{
SIMpoisson_glmm uses \code{glm} to obtain the summary statistics. This probably has a lot of overhead and could maybe be replaced by a faster alternative.
}
\section{Functions}{
\itemize{
\item \code{SIMpoisson_glmm_beta}: SIMpoisson_glmm_beta performs the same simulations as SIMpoisson_glmm, but it is used to estimate beta only when sigma is known.

\item \code{SIMpoisson_glmm_sigma}: SIMpoisson_glmm_sigma performs the same simulations as SIMpoisson_glmm, but it is used to estimate sigma only when beta is known.
}}

\examples{

set.seed(1234)
x = matrix(runif(100))

testsim = SIMpoisson_glmm_regr(200, c(2,1), x)
plot(testsim[,1], testsim[,2])
apply(testsim, 2, mean)

testsim2 = SIMpoisson_glmm_regr(200, c(5,1), x)
plot(testsim2[,1], testsim2[,2])
apply(testsim2, 2, mean)

x = matrix(c(rep(1, 50), rep(0, 100), rep(1, 50)), ncol=2 )
testsim3 = SIMpoisson_glmm_regr(100, c(3,5,1), x)
apply(testsim3, 2, mean)

testsim_beta = SIMpoisson_glmm_beta(200, 2, sigma=1, x)

testsim_sigma = SIMpoisson_glmm_sigma(20, 1, beta = 2, x)

}
\seealso{
\link{\code{glm}}
}
\author{
Johanna Bertl
}
