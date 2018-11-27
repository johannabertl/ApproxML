#' Simulation function for a Poisson model with one random component
#'
#' The function simulates data under a Poisson model with given explanatory variables \code{x} and parameter vector \code{theta}, consisting of the coefficients of the fixed effects, beta, and sigma, the standard deviation of a single random effect with distribution \code{N(0, sigma^2)}. For each of the \code{nk} simulated datasets, a vector of summary statistics is computed, consisting of the regression coefficients estimated in a Poisson model without random effects, and an estimate of the dispersion based on the Pearson statistic (sum of the squared Pearson residuals) of this model. The function is designed to estimate theta by the approximate maximum likelihood algorithm in \code{KDKW.FD} or \code{KDKW.SP}.
#'
#' The simulations are obtained with the given seed (designed for the use of Common Random Numbers in the Approximate Maximum Likelihood Algorithm).
#'
#' SIMpoisson_glmm uses \code{glm} to obtain the summary statistics. This probably has a lot of overhead and could maybe be replaced by a faster alternative.
#'
#' @param nk integer, number of datasets to be simulated
#' @param theta numeric, parameter vector. theta = c(beta, sigma)
#' @param seed integer. Seed to simulate the random effect and poisson distributed response.
#' @param x numeric matrix, explanatory variables
#'
#' @author Johanna Bertl
#'
#' @seealso \link{\code{glm}}
#'
#' @examples
#'
#'
#' x = matrix(c(rep(1, 50), rep(0, 100), rep(1, 50)), ncol=2 )
#' testsim = SIMCRNpoisson_glmm_regr(100, c(3,5,1), seed = 1234, x)
#' apply(testsim, 2, mean)
#' testsim2 = SIMCRNpoisson_glmm_regr(100, c(3,5,1), seed = 1234, x)
#' apply(testsim2, 2, mean)
#'
#' @export

SIMCRNpoisson_glmm_regr = function(nk, theta, seed, x){

  x = as.matrix(x)
  n = nrow(x)
  k = ncol(x)
  beta = theta[1:k]
  sigma = theta[length(theta)]
  sumstat = matrix(NA, ncol=k+1, nrow=nk)

  set.seed(seed)

  for(i in 1:nk){
    # simulate random effect
    u = rnorm(n, 0, sigma)
    # simulate Poisson data
    eta = x%*%beta + u
    y = rpois(n, exp(eta))

    # summary statistics
    fit_glm = glm(y ~ 0+x, family=poisson)
    sumstat[i,k+1] = sqrt(sum(residuals(fit_glm, type="pearson")^2))
    sumstat[i,1:k] = fit_glm$coefficients

  }

  sumstat

}


#' @describeIn SIMCRNpoisson_glmm_regr SIMCRNpoisson_glmm_beta performs the same simulations as SIMpoisson_glmm, but it is used to estimate beta only when sigma is known.
#' @export

SIMCRNpoisson_glmm_beta = function(nk, theta, seed, sigma, x){

  x = as.matrix(x)
  n = nrow(x)
  k = ncol(x)
  beta = theta
  sumstat = matrix(NA, ncol=1, nrow=nk)

  set.seed(seed)

  for(i in 1:nk){
    # simulate random effect
    u = rnorm(n, 0, sigma)
    # simulate Poisson data
    eta = x%*%beta + u
    y = rpois(n, exp(eta))

    # summary statistics
    fit_glm = glm(y ~ 0+x, family=poisson)
    sumstat[i,1] = fit_glm$coefficients

  }

  sumstat

}


#' @describeIn SIMCRNpoisson_glmm_regr SIMCRNpoisson_glmm_sigma performs the same simulations as SIMpoisson_glmm, but it is used to estimate sigma only when beta is known.
#' @export

SIMCRNpoisson_glmm_sigma = function(nk, theta, seed, beta, x){

  x = as.matrix(x)
  n = nrow(x)
  k = ncol(x)
  sigma = theta
  sumstat = matrix(NA, ncol=1, nrow=nk)

  set.seed(seed)

  for(i in 1:nk){
    # simulate random effect
    u = rnorm(n, 0, sigma)
    # simulate Poisson data
    eta = x%*%beta + u
    y = rpois(n, exp(eta))

    # summary statistics

    r = ( y - exp(x%*%beta) ) / sqrt(exp(x%*%beta))
    sumstat[i,1] = sqrt(sum(r^2))

  }

  sumstat

}



