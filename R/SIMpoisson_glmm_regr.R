#' Simulation function for a Poisson model with one random component
#'
#' The function simulates data under a Poisson model with given explanatory variables \code{x} and parameter vector \code{theta}, consisting of the coefficients of the fixed effects, beta, and sigma, the standard deviation of a single random effect with distribution \code{N(0, sigma^2)}. For each of the \code{nk} simulated datasets, a vector of summary statistics is computed, consisting of the regression coefficients estimated in a Poisson model without random effects, and an estimate of the dispersion based on the Pearson statistic (sum of the squared Pearson residuals) of this model. The function is designed to estimate theta by the approximate maximum likelihood algorithm in \code{KDKW.FD} or \code{KDKW.SP}.
#'
#' SIMpoisson_glmm uses \code{glm} to obtain the summary statistics. This probably has a lot of overhead and could maybe be replaced by a faster alternative.
#'
#' @param nk integer, number of datasets to be simulated
#' @param theta numeric, parameter vector. theta = c(beta, sigma)
#' @param x numeric matrix, explanatory variables
#'
#' @author Johanna Bertl
#'
#' @seealso \link{\code{glm}}
#'
#' @examples
#'
#' set.seed(1234)
#' x = matrix(runif(100))
#'
#' testsim = SIMpoisson_glmm_regr(200, c(2,1), x)
#' plot(testsim[,1], testsim[,2])
#' apply(testsim, 2, mean)
#'
#' testsim2 = SIMpoisson_glmm_regr(200, c(5,1), x)
#' plot(testsim2[,1], testsim2[,2])
#' apply(testsim2, 2, mean)
#'
#' x = matrix(c(rep(1, 50), rep(0, 100), rep(1, 50)), ncol=2 )
#' testsim3 = SIMpoisson_glmm_regr(100, c(3,5,1), x)
#' apply(testsim3, 2, mean)
#'
#' testsim_beta = SIMpoisson_glmm_beta(200, 2, sigma=1, x)
#'
#' testsim_sigma = SIMpoisson_glmm_sigma(20, 1, beta = 2, x)
#'
#' @export

SIMpoisson_glmm_regr = function(nk, theta, x){

  x = as.matrix(x)
  n = nrow(x)
  k = ncol(x)
  beta = theta[1:k]
  sigma = theta[length(theta)]
  sumstat = matrix(NA, ncol=k+1, nrow=nk)

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


#' @describeIn SIMpoisson_glmm_regr SIMpoisson_glmm_beta performs the same simulations as SIMpoisson_glmm, but it is used to estimate beta only when sigma is known.
#' @export

SIMpoisson_glmm_beta = function(nk, theta, sigma, x){

  x = as.matrix(x)
  n = nrow(x)
  k = ncol(x)
  beta = theta
  sumstat = matrix(NA, ncol=1, nrow=nk)

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


#' @describeIn SIMpoisson_glmm_regr SIMpoisson_glmm_sigma performs the same simulations as SIMpoisson_glmm, but it is used to estimate sigma only when beta is known.
#' @export

SIMpoisson_glmm_sigma = function(nk, theta, beta, x){

  x = as.matrix(x)
  n = nrow(x)
  k = ncol(x)
  sigma = theta
  sumstat = matrix(NA, ncol=1, nrow=nk)

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



