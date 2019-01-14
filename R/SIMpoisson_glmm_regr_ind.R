#' Simulation function for a Poisson model with one random component
#'
#' The function simulates data under a Poisson model with given explanatory variables \code{x} and parameter vector \code{theta}, consisting of the coefficients of the fixed effects, beta, and sigma, the standard deviation of a single random effect with distribution \code{N(0, sigma^2)}. For each of the \code{nk} simulated datasets, a vector of summary statistics is computed, consisting of the regression coefficients estimated in a Poisson model without random effects, and an estimate of the dispersion based on the Pearson statistic (sum of the squared Pearson residuals) of this model. The function is designed to estimate theta by the approximate maximum likelihood algorithm in \code{KDKW.FD} or \code{KDKW.SP}.
#'
#'
#' SIMpoisson_glmm uses \code{glm} to obtain the summary statistics. This probably has a lot of overhead and could maybe be replaced by a faster alternative.
#'
#' Note that this is an experimental function that is only for use with x consisting of a set of dummy variables that is equivalent to a single factor that defines groups. The summary statistics are computed separately for each independent group.
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
#'
#' x = matrix(c(rep(1, 50), rep(0, 100), rep(1, 50)), ncol=2 )
#' testsim = SIMpoisson_glmm_regr_ind(100, c(3,5,1), x)
#' apply(testsim, 2, mean)
#'
#' @export

SIMpoisson_glmm_regr_ind = function(nk, theta, x){

  x = as.matrix(x)
  n = nrow(x)
  k = ncol(x)
  beta = theta[1:k]
  sigma = theta[length(theta)]
  sumstat = matrix(NA, ncol=k*2, nrow=nk)

  x.logical = matrix(as.logical(x), nrow=nrow(x), ncol = ncol(x))

  for(i in 1:nk){
    # simulate random effect
    u = rnorm(n, 0, sigma)
    # simulate Poisson data
    eta = x%*%beta + u
    y = rpois(n, exp(eta))

    # summary statistics

    for(j in 1:k){

      fit_glm = glm(y[x.logical[,j]] ~ 1, family=poisson)
      sumstat[i,2*j] = sqrt(sum(residuals(fit_glm, type="pearson")^2))
      sumstat[i,2*j-1] = fit_glm$coefficients

    }

  }

  sumstat

}


