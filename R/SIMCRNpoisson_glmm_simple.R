#' Simulation function for a Poisson model with fixed mean and random intercept
#'
#' The function simulates data under a Poisson model with fixed mean lambda and a random intercept with with distribution \code{N(0, sigma^2)}. For each of the \code{nk} simulated datasets, a vector of summary statistics is computed, consisting of the mean and an estimate of the dispersion based on the Pearson statistic (sum of the squared Pearson residuals) of this model. The function is designed to estimate theta by the approximate maximum likelihood algorithm in \code{KDKW.FD} or \code{KDKW.SP}.
#'
#' The simulations are obtained with the given seed (designed for the use of Common Random Numbers in the Approximate Maximum Likelihood Algorithm).
#'
#' SIMpoisson_glmm_simple uses simGLMM::simDataGLMM.
#'
#' This function might be replaced by a more flexible function based on simGLMM::simDataGLMMdesign that uses a design matrix.
#'
#' @param nk integer, number of datasets to be simulated
#' @param theta numeric, parameter vector. theta = c(lambda, sigma)
#' @param seed integer. Seed to simulate the random effect and poisson distributed response.
#' @param n integer, number of observations
#'
#' @author Johanna Bertl
#'
#' @examples
#'
#' SIMCRNpoisson_glmm_simple(1, c(1,0.5), 1234, 10)
#'
#' set.seed(1234)
#' dat = simDataGLMM(clus=10, rep=1, fixedMean=1, fixedCoef=NULL, fixedMat=NULL, covM = 0.5, disFam=poisson())
#' mean(dat$res)
#'
#' @export


SIMCRNpoisson_glmm_simple = function(nk, theta, seed, n){

  if(theta[2]<=0) stop("SIMCRNpoisson_glmm_simple: sigma (theta[2]) must be positive.")

  require(simGLMM)

  set.seed(seed)

  data = simDataGLMM(clus=n*nk, rep=1, fixedMean=theta[1], covM = as.matrix(theta[2]), disFam=poisson())

  h = factor(rep(1:nk, each=n))

  mean.vec = unlist(by(data$res, INDICES=h, mean))
  mean.vec.ext = rep(mean.vec, each=n)

  pearson = (data$res - mean.vec.ext)/mean.vec.ext
  pearson.vec = unlist(by(pearson, INDICES=h, FUN = function(x) sqrt(sum(x^2))))

  return(cbind(mean.vec, pearson.vec))

}

