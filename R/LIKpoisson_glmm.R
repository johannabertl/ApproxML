#' Likelihood of the Poisson GLMM
#'
#' The likelihood of the parameters beta and sigma in a Poisson GLMM is computed. The integral is calculated with Gauss-Hermite quadrature.
#'
#' The function dpois is used to evaluate the conditional likelihood, as it is fast (implemented in C) and can handle large values of y, that can otherwise be problematic in the evaluation of factorial(y).
#'
#' The loglikelihood can sometimes be obtained even when the likelihood evaluates to zero (numerically).
#'
#' @param beta scalar
#' @param sigma scalar
#' @param x vector, regressor
#' @param y vector, regressand
#'
#' @examples
#'
#' x = 1:3
#' y = c(2,4,6)
#' beta = 1
#' sigma = 1
#'
#' betavec = seq(0, 5, by =0.1)
#' likvec = numeric(length(betavec))
#' for(i in 1:length(likvec)) likvec[i] = LIKpoisson_glmm(betavec[i], sigma, x, y)
#' plot(betavec, likvec)
#'
#'
#' @export

LIKpoisson_glmm = function(beta, sigma, x, y){
  prod(mapply(LIKpoisson_glmm_single, x, y, MoreArgs=list(beta = beta, sigma = sigma)))
}

#' @describeIn LIKpoisson_glmm LOGLIKpoisson_glmm directly sums the log of the likelihood terms that are computed in the same way as by the function LIKpoisson_glmm.
#' @export

LOGLIKpoisson_glmm = function(beta, sigma, x, y){
  sum(log(mapply(LIKpoisson_glmm_single, x, y, MoreArgs=list(beta = beta, sigma = sigma))))
}

LIKpoisson_glmm_single = function(x, y, beta, sigma){
  HermiteIntegration(CONDLIKpoisson_glmm_single, sig2 = sigma^2, beta = beta, x = x, y = y)
}

CONDLIKpoisson_glmm_single = function(u, beta, x, y) {
  dpois(y, exp( x %*% beta + u))
}

