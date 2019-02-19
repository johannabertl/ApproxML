#' Simulation function for a multivariate normal distribution 
#' 
#' The function SIMnormalmean.anydim simulates from a multivariate normal distribution and computes the arithmetic means per dimension as summary statistics with the goal to estimate the mean vector. The full variance covariance matrix can be specified (but not estimated).
#' 
#' SIMnormalmean.anydim uses the library mvtnorm. 
#' 
#' @param nk number of datasets
#' @param mu mean vector
#' @param sigma VC matrix
#' @param n number of observations
#' 
#' @export

SIMnormalmean.anydim = function(nk, mu, sigma, n){
  rmvnorm(nk, mu, sigma/n)
}