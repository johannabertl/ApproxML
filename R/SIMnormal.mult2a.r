#' Simulation function for bivariate normal distribution with experimental summary statistics
#' 
#' SIMnormal.mult2a simulates data from a bivariate normal distributions with a fixed diagonal variance covariance matrix and computes summary statistics to estimate the mean vector. 
#' 
#' Summary statistics: mean of ALL values, product of mean1 and mean2
#' 
#' @param nk number of simulations
#' @param mu 2-vector of the mean
#' @param sigma 2-vector of the variances
#' @param n number of observations

#' @export

SIMnormal.mult2a = function(nk, mu, sigma = c(1,1), n){
	
	x1 = matrix(rnorm(nk*n, mu[1], sqrt(sigma[1])), ncol = nk, nrow = n)
	x2 = matrix(rnorm(nk*n, mu[2], sqrt(sigma[2])), ncol = nk, nrow = n)
	
	mean1 = apply(x1, 2, mean)
	mean2 = apply(x2, 2, mean)
	
	sumstat = cbind(mean1+mean2, mean1*mean2)
	sumstat		
}