# function to simulate data from 2 univariate normal distributions with fixed variances and compute summary statistics
# (goal: estimate the 2 means)

# summary statistics: mean of ALL values, product of mean1 and mean2

# nk ... number of simulations
# mu ... 2-vector of the mean
# sigma ... 2-vector of the variances
# n ... number of observations

SIMnormal.mult2a = function(nk, mu, sigma = c(1,1), n){
	
	x1 = matrix(rnorm(nk*n, mu[1], sqrt(sigma[1])), ncol = nk, nrow = n)
	x2 = matrix(rnorm(nk*n, mu[2], sqrt(sigma[2])), ncol = nk, nrow = n)
	
	mean1 = apply(x1, 2, mean)
	mean2 = apply(x2, 2, mean)
	
	sumstat = cbind(mean1+mean2, mean1*mean2)
	sumstat		
}