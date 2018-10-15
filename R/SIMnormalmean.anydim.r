# function to simulate data from a univariate normal distribution with 
# parameter mu and fixed Sigma and compute summary statistics from it. 
# summary statistics: arithmetic means.
# uses the library mvtnorm

SIMnormalmean.anydim = function(nk, mu, sigma, n){
  rmvnorm(nk, mu, sigma/n)
}