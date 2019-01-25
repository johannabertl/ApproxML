

# This file contains multivariate kernel density estimation functions, that evaluate the kde in a single point.


#### Kernel density estimation ####

# The function is faster than kde from ks (see ../Version5/testkde.r for the comparison for the Gaussian kernel).

# The kde is only evaluated at a single value. A bandwidth matrix needs to be supplied.

# dat ... data matrix
# x ... point where the kde is evaluated
# bw ... bandwidth matrix
# kernel ... kernel function (a few kernel functions are given below)

#' @export

kde = function(dat, x, bw, kernel){
	sum(apply(as.matrix(dat), 1, FUN=kernel, x=x, bw=bw), na.rm=T)/sum(apply(!is.na(as.matrix(dat)), 1, all))
}


#### kernel functions ####

# multivariate Gaussian kernel #
# only works for more than 1 dimension! #
# require(mvtnorm)

#' @export

dmvnorm2 = function(mean, x, bw){dmvnorm(x, mean, sigma=bw)}

# univariate Gaussian kernel #

#' @export
dnorm2 = function(mean, x, bw){dnorm(x, mean, bw)}

# Alternative kernel for robustification #
# unscaled, i. e. not scaled to total weight 1 (for speed)

#' @export
robust.unscaled = function(mean, x, bw){
  squared = t(x-mean) %*% solve(bw) %*% (x-mean)
  if(squared<=1) {exp(-0.5*squared)} else {exp(-0.5*sqrt(squared))}
}

# Alternative kernel for robustification #
# unscaled, i. e. not scaled to total weight 1 (for speed)
# requires the bandwidth matrix to be diagonal! (that's no problem, only diagonal bandwidth matrices are implemented in bw.r) probably faster than the above version that uses "solve" (test this!).

#' @export
robust.unscaled.diagonal = function(mean, x, bw){
  squared = t(x-mean) %*% diag(1/diag(bw)) %*% (x-mean)
  if(squared<=1) {exp(-0.5*squared)} else {exp(-0.5*sqrt(squared))}
}


# t-distribution kernel with 5 df#
# (The scale matrix sigma is set equal to the bandwidth function - this is correct, I checked with the density function of the multivariate t-distribution if we get K_H(x) = |H|^(-1/2)*K(H^(-1/2)x).)
# require(mvtnorm)

#' @export
dmvt2.df5 = function(mean, x, bw){dmvt(x, delta=mean, sigma=bw, df=5, log=F)}

# t-distribution kernel with 3 df#
# (see comment above)
# require(mvtnorm)

#' @export
dmvt2.df3 = function(mean, x, bw){dmvt(x, delta=mean, sigma=bw, df=3, log=F)}

