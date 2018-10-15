

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


#### Complex step derivative ####

# The function csd computes the complex step derivative of a kernel density estimator with Gaussian kernel and diagonal bandwidth matrix

# s ... simulated sum stats matrix (nk x d)
# s.obs ... observed sum stats
# theta ... simulated parameter matrix (nk x p)
# theta.k ... current theta value
# csd.const ... precision parameter of the complex step derivative
# H ... bandwidth matrix (diagonal!)

csd = function(s, s.obs, theta, theta.k, csd.const, H){

  d = length(s.obs)
  p = length(theta.k)
  nk = nrow(s)

  # part 1: normal kernel kde (nk x 1)
  kde.vec = apply(cbind(s, theta), 1, dmvnorm2, x=c(s.obs, theta.k), bw=H)

  # part 2: exponential part (p x 1)
  exp.vec = exp(csd.const/(2*diag(H)[(d+1):(d+p)]))

  # part 3: sinus (nk x p)
  sin.mat = sin(csd.const*(theta - matrix(theta.k, byrow=T, ncol=p, nrow=nk))/matrix(diag(H)[(d+1):(d+p)], byrow=T, ncol=p, nrow=nk))

  # sum and scaling (with nk and csd.const)
  colSums(kde.vec * sin.mat)*exp.vec/(nk*csd.const)

}


#### Kernel density derivative estimator ####

# The function kdde computes the derivative of a a kernel density estimator by the kernel density derivative method, i. e. the kernel is replaced by its derivative. Here, the normal kernel is used.
# NOTE: Only diagonal bandwidth matrices should be used! The implementation is incorrect for general bw matrices.

kdde = function(s, s.obs, theta, theta.k, H){

  d = length(s.obs)
  p = length(theta.k)
  nk = nrow(s)

  # part1: normal kernel of all dimensions on each point: (1 x nk)
  kde.vec = apply(cbind(s, theta), 1, dmvnorm2, x=c(s.obs, theta.k), bw=H)

  # part2: matrix of (internal) derivative vectors for the theta components (nk x p)
  theta.k.mat = matrix(theta.k, ncol=p, nrow=nk, byrow=T)
  H.mat = matrix(diag(H)[(d+1):(d+p)], ncol=p, nrow=nk, byrow=T)
  der.mat = (theta - theta.k.mat)/H.mat

  # gradient:
  (kde.vec %*% der.mat)/nk

}
