#' Likelihood estimation for starting points
#'
#' This function draws starting points from a given matrix of intervals (randomly samples from a uniform distribution or on a grid). For each starting point, the likelihood is simulated.
#'
#' For points on a grid, the number of starting points N1 is approximate. For each dimension, round(N1^(1/p)) equally spaced points in the interval are chosen, with p being the number of dimensions.
#'
#' Note that this could easily be parallelized. This is currently not implemented, since it depends on the system the program is run on.
#'
#' @param int matrix of intervals (2 columns, p rows)
#' @param sampling character. One of "random" or "grid"
#' @param N1 number of starting points
#' @param s.obs observed summary statistics
#' @param simfun name of the function which is used for the simulation and computes the summary statistics
#' @param Hfun function which computes the bandwidth matrix. Can be one of the functions in bw.R
#' @param kernel kernel function for KDE
#' @param nk number of simulations for the estimation of the likelihood with kernel density estimation
#' @param ... further parameters for simfun
#'
#' @examples
#'
#' int = matrix(c(-2, 5, 0, 2), ncol=2, nrow=2, byrow=T)
#'
#' test_start_grid = STARTrandomlik(int, sampling="grid", N1=20, s.obs = c(2.2, 9), simfun = SIMpoisson_glmm, kernel = robust.unscaled.diagonal, nk=20, x = 1:3)
#'
#' test_start_random = STARTrandomlik(int, sampling="random", N1=20, s.obs = c(2.2, 9), simfun = SIMpoisson_glmm, kernel = robust.unscaled.diagonal, nk=20, x = 1:3)
#'
#' @export


STARTlik = function(int, sampling, N1, s.obs, simfun, Hfun=bw.nrd0.mult, kernel, nk, ...){

  t1=proc.time()[3]

  #### Data structures ####

  # dimension of parameterspace
  p = nrow(int)

  # matrix of starting values (1 column per parameter)

  int.list = apply(int, 1, list)
  int.list = lapply(int.list, unlist)

  if(sampling=="random"){
    val.list = lapply(int.list, runif2, N1)
    val.mat = matrix(unlist(val.list), ncol=p, nrow=N1)

  } else {

    val.list = lapply(int.list, seq2, round(N1^(1/p)))
    val.mat = as.matrix(expand.grid(val.list))

  }

  start.list = apply(val.mat, 1, list)
  start.list = lapply(start.list, unlist)


  #### Likelihood computation ####

  lik.list = lapply(start.list, liksim, simfun=simfun, nk=nk, Hfun=Hfun, kernel=kernel, s.obs=s.obs, ...)
  lik = unlist(lik.list)


  #### Output ####

  t2 = proc.time()[3]
  t.total = t2 - t1

  list(prtime = as.numeric(t.total), order = order(lik, decreasing=T), lik = lik, start = val.mat)

}

# functions for apply:

runif2 = function(m, n){runif(n, m[1], m[2])}
seq2 = function(int, length.out) {seq(from=int[1], to = int[2], length.out = length.out)}

liksim = function(theta, simfun, nk, Hfun, kernel, s.obs, ...){
  sk = simfun(nk, theta, ...)

  Hk = Hfun(sk)
  kde(sk, s.obs, Hk, kernel)
}
