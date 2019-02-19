#' Simulation function for a simple coalescent model
#' 
#' The function simulates a coalescent tree and mutations and calculates the number of segregating sites under the infinite sites assumption. The number of segregating sites are used as summary statistics. 
#' The only parameters of the coalescent model are the scaled mutation rate theta and the number of samples, n. 
#' 
#' @param nk Scalar. Number of simulated data sets (trees). 
#' @param theta Scalar. Scaled mutation rate.
#' @param n Scalar. Sample size.
#' 
#' @author Johanna Bertl
#' 
#' @examples 
#' 
#' # Simulating 10 coalescent trees with scaled mutation rate 5 and 20 branches each. The number of segregating sites are: 
#' 
#' SIMcoal.theta(10, 5, 20)
#' 
#' # Estimation with KDKW.FD: 
#' test2 = KDKW.FD(s.obs = 3, theta.0 = 1, rest = matrix(c(0, 10), nrow=1), simfun = SIMcoal.theta, K = 500, a = 10, ce = 1, nk = 10, n = 5, Hnum = 1, kernel = dnorm2)
#' plot(test2$theta, t="l")
#' 
#' @export

SIMcoal.theta = function(nk, theta, n){
  
  # tree length: 
  # simulate overall tree length 
  h = 2:n
  lambda = choose(h, 2) 
  lambda.mat = matrix(lambda, nrow = n-1, ncol=nk, byrow = F)
  l.mat = mapply(rexp2, lambda.mat)  
  l.mat = matrix(l.mat,nrow = n-1, ncol = nk, byrow = F)
  l.mat2 = l.mat*h
  l.vec = apply(l.mat2, 2, sum) # nk-vector of overall tree lengths
  
  # number of segregating sites:
  lambda = (theta/2)*l.vec
  seg = apply(as.matrix(lambda), 1, rpois2)
  
  return(seg)
}


### Functions used in the above function:

# function to simulate exponentially distributed random variables
# (for easier use of rexp in mapply):
rexp2 = function(lambda){rexp(n = 1,rate = lambda)}

# for easier use in apply
rpois2 = function(lambda){rpois(n = 1,lambda)}