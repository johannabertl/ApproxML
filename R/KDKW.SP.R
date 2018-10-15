#' Approximate maximum likelihood algorithm based on stochastic approximation with simultaneous perturbations
#'
#' This function approximates the maximum likelihood estimate of the multivariate parameter theta. It requires a function to simulate data and compute summary statistics at each position of the parameter space. These simulations are used to obtain estimates the likelihood in a stochastic approximation algorithm based where the ascent direction is given by simultaneous perturbations.
#'
#' @inheritParams KDKW.FD
#'
#' @references
#'
#' \itemize{
#' \item Bertl, J.; Ewing, G.; Kosiol, C. & Futschik, A. Approximate Maximum Likelihood Estimation for Population Genetic Inference. Statistical Applications in Genetics and Molecular Biology, 2017, 16, p. 291-312
#' \item Spall, J. C. Multivariate stochastic approximation using a simultaneous perturbation gradient approximation. IEEE Transactions on Automatic Control, 1992, 37, 352-355
#' \item Spall, J. C. Introduction to Stochastic Search and Optimization: Estimation, Simulation and Control. Wiley, 2003
#' }
#'
#' @author Johanna Bertl
#'
#' @examples
#'
#' # bivariate normal distribution
#'
#' test1 = KDKW.SP(s.obs=c(0,0), theta.0=c(2,2), simfun=SIMnormalmean.anydim, rest=matrix(c(-100, 100), ncol=2, nrow=2, byrow=T), ce=2, gamma = 1/6, C=0, a = 10, alpha = 1, A=10, K=500, nk=50, Hfun = bw.nrd0.mult, Hnum="once", dmvnorm2, lg = T, sigma=diag(c(1,1)), n=1)
#' matplot(test1$theta, t="l")
#'
#' @export


KDKW.SP = function(s.obs, theta.0, simfun, rest=matrix(c(-Inf, Inf), ncol=2, nrow=length(theta.0), byrow=T),
  ce, gamma = 1/6, C=0, a, alpha = 1, A = 0, K, nk,
  Hfun = bw.nrd0.flex, Hnum, kernel,
  lg = T,
  lib = "~/Rpack",
  ...){

  t1=proc.time()[3]

  #### data structures ####

  p = length(theta.0)

  # empty matrix for the theta values:
  theta = matrix(NA, nrow = K, ncol = p)
  theta[1,] = theta.0

  # data structures to record shifts into the parameterspace
  adaptcount=0
  adapt.l = matrix(0, nrow=K, ncol=p)
  adapt.u = matrix(0, nrow=K, ncol=p)

  # gain sequences
 	kvec = 1:K
	ck = ce/((kvec+1+C)^gamma)
	ak = a/((kvec+1+A)^alpha)


  #### Algorithm ####

  for (k in 2:K){

    gradient = numeric(p)

    # random direction:
    delta = rbinom(p, 1, 0.5)*2 - 1

    theta.plus = theta[k-1,] + ck[k]*delta
    theta.minus =theta[k-1,] - ck[k]*delta

    # corrections of theta.plus, theta.minus and theta[k-1,]according to restrictions.

    for(j in 1:p){

        theta.max = max(theta.plus[j], theta.minus[j])
        theta.min = min(theta.minus[j], theta.plus[j])

        if ((theta.min < rest[j,1]) & (theta.max > rest[j,2])){
            stop("ERROR 1: Specification of boundaries and tuning sequences does not match.")
        }

        if (theta.min < rest[j,1]){
            d = rest[j,1] - theta.min
            theta.minus[j] = theta.minus[j] + d
            theta.plus[j] = theta.plus[j] + d
            theta[k-1,j] = theta[k-1,j] + d
            adaptcount = adaptcount+1
            adapt.l[k,j] = 1
            theta.max = theta.max + d
            if (theta.max > rest[j,2]){
                stop("ERROR 2: Specification of boundaries and tuning sequences does not match.")
            }

        }

        if (theta.max > rest[j,2]){
            d = theta.max - rest[j,2]
            theta.plus[j] = theta.plus[j] - d
            theta.minus[j] = theta.minus[j] - d
            theta[k-1,j] = theta[k-1,j] - d
            adaptcount = adaptcount+1
            adapt.u[k,j] = 1
            theta.min = theta.min - d
            if (theta.min < rest[j,1]){
                stop("ERROR 3: Specification of boundaries and tuning sequences does not match.")
            }
        }
    }

    # estimation of the gradient:

    # simulate summary statistics:

    s.plus = simfun(nk, theta.plus, ...)
    s.minus = simfun(nk, theta.minus, ...)

    # compute bandwidth matrix

    if (Hnum == "once"){
    	if(k == 2){

 		    H.plus = Hfun(s.plus)
        H.minus = Hfun(s.minus)

    	}
    } else {
      if(((k-1) %% Hnum) == 1 | Hnum == 1){

 		    H.plus = Hfun(s.plus)
        H.minus = Hfun(s.minus)


      }
   }


    # kernel density estimation (the function kdenorm from kde.r is used)

    dens.plus = kde(s.plus, s.obs, H.plus, kernel)
    dens.minus = kde(s.minus, s.obs, H.minus, kernel)


    if (lg==F){
      gradient = (dens.plus - dens.minus)/(2*ck[k]*delta)
   	} else {
    	if (dens.plus==0 | dens.minus==0)
      		{
      			stop(paste("ERROR 4: log(0) gives infinite gradient in step",k))
      		}
        gradient = (log(dens.plus) - log(dens.minus))/(2*ck[k]*delta)
      }


    # UPDATE #

    theta[k,] = theta[k-1,] + ak[k]*gradient

  } # end of k in 1:K loop


  #### Output ####

  t2 = proc.time()[3]
  t.total = as.numeric(t2 - t1)

  # Output list:
  list(prtime = t.total,
       a = a, ce = ce, type="SP",
      adapt = adaptcount, adapt.l = adapt.l, adapt.u = adapt.u,
      theta = theta)
}


###########################################################################

### Functions for apply (and mclapply):

#' @describeIn KDKW.SP The same function for a set of different starting values, use e. g. with apply or mclapply.
#' @export

KDKW.SP.theta.0 = function(theta.0, s.obs, simfun, rest,
  ce, gamma, C, a, alpha, A, K, nk,
  Hfun, Hnum,
  lg,
  lib,
  ...){

	try(KDKW.SP(s.obs=s.obs, theta.0=theta.0, simfun=simfun, rest=rest,
  ce=ce, gamma=gamma, C=C, a=a, alpha=alpha, A=A, K=K, nk=nk,
  Hfun = Hfun, Hnum=Hnum,
  lg = lg,
  lib = lib,
  ...), silent=TRUE)

}

#' @describeIn KDKW.SP The same function for a set of different summary statistics, use e. g. with apply or mclapply.
#' @export

KDKW.SP.s.obs = function(s.obs, theta.0, simfun, rest,
  ce, gamma, C, a, alpha, A, K, nk,
  Hfun, Hnum,
  lg,
  lib,
  ...){

	try(KDKW.SP(s.obs=s.obs, theta.0=theta.0, simfun=simfun, rest=rest,
  ce=ce, gamma=gamma, C=C, a=a, alpha=alpha, A=A, K=K, nk=nk,
  Hfun = Hfun, Hnum=Hnum,
  lg = lg,
  lib = lib,
  ...), silent=TRUE)

}

#' @describeIn KDKW.SP The same function for a set of different starting values and different summary statistics, use e. g. with apply or mclapply.
#' @export

KDKW.SP.mult1 = function(mult, simfun, rest,
                         ce, gamma, C, a, alpha, A, K, nk,
                         Hfun, Hnum,
                         lg,
                         lib,
                         ...){

  try(KDKW.SP(s.obs=mult[[1]], theta.0=mult[[2]], simfun=simfun, rest=rest,
              ce=ce, gamma=gamma, C=C, a=a, alpha=alpha, A=A, K=K, nk=nk,
              Hfun = Hfun, Hnum=Hnum,
              lg = lg,
              lib = lib,
              ...), silent=TRUE)

}
