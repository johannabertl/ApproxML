#' Approximate maximum likelihood algorithm based on Kiefer-Wolfowitz stochastic approximation with common random numbers
#'
#' This function approximates the maximum likelihood estimate of the multivariate parameter theta. It requires a function to simulate data given a random seed and compute summary statistics at each position of the parameter space. These simulations are used to obtain estimates the likelihood in a stochastic approximation algorithm based on finite differences approximation of the gradient (Kiefer-Wolfowitz algorithm). The random seed is used to generate common random numbers (CRNs) to reduce the variance of the gradient estimates and thereby the final parameter estimate.
#'
#' Note that currently in each iteration, a single integer valued seed is passed on to the simulation function for simulation at theta+ and theta-. This might not be enough for complex simulations and should probably be replaced by an implementation e. g. using setRNG::setRNG that sets a vector of 264 integers as seed, and maybe a flexible number of seeds, if multiple random samples need to be generated.
#'
#' The simulation function needs to have an argument "seed". At the moment, this is not implemented for most simulation functions.
#'
#' @section Gain sequences:
#' The gain sequences are defined following Spall (2013):
#' \deqn{ a_k = \frac{a}{(k + A)^{\alpha}} \text{ and }
#' c_k = \frac{c}{(k + C)^{\gamma}}. }{ a_k = a/(k + A)^alpha and c_k = c/(k + C)^gamma. }
#' Usually, C>0 is only used if the algorithm is resumed. The default values alpha=1 and gamma=1/6 are given in Spall (2013).
#'
#' @param s.obs Vector of observed summary statistics
#' @param theta.0 Vector of starting point
#' @param simfun name of the function which is used to simulate data (further arguments to simfun are given with ...). Note that simfun needs to have an argument "seed" for common random number generation.
#' @param rest Matrix with restrictions on the parameterspace (matrix with 2 columns that contain the lower and upper bounds for the elements of theta). Use -Inf and Inf for unrestricted parameters. The default value is an unrestricted parameter space.
#'
#' @param ce numeric
#' @param gamma numeric
#' @param C integer
#' @param a numeric
#' @param alpha numeric
#' @param A integer
#'
#' @param K number of steps
#' @param nk number of simulations per step
#' @param Hfun function which computes the bandwidth matrix. Can be one of the functions in bw.r
#' @param Hnum how often shall the bandwidth matrix be computed? Possible values: "once" or an integer indicating after how many iterations the bandwidth selection is to be renewed (e. g. 1 => bandwidth selection in each step)
#' @param kernel kernel function
#' @param lg logical. Use the log-likelihood instead of the likelihood?
#'
#' @param ... further arguments which are passed on to simulate the data with simfun
#'
#' @references
#'
#' \itemize{
#' \item Bertl, J.; Ewing, G.; Kosiol, C. & Futschik, A. Approximate Maximum Likelihood Estimation for Population Genetic Inference. Statistical Applications in Genetics and Molecular Biology, 2017, 16, p. 291-312
#' \item Kiefer, J. & Wolfowitz, J. Stochastic estimation of the maximum of a regression function. The Annals of Mathematical Statistics, 1952, 23, p. 462-466
#' \item Blum, J. R. Multidimensional stochastic approximation methods. The Annals of Mathematical Statistics, 1954, 25, 737-744
#' \item Spall, J. C. Introduction to Stochastic Search and Optimization: Estimation, Simulation and Control. Wiley, 2003
#' }
#'
#' @author Johanna Bertl
#'
#' @examples
#'
#' # Comparison between the algorithm with common random numbers and the standard version with independent random numbers on a simple Poisson GLMM with a single fixed mean and a random intercept.
#'
#' test1 = KDKW.FD.CRN(s.obs=c(3.5, 10), theta.0=c(1,0.5), simfun=SIMCRNpoisson_glmm_simple, rest=matrix(c(-1, 5, 0.1, 3), ncol=2, nrow=2, byrow=T), ce=1, gamma = 1/6, C=0, a = 2, alpha = 1, A=0, K=200, nk=20, Hfun = bw.nrd0.flex, Hnum="once", kernel=robust.unscaled.diagonal, lg = T, n=100)
#'
#' test_noCRN = KDKW.FD(s.obs=c(3.5, 10), theta.0=c(1,0.5), simfun=SIMpoisson_glmm_simple, rest=matrix(c(-1, 5, 0.1, 3), ncol=2, nrow=2, byrow=T), ce=1, gamma = 1/6, C=0, a = 2, alpha = 1, A=0, K=200, nk=20, Hfun = bw.nrd0.flex, Hnum="once", kernel=robust.unscaled.diagonal, lg = T, n=100)
#'
#' par(mfrow=c(2,1))
#' matplot(test1$theta, t="l", ylim=c(0.5, 1.5))
#' matplot(test_noCRN$theta, t="l", ylim=c(0.5, 1.5))
#' par(mfrow=c(1,1))
#'
#' @export

KDKW.FD.CRN = function(s.obs, theta.0, simfun, rest=matrix(c(-Inf, Inf), ncol=2, nrow=length(theta.0), byrow=T),
  ce, gamma = 1/6, C = 0, a, alpha = 1, A = 0, K, nk,
  Hfun = bw.nrd0.mult, Hnum, kernel,
  lg = T,
  ...){

  t1=proc.time()[3]

  #### data structures ####

  p = length(theta.0)

  # matrix of p unit vectors
  E = diag(1, p)

  # empty matrix for the theta values:
  theta = matrix(NA, nrow = K, ncol = p)
  theta[1,] = theta.0

  # data structures to record shifts into parameter space
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

    # matrices with theta.plus and theta.minus values in the diagonals,
    # the other values are copied from theta_k-1

    theta.plus = matrix(theta[k-1,], ncol=p, nrow=p, byrow=T) + ck[k]*E
    theta.minus = matrix(theta[k-1,], ncol=p, nrow=p, byrow=T) - ck[k]*E

    # corrections of theta.plus, theta.minus and theta[k-1,]according to restrictions:

    for(j in 1:p){

        if (theta.minus[j,j]<rest[j,1] & theta.plus[j,j]>rest[j,2]){
            stop("ERROR 1: Specification of boundaries and tuning sequences does not match.")
        }

        if (theta.minus[j,j]<rest[j,1]){
            d = rest[j,1] - theta.minus[j,j]
            theta.minus[,j] = theta.minus[,j] + d
            theta.plus[,j] = theta.plus[,j] + d
            theta[k-1,j] = theta[k-1,j] + d
            adaptcount = adaptcount+1
            adapt.l[k,j] = 1
            if (theta.plus[j,j]>rest[j,2]){
                stop("ERROR 2: Specification of boundaries and tuning sequences does not match.")
            }

        }

        if (theta.plus[j,j]>rest[j,2]){
            d = theta.plus[j,j] - rest[j,2]
            theta.plus[,j] = theta.plus[,j] - d
            theta.minus[,j] = theta.minus[,j] - d
            theta[k-1,j] = theta[k-1,j] - d
            adaptcount = adaptcount+1
            adapt.u[k,j] = 1
            if (theta.minus[j,j] < rest[j,1]){
                stop("ERROR 3: Specification of boundaries and tuning sequences does not match.")
            }
        }
    }

    # estimation of the gradient:

    for(j in 1:p){

        # seed for CRNs:
        seed = sample.int(2^30, 1)

        # simulate summary statistics:

        s.plus = simfun(nk, theta.plus[j,], seed, ...)
        s.minus = simfun(nk, theta.minus[j,], seed, ...)

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


        # kernel density estimation

        dens.plus = kde(s.plus, s.obs, H.plus, kernel)
        dens.minus = kde(s.minus, s.obs, H.minus, kernel)


      	if (lg==F){
        	gradient[j] = (dens.plus - dens.minus)/(2*ck[k])
      	} else {
      		if (dens.plus==0 | dens.minus==0)
      			{
      				stop(paste("ERROR 4: log(0) gives infinite gradient in step",k))
      			}
      		gradient[j] = (log(dens.plus) - log(dens.minus))/(2*ck[k])
      	}


    } # end of j in 1:p loop


    # UPDATE #

    theta[k,] = theta[k-1,] + ak[k]*gradient

  } # end of k in 1:K loop

  #### Output ####

  t2 = proc.time()[3]
  t.total = as.numeric(t2 - t1)

  # Output list:
  list(prtime = t.total,
      adapt = adaptcount, adapt.l = adapt.l, adapt.u = adapt.u,
      theta = theta,
      ce = ce, a = a, type = "FD")
}


###########################################################################

### Functions for apply (and mclapply):

#' @describeIn KDKW.FD.CRN The same function for a set of different starting values, use e. g. with apply or mclapply. It uses \code{try} catch errors.
#' @export

KDKW.FD.CRN.theta.0 = function(theta.0, s.obs, simfun, rest=matrix(c(-Inf, Inf), ncol=2, nrow=length(theta.0), byrow=T),
  ce, gamma=1/6, C=0, a, alpha=1, A=0, K, nk,
  Hfun = bw.nrd0.mult, Hnum, kernel,
  lg=T,
  ...){

	try(KDKW.FD.CRN(s.obs=s.obs, theta.0=theta.0, simfun=simfun, rest=rest,
  ce=ce, gamma=gamma, C=C, a=a, alpha=alpha, A=A, K=K, nk=nk,
  Hfun = Hfun, Hnum=Hnum, kernel=kernel,
  lg = lg,
  ...), silent=TRUE)

}


#' @describeIn KDKW.FD.CRN The same function for a set of different summary statistics, use e. g. with apply or mclapply. It uses \code{try} catch errors.
#' @export

KDKW.FD.CRN.s.obs = function(s.obs, theta.0, simfun, rest=matrix(c(-Inf, Inf), ncol=2, nrow=length(theta.0), byrow=T),
  ce, gamma=1/6, C=0, a, alpha=1, A=0, K, nk,
  Hfun = bw.nrd0.mult, Hnum, kernel,
  lg=T,
  ...){

	try(KDKW.FD.CRN(s.obs=s.obs, theta.0=theta.0, simfun=simfun, rest=rest,
  ce=ce, gamma=gamma, C=C, a=a, alpha=alpha, A=A, K=K, nk=nk,
  Hfun = Hfun, Hnum=Hnum, kernel=kernel,
  lg = lg,
  ...), silent=TRUE)

}
