#include <RcppArmadillo.h>
#include "bw.h"
#include "kernels.h"
#include "Model_Factory.h"
#include "Model.h"
#include "Normal.h"
#include "Coalescent_theta.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Approximate maximum likelihood algorithm based on Kiefer-Wolfowitz stochastic approximation
//'
//' This function approximates the maximum likelihood estimate of the multivariate parameter theta.
//' It requires a function to simulate data and compute summary statistics at each position of 
//' the parameter space. These simulations are used to obtain estimates of the likelihood in a 
//' stochastic approximation algorithm based on finite differences approximation of the gradient 
//' (Kiefer-Wolfowitz algorithm).
//'
//' @section Gain sequences:
//' The gain sequences are defined following Spall (2013):
//' \deqn{ a_k = \frac{a}{(k + A)^{\alpha}} \text{ and }
//' c_k = \frac{c}{(k + C)^{\gamma}}. }{ a_k = a/(k + A)^alpha and c_k = c/(k + C)^gamma. }
//' Usually, C>0 is only used if the algorithm is resumed. The default values alpha=1 and 
//' gamma=1/6 are given in Spall (2013).
//' 
//' @section Implementation:
//' The function is implemented using Rcpp and RcppArmadillo.
//'
//' @param s_obs Vector of observed summary statistics
//' @param theta_0 Vector of starting point
//' @param theta_min Vector of minimal values for theta
//' @param theta_max Vector of maximal values for theta
//' @param simfun name of the function which is used to simulate data (further arguments 
//' to simfun are given with ...)
//'
//' @param ce numeric
//' @param gamma numeric
//' @param C integer
//' @param a numeric
//' @param alpha numeric
//' @param A integer
//'
//' @param K number of steps
//' @param nk number of simulations per step
//' @param Hfun function which computes the bandwidth matrix. Can be one of the functions in bw.r
//' @param Hnum how often shall the bandwidth matrix be computed? Possible values: "once" or an integer indicating after how many iterations the bandwidth selection is to be renewed (e. g. 1 => bandwidth selection in each step)
//' @param kernel kernel function
//' @param lg logical. Use the log-likelihood instead of the likelihood?
//'
//'
//' @references
//'
//' \itemize{
//' \item Bertl, J.; Ewing, G.; Kosiol, C. & Futschik, A. Approximate Maximum Likelihood Estimation for Population Genetic Inference. Statistical Applications in Genetics and Molecular Biology, 2017, 16, p. 291-312
//' \item Kiefer, J. & Wolfowitz, J. Stochastic estimation of the maximum of a regression function. The Annals of Mathematical Statistics, 1952, 23, p. 462-466
//' \item Blum, J. R. Multidimensional stochastic approximation methods. The Annals of Mathematical Statistics, 1954, 25, 737-744
//' \item Spall, J. C. Introduction to Stochastic Search and Optimization: Estimation, Simulation and Control. Wiley, 2003
//' }
//'
//' @author Johanna Bertl
//'
//' @examples
//'
//' # bivariate normal distribution
//'
//' test = KDKW_FD_Rcpp(s_obs = c(1,2), theta_0 = c(0.5,1.5), theta_min = c(-5, -5), 
//' theta_max = c(5, 5), K = 3000, a = 10, ce = 2, nk = 10, simfun = "Normal",
//' fixed_parameters = list(VC = diag(2))) 
//' matplot(test$theta, t="l")
//'
//' @export
// [[Rcpp::export]]
List KDKW_FD_Rcpp(NumericVector s_obs, NumericVector theta_0, NumericVector theta_min, 
                           NumericVector theta_max, String simfun, List fixed_parameters,
                           int K, int nk,
                           double a, double ce, double alpha = 1, double gamma = 1/6,
                           int A = 0, int C = 0, bool use_log = true) {
  
  // Test if something is happening at all
  //Rcout << "hello";
  
  // Number of paramters
  int p = theta_0.size();
  
  // initialize Model objects
  Model_Factory* myFactory = new Model_Factory();
  Model* model_plus = myFactory->Create_Model(simfun);
  Model* model_minus = myFactory->Create_Model(simfun);
  model_plus->fixed_parameters = fixed_parameters;
  model_minus->fixed_parameters = fixed_parameters;
  
  // Empty matrix for the iterates
  NumericMatrix theta(K, p);
  theta(0,_) = theta_0;
  
  // data structures to record shifts into parameter space
  IntegerVector min_shift(K);
  IntegerVector max_shift(K);
  
  // gain sequences
  IntegerVector kvec = seq(1,K);
  DoubleVector ck = ce/(pow(kvec+1+C, gamma));
  DoubleVector ak = a/(pow(kvec+1+A, alpha));

  
  // Algorithm
  for(int k = 1; k < K; k++){
    
    
    // TESTING:
    // Rcout << k;
    // Rcout << "\n";
    
    // empty gradient vector
    DoubleVector gradient(p);
    
    
    // correction of theta[k-1,j] according to restrictions and ck:
    
    for(int j = 0; j < p; j++){
      if(theta(k-1,j) < theta_min[j] + ck[k]){
        theta(k-1,j) = theta_min[j] + ck[k];
        min_shift[k] = min_shift[k] + 1;
        // Rcout << "min corr \n";
      }
      if(theta(k-1,j) > theta_max[j] - ck[k]){
        theta(k-1,j) = theta_max[j] - ck[k];
        max_shift[k] = max_shift[k] + 1;
        // Rcout << "max corr \n";
      }
    }
    
    
    // computing gradient (by element)
    for(int j = 0; j < p; j++){
      
      
      arma::vec theta_plus = theta(k-1,_);
      theta_plus[j] = theta_plus[j] + ck[k];
      arma::vec theta_minus = theta(k-1,_);
      theta_minus[j] = theta_minus[j] - ck[k];
      // Is this stupid? In R, I used matrix manipulation outside this loop.
      // Is overwriting a vector problematic?
      
      model_plus->parameters = theta_plus;
      model_minus->parameters = theta_minus;
      arma::mat s_plus = model_plus->simulate(nk);
      arma::mat s_minus = model_minus->simulate(nk);    
      
      
      // compute bandwidth matrix
      // LATER: only every xth iteration, etc. 
          
      arma::vec H_plus = bw_nrd0(s_plus);
      arma::vec H_minus = bw_nrd0(s_minus);
      
      // H_plus.print("H_plus: ");
      // H_minus.print("H_minus: ");
        
      
      // compute KDE
      
      double kde_plus = mean(normal_diag(s_plus, s_obs, H_plus));
      double kde_minus = mean(normal_diag(s_minus, s_obs, H_minus));
      
      // Rcout << "\n";
      // Rcout << kde_plus;
      // Rcout << " - ";
      // Rcout << kde_minus;
      // Rcout << "\n";
      
      if(use_log){
        if(kde_plus == 0 | kde_minus == 0){
          stop("Kernel density estimate is 0, cannot compute logarithm.");
        } else {
          gradient[j] = (log(kde_plus) - log(kde_minus))/(2*ck[k]);
        }
      } else {
        gradient[j] = (kde_plus - kde_minus)/(2*ck[k]);
      }
    }
    
    // UPDATE
    theta(k,_) = theta(k-1,_) + ak[k]*gradient;
    
  }
  
  return List::create(_["theta"] = theta, _["min_shift"] = min_shift, _["max_shift"] = max_shift);
}


/*** R



# 5-dim normal distribution: estimating the mean vector
# Comparison between Rcpp and R function

t1 = proc.time()[3]
test = KDKW_FD_Rcpp(s_obs = 1:5, theta_0 = c(0.5, 1.5, 2.5, 3.5, 4.5), 
                    theta_min = rep(-5, 5), theta_max = rep(10, 5), simfun = "Normal",
                    fixed_parameters = list(VC = diag(5)),
                    K = 500, a = 500, ce = 2, nk = 10, use_log = F) 
proc.time()[3] - t1
matplot(test$theta, t="l")

t1 = proc.time()[3]
test = KDKW.FD(s.obs = 1:5, theta.0 = c(0.5, 1.5, 2.5, 3.5, 4.5), 
                    rest = matrix(c(-5, 10), nrow=5, ncol=2, byrow=T), simfun = SIMnormalmean.anydim,
                    sigma = diag(5), n = 1, Hnum = 1, kernel = dmvnorm2,
                    K = 500, a = 500, ce = 2, nk = 10, lg = F) 
proc.time()[3] - t1
matplot(test$theta, t="l")

# 5-dim normal distribution: 
# Starting too far from the maximum, where the likelihood is (numerically) zero.

t1 = proc.time()[3]
test = KDKW_FD_Rcpp(s_obs = 1:5, theta_0 = c(0.5,1.5, 2.5, 3.5, 4.5), 
                    theta_min = rep(-5, 5), theta_max = rep(10, 5), simfun = "Normal",
                    fixed_parameters = list(VC = diag(5)),
                    K = 3000, a = 500, ce = 2, nk = 10, use_log = T) 
  proc.time()[3] - t1
  matplot(test$theta, t="l")


# Coalescent: estimating theta
# Comparison between KDKW_FD_Rcpp and KDKW.FD.

t1 = proc.time()[3]
test2 = KDKW_FD_Rcpp(s_obs = 3, theta_0 = 1, 
                     theta_min = 0, theta_max = 10, simfun = "Coalescent_theta",
                     fixed_parameters = list(n = 5),
                     K = 1500, a = 10, ce = 1, nk = 10)
proc.time()[3] - t1
plot(test2$theta, t="l")

t1 = proc.time()[3]
test2 = KDKW.FD(s.obs = 3, theta.0 = 1, 
                  rest = matrix(c(0, 10), nrow=1), simfun = SIMcoal.theta,
                     n = 5, Hnum = 1, kernel = dnorm2,
                     K = 1500, a = 10, ce = 1, nk = 10)
proc.time()[3] - t1
plot(test2$theta, t="l")

*/

/* missing features: 
 * - log likelihood
 * - different bandwidth functions
 */

/* improvements: 
 * - use .at() instead of [] (maybe not for Armadillo structures? Find out!)
 * - change the fixed_parameters List to a better c++ datatype, to avoid using as<..>(..) each 
 * time simulate() is called.
 * - allow infinity in the constraints?
 * - handle missing values?
 * - use Rcpp::checkUserInterrupt() 
 */

/* Try to catch the following errors: 
 * - nk = 1 needs a fixed bw. Doesn't work right now. Throw an error.
 * - make sure that theta_0, theta_min and theta_max have the same length. 
 * - make sure that s_obs and the output of the simulation function have the 
 * same dimension.
 * - make sure that the parameters and fixed_parameters are correct (e. g. for normal 
 * distribution, check if dimensions fit -- implement check_parameters method)
 * - check constraints and ck[1] before starting the algorithm
 */

/* Before uploading to CRAN, read: https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Plummer.pdf 
 * 
 */