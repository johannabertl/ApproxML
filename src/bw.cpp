#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Estimating the bandwidth
//' 
//' The function bw_nrd0 estimates the entries of a diagonal bandwidth matrix, based on
//' the multivariate extension of "Silverman's rule of thumb". 
//' 
//' Note that the implementation is not the same as bw.nrd0 in R, which is more robust. 
//' In some cases, bw.nrd0 uses the IQR or sets the bandwidth to 1. This is not 
//' implemented here to save runtime.

// [[Rcpp::export]]
arma::vec bw_nrd0(arma::mat x){
  
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec hi(k);
  //arma::vec lo(k);
  
  for(int i = 0; i<k; i++){
    hi[i] = stddev(x.col(i));
    //lo[i] = hi[i]; // extend!
  }
  
  //return 0.9 * lo * pow(n*k, - 1/(k+4));
  return 0.9 * hi * pow(n*k, - 1/(k+4));
}

/* TO DO:
 * - handle incomplete cases?
 * - implement the complete formula (including IQR)? Comments on this: 
 *    - quantile() is not implemented in Rcpp. In armadillo, there is no 
 *    implementation either, partial_sort could be used (I've used that in one of the 
 *    exercises from Hadley's book, see /home/au527055/Documents/ThemenAarhusMath/AML/learn_Rcpp/5_The_STL3_exercises.cpp). 
 *    In boost, there is an implementation (https://stats.stackexchange.com/questions/7358/c-libraries-for-statistical-computing),
 *    but I don't want to introduce further dependencies. 
 *    - Maybe it is faster to use the R quantile function instead of implementing
 *    something myself. Try this!
 */

/*** R
set.seed(1234)
x = matrix(rnorm(20), ncol=2, nrow=10)
bw_nrd0(x)

# bw.nrd0
# function (x) 
# {
#   hi <- sd(x)
#   if (!(lo <- min(hi, IQR(x)/1.34))) 
#     (lo <- hi) || (lo <- abs(x[1L])) || (lo <- 1)
#   0.9 * lo * length(x)^(-0.2)
# }

*/
