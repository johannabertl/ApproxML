#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec bw_nrd0(arma::mat x){
  
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec hi(k);
  arma::vec lo(k);
  
  for(int i = 0; i<k; i++){
    hi[i] = stddev(x.col(i));
    lo[i] = hi[i]; // extend!
  }
  
  return 0.9 * lo * pow(n*k, - 1/(k+4));
}

/* TO DO:
 * - handle incomplete cases?
 * - implement the complete formula (including IQR)
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
