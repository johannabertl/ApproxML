#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// NOTE: The values in H are the standard deviations of the kernel. 
// Take care that it is not the variance!

// [[Rcpp::export]]
arma::vec normal_diag(arma::mat dat, arma::vec x, arma::vec H) {
  int p = x.size();
  int n = dat.n_rows;
  arma::mat univar_dens(n,p);
  
  for(int i = 0; i < p; i++){
    //double h = sqrt(H[i]);
    univar_dens.col(i) = normpdf(dat.col(i), x[i], H[i]);
  }
  
  arma::vec res = prod(univar_dens, 1);
  
  return res;
}

/*** R

dat = matrix(rnorm(10), ncol=2, nrow=5)
x = c(0,0)
H = c(1,1)

normal_diag(dat, x, H)
apply(cbind(dnorm(x[1], dat[,1], sqrt(H[1])), dnorm(x[2], dat[,2], sqrt(H[2]))), 1, prod)
*/
