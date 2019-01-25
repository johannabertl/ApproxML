#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat SIMtestC(int nk, arma::vec theta, arma::mat sigma) {
  return mvnrnd(theta, sigma, nk).t();
}


// Test

/*** R
SIMtestC(10, c(0, 0), diag(2))
*/
