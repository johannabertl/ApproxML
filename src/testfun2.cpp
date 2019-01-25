#include <RcppArmadillo.h>
#include "Model.h"
#include "Normal.h"
#include "Model_Factory.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat testfun2(arma::vec parameters, arma::mat fixed_parameters){
  Normal std_norm;
  std_norm.parameters = parameters;
  std_norm.fixed_parameters = fixed_parameters;
  return std_norm.simulate(10);
}




/*** R
testfun2(c(0,0), diag(2))
*/
