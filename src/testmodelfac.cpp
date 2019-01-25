#include <RcppArmadillo.h>
#include "Model_Factory.h"
#include "Model.h"
#include "Normal.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat testmodelfac(String type, arma::vec parameters, arma::mat fixed_parameters)
{
  Model_Factory* myFactory = new Model_Factory();
  Model* std_norm = myFactory->Create_Model(type);

  std_norm->parameters = parameters;
  std_norm->fixed_parameters = fixed_parameters;
  arma::mat simres = std_norm->simulate(10);
  return simres;
  //return parameters;
}


/*** R
testmodelfac("Normal", c(0,0), diag(2))
*/
