#ifndef MODEL_H
#define MODEL_H

#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Model {
public:
  arma::vec parameters;
  List fixed_parameters;
  //arma::mat data;
  static Model* Create_Model(String type);
  // virtual void set_fixed_parameters(List fixed_parameters);
  virtual arma::mat simulate(int nk) = 0;
};

#endif
