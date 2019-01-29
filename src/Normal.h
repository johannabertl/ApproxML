#ifndef NORMAL_H
#define NORMAL_H

#include <RcppArmadillo.h>
#include "Model.h"

using namespace arma;

class Normal: public Model{
  // child of Model
public:
  // Simulation function
  arma::mat simulate(int nk) {
    arma::mat VC = as<arma::mat>(Model::fixed_parameters["VC"]);
    return mvnrnd(Model::parameters, VC, nk).t();
  }
};

#endif

/* COMMENT
 * 
 */