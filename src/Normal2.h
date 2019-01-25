#ifndef NORMAL2_H
#define NORMAL2_H

#include <RcppArmadillo.h>
#include "Model.h"

using namespace arma;

class Normal2: public Model{
  // child of Model
public:
  arma::mat VC;
  // Simulation function
  arma::mat simulate(int nk) {
    return mvnrnd(Model::parameters, Model::fixed_parameters, nk).t();
  };
};

#endif

/* COMMENT
 * 
 */