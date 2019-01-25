#ifndef NORMAL_H
#define NORMAL_H

#include <RcppArmadillo.h>
#include "Model.h"

using namespace arma;

class Normal: public Model{
  // child of Model
public:
  //arma::mat VC;
  // Function to extract the parameters of the fixed_parameters list
  // void set_fixed_parameters(List fixed_parameters){
  //   // arma::vec VC_vec = fixed_parameters("VC");
  //   // int dim = sqrt(VC_vec.n_elem);
  //   // VC = mat(&VC_vec, dim, dim, true);
  //   //arma::mat VC = as<arma::mat>(fixed_parameters["VC"]);
  //   VC = { {1, 0}, {0, 1} };
  // }
  // Simulation function
  arma::mat simulate(int nk) {
    return mvnrnd(Model::parameters, Model::fixed_parameters, nk).t();
  }
};

#endif

/* COMMENT
 * 
 */