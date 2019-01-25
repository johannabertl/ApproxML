#ifndef MODEL_H
#define MODEL_H

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

Model* Model_Factory::Create_Model(String type);

#endif