#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H


#include <RcppArmadillo.h>
#include "Model.h"
#include "Normal.h"
#include "Coalescent_theta.h"

using namespace arma;


class Model_Factory {
public:
  Model* Create_Model(String type) {
    if ( type == "Normal" ) return new Normal();
    if ( type == "Coalescent_theta" ) return new Coalescent_theta();
    return NULL;
  }
};

#endif 