#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include <RcppArmadillo.h>
#include <Model.h>
#include <Normal.h>

using namespace arma;


class Model_Factory {
public:
  Model* Create_Model(String type) {
    if ( type == "Normal" ) return new Normal();
    //  if ( type == "square" ) return new Square();
    return NULL;
  }
};

#endif 