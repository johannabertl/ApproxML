#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include "Model.h"
#include <RcppArmadillo.h>

using namespace arma;

class Model_Factory {
public:
  Model* Create_Model(String type);
};

#endif 