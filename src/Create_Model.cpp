#include "Model_Factory.h"
#include "Model.h"
#include "Normal.h"
#include <RcppArmadillo.h>
//#include <memory>

using namespace Rcpp;
using namespace arma;
using namespace std;

Model* Model_Factory::Create_Model(String type) {
  if ( type == "Normal" ) return new Normal();
  //  if ( type == "square" ) return new Square();
  return NULL;
}
