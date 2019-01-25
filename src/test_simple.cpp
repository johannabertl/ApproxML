#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

  arma::mat set_fixed_parameters(List fixed_parameters){
    arma::mat VC = fixed_parameters("VC");
    // int dim = sqrt(VC.size());
    // return VC.reshape(dim); 
    return VC;
  }



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set_fixed_parameters(list(VC=diag(2)))
*/
