#ifndef COALESCENT_THETA_H
#define COALESCENT_THETA_H

#include <RcppArmadillo.h>
#include "Model.h"

using namespace arma;

class Coalescent_theta: public Model{
  // child of Model
public:
  
  // Simulation function
  arma::mat simulate(int nk) {
    
    int n = as<int>(Model::fixed_parameters["n"]);
    int k = 2;
    
    // "choose from" vector (choose (n, 2), (n-1, 2), etc.)
    IntegerVector choose_from = seq(2, n);
    // make a numeric vector, because choose (below) can't use integer vectors
    NumericVector choose_from_num = as<NumericVector>(choose_from);
    // rate vector
    NumericVector rates = choose(choose_from_num, 2);
    
    // matrix of random branch lengths:
    NumericMatrix branches(n-1, nk);
    
    for(int i = 0; i < n-1; i++){
      branches(i,_) = rexp(nk, rates[i]);
    }
    
    // total branchlengths
    NumericVector total_branches = colSums(branches);
    
    // number of segregating sites
    NumericVector seg_sites(nk);
    
    for(int j = 0; j < nk; j++){
      seg_sites[j] = rpois(1, total_branches[j])[0];
    }
    
    return seg_sites;
    
    
  }
};

#endif

