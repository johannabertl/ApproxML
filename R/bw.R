# Functions to compute the bandwidth matrix


### 1a. Flexible Hpi.diag function with 1 stage estimation

# Warning: This function requires the package ks.
# It can only be used for up to 6 dimensions.

#' @export

Hpi.diag.flex1 = function(x){
  if(is.matrix(x)){Hpi.diag(x, nstage = 1)}
  else {hpi(x)}
}

### 1b. Flexible Hpi.diag function with 2 stage estimation

# Warning: This function requires the package ks.
# It can only be used for up to 6 dimensions.

#' @export

Hpi.diag.flex2 = function(x){
  if(is.matrix(x)){Hpi.diag(x, nstage = 2)}
  else {hpi(x)}
}


### 2. Flexible bw.nrd function (Silverman)

#' @export

bw.nrd.flex = function(x){
  if(is.matrix(x)){
    bw = apply(x, 2, bw.nrd)
    diag(bw, ncol(x))}
  else {bw.nrd(x)}
}

### 3. Flexible bw.nrd0 function (Scott)
# (using the modified bw.nrd0 function that allows NA values)

#' @export

bw.nrd0.mod = function(x){bw.nrd0(x[complete.cases(x)])}

#' @export

bw.nrd0.flex = function(x){
  if(is.matrix(x)){
    bw = apply(x, 2, bw.nrd0.mod)
    diag(bw, ncol(x))}
  else {bw.nrd0.mod(x)}
}

### 3b. Generalization of Scott's rule (Haerdle et al 2004, p. 73)
# For multivariate bandwidth estimation, this rule should be used instead of bw.nrd0.flex (or bw.nrd.flex), because it takes the number of dimensions into account.

#' @export

bw.nrd0.mult = function(x){
  if(is.matrix(x)){
    n = nrow(x)
    bw = apply(x, 2, bw.nrd0.mod)*(n^(1/5-1/(ncol(x)+4)))
    diag(bw, ncol(x))}
  else {bw.nrd0.mod(x)}
}


#### Functions for kernel density derivative estimation ####
# (only first derivative)

# normal scale diagonal bandwidth matrix by Chacon et al, Stat. Sin. 2011, eq. 3.2 (extension of Silverman/Scott's rule for derivatives)

# univariate:

#' @export

bw.chacon.uni = function(x){
  var(x, na.rm=T)*(4/(length(x)*(5)))^(2/(7))
}

# multivariate:

#' @export
bw.chacon = function(x){
  if(is.matrix(x)){
    d = ncol(x)
    n = nrow(x)
    h = diag(var(x, na.rm=T))*(4/(n*(d + 4)))^(2/(d+6))
    diag(h, d)
  } else {
    bw.chacon.uni(x)
  }
}
