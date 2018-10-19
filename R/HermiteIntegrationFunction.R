#' Gauss-Hermite quadrature
#'
#' The one-dimensional integral of a function is numerically approximated with Gauss-Hermite quadrature with 20 nodes. See Abramowitz and Stegun 25.4.46 pp 890 and Table 25.10 pp 924.
#'
#' Integrals of the form int f(x) phi(x, sigma^2) dx are computed, where phi is the normal distribution density with mean 0 and variance sigma^2. FUN is the function f() and sig2 is sigma^2.
#'
#' @references Abramowitz, M. & Stegun, I. A. (Eds.) Handbook of mathematical functions: with formulas, graphs, and mathematical tables. Dover Publications New York, 1972
#'
#' @param FUN function to be integrated
#' @param sig2 scalar
#' @param ... Additional arguments to the function given in FUN
#'
#' @export

HermiteIntegration <- function(FUN, sig2 = 1, ...){

  # Setting the weights
  W <- c(2.229393645534e-13, 4.399340992273e-10, 1.086069370769e-07, 7.802556478532e-06, 2.283386360163e-04,
         3.243773342238e-03, 2.481052088746e-02, 1.090172060200e-01, 2.866755053628e-01, 4.622436696006e-01,
         4.622436696006e-01, 2.866755053628e-01, 1.090172060200e-01, 2.481052088746e-02, 3.243773342238e-03,
         2.283386360163e-04, 7.802556478532e-06, 1.086069370769e-07, 4.399340992273e-10, 2.229393645534e-13)
  # Seting the nodes
  X <- c(-5.387480900112, -4.6036824495507, -3.9447640401156, -3.3478545673832, -2.7888060584281,
         -2.2549740020893, -1.7385377121166, -1.2340762153953, -0.7374737285454, -0.2453407083009,
         0.2453407083009,  0.7374737285454,  1.2340762153953,  1.7385377121166,  2.2549740020893,
         2.7888060584281,  3.3478545673832,  3.9447640401156,  4.6036824495507,  5.387480900112)

  FUNv <- Vectorize(FUN = FUN)

  # Integration adapted to use the normal density with mean = 0 and var = sigma2
  Integral <- sum(FUNv(X * sqrt(sig2*2), ...) * W / sqrt(pi))

  return(Integral)
}
