// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// KDKW_FD_Rcpp
List KDKW_FD_Rcpp(NumericVector s_obs, NumericVector theta_0, NumericVector theta_min, NumericVector theta_max, String simfun, List fixed_parameters, int K, int nk, double a, double ce, double alpha, double gamma, int A, int C);
RcppExport SEXP _ApproxML_KDKW_FD_Rcpp(SEXP s_obsSEXP, SEXP theta_0SEXP, SEXP theta_minSEXP, SEXP theta_maxSEXP, SEXP simfunSEXP, SEXP fixed_parametersSEXP, SEXP KSEXP, SEXP nkSEXP, SEXP aSEXP, SEXP ceSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP ASEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type s_obs(s_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_0(theta_0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_min(theta_minSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_max(theta_maxSEXP);
    Rcpp::traits::input_parameter< String >::type simfun(simfunSEXP);
    Rcpp::traits::input_parameter< List >::type fixed_parameters(fixed_parametersSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type ce(ceSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(KDKW_FD_Rcpp(s_obs, theta_0, theta_min, theta_max, simfun, fixed_parameters, K, nk, a, ce, alpha, gamma, A, C));
    return rcpp_result_gen;
END_RCPP
}
// SIMtestC
arma::mat SIMtestC(int nk, arma::vec theta, arma::mat sigma);
RcppExport SEXP _ApproxML_SIMtestC(SEXP nkSEXP, SEXP thetaSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(SIMtestC(nk, theta, sigma));
    return rcpp_result_gen;
END_RCPP
}
// bw_nrd0
arma::vec bw_nrd0(arma::mat x);
RcppExport SEXP _ApproxML_bw_nrd0(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(bw_nrd0(x));
    return rcpp_result_gen;
END_RCPP
}
// normal_diag
arma::vec normal_diag(arma::mat dat, arma::vec x, arma::vec H);
RcppExport SEXP _ApproxML_normal_diag(SEXP datSEXP, SEXP xSEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_diag(dat, x, H));
    return rcpp_result_gen;
END_RCPP
}
// chooseRcpp
arma::mat chooseRcpp(int n, double theta, int nk);
RcppExport SEXP _ApproxML_chooseRcpp(SEXP nSEXP, SEXP thetaSEXP, SEXP nkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    rcpp_result_gen = Rcpp::wrap(chooseRcpp(n, theta, nk));
    return rcpp_result_gen;
END_RCPP
}
// testmodelfac
arma::mat testmodelfac(String type, arma::vec parameters, List fixed_parameters);
RcppExport SEXP _ApproxML_testmodelfac(SEXP typeSEXP, SEXP parametersSEXP, SEXP fixed_parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< List >::type fixed_parameters(fixed_parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(testmodelfac(type, parameters, fixed_parameters));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ApproxML_KDKW_FD_Rcpp", (DL_FUNC) &_ApproxML_KDKW_FD_Rcpp, 14},
    {"_ApproxML_SIMtestC", (DL_FUNC) &_ApproxML_SIMtestC, 3},
    {"_ApproxML_bw_nrd0", (DL_FUNC) &_ApproxML_bw_nrd0, 1},
    {"_ApproxML_normal_diag", (DL_FUNC) &_ApproxML_normal_diag, 3},
    {"_ApproxML_chooseRcpp", (DL_FUNC) &_ApproxML_chooseRcpp, 3},
    {"_ApproxML_testmodelfac", (DL_FUNC) &_ApproxML_testmodelfac, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ApproxML(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}