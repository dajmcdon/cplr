// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// r_sign
NumericVector r_sign(const int& n);
RcppExport SEXP _cplr_r_sign(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(r_sign(n));
    return rcpp_result_gen;
END_RCPP
}
// compressCpp
List compressCpp(arma::mat const& X, int q, arma::colvec const& y, double s);
RcppExport SEXP _cplr_compressCpp(SEXP XSEXP, SEXP qSEXP, SEXP ySEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< arma::colvec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(compressCpp(X, q, y, s));
    return rcpp_result_gen;
END_RCPP
}
// getDesign
arma::mat getDesign(arma::cube const& arr, arma::umat mask, int rad);
RcppExport SEXP _cplr_getDesign(SEXP arrSEXP, SEXP maskSEXP, SEXP radSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube const& >::type arr(arrSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< int >::type rad(radSEXP);
    rcpp_result_gen = Rcpp::wrap(getDesign(arr, mask, rad));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cplr_r_sign", (DL_FUNC) &_cplr_r_sign, 1},
    {"_cplr_compressCpp", (DL_FUNC) &_cplr_compressCpp, 4},
    {"_cplr_getDesign", (DL_FUNC) &_cplr_getDesign, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_cplr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
