// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
    {"_cplr_getDesign", (DL_FUNC) &_cplr_getDesign, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_cplr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
