// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// critcpp
SEXP critcpp(SEXP a, SEXP b);
RcppExport SEXP _ClustVarLV_critcpp(SEXP aSEXP, SEXP bSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type a(aSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(critcpp(a, b));
    return rcpp_result_gen;
 END_RCPP
 }
// mincpp
SEXP mincpp(SEXP a);
RcppExport SEXP _ClustVarLV_mincpp(SEXP aSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(mincpp(a));
    return rcpp_result_gen;
 END_RCPP
 }
 // powerEigen
 Rcpp::List powerEigen(Eigen::MatrixXd& X);
 RcppExport SEXP _ClustVarLV_powerEigen(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(powerEigen(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ClustVarLV_critcpp", (DL_FUNC) &_ClustVarLV_critcpp, 2},
    {"_ClustVarLV_mincpp", (DL_FUNC) &_ClustVarLV_mincpp, 1},
    {"_ClustVarLV_powerEigen", (DL_FUNC) &_ClustVarLV_powerEigen, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_ClustVarLV(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
