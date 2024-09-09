// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_wildfire_type
int get_wildfire_type();
RcppExport SEXP _RCAFE_get_wildfire_type() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_wildfire_type());
    return rcpp_result_gen;
END_RCPP
}
// TESTER
IntegerVector TESTER(int n);
RcppExport SEXP _RCAFE_TESTER(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(TESTER(n));
    return rcpp_result_gen;
END_RCPP
}
// doFire
List doFire(const IntegerMatrix& tsf, const List& regime);
RcppExport SEXP _RCAFE_doFire(SEXP tsfSEXP, SEXP regimeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type tsf(tsfSEXP);
    Rcpp::traits::input_parameter< const List& >::type regime(regimeSEXP);
    rcpp_result_gen = Rcpp::wrap(doFire(tsf, regime));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RCAFE_get_wildfire_type", (DL_FUNC) &_RCAFE_get_wildfire_type, 0},
    {"_RCAFE_TESTER", (DL_FUNC) &_RCAFE_TESTER, 1},
    {"_RCAFE_doFire", (DL_FUNC) &_RCAFE_doFire, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_RCAFE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
