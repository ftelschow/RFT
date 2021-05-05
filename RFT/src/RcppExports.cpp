// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// tupel2index
int tupel2index(int x0, int x1, IntegerVector Dim);
RcppExport SEXP _RFT_tupel2index(SEXP x0SEXP, SEXP x1SEXP, SEXP DimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Dim(DimSEXP);
    rcpp_result_gen = Rcpp::wrap(tupel2index(x0, x1, Dim));
    return rcpp_result_gen;
END_RCPP
}
// triple2index
int triple2index(int x0, int x1, int x2, IntegerVector Dim);
RcppExport SEXP _RFT_triple2index(SEXP x0SEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP DimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< int >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Dim(DimSEXP);
    rcpp_result_gen = Rcpp::wrap(triple2index(x0, x1, x2, Dim));
    return rcpp_result_gen;
END_RCPP
}
// ECcrit1D_C
NumericMatrix ECcrit1D_C(NumericMatrix X);
RcppExport SEXP _RFT_ECcrit1D_C(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(ECcrit1D_C(X));
    return rcpp_result_gen;
END_RCPP
}
// ECcrit2D_C
NumericVector ECcrit2D_C(NumericVector X, IntegerVector dimX);
RcppExport SEXP _RFT_ECcrit2D_C(SEXP XSEXP, SEXP dimXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dimX(dimXSEXP);
    rcpp_result_gen = Rcpp::wrap(ECcrit2D_C(X, dimX));
    return rcpp_result_gen;
END_RCPP
}
// ECcurve1D_C
NumericMatrix ECcurve1D_C(NumericMatrix Y, NumericMatrix dEC, NumericVector u);
RcppExport SEXP _RFT_ECcurve1D_C(SEXP YSEXP, SEXP dECSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dEC(dECSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(ECcurve1D_C(Y, dEC, u));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _RFT_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// meanC
double meanC(NumericVector x);
RcppExport SEXP _RFT_meanC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(meanC(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RFT_tupel2index", (DL_FUNC) &_RFT_tupel2index, 3},
    {"_RFT_triple2index", (DL_FUNC) &_RFT_triple2index, 4},
    {"_RFT_ECcrit1D_C", (DL_FUNC) &_RFT_ECcrit1D_C, 1},
    {"_RFT_ECcrit2D_C", (DL_FUNC) &_RFT_ECcrit2D_C, 2},
    {"_RFT_ECcurve1D_C", (DL_FUNC) &_RFT_ECcurve1D_C, 3},
    {"_RFT_rcpp_hello_world", (DL_FUNC) &_RFT_rcpp_hello_world, 0},
    {"_RFT_meanC", (DL_FUNC) &_RFT_meanC, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_RFT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
