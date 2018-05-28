// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// parallel_pairwise_score
NumericMatrix parallel_pairwise_score(List seqs, NumericMatrix mtx, NumericVector weight);
RcppExport SEXP _alndist_parallel_pairwise_score(SEXP seqsSEXP, SEXP mtxSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mtx(mtxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(parallel_pairwise_score(seqs, mtx, weight));
    return rcpp_result_gen;
END_RCPP
}
// seq_to_idx
IntegerVector seq_to_idx(CharacterVector x, NumericMatrix mtx);
RcppExport SEXP _alndist_seq_to_idx(SEXP xSEXP, SEXP mtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mtx(mtxSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_to_idx(x, mtx));
    return rcpp_result_gen;
END_RCPP
}
// two_seq_score
double two_seq_score(CharacterVector x, CharacterVector y, NumericMatrix mtx, NumericVector weight);
RcppExport SEXP _alndist_two_seq_score(SEXP xSEXP, SEXP ySEXP, SEXP mtxSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mtx(mtxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(two_seq_score(x, y, mtx, weight));
    return rcpp_result_gen;
END_RCPP
}
// pairwise_score
NumericMatrix pairwise_score(List seqs, NumericMatrix mtx, NumericVector weight);
RcppExport SEXP _alndist_pairwise_score(SEXP seqsSEXP, SEXP mtxSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mtx(mtxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwise_score(seqs, mtx, weight));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_alndist_parallel_pairwise_score", (DL_FUNC) &_alndist_parallel_pairwise_score, 3},
    {"_alndist_seq_to_idx", (DL_FUNC) &_alndist_seq_to_idx, 2},
    {"_alndist_two_seq_score", (DL_FUNC) &_alndist_two_seq_score, 4},
    {"_alndist_pairwise_score", (DL_FUNC) &_alndist_pairwise_score, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_alndist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
