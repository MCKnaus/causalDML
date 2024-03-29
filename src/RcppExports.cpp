// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// norm_drl_rcpp
arma::vec norm_drl_rcpp(arma::sp_mat alpha, arma::mat y_mat, arma::colvec y, arma::mat w_mat, arma::mat e_mat);
RcppExport SEXP _causalDML_norm_drl_rcpp(SEXP alphaSEXP, SEXP y_matSEXP, SEXP ySEXP, SEXP w_matSEXP, SEXP e_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_mat(y_matSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w_mat(w_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type e_mat(e_matSEXP);
    rcpp_result_gen = Rcpp::wrap(norm_drl_rcpp(alpha, y_mat, y, w_mat, e_mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_causalDML_norm_drl_rcpp", (DL_FUNC) &_causalDML_norm_drl_rcpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_causalDML(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
