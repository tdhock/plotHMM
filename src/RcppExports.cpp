// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// forward_interface
Rcpp::List forward_interface(Rcpp::NumericMatrix log_emission_mat, Rcpp::NumericMatrix transition_mat, Rcpp::NumericVector initial_prob_vec);
RcppExport SEXP _plotHMM_forward_interface(SEXP log_emission_matSEXP, SEXP transition_matSEXP, SEXP initial_prob_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type log_emission_mat(log_emission_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type transition_mat(transition_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type initial_prob_vec(initial_prob_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_interface(log_emission_mat, transition_mat, initial_prob_vec));
    return rcpp_result_gen;
END_RCPP
}
// backward_interface
Rcpp::NumericMatrix backward_interface(Rcpp::NumericMatrix log_emission_mat, Rcpp::NumericMatrix transition_mat);
RcppExport SEXP _plotHMM_backward_interface(SEXP log_emission_matSEXP, SEXP transition_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type log_emission_mat(log_emission_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type transition_mat(transition_matSEXP);
    rcpp_result_gen = Rcpp::wrap(backward_interface(log_emission_mat, transition_mat));
    return rcpp_result_gen;
END_RCPP
}
// multiply_interface
Rcpp::NumericMatrix multiply_interface(Rcpp::NumericMatrix log_alpha_mat, Rcpp::NumericMatrix log_beta_mat);
RcppExport SEXP _plotHMM_multiply_interface(SEXP log_alpha_matSEXP, SEXP log_beta_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type log_alpha_mat(log_alpha_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type log_beta_mat(log_beta_matSEXP);
    rcpp_result_gen = Rcpp::wrap(multiply_interface(log_alpha_mat, log_beta_mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_plotHMM_forward_interface", (DL_FUNC) &_plotHMM_forward_interface, 3},
    {"_plotHMM_backward_interface", (DL_FUNC) &_plotHMM_backward_interface, 2},
    {"_plotHMM_multiply_interface", (DL_FUNC) &_plotHMM_multiply_interface, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_plotHMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
