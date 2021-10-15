#include <Rcpp.h>
#include "forward.h"

// [[Rcpp::export]]
Rcpp::List forward_interface
(Rcpp::NumericMatrix log_emission_mat,
 Rcpp::NumericVector initial_prob_vec,
 Rcpp::NumericMatrix transition_mat
 ) {
  int N_data = log_emission_mat.nrow();
  int N_states = log_emission_mat.ncol();
  if(N_data < 1){
    Rcpp::stop("log_emission_mat must have at least one row");
  }
  if(N_states < 1){
    Rcpp::stop("log_emission_mat must have at least one col");
  }
  if(initial_prob_vec.length() != N_states){
    Rcpp::stop("length of initial_prob_vec must be same as number of columns of log_emission_mat");
  }
  if(transition_mat.nrow() != N_states){
    Rcpp::stop("nrow(transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(transition_mat.ncol() != N_states){
    Rcpp::stop("ncol(transition_mat) must be same as ncol(log_emission_mat)");
  }
  Rcpp::NumericMatrix log_alpha_mat(N_data, N_states);
  Rcpp::NumericVector log_lik_vec(1);
  int status = forward
    (N_data,
     N_states,
     &log_emission_mat[0],
     &initial_prob_vec[0],
     &transition_mat[0],
     //inputs above, outputs below.
     &log_alpha_mat[0],
     &log_lik_vec[0]
     );
  if(status == ERROR_INITIAL_PROB_VEC_ENTRIES_MUST_BE_BETWEEN_ZERO_AND_ONE){
    Rcpp::stop("initial_prob_vec entries must be between zero and one");
  }
  return Rcpp::List::create
    (Rcpp::Named("log_alpha", log_alpha_mat),
     Rcpp::Named("log_lik", log_lik_vec));
}

// TODO use arma::cube which is 3d tensor. cube(ptr_aux_mem, n_rows, n_cols, n_slices, copy_aux_mem = true, strict = false)
