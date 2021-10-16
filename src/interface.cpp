#include <Rcpp.h>
#include "forward.h"
#include "backward.h"

// [[Rcpp::export]]
Rcpp::List forward_interface
(Rcpp::NumericMatrix log_emission_mat,
 Rcpp::NumericMatrix transition_mat,
 Rcpp::NumericVector initial_prob_vec
 ) {
  int N_data = log_emission_mat.nrow();
  int N_states = log_emission_mat.ncol();
  if(N_data < 1){
    Rcpp::stop("log_emission_mat must have at least one row");
  }
  if(N_states < 1){
    Rcpp::stop("log_emission_mat must have at least one col");
  }
  if(transition_mat.nrow() != N_states){
    Rcpp::stop("nrow(transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(transition_mat.ncol() != N_states){
    Rcpp::stop("ncol(transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(initial_prob_vec.length() != N_states){
    Rcpp::stop("length of initial_prob_vec must be same as number of columns of log_emission_mat");
  }
  Rcpp::NumericMatrix log_alpha_mat(N_data, N_states);
  Rcpp::NumericVector log_lik_vec(1);
  int status = forward
    (N_data,
     N_states,
     &log_emission_mat[0],
     &transition_mat[0],
     &initial_prob_vec[0],
     //inputs above, outputs below.
     &log_alpha_mat[0],
     &log_lik_vec[0]
     );
  if(status == ERROR_FORWARD_INITIAL_PROB_VEC_ENTRIES_MUST_BE_BETWEEN_ZERO_AND_ONE){
    Rcpp::stop("initial_prob_vec entries must be between zero and one");
  }
  return Rcpp::List::create
    (Rcpp::Named("log_alpha", log_alpha_mat),
     Rcpp::Named("log_lik", log_lik_vec));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix backward_interface
(Rcpp::NumericMatrix log_emission_mat,
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
  if(transition_mat.nrow() != N_states){
    Rcpp::stop("nrow(transition_mat) must be same as ncol(log_emission_mat)");
  }
  if(transition_mat.ncol() != N_states){
    Rcpp::stop("ncol(transition_mat) must be same as ncol(log_emission_mat)");
  }
  Rcpp::NumericMatrix log_beta_mat(N_data, N_states);
  int status = backward
    (N_data,
     N_states,
     &log_emission_mat[0],
     &transition_mat[0],
     //inputs above, outputs below.
     &log_beta_mat[0]
     );
  if(status == ERROR_BACKWARD_TRANSITION_MAT_ENTRIES_MUST_BE_BETWEEN_ZERO_AND_ONE){
    Rcpp::stop("transition_mat entries must be between zero and one");
  }
  return log_beta_mat;
}

// TODO use arma::cube which is 3d tensor. cube(ptr_aux_mem, n_rows, n_cols, n_slices, copy_aux_mem = true, strict = false)
