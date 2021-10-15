#include <armadillo>
#include <math.h>
#include "forward.h"
#include "eln.h"

int forward
(int N_data,
 int N_states,
 double *log_emission_ptr,
 double *initial_prob_ptr,
 double *transition_ptr,
 //inputs above, outputs below.
 double *log_alpha_ptr,
 double *log_lik_ptr
 ){
  arma::mat log_emission_mat    //copy_aux_mem, strict(no size change)
    (log_emission_ptr, N_data, N_states, false, true);
  arma::vec initial_prob_vec
    (initial_prob_ptr, N_states, false, true);
  arma::mat transition_mat
    (transition_ptr, N_states, N_states, false, true);
  arma::mat log_alpha_mat
    (log_alpha_ptr, N_data, N_states, false, true);
  for(int state_i=0; state_i<N_states; state_i++){
    double prob = initial_prob_vec(state_i);
    if(!(0 <= prob && prob <= 1)){
      return ERROR_INITIAL_PROB_VEC_ENTRIES_MUST_BE_BETWEEN_ZERO_AND_ONE;
    }
    log_alpha_mat(0,state_i) = elnproduct
      (log(prob), log_emission_mat(0,state_i));
  }
  for(int data_t=1; data_t<N_data; data_t++){
    for(int state_j=0; state_j<N_states; state_j++){
      double log_alpha = -INFINITY;
      for(int state_i=0; state_i<N_states; state_i++){
	double alpha_trans = elnproduct
	  (log_alpha_mat(data_t-1,state_i),
	   log(transition_mat(state_i,state_j)));
	log_alpha = elnsum(log_alpha, alpha_trans);
      }
      log_alpha_mat(data_t,state_j) = elnproduct
	(log_alpha, log_emission_mat(data_t, state_j));
    }
  }
  *log_lik_ptr = -INFINITY;
  for(int state_j=0; state_j<N_states; state_j++){
    *log_lik_ptr = elnsum(*log_lik_ptr, log_alpha_mat(N_data-1, state_j));
  }
  return 0;
}
