#define ERROR_FORWARD_INITIAL_PROB_VEC_ENTRIES_MUST_BE_BETWEEN_ZERO_AND_ONE 1

int forward
(int N_data,
 int N_states,
 double *log_emission_ptr,
 double *transition_ptr,
 double *initial_prob_ptr,
 //inputs above, outputs below.
 double *log_alpha_ptr,
 double *log_lik_ptr
 );
