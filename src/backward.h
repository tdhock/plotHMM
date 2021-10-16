#define ERROR_BACKWARD_TRANSITION_MAT_ENTRIES_MUST_BE_BETWEEN_ZERO_AND_ONE 1

int backward
(int N_data,
 int N_states,
 double *log_emission_ptr,
 double *transition_ptr,
 //inputs above, outputs below.
 double *log_beta_ptr
 );
