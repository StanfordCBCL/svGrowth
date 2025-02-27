#ifndef FUNCTIONS
#define FUNCTIONS
using json = nlohmann::json;
#include "json.hpp"

/*
void update_time_step_explicit(vessel& curr_vessel);
void update_sigma_handshake(void* curr_vessel);
void update_time_step_handshake(vessel& curr_vessel, int iter_arg = 0);
*/

int print_state_mr(size_t iter, gsl_multiroot_fsolver* s);
int ramp_pressure_test(void* curr_vessel, double P_low, double P_high);
int run_pd_test(vessel& curr_vessel, double P_low, double P_high, double lambda_z_test);
int find_equil_geom(void* curr_vessel);
int equil_obj_f(const gsl_vector* x, void* curr_vessel, gsl_vector* f);
int find_tf_geom(void* vessel_in);
int tf_obj_f(const gsl_vector* x, void* vessel_in, gsl_vector* f);
void update_time_step(vessel& curr_vessel);
void update_kinetics(layer& curr_layer);
int find_iv_geom(void* vessel_in);
double iv_obj_f(double a_mid_guess_inner, void *input_vessel);
void update_sigma(void* curr_layer);
vector<double> constitutive(void* curr_const, double lambda_alpha_s, int ts);
double get_app_visc(void* curr_vessel, int sn);
#endif /* GNR_FUNCTIONS */
