#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>


#include "vessel.h"
#include "layer.h"
#include "json.hpp"
#include "constituent.h"
#include "functions.h"
#include "matfun.h"
#include "exprtk.hpp"
// #define COMBINED_DELTA_SIGMA 
#define ZERO_LAGRANGE_BC
using std::string;
using std::vector;
using std::cout;



int print_state_mr(size_t iter, gsl_multiroot_fsolver* s) {

    printf("iter = %3u", iter);
    printf("x = ");
    for (int i = 0; i < s->x->size; i++) {
        printf("%.3f ", gsl_vector_get(s->x, i));
    }
    printf("f(x) = ");
    for (int i = 0; i < s->x->size; i++) {
        printf("%.3f ", gsl_vector_get(s->f, i));
    }
    printf("\n");

    return 0;
}

int ramp_pressure_test(void* curr_vessel, double P_low, double P_high) {
    int sn =((struct vessel*) curr_vessel)->sn;
    int equil_check = 0;
    int num_P = 100;
    double P_incr = (P_high - P_low) / num_P;


    ((struct vessel*) curr_vessel)->P[sn] = P_low;
    //std::cout << "Start ramp..." << std::endl;
    for (int i = 0; i < num_P; i++) {
        equil_check = find_iv_geom(curr_vessel);
        ((struct vessel*) curr_vessel)->P[sn] += P_incr;
    }
    return equil_check;
}

int run_pd_test(vessel& curr_vessel, double P_low, double P_high, double lambda_z_test) {

    int sn = curr_vessel.sn;
    double lambda_z_store = curr_vessel.lambda_z_tau[sn];
    double P_store = curr_vessel.P[sn];

    int equil_check = 0;
    int num_P = 100;
    double P_incr = (P_high - P_low) / num_P;
    curr_vessel.P[sn] = P_low;
    curr_vessel.lambda_z_tau[sn] = lambda_z_test;

    for (int i = 0; i < num_P; i++) {

        equil_check = find_iv_geom(&curr_vessel);
        curr_vessel.Exp_out << curr_vessel.P[sn];
        for (int l = 0; l < curr_vessel.layers.size(); l++) {
            curr_vessel.Exp_out << "\t" << curr_vessel.layers[l].a_mid[sn];
        }
        curr_vessel.Exp_out << "\n";
        curr_vessel.P[sn] += P_incr;

    }

    curr_vessel.lambda_z_tau[sn] = lambda_z_store;
    curr_vessel.P[sn] = P_store;

    return 0;
}

int find_equil_geom(void* vessel_in) {
    //Finds the mechanobiologically equilibrated geometry for a given set of loads inclduing
    //pressure, flow, and axial stretch with a set of G&R parameter values from the original
    //homeostatic state
    struct vessel *curr_vessel = (struct vessel *) vessel_in;
    if (curr_vessel->layers.size() > 1) throw std::invalid_argument("Input to find_equil_geom must be single layer vessel");
    int sn = curr_vessel->sn;

    //Loading changes
    double gamma = curr_vessel->P[sn] / curr_vessel->P_h; //Fold change in pressure from homeostatic
    double epsilon = curr_vessel->Q[sn] / curr_vessel->Q_h; //Fold chance in flow from homeostatic
    double lambda = curr_vessel->lambda_z_tau[sn] / curr_vessel->lambda_z_h; //Fold change in axial stretch from homeostatic

    //Homeostatic geometry
    double a_h = curr_vessel->layers[0].a_h;
    double h_h = curr_vessel->layers[0].h_h;

    //Initial guesses based on the loading changes
    double a_e_guess = pow(epsilon, 1.0 / 3.0) * a_h;
    double h_e_guess = gamma * pow(epsilon, 1.0 / 3.0) * h_h;
    double C_guess = 1.0;
    double f_z_e_guess = curr_vessel->f_h* (h_e_guess * (2 * a_e_guess + h_e_guess)) / (h_h * (2 * a_h + h_h));

    //For psuedo-time dependent evolutions
    //if (0.98 * ((struct vessel*) curr_vessel)->a_h > ((struct vessel*) curr_vessel)->a_e || ((struct vessel*) curr_vessel)->a_e > 1.02 * ((struct vessel*) curr_vessel)->a_h){
        //a_e_guess = ((struct vessel*) curr_vessel)-> a_e;
        //h_e_guess = ((struct vessel*) curr_vessel)-> h_e;
        //rho_c_e_guess = ((struct vessel*) curr_vessel)-> rho_c_e;
        //f_z_e_guess = ((struct vessel*) curr_vessel)-> f_z_e;
    //}
    
    //printf("%s %f %s %f %s %f\n", "Time:", ((struct vessel*) curr_vessel)->s, "a_e_guess: ", ((struct vessel*) curr_vessel)->a_e, "h_e_guess:", ((struct vessel*) curr_vessel)->h_e);
    //fflush(stdout);

    const gsl_multiroot_fsolver_type* T;
    gsl_multiroot_fsolver* s;

    int status;
    size_t iter = 0;

    const size_t n = 4;

    gsl_multiroot_function f = { &equil_obj_f, n, curr_vessel };
    double x_init[4] = {a_e_guess, h_e_guess, C_guess, f_z_e_guess};
    gsl_vector* x = gsl_vector_alloc(n);

    for (int i = 0; i < n; i++) {
        gsl_vector_set(x, i, x_init[i]);
    }

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, n);

    gsl_multiroot_fsolver_set(s, &f, x);

    //print_state_mr(iter, s);

    double epsabs = 1e-7;
    double epsrel = 1e-3;
    int check_tol = 0;
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        //print_state_mr(iter, s);

        if (status)
            break;

        //status = gsl_multiroot_test_delta(s->dx, s->x, epsabs, epsrel);
        status = gsl_multiroot_test_residual(s->f, epsabs);
        //for (int i = 0; i < n; i++)
            //check_tol += (gsl_vector_get(s->f, i) > epsabs);

        //if (check_tol >= 0)

    } while (status == GSL_CONTINUE && iter < 1000);
    
    //printf("status = %s \n", gsl_strerror(status));

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;

}

int equil_obj_f(const gsl_vector* x, void* vessel_in, gsl_vector* f) {
    struct vessel *curr_vessel = (struct vessel *) vessel_in;
    int sn = curr_vessel->sn;
    //Mechanobiologically equilibrated objective function
    //Unknown input variables
    const double a_e_guess = gsl_vector_get(x, 0);
    const double h_e_guess = gsl_vector_get(x, 1);
    const double C_guess = gsl_vector_get(x, 2);
    const double f_z_e_guess = gsl_vector_get(x, 3);

    int n_alpha = curr_vessel->layers[0].constituents.size();
    int first_gnr_idx = 0;
    for (int alpha = 0; alpha < n_alpha; alpha ++) {
        if (curr_vessel->layers[0].constituents[alpha].K_tauw_p_alpha_h > 0) {
            first_gnr_idx = alpha;
            break;
        }
    }

    //Equations for J1
    //Stress values from equilibrium equations
    double sigma_e_th_lmb = curr_vessel->P[sn] * a_e_guess / h_e_guess;
    double sigma_e_z_lmb = f_z_e_guess / (M_PI * h_e_guess * (2 * a_e_guess + h_e_guess));

    //WSS from Pousielle flow, constant viscosity
    int nts = curr_vessel->nts;
    double a_store = curr_vessel->layers[0].a[sn];
    curr_vessel->layers[0].a[sn] = a_e_guess;
    double mu = get_app_visc(curr_vessel, sn);
    curr_vessel->layers[0].a[sn] = a_store;
    double bar_tauw_e = 4*mu*curr_vessel->Q[sn]/(3.14159265*pow(a_e_guess*100, 3)); 

    //Ratio of stress:WSS mediated matrix production
    double eta_K = curr_vessel->layers[0].constituents[first_gnr_idx].K_sigma_p_alpha_h /
    curr_vessel->layers[0].constituents[first_gnr_idx].K_tauw_p_alpha_h;

    //Stress and WSS deviations from HS state
    double delta_sigma = (sigma_e_th_lmb + sigma_e_z_lmb) /
    (curr_vessel->layers[0].sigma_h[1] + curr_vessel->layers[0].sigma_h[2]) - 1;
    double delta_tauw = bar_tauw_e / curr_vessel->bar_tauw_h - 1;

    //Equations for J2
    //Homeostatic geometry
    double a_h = curr_vessel->layers[0].a_h;
    double h_h = curr_vessel->layers[0].h_h;

    //Equilibrated stretches
    double lambda_r_e = h_e_guess / h_h;
    double lambda_th_e = (a_e_guess + h_e_guess / 2) / (a_h + h_h / 2);
    double lambda_z_e = curr_vessel->lambda_z_tau[sn];
    double F_e[3] = { lambda_r_e, lambda_th_e, lambda_z_e };

    //Equilibrated volume change
    double J_e = lambda_r_e * lambda_th_e * lambda_z_e; 






    //Equilibrated mass densities
    vector<double> rho_alpha = {0.0};
    rho_alpha.resize(n_alpha);
    
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        if (curr_vessel->layers[0].constituents[alpha].k_alpha_h == 0) {
            //Equilibrated elastin density, account for possible elastin degradation
            rho_alpha[alpha] = curr_vessel->layers[0].constituents[alpha].rhoR_alpha[sn] / J_e; 
        } else {
            double eta = curr_vessel->layers[0].constituents[alpha].k_alpha_h * curr_vessel->layers[0].constituents[alpha].K_sigma_p_alpha_h;
            rho_alpha[alpha] = curr_vessel->layers[0].constituents[alpha].rhoR_alpha_h / J_e * pow(C_guess, eta);
            curr_vessel->layers[0].constituents[alpha].rho_alpha_e = rho_alpha[alpha] * J_e;

        }
    }

    //Array of equilibrated densities
    //Add individual collagen densities
    double rho_h = curr_vessel->layers[0].rhoR_h;

    //cauchy stress hat for each const. in each direction
    vector<double> hat_sigma_alpha_dir(3 * n_alpha, 0);
    vector<double> sigma_e_dir(3, 0);

    //equilibrated stresses
    vector<double> sigma_e(3, 0);



    //equilibrated constituent strech
    double lambda_alpha_ntau_s = 0.0;

    //equilbirated 2nd PK stress hat
    double hat_S_alpha = 0.0;

    double C = 0;
    double lambda_act = 1.0;
    double parab_act = 0;
    double hat_sigma_act_e = 0;
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        for (int dir = 0; dir < 3; dir++) {
            //Check if constituent is anisotropic
            if (curr_vessel->layers[0].constituents[alpha].eta_alpha_h >= 0) {
                //Constituent stretch is equal to deposition stretch
                lambda_alpha_ntau_s = curr_vessel->layers[0].constituents[alpha].g_alpha_h[sn];

                //2nd PK stress hat at equilibrium
                hat_S_alpha = curr_vessel->layers[0].constituents[alpha].c_alpha_h[0 * nts + sn] * (pow(lambda_alpha_ntau_s, 2) - 1) *
                exp(curr_vessel->layers[0].constituents[alpha].c_alpha_h[1 * nts + sn] * pow(pow(lambda_alpha_ntau_s, 2) - 1, 2));

                //Cauchy stress hat at equilibrium
                hat_sigma_alpha_dir[alpha * dir + dir] = curr_vessel->layers[0].constituents[alpha].G_alpha_h[dir * nts + sn] * 
                hat_S_alpha * curr_vessel->layers[0].constituents[alpha].G_alpha_h[dir * nts + sn];

            }
            else {
                //2nd PK stress hat at equilibrium
                hat_S_alpha = curr_vessel->layers[0].constituents[alpha].c_alpha_h[0 * nts + sn];

                //Cauchy stress hat at equilibrium
                hat_sigma_alpha_dir[alpha * dir + dir] = curr_vessel->layers[0].constituents[alpha].G_alpha_h[dir * nts + sn] * hat_S_alpha *
                curr_vessel->layers[0].constituents[alpha].G_alpha_h[dir * nts + sn];

                //Check if the consituent is present from the initial time point
                //Account for volume change and mixture deformation
                if (curr_vessel->layers[0].constituents[alpha].k_alpha_h == 0) {
                    hat_sigma_alpha_dir[alpha * dir + dir] = F_e[dir] * hat_sigma_alpha_dir[alpha * dir + dir] * F_e[dir];
                }

            }

            sigma_e_dir[dir] += rho_alpha[alpha] / rho_h * hat_sigma_alpha_dir[alpha * dir + dir];

            if (curr_vessel->layers[0].constituents[alpha].alpha_active == 1 && dir == 1) {
                //equilibrated active stress
                C = curr_vessel->layers[0].constituents[alpha].CB - curr_vessel->layers[0].constituents[alpha].CS * delta_tauw;
                parab_act = 1 - pow((curr_vessel->layers[0].constituents[alpha].lambda_m - lambda_act) /
                    (curr_vessel->layers[0].constituents[alpha].lambda_m - curr_vessel->layers[0].constituents[alpha].lambda_0), 2);
                hat_sigma_act_e = curr_vessel->layers[0].constituents[alpha].T_act * (1 - exp(-pow(C, 2))) * lambda_act * parab_act;

                sigma_e_dir[dir] += rho_alpha[alpha] / rho_h * hat_sigma_act_e;
            }
        }
    }

    //Four objective equations
    double J1 = eta_K*delta_sigma - delta_tauw; //mechano-mediated matrix production equation
    double J2 = - rho_h; //mixture mass balance equation
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        J2 += rho_alpha[alpha];
    }
    double J3 = sigma_e_dir[1] - sigma_e_dir[0] - sigma_e_th_lmb; //circumferential const - lmb
    double J4 = sigma_e_dir[2] - sigma_e_dir[0] - sigma_e_z_lmb; //axial const - lmb

    gsl_vector_set(f, 0, J1);
    gsl_vector_set(f, 1, J2);
    gsl_vector_set(f, 2, J3);
    gsl_vector_set(f, 3, J4);

    //Store equilibrated results  
    //if (a_e_guess == a_e_guess) {
    curr_vessel->layers[0].a_e = a_e_guess;
    curr_vessel->layers[0].h_e = h_e_guess;
    curr_vessel->f_z_e = f_z_e_guess;
    curr_vessel->mb_equil_e = 1 + curr_vessel->layers[0].constituents[first_gnr_idx].K_sigma_p_alpha_h * delta_sigma -
    curr_vessel->layers[0].constituents[first_gnr_idx].K_tauw_p_alpha_h * delta_tauw;
    //}     

    return GSL_SUCCESS;
}

int find_tf_geom(void* vessel_in) {
    struct vessel *curr_vessel = (struct vessel *) vessel_in;
    int sn = curr_vessel->sn;


    //For geometries defined in the loaded configuration this code solves for the unloaded geometric 
    //variables at the current time point and stores them in the current v
    //vessel structure. Unloaded is also
    //referred to as traction free (but not stress free).

    //Local vars to store current pressure, force, and stretch
    double P_temp = curr_vessel->P[sn];
    double f_temp = curr_vessel->f;
    double lambda_z_temp = curr_vessel->lambda_z_tau[sn];

    //Get initial guesses from the loaded geometry
    double lambda_th_ul = 0.95 * curr_vessel->layers[0].a_mid[sn] / curr_vessel->layers[0].a_mid[0];
    double lambda_z_ul = 0.95 * lambda_z_temp;

    //Update vesseel loads to zero for traction free
    curr_vessel->P[sn] = 0.0;
    curr_vessel->f = 0.0;
    //((struct vessel*) curr_vessel)->T_act = 0.0;

    const gsl_multiroot_fsolver_type* T;
    gsl_multiroot_fsolver* s;

    int status;
    size_t iter = 0;

    const size_t n = 2;

    gsl_multiroot_function f = { &tf_obj_f, n, curr_vessel };
    double x_init[2] = { lambda_th_ul, lambda_z_ul };
    gsl_vector* x = gsl_vector_alloc(n);

    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, 2);

    gsl_multiroot_fsolver_set(s, &f, x);

    //print_state_mr(iter, s);

    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        //print_state_mr(iter, s);

        if (status)   // check if solver is stuck
            break;

        status = gsl_multiroot_test_residual(s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 100);

    //printf("status = %s\n", gsl_strerror(status));

    //Store the traction free results
    for (int layer = 0; layer < curr_vessel->layers.size(); layer++) {
        curr_vessel->layers[layer].A_mid[sn] = gsl_vector_get(s->x, 0) * curr_vessel->layers[layer].a_mid[0];
        curr_vessel->layers[layer].lambda_z_pre[sn] = 1 / gsl_vector_get(s->x, 1);
        curr_vessel->layers[layer].H[sn] = 1.0 / (gsl_vector_get(s->x, 0) * gsl_vector_get(s->x, 1)) * curr_vessel->layers[layer].h[0];
        curr_vessel->layers[layer].A[sn] = curr_vessel->layers[layer].A_mid[sn] - curr_vessel->layers[layer].H[sn] / 2;

        //Return current vars to the loaded conditions
    }


    curr_vessel->lambda_z_tau[sn] = lambda_z_temp;
    curr_vessel->P[sn] = P_temp;
    curr_vessel->f = f_temp;

    status = find_iv_geom(curr_vessel);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;

}

int tf_obj_f(const gsl_vector* x, void* vessel_in, gsl_vector* f) {

    struct vessel *curr_vessel = (struct vessel *) vessel_in;
    //Seperate out inputs, both are absolute strains (not relative to current values)
    const double lambda_th_ul_guess = gsl_vector_get(x, 0);
    const double lambda_z_ul_guess = gsl_vector_get(x, 1);

    //Finds the difference in the theoretical stress from Laplace for deformed mixture
    //from the stress calculated from the mixture equations
    int n_layers = curr_vessel->layers.size();
    int sn = curr_vessel->sn;
    double a = 0.0, h = 0.0, lambda_t = 0.0, J_s = 1.0, a_mid = 0.0;
    double lambda_z = lambda_z_ul_guess;
    double pa_calc = 0.0, fz_total = 0.0;

    curr_vessel->lambda_z_tau[sn] = lambda_z;

    for (int layer = 0; layer < n_layers; layer ++) {

        J_s = curr_vessel->layers[layer].rhoR[sn] / curr_vessel->layers[layer].rho[sn];

        if (layer == 0) {
            lambda_t = lambda_th_ul_guess;
            a_mid = lambda_t * curr_vessel->layers[layer].a_mid[0];
        } else {
            double h_h = curr_vessel->layers[layer].h[0];
            double a_h = curr_vessel->layers[layer].a[0];
            // want to solve the equation A lambdat^2 + B lambdat + C = 0
            double A = a_h + h_h / 2;
            double B = -(a + h);
            double C = -h_h / 2 * J_s / lambda_z;
            lambda_t = -B / 2 + sqrt(B * B / 4 - C);
            a_mid = lambda_t * (a_h + h_h / 2);
        }

        h = J_s / (lambda_t * lambda_z) * curr_vessel->layers[layer].h[0];
        a = a_mid - h / 2;

        //Update current total stretches
        curr_vessel->layers[layer].lambda_th_curr = lambda_t;

        update_sigma(&curr_vessel->layers[layer]);

        //Should be 0 for the traction-free configuration
        pa_calc += h * curr_vessel->layers[layer].sigma[1];
        fz_total += M_PI * h * (2 * a + h) * curr_vessel->layers[layer].sigma[2];
    }
    gsl_vector_set(f, 0, pa_calc);
    gsl_vector_set(f, 1, fz_total);

    return GSL_SUCCESS;

}

void update_time_step(vessel& curr_vessel) {


    int n_layers = curr_vessel.layers.size();
    int nts = curr_vessel.nts;
    int sn = curr_vessel.sn;
    double s = sn * curr_vessel.dt;
    vector<double> rhoR_s0 = {}; //Total mass at current step
    double rhoR_s1 = 0;
    double tol = 1E-14; //Convergence tolerance
    int iter = 0, equil_check = 0, unloaded_check = 0;
    double delta_sigma, delta_tauw, mb_equil;
    //Get initial guess for radius
    for (int layer = 0; layer < n_layers; layer++) {
        curr_vessel.layers[layer].a_mid[sn] = curr_vessel.layers[layer].a_mid[sn - 1];
        for (int alpha = 0; alpha < curr_vessel.layers[layer].constituents.size(); alpha++) {
            curr_vessel.layers[layer].constituents[alpha].a_act[sn] = curr_vessel.layers[layer].constituents[alpha].a_act[sn - 1];
        }

        update_kinetics(curr_vessel.layers[layer]);
    }


    //Find new equilibrium geometry.
    //std::cout << curr_vessel.a_mid[sn]<< std::endl;
    // step up gradually if needed
    double delta_P_th = 0.05;
    if (abs((curr_vessel.P[sn] - curr_vessel.P[sn - 1]) / curr_vessel.P[sn - 1]) > delta_P_th) {
        ramp_pressure_test(&curr_vessel, curr_vessel.P[sn - 1], curr_vessel.P[sn]);
    }

    equil_check = find_iv_geom(&curr_vessel);

    //Find real mass density production at the current time step iteratively
    double mass_check = 0.0;
    do {
        iter++;
        //Value from previous prediction
        rhoR_s0 = {};
        mass_check = 0;
        //Update prediction
        for (int layer = 0; layer < curr_vessel.layers.size(); layer++) {
            rhoR_s0.push_back(curr_vessel.layers[layer].rhoR[sn]);
            update_kinetics(curr_vessel.layers[layer]);
        }   

        equil_check = find_iv_geom(&curr_vessel);

        for (int layer = 0; layer < curr_vessel.layers.size(); layer++) {
            rhoR_s1 = curr_vessel.layers[layer].rhoR[sn];
            mass_check = fmax(abs((rhoR_s1 - rhoR_s0[layer]) / rhoR_s0[layer]), mass_check);
        }   
    } while (mass_check > tol&& iter < 100);

    ////Find current mechanobiological perturbation

    // Update is_contacting variable
    for (int layer = 1; layer < curr_vessel.layers.size(); layer++) {

        double b_inner = curr_vessel.layers[layer - 1].a[sn] + curr_vessel.layers[layer - 1].h[sn];
        double a_h_outer = curr_vessel.layers[layer].a_h;
        curr_vessel.layers[layer].is_contacting = curr_vessel.layers[layer].is_contacting || (b_inner > a_h_outer);
    }


    delta_sigma = (curr_vessel.layers[0].sigma_inv / curr_vessel.layers[0].sigma_inv_h) - 1;
    delta_tauw = (curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1;

    // TODO why is this hard coded ?
    mb_equil = 1 + curr_vessel.layers[0].constituents[2].K_sigma_p_alpha_h * delta_sigma -
    curr_vessel.layers[0].constituents[2].K_tauw_p_alpha_h * delta_tauw;

}


void update_kinetics(layer& curr_layer) {
//    std cout << "Begin update kinetics ----------------------------" << std::endl;

    //This function updates the kinetics for G&R.
    int n_alpha = curr_layer.constituents.size();
    int nts = curr_layer.parent_vessel->nts;
    double dt = curr_layer.parent_vessel->dt;
    double s = curr_layer.parent_vessel->s;
    int sn = curr_layer.parent_vessel->sn;
    int taun_min = 0;

    double tau_max = nts + 1; 
    // for (int alpha = 0; alpha < n_alpha; alpha++) {
    //     double k_alpha_h = curr_layer.constituents[alpha].k_alpha_h;
    //     if (k_alpha_h > 0 && 10 / k_alpha_h > tau_max) {
    //         tau_max = 10 / k_alpha_h;
    //     }
    // }
    //Differences in current mechanical state from the reference state
    //Stress invariant
    double delta_sigma = (curr_layer.sigma_inv / curr_layer.sigma_inv_h) - 1;
    //Wall shear stress
    double delta_tauw = (curr_layer.parent_vessel->bar_tauw / curr_layer.parent_vessel->bar_tauw_h) - 1;

    //Initialize pars for looping later
    double K_sigma_p = 0, K_tauw_p = 0, K_sigma_d = 0, K_tauw_d = 0, K_infl_p = 0;
    double upsilon_p = 0, upsilon_d = 0;

    double k_alpha_s = 0;
    double mR_alpha_s = 0;
    double rhoR_alpha_calc = 0;

    double mq_0 = 0, mq_1 = 0, mq_2;
    double q_0 = 0, q_1 = 0, q_2;
    double k_0 = 0, k_1 = 0, k_2;
    double rhoR_s = 0, rhoR_alpha_s = 0;
    double J_s = 0;

    int n = 0; //number of points in integration interval

    bool deg_check = 0;

    //Check if we've exceeded the initial history
    if (s > tau_max) {
        taun_min = sn - int(tau_max / dt);
    }
    else {
        taun_min = 0;
    }

    n = (sn - taun_min) + 1; //find number of integration pts
    bool even_n = n % 2 == 0; //check if even # of int points

    //Loop through each constituent to update its mass density
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Find if constituent degrades
        deg_check = curr_layer.constituents[alpha].k_alpha_h > 0;

        if ((sn > 0 && deg_check) && curr_layer.parent_vessel->pol_only_flag == 0) {

            if (curr_layer.constituents[alpha].alpha_infl == 0) {

                //Get the gains for the current constituent
                K_sigma_p = curr_layer.constituents[alpha].K_sigma_p_alpha_h;
                K_tauw_p = curr_layer.constituents[alpha].K_tauw_p_alpha_h;

                K_sigma_d = curr_layer.constituents[alpha].K_sigma_d_alpha_h;
                K_tauw_d = curr_layer.constituents[alpha].K_tauw_d_alpha_h;

                K_infl_p = curr_layer.constituents[alpha].K_infl_p_alpha[sn];


                //Update the stimulus functions for each constituent
                upsilon_p = 1 + K_sigma_p * delta_sigma - K_tauw_p * delta_tauw + K_infl_p;
                if (upsilon_p < 0.1) {
                    upsilon_p = 0.1;
                }

                upsilon_d = 1 + K_sigma_d * pow(delta_sigma, 2) + K_tauw_d * pow(delta_tauw, 2);

                if (curr_layer.constituents[alpha].rhoR_alpha[sn] <= curr_layer.constituents[alpha].rhoR_alpha_h) {
                    rhoR_alpha_calc = curr_layer.constituents[alpha].rhoR_alpha_h;
                }
                else {
                    rhoR_alpha_calc = curr_layer.constituents[alpha].rhoR_alpha[sn];
                }
            }
            else {

                upsilon_p = curr_layer.ups_infl_p[sn];
                upsilon_d = 1 + curr_layer.ups_infl_d[sn];

                rhoR_alpha_calc = curr_layer.constituents[alpha].rhoR_alpha_h;

            }

            //Make sure productions don't become negative
            upsilon_p = (upsilon_p > 0)* upsilon_p;
            upsilon_d = (upsilon_d > 0)* upsilon_d;

            curr_layer.constituents[alpha].ups_saved = upsilon_p;

            //Reset these to zero for each constituent
            k_alpha_s = 0;
            mR_alpha_s = 0;
            rhoR_alpha_s = 0;

            //Update kinetic values for the current time
            k_alpha_s = curr_layer.constituents[alpha].k_alpha_h * upsilon_d;
            
            mR_alpha_s = curr_layer.constituents[alpha].k_alpha_h * rhoR_alpha_calc * upsilon_p; //curr_layer.rhoR_alpha_h[alpha];

            curr_layer.constituents[alpha].k_alpha[sn] = k_alpha_s;
            curr_layer.constituents[alpha].mR_alpha[sn] = mR_alpha_s;

            k_2 = curr_layer.constituents[alpha].k_alpha[sn];
            q_2 = 1.0;
            mq_2 = curr_layer.constituents[alpha].mR_alpha[sn] * q_2;

            //loop through and update constituent densities from previous time points
            //starting from the current time point and counting down is more efficient
            for (int taun = sn - 1; taun >= taun_min; taun = taun - 1) { // For simpson's rule: taun = taun - 2 and taun >= taun_min + 1

                //Trapezoidal rule   
                k_1 = curr_layer.constituents[alpha].k_alpha[taun];
                q_1 = exp(-(k_2 + k_1) * dt / 2) * q_2;
                mq_1 = curr_layer.constituents[alpha].mR_alpha[taun] * q_1;

                // If using Simpson's rule:
                // k_0 = curr_layer.constituents[alpha].k_alpha[taun - 1];
                // q_0 = exp(-(k_2 + 4.0 * k_1 + k_0) * dt / 3) * q_2;
                // mq_0 = curr_layer.constituents[alpha].mR_alpha[taun - 1] * q_0;

                rhoR_alpha_s += (mq_2 + mq_1) * dt / 2;
                k_2 = k_1;
                q_2 = q_1;
                mq_2 = mq_1;

                // If using Simpson's rule:
                // rhoR_alpha_s += (mq_2 + 4.0 * mq_1 + mq_0) * dt / 3;
                // k_2 = k_0;
                // q_2 = q_0;
                // mq_2 = mq_0;

            }

            //At last time step, doing trapezoidal integration if even integration pts (only if using Simpson's rule)
            //  if (even_n) {
            //     k_0 = curr_layer.constituents[alpha].k_alpha[taun_min];
            //     q_0 = exp(-(k_2 + k_0) * dt / 2) * q_2;
            //     mq_0 = curr_layer.constituents[alpha].mR_alpha[taun_min] * q_0;

            //     rhoR_alpha_s += (mq_2 + mq_0) * dt / 2;

            // } 

            //Account for the cohort of material present initially
            if (taun_min == 0) {
                rhoR_alpha_s += curr_layer.constituents[alpha].rhoR_alpha[0] * q_1; // q_0 if using Simpson's rule
            }

            //Update referential volume fraction
            curr_layer.constituents[alpha].epsilonR_alpha[sn] = rhoR_alpha_s / curr_layer.constituents[alpha].rho_hat_alpha_h;
            J_s += curr_layer.constituents[alpha].epsilonR_alpha[sn];
        }
        else {
            //Precalculate maintenance or loss of cosntituents not produced
            rhoR_alpha_s = curr_layer.constituents[alpha].rhoR_alpha[sn];
            J_s += curr_layer.constituents[alpha].epsilonR_alpha[sn];
        }

        curr_layer.constituents[alpha].rhoR_alpha[sn] = rhoR_alpha_s;
        rhoR_s += rhoR_alpha_s;

    }
    //Update spatial volume fractions
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        curr_layer.constituents[alpha].epsilon_alpha[sn] = curr_layer.constituents[alpha].epsilonR_alpha[sn] / J_s;
    }
    curr_layer.rhoR[sn] = rhoR_s;
    curr_layer.rho[sn] = rhoR_s / J_s;
}

int find_iv_geom(void* vessel_in) {
    //Finds the loaded conifiguration for a given pressure and axial stretch


    struct vessel *curr_vessel = (struct vessel *) vessel_in;
    int status;
    int iter = 0;
    int max_iter = 10000;

    int sn = curr_vessel->sn;

    const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
    gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
    gsl_function f = { &iv_obj_f, vessel_in};

    double q = iv_obj_f(curr_vessel->layers[0].a_mid[sn], curr_vessel);

    //Set search range for new mid radius
    double a_mid_act, a_mid_high, a_mid_low;
    a_mid_low = 0.90 * curr_vessel->layers[0].a_mid[sn];
    a_mid_high = 1.1 * curr_vessel->layers[0].a_mid[sn];
    gsl_root_fsolver_set(s, &f, a_mid_low, a_mid_high);

    //printf("Using %s method \n", gsl_root_fsolver_name(s));
    //printf("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err (est)");

    a_mid_act = 0.0;
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        a_mid_act = gsl_root_fsolver_root(s);
        a_mid_low = gsl_root_fsolver_x_lower(s);
        a_mid_high = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(a_mid_low, a_mid_high, 0, pow(10, -12));

        //if (status == GSL_SUCCESS) {
        //    printf("Loaded Config Converged:\n");
        //    printf("%5d [%.7f, %.7f] %.7f %.7fs\n", iter, a_mid_low, a_mid_high, a_mid_act, a_mid_high - a_mid_low);
        //}

        //printf("%5d [%.7f, %.7f] %.7f %.7fs\n", iter, a_mid_low, a_mid_high, a_mid_act, a_mid_high - a_mid_low);
    } while (status == GSL_CONTINUE && iter < max_iter);




    if (iter > max_iter-2){
        printf("status = NOT CONVERGED\n");
    }
    //printf("status = %s\n", gsl_strerror(status));

    //int set_iv = iv_obj_f(a_mid_act, curr_vessel);
    //((struct vessel*) curr_vessel)->a_mid[sn] = a_mid_act;

    return status;
}
// Pass a guess for the middle radius of the inner layer, a list of vessels representing the different layers, and the number of layers 
double iv_obj_f(double a_mid_guess_inner, void *input_vessel) {

    struct vessel *curr_vessel = (struct vessel *) input_vessel;

    //Finds the difference in the theoretical stress from Laplace for deformed mixture
    //from the stress calculated from the mixture equations
    int n_layers = curr_vessel->layers.size();
    int sn = curr_vessel->sn;
    double a = 0.0, h = 0.0, lambda_t = 0.0, lambda_z = 0.0, J_s = 1.0, a_mid = 0.0;
    double mu = 0.0;
    double pa_calc = 0.0;
    double fz_total = 0.0, h_total = 0.0; 
    bool equil_check = !(sn > 0 || curr_vessel->num_exp_flag == 1); // True if trying to determine the homeostatic equilibrium state

    lambda_z = curr_vessel->lambda_z_tau[sn] / curr_vessel->lambda_z_tau[0];  // Read z stretch from layer

    for (int layer = 0; layer < n_layers; layer ++) {
        if (!equil_check) J_s = curr_vessel->layers[layer].rhoR[sn] / curr_vessel->layers[layer].rho[sn]; // Compute volumetric strain from density

        if (layer == 0) {
            // For first layer, set hoop stretch from a_mid_guess
            lambda_t = a_mid_guess_inner / curr_vessel->layers[layer].a_mid_h; // L: Circumferential stretch a_mid/a_mid_h
            a_mid = a_mid_guess_inner;
        } else {
            // For subsequent layers, compute hoop stretch by solving the equation  A lambdat^2 + B lambdat + C = 0
            // This equation is derived as follows: 
            // a + h = a_new
            // a + h = a_mid_new - h_new / 2
            // a + h = (a_h + h_h / 2) * lambda_t - h_h * lambda_r / 2
            // a + h = (a_h + h_h / 2) * lambda_t - h_h * J_s / lambda_t / lambda_z / 2
            // 0 = (a_h + h_h / 2) * lambda_t^2 - (a + h) * lambda_t - h_h * J_s / 2 / lambda_z
            double h_h = curr_vessel->layers[layer].h_h;
            double a_h = curr_vessel->layers[layer].a_h;
            bool now_contacting = a + h >= a_h && sn > 0; // not allowed to begin contact at t = 0
            if (!now_contacting && !curr_vessel->layers[layer].is_contacting) {
                if (sn > 0) {
                    curr_vessel->layers[layer].a_mid[sn] = curr_vessel->layers[layer].a_mid_h; // set to 0 for debugging purposes
                    curr_vessel->layers[layer].a[sn] = curr_vessel->layers[layer].a_h;
                    curr_vessel->layers[layer].h[sn] = curr_vessel->layers[layer].h_h;
                }
                break; // if inner layer does not contact current layer, then ignore it
            }
            double A = a_h + h_h / 2;
            double B = -(a + h);
            double C = -h_h / 2 * J_s / lambda_z;
            lambda_t = (-B  + sqrt(B * B - 4 * A * C)) / (2 * A);
            a_mid = lambda_t * (a_h + h_h / 2);
            // printf("a_h:%f, h_h: %f, lambda_t:%f\n", a_h, h_h, lambda_t);

        }
        // Compute a and h from stretches
        h = J_s / (lambda_t * lambda_z) * curr_vessel->layers[layer].h_h;
        a = a_mid - h / 2;
        // printf("a:%f, h:%f\n", a, h);


        //Update vessel geometry before running update_sigma
        curr_vessel->layers[layer].a_mid[sn] = a_mid;
        curr_vessel->layers[layer].a[sn] = a;
        curr_vessel->layers[layer].h[sn] = h;
        if (equil_check) {
            for (int alpha = 0; alpha < curr_vessel->layers[layer].constituents.size(); alpha ++) {
                curr_vessel->layers[layer].constituents[alpha].a_act[sn] = a;
            }
            curr_vessel->layers[layer].lambda_th_curr = 1.0; // If homeostatic, then hoop stretch is zero
        } else {
           curr_vessel->layers[layer].lambda_th_curr = lambda_t; // otherwise saze stretch
        }
        //Update WSS from Q Flow
        // If not the inner most layer, copy wall shear stress from inner most layer
        if (!equil_check && curr_vessel->wss_calc_flag > 0) {
            if (layer == 0) {
                mu = get_app_visc(curr_vessel, sn);
                curr_vessel->bar_tauw = 4 * mu * curr_vessel->Q[sn] / (3.14159265 * pow(a * 100, 3));
            } 
        }
        update_sigma(&curr_vessel->layers[layer]);
        pa_calc += h * curr_vessel->layers[layer].sigma[1];
        fz_total += M_PI * h * (2 * a + h) * curr_vessel->layers[layer].sigma[2];
        h_total += h;

    }

    //Calculating sigma_t_th from pressure P
    double pa_th = curr_vessel->P[sn] * curr_vessel->layers[0].a[sn];
        // Compute sigma invariant
    curr_vessel->f = fz_total;

    double J = pa_calc - pa_th;

    // printf("pa: %f, %f", pa_calc, pa_th);
#ifdef COMBINED_DELTA_SIGMA

    double sigma_inv = pa_th / h_total + fz_total / (M_PI * h_total * (2 * (a + h) - h_total));
    for (int layer = 0; layer < n_layers; layer++) {
        curr_vessel->layers[layer].sigma_inv = sigma_inv;
    }

#endif

    return J;
}

void update_sigma(void* curr_layer) {

    struct layer *curr_lay = (struct layer*) curr_layer;
    int n_alpha = curr_lay->constituents.size();
    //Get current time index
    double s = curr_lay->parent_vessel->s;
    int sn = curr_lay->parent_vessel->sn;
    int nts = curr_lay->parent_vessel->nts;
    double dt = curr_lay->parent_vessel->dt;
    int taun_min = 0;
    double tau_max = nts + 1; 
    // for (int alpha = 0; alpha < n_alpha; alpha++) {
    //     double k_alpha_h = curr_lay->constituents[alpha].k_alpha_h;
    //     if (k_alpha_h > 0 && 10 / k_alpha_h > tau_max) {
    //         tau_max = 10 / k_alpha_h;
    //     }
    // }
    //Specify vessel geometry
    double a0 = curr_lay->a[0];
    double h0 = curr_lay->h[0];

    //Calculate vessel stretches
    double lambda_th_s = curr_lay->lambda_th_curr;
    double lambda_z_s = curr_lay->parent_vessel->lambda_z_tau[sn] / curr_lay->parent_vessel->lambda_z_tau[0];


    //Calculate constituent specific stretches for evolving constituents at the current time
    vector<double> lambda_alpha_s(n_alpha, 0);
    double eta_alpha = 0;
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Check to see if constituent is isotropic
        eta_alpha = curr_lay->constituents[alpha].eta_alpha_h;

        if (eta_alpha >= 0) {

            //Stretch is equal to the sqrt of I4
            lambda_alpha_s[alpha] = sqrt(pow(lambda_z_s * cos(eta_alpha), 2)
                + pow(lambda_th_s * sin(eta_alpha), 2));

            //Update stored current stretch if not numerical experiment

            if (curr_lay->parent_vessel->num_exp_flag == 0) {

                curr_lay->constituents[alpha].lambda_alpha_tau[sn] = lambda_alpha_s[alpha];
            }
        }

    }  
    //Find the current deformation gradient
    double J_s = curr_lay->rhoR[sn] / curr_lay->rho[sn];
    double F_s[3] = { J_s / (lambda_th_s * lambda_z_s), lambda_th_s, lambda_z_s };

    //Find the mechanical contributions of each constituent for each direction
    double a, h;
    double lambda_th_tau = 0;
    double lambda_z_tau = 0; //Assume constant axial stretch
    double J_tau = 1;
    double F_tau[3] = { 1, 1, 1 };
    double lambda_alpha_ntau_s = 0;
    double Q1 = 0, Q2 = 0;
    double F_alpha_ntau_s = 0;
    double hat_S_alpha = 0;
    double sigma[3] = { 0 };
    double lagrange = 0;
    double pol_mod = 0;

    //Local active variables
    double C = 0;
    double lambda_act = 0;
    double parab_act = 0;
    double hat_sigma_act = 0, sigma_act = 0;

    //Stiffness variables
    double hat_dSdC_alpha = 0;
    double hat_dSdC_act = 0;
    double Cbar[3] = { 0 };
    double Cbar_act = 0;
    vector<double> constitutive_return = { 0, 0 };

    //Integration variables
    //For mass
    double mq_0 = 0, mq_1 = 0, mq_2 = 0;
    double q_0 = 1.0, q_1 = 1.0, q_2 = 1.0;
    double k_0 = 0, k_1 = 0, k_2 = 0;

    int n = 0; //number of pts in integration interval

    //For stress
    vector<double> hat_sigma_0 = { 0, 0, 0 }, hat_sigma_1 = { 0, 0, 0 }, hat_sigma_2 = { 0, 0, 0 };
    //For active stress
    double a_act = 0;
    double q_act_0 = 0, q_act_1 = 0, q_act_2 = 0;
    double a_0 = 0, a_1 = 0, a_2 = 0;

    //For stiffness
    vector<double> hat_Cbar_0 = { 0, 0, 0 }, hat_Cbar_1 = { 0, 0, 0 }, hat_Cbar_2 = { 0, 0, 0 };

    //Boolean for checks
    bool deg_check = 0;

    //Determine if beyond initial time history
    if (s > tau_max) {
        taun_min = sn - int(tau_max / dt);
    }
    else {
        taun_min = 0;
    }

    n = (sn - taun_min) + 1;; //number of integration pts
    bool even_n = n % 2 == 0;

    //Similar integration to that used for kinematics
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        double k_act = curr_lay->constituents[alpha].k_act;
        //Trapz rule allows for fast heredity integral evaluation
        k_2 = curr_lay->constituents[alpha].k_alpha[sn];
        q_2 = 1.0;
        mq_2 = curr_lay->constituents[alpha].mR_alpha[sn];

        //Find active radius from current cohort
        if (curr_lay->constituents[alpha].alpha_active == 1) {
            a_2 = curr_lay->a[sn];
            q_act_2 = 1.0;
        }

        //Kinematics
        F_tau[0] = F_s[0], F_tau[1] = F_s[1], F_tau[2] = F_s[2];
        J_tau = J_s;
        //Find stress from current cohort
        for (int dir = 0; dir < 3; dir++) {

            constitutive_return = constitutive(&curr_lay->constituents[alpha], lambda_alpha_s[alpha], sn);
            hat_S_alpha = constitutive_return[0];
            hat_dSdC_alpha = constitutive_return[1];
            F_alpha_ntau_s = F_s[dir] / F_tau[dir] * curr_lay->constituents[alpha].G_alpha_h[dir * nts + sn];
            hat_sigma_2[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;

            hat_Cbar_2[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;

        }

        //Boolean for whether the constituent increases ref mass density
        deg_check = curr_lay->constituents[alpha].mR_alpha_h > 0;

        //Check if during G&R or at initial time point
        if (sn > 0 && deg_check) {

            for (int taun = sn - 1; taun >= taun_min; taun = taun - 1) { // For simpson's rule: taun = taun - 2 and taun >= taun_min + 1

                //Find the 1st intermediate deformation gradient
                a = curr_lay->a[taun];
                h = curr_lay->h[taun];
                lambda_th_tau = (a + h / 2) / (a0 + h0 / 2);
                lambda_z_tau = curr_lay->parent_vessel->lambda_z_tau[taun];
                J_tau = curr_lay->rhoR[taun] / curr_lay->rho[taun];
                F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau);
                F_tau[1] = lambda_th_tau;
                F_tau[2] = lambda_z_tau;

                //Find 1st intermediate kinetics
                k_1 = curr_lay->constituents[alpha].k_alpha[taun];
                q_1 = exp(-(k_2 + k_1) * dt / 2) * q_2;
                //mRq -> referential
                mq_1 = curr_lay->constituents[alpha].mR_alpha[taun] * q_1;

                //Find intermediate active state
                if (curr_lay->constituents[alpha].alpha_active == 1) {
                    a_1 = a;
                    q_act_1 = exp(-k_act * dt) * q_act_2;
                }

                //Find 1st intermeidate stress component in each direction
                for (int dir = 0; dir < 3; dir++) {
                    constitutive_return = constitutive(&curr_lay->constituents[alpha], lambda_alpha_s[alpha], taun);
                    hat_S_alpha = constitutive_return[0];
                    hat_dSdC_alpha = constitutive_return[1];
                    F_alpha_ntau_s = F_s[dir] / F_tau[dir] * curr_lay->constituents[alpha].G_alpha_h[dir * nts + taun];
                    hat_sigma_1[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                    hat_Cbar_1[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;
                }

                //Find the 2nd intermediate deformation gradient (for Simpson's rule only)
                // a = curr_lay->a[taun - 1];
                // h = curr_lay->h[taun - 1];
                // lambda_th_tau = (a + h / 2) / (a0 + h0 / 2);
                // lambda_z_tau = curr_lay->parent_vessel->lambda_z_tau[taun - 1];
                // J_tau = curr_lay->rhoR[taun - 1] / curr_lay->rho[taun - 1];;
                // F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau);
                // F_tau[1] = lambda_th_tau;
                // F_tau[2] = lambda_z_tau;

                // //Find 2nd intermediate kinetics
                // k_0 = curr_lay->constituents[alpha].k_alpha[taun - 1];
                // q_0 = exp(-(k_2 + 4 * k_1 + k_0) * dt / 3) * q_2;
                // mq_0 = curr_lay->constituents[alpha].mR_alpha[taun - 1] * q_0;

                // //Find intermediate active state
                // if (curr_lay->constituents[alpha].alpha_active == 1) {
                //     a_0 = a;
                //     q_act_0 = exp(-k_act * dt) * q_act_1;
                // }

                //Find component in each direction
                for (int dir = 0; dir < 3; dir++) {
                    
                    // If using Simpson's rule:
                    // constitutive_return = constitutive(&curr_lay->constituents[alpha], lambda_alpha_s[alpha], taun - 1);
                    // hat_S_alpha = constitutive_return[0];
                    // hat_dSdC_alpha = constitutive_return[1];
                    // F_alpha_ntau_s = F_s[dir] / F_tau[dir] * curr_lay->constituents[alpha].G_alpha_h[dir * nts + taun - 1];
                    // hat_sigma_0[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                    // hat_Cbar_0[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;

                    // //Add to the stress and stiffness contribution in the given direction
                    // sigma[dir] += (mq_2 * hat_sigma_2[dir] + 4 * mq_1 * hat_sigma_1[dir] + mq_0 * hat_sigma_0[dir])
                    // / curr_lay->constituents[alpha].rho_hat_alpha_h * dt / 3;
                    // Cbar[dir] += (mq_2 * hat_Cbar_2[dir] + 4 * mq_1 * hat_Cbar_1[dir] + mq_0 * hat_Cbar_0[dir])
                    // / curr_lay->constituents[alpha].rho_hat_alpha_h * dt / 3;
                    
                    //Add to the stress and stiffness contribution in the given direction
                    sigma[dir] += (mq_2 * hat_sigma_2[dir] + mq_1 * hat_sigma_1[dir])
                    / curr_lay->constituents[alpha].rho_hat_alpha_h * dt / 2;
                    Cbar[dir] += (mq_2 * hat_Cbar_2[dir] + mq_1 * hat_Cbar_1[dir])
                    / curr_lay->constituents[alpha].rho_hat_alpha_h * dt / 2;
                }

                //Store active vars for next iteration
                //Find intermediate active state
                if (curr_lay->constituents[alpha].alpha_active == 1) {
                    // a_act += k_act * (q_act_2 * a_2 + 4 * q_act_1 * a_1 + q_act_0 * a_0) * dt / 3;
                    a_act += k_act * (q_act_2 * a_2 + q_act_1 * a_1) * dt / 2;
                    a_2 = a_1;
                    q_act_2 = q_act_1;
                }

                //Store intermediate kinetics for next iteration (for Simpson's rule only)
                // k_2 = k_0;
                // q_2 = q_0;
                // mq_2 = mq_0;

                //Store intermediate kinetics for next iteration
                k_2 = k_1;
                q_2 = q_1;
                mq_2 = mq_1;

                //Store intermediate stress and stiffness for next iteration (for Simpson's rule only)
                // hat_sigma_2 = hat_sigma_0;
                // hat_Cbar_2 = hat_Cbar_0;

                //Store intermediate stress and stiffness for next iteration
                hat_sigma_2 = hat_sigma_1;
                hat_Cbar_2 = hat_Cbar_2;
            }

            //  if (even_n) { // for Simpson's rule only


            //     //Find the 2nd intermediate deformation gradient
            //     a = curr_lay->a[taun_min];
            //     h = curr_lay->h[taun_min];
            //     lambda_th_tau = (a + h / 2) / (a0 + h0 / 2);
            //     lambda_z_tau = curr_lay->parent_vessel->lambda_z_tau[taun_min];
            //     J_tau = curr_lay->rhoR[taun_min] / curr_lay->constituents[alpha].rho_hat_alpha_h;
            //     F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau);
            //     F_tau[1] = lambda_th_tau;
            //     F_tau[2] = lambda_z_tau;

            //     //Find 2nd intermediate kinetics
            //     k_0 = curr_lay->constituents[alpha].k_alpha[taun_min];
            //     q_0 = exp(-(k_2 + k_0) * dt / 2) * q_2;
            //     mq_0 = curr_lay->constituents[alpha].mR_alpha[taun_min] * q_0;

            //     //Find intermediate active state
            //     if (curr_lay->constituents[alpha].alpha_active == 1) {
            //         a_0 = a;
            //         q_act_0 = exp(-k_act * dt) * q_act_2;
            //     }

            //     //Find component in each direction
            //     for (int dir = 0; dir < 3; dir++) {

            //         constitutive_return = constitutive(&curr_lay->constituents[alpha], lambda_alpha_s[alpha], taun_min);
            //         hat_S_alpha = constitutive_return[0];
            //         hat_dSdC_alpha = constitutive_return[1];
            //         F_alpha_ntau_s = F_s[dir] / F_tau[dir] * curr_lay->constituents[alpha].G_alpha_h[dir * nts + taun_min];
            //         hat_sigma_0[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
            //         hat_Cbar_0[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;

            //         //Add to the stress and stiffness contribution in the given direction
            //         sigma[dir] += (mq_2 * hat_sigma_2[dir] + mq_0 * hat_sigma_0[dir])
            //         / curr_lay->constituents[alpha].rho_hat_alpha_h * dt / 2;
            //         Cbar[dir] += (mq_2 * hat_Cbar_2[dir] + mq_0 * hat_Cbar_0[dir])
            //         / curr_lay->constituents[alpha].rho_hat_alpha_h * dt / 2;

            //     }

            //     if (curr_lay->constituents[alpha].alpha_active == 1) {
            //         a_act += k_act * (q_act_2 * a_2 + q_act_0 * a_0) * dt / 2;
            //     }
            // }
 
            //Add in the stress and stiffness contributions of the initial material
            if (taun_min == 0) {
                for (int dir = 0; dir < 3; dir++) {
                    sigma[dir] += curr_lay->constituents[alpha].rhoR_alpha[0]
                    / curr_lay->constituents[alpha].rho_hat_alpha_h * q_1 * hat_sigma_1[dir]; // q_0 * hat_sigma_0[dir] for Simpson's rule
                    Cbar[dir] += curr_lay->constituents[alpha].rhoR_alpha[0]
                    / curr_lay->constituents[alpha].rho_hat_alpha_h * q_1 * hat_Cbar_1[dir]; // q_0 * hat_Cbar_0[dir] for Simpson's rule
                }
            }

        }
        //Initial time point and constituents with prescribed degradation profiles
        else {
            //Find stress from initial cohort          
            for (int dir = 0; dir < 3; dir++) {

                constitutive_return = constitutive(&curr_lay->constituents[alpha], lambda_alpha_s[alpha], 0);
                hat_S_alpha = constitutive_return[0];
                hat_dSdC_alpha = constitutive_return[1];
                F_alpha_ntau_s = F_s[dir] * curr_lay->constituents[alpha].G_alpha_h[dir * nts + sn];
                hat_sigma_2[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                hat_Cbar_2[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;
                // Print here: curr_lay->constituents[alpha].rhoR_alpha[sn] / curr_lay->constituents[alpha].rho_hat_alpha_h * hat_sigma_2[dir];
                sigma[dir] += curr_lay->constituents[alpha].rhoR_alpha[sn] /
                curr_lay->constituents[alpha].rho_hat_alpha_h * hat_sigma_2[dir];
                Cbar[dir] += curr_lay->constituents[alpha].rhoR_alpha[sn] /
                curr_lay->constituents[alpha].rho_hat_alpha_h * hat_Cbar_2[dir];
            }

        }

        if (taun_min == 0 && curr_lay->constituents[alpha].alpha_active == 1) {
            a_act += curr_lay->constituents[alpha].a_act[0] * q_act_1; // q_act_0 for Simpson's rule
        }


    }
    //Find active stress contribtion
    //add in initial active stress radius contribution
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        if (curr_lay->constituents[alpha].alpha_active == 1) {
            C = curr_lay->constituents[alpha].CB -
            curr_lay->constituents[alpha].CS * (curr_lay->parent_vessel->bar_tauw /
                curr_lay->parent_vessel->bar_tauw_h - 1);

            lambda_act = curr_lay->a[sn] / curr_lay->constituents[alpha].a_act[sn];

            if (sn == 0) {
                a_act = curr_lay->constituents[alpha].a_act[0];
                C = curr_lay->constituents[alpha].CB;
                lambda_act = 1.0;
            }

            parab_act = 1 - pow((curr_lay->constituents[alpha].lambda_m - lambda_act) /
                (curr_lay->constituents[alpha].lambda_m - curr_lay->constituents[alpha].lambda_0), 2);

            hat_sigma_act = curr_lay->constituents[alpha].T_act * (1 - exp(-pow(C, 2))) * lambda_act * parab_act;

            sigma_act += curr_lay->constituents[alpha].rhoR_alpha[sn] / J_s / 
            curr_lay->rhoR_h * hat_sigma_act;

            hat_dSdC_act = curr_lay->constituents[alpha].T_act * (pow(lambda_act, -2) / 2 * 
             ((curr_lay->constituents[alpha].lambda_m - lambda_act) / 
                 pow(curr_lay->constituents[alpha].lambda_m - curr_lay->constituents[alpha].lambda_0, 2)) 
             - pow(lambda_act, -3) / 4 * (parab_act));

            Cbar_act += curr_lay->constituents[alpha].rhoR_alpha[sn] / J_s / curr_lay->rhoR_h * 
            lambda_act * lambda_act * lambda_act * lambda_act * hat_dSdC_act;
            curr_lay->constituents[alpha].a_act[sn] = a_act;
        }
    }
    //The Lagrange multiplier is the radial stress component
    //subtract from each direction 
    // use boundary condition sigmarr = P * (amidh - ainner) / (aouter - ainner), same as P/2 for single layer
    int n_layers = curr_lay->parent_vessel->layers.size();
    double inner_radius_h = curr_lay->parent_vessel->layers[0].a_h;
    double outer_radius_h = curr_lay->parent_vessel->layers[n_layers - 1].a_h + curr_lay->parent_vessel->layers[n_layers - 1].h_h;
    double a_mid_h = curr_lay->a_mid_h;
    lagrange = sigma[0];
    #ifndef ZERO_LAGRANGE_BC
        lagrange += (a_mid_h - inner_radius_h) / (outer_radius_h - inner_radius_h) * curr_lay->parent_vessel->P[sn];
    #endif
    for (int dir = 0; dir < 3; dir++) {

        //Accounting for active stress
        if (dir == 1) {
            sigma[dir] += sigma_act;
            Cbar[dir] += Cbar_act;
        }
        //Calculating stiffness using extra stresses
        Cbar[dir] = 2 * sigma[dir] + 2 * Cbar[dir];
        //Calculating full cauchy stress
        sigma[dir] = sigma[dir] - lagrange;
        
        curr_lay->sigma[dir] = sigma[dir];

        curr_lay->Cbar[dir] = Cbar[dir];
    }
    //Save updated active radius (L: it's not active radius anymore?)
    curr_lay->sigma_inv = sigma[0]+sigma[1]+sigma[2];



}

vector<double> constitutive(void* curr_constituent, double lambda_alpha_s, int ts) {
    struct constituent *curr_const = (struct constituent *) curr_constituent;

    double lambda_alpha_ntau_s = 0;
    double Q1 = 0;
    double Q2 = 0;
    double hat_S_alpha = 0;
    double hat_dS_dlambda2_alpha = 0;
    double pol_mod = 0;
    double epsilon_curr = 0;
    
    int nts = curr_const->parent_layer->parent_vessel->nts;
    int sn = curr_const->parent_layer->parent_vessel->sn;
    
    vector<double> return_constitutive = { 0, 0 };

    //Check if anisotropic
    if (curr_const->eta_alpha_h >= 0) {

        lambda_alpha_ntau_s = curr_const->g_alpha_h[ts] * lambda_alpha_s / curr_const->lambda_alpha_tau[ts];

        if (lambda_alpha_ntau_s < 1) {
            lambda_alpha_ntau_s = 1;
        }

        Q1 = (pow(lambda_alpha_ntau_s, 2) - 1);
        Q2 = curr_const->c_alpha_h[1 * nts + ts] * pow(Q1, 2);
        hat_S_alpha = curr_const->c_alpha_h[0 * nts + ts] * Q1 * exp(Q2);
        hat_dS_dlambda2_alpha = curr_const->c_alpha_h[0 * nts + ts] * exp(Q2) * (1 + 2 * Q2);
        if (sn == 0) {
        }

    }
    else {

        if (curr_const->is_pol > 0) {

            if (curr_const->epsilon_alpha[sn] < curr_const->epsilon_pol_min) {
                curr_const->epsilon_pol_min = curr_const->epsilon_alpha[sn];
            }
            epsilon_curr = curr_const->epsilon_pol_min;
            pol_mod = 0.03 * pow(epsilon_curr, 2);
        }
        else {
            pol_mod = 1;
        }

        hat_S_alpha = pol_mod * curr_const->c_alpha_h[0 * nts + ts];
    }

    return_constitutive = { hat_S_alpha , hat_dS_dlambda2_alpha };
    return return_constitutive;

}

double get_app_visc(void* curr_vessel, int sn){
    //Returns apparent viscosity if diameter dependent factors from Secomb 2017 are in effect
    //Otherwise returns default for blood
    double d = 0.0;
    double mu = 0.00;

    if (((struct vessel*)curr_vessel)->app_visc_flag == 1){
        d = ((struct vessel*)curr_vessel)->layers[0].a[sn] * 2 * 1000000;
        mu = (1+(6*exp(-0.0858*d)+3.2-2.44*exp(-0.06*pow(d,0.645))-1)*pow(d/(d-1.1),2)*pow(d/(d-1.1),2)) * 0.0124;
    }
    else{
        mu = 0.04;
    }

    return mu;
}
