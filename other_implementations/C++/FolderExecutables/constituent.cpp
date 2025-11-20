// constituent.cpp
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include "json.hpp"
#include "vessel.h"
#include "layer.h"
#include "constituent.h"
#include "functions.h"

using std::string;
using std::vector;
using std::cout;

using json = nlohmann::json;

constituent::constituent() { } //Default constructor
void constituent::initializeJSON(json& json_in, layer* parent_l) { 
    ups_saved = 0.0;

    // Do not save references to parent yet, do that from vessel.cpp once everything is in the appropriate list
    parent_layer = parent_l;
    int nts = parent_l->parent_vessel->nts;
    double dt = parent_l->parent_vessel->dt;
    phi_alpha_h = json_in["mass_ratio"];
    string model = json_in["constituent_model"];
    c_alpha_h = {};
    g_alpha_h = {};
    G_alpha_h = {};
    if (model.compare("neohookean") == 0) {
        parent_l->parent_vessel->evaluate_expr(json_in["c1"], c_alpha_h, dt, nts);
        g_alpha_h = {0.0};
        g_alpha_h.resize(nts);
        eta_alpha_h = -1.0;
        // Appending prestretches to G_alpha_h for each direction
        parent_l->parent_vessel->evaluate_expr(json_in["prestretch_r"], G_alpha_h, dt, nts);
        parent_l->parent_vessel->evaluate_expr(json_in["prestretch_th"], G_alpha_h, dt, nts);   
        parent_l->parent_vessel->evaluate_expr(json_in["prestretch_z"], G_alpha_h, dt, nts);      
    } else if (model.compare("fung") == 0) {
        parent_l->parent_vessel->evaluate_expr(json_in["c1"], c_alpha_h, dt, nts);
        parent_l->parent_vessel->evaluate_expr(json_in["c2"], c_alpha_h, dt, nts);
        parent_l->parent_vessel->evaluate_expr(json_in["prestretch"], g_alpha_h, dt, nts);
        eta_alpha_h = json_in["orientation"];
        eta_alpha_h = eta_alpha_h * M_PI / 180.0;
        G_alpha_h = {0.0};
        G_alpha_h.resize(nts * 3);
        for (int sn = 0; sn < nts; sn ++) {
            G_alpha_h[0 * nts + sn] = 0.0;
            G_alpha_h[1 * nts + sn] = g_alpha_h[sn] * sin(eta_alpha_h);
            G_alpha_h[2 * nts + sn] = g_alpha_h[sn] * cos(eta_alpha_h);

        }
    }

    K_sigma_p_alpha_h = 0; 
    K_sigma_d_alpha_h = 0;
    K_tauw_p_alpha_h = 0;
    K_tauw_d_alpha_h = 0;
    k_alpha_h = 0;
    if (json_in.contains("degradation")) k_alpha_h = json_in["degradation"];
    if (json_in.contains("stress_mediated_production")) K_sigma_p_alpha_h = json_in["stress_mediated_production"];
    if (json_in.contains("stress_mediated_degradation")) K_sigma_d_alpha_h = json_in["stress_mediated_degradation"];
    if (json_in.contains("wss_mediated_production")) K_tauw_p_alpha_h = json_in["wss_mediated_production"];
    if (json_in.contains("wss_mediated_degradation")) K_tauw_d_alpha_h = json_in["wss_mediated_degradation"];

    if (json_in.contains("degradation_profile")) {
        degradation_profile = {};
        parent_l->parent_vessel->evaluate_expr(json_in["degradation_profile"], degradation_profile, dt, nts);
    } else {
        degradation_profile = {0.0};
        degradation_profile.resize(nts);
        for (int i = 0; i < nts; i ++) degradation_profile[i] = 1.0;
    }


    if (json_in.contains("inflammation_production")) {
        K_infl_p_alpha = {};
        parent_l->parent_vessel->evaluate_expr(json_in["inflammation_production"], K_infl_p_alpha, dt, nts);
    } else {
        K_infl_p_alpha = {0.0};
        K_infl_p_alpha.resize(nts);
    }

    if (json_in["is_active"]) {
        alpha_active = 1;
        T_act_h = json_in["active_stress_max_homeostatic"];
        k_act = json_in["active_remodeling_rate"];
        lambda_0 = json_in["active_stretch_min"];
        lambda_m = json_in["active_stretch_max"];
        CB = json_in["active_cb"];
        CS = json_in["active_cs"];
    } else {
        alpha_active = 0;
        T_act_h = 0;
        k_act = 0;
        lambda_0 = 0;
        lambda_m = 1;
        CB = 0;
        CS = 0;
    }

    // Set TEVG parameters to 0
    alpha_infl = 0;
    epsilon_pol_min = 0; 
    epsilonR_alpha_0 = 0; 


    // Computed values
    rhoR_alpha_h = phi_alpha_h * parent_l->rho_hat_h;
    rho_hat_alpha_h = parent_l->rho_hat_h;
    mR_alpha_h = rhoR_alpha_h * k_alpha_h;
    T_act = T_act_h;
    rho_alpha_e = 0;
    is_pol = 0;

    //Histories
    rhoR_alpha = { 0 }, mR_alpha = { 0 }, k_alpha = { 0 }, lambda_act = {0}, a_act = {0};
    epsilonR_alpha = { 0 }, epsilon_alpha = { 0 }, lambda_alpha_tau = { 0 };

    rhoR_alpha.resize(nts); //referential apparent mass densities (time history)
    epsilon_alpha.resize(nts); //current volume fractions
    epsilonR_alpha.resize(nts); //referential volume fractions
    mR_alpha.resize(nts); //referential mass production rate (time history)
    k_alpha.resize(nts); //mass removal decay (time history)
    lambda_alpha_tau.resize(nts);
    lambda_act.resize(nts);
    a_act.resize(nts);

    a_act[0] = parent_l->a_h;
    lambda_act[0] = 1.0;
    mR_alpha[0] = mR_alpha_h;
    k_alpha[0] = k_alpha_h;
    epsilonR_alpha[0] = phi_alpha_h;
    epsilon_alpha[0] = phi_alpha_h;
    rhoR_alpha[0] = rhoR_alpha_h;
    lambda_alpha_tau[0] = 1.0;

} //Default constructor
