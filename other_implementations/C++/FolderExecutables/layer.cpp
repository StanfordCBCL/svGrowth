// layer.cpp
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include "json.hpp"
#include "layer.h"
#include "vessel.h"
#include "constituent.h"
#include "functions.h"

using std::string;
using std::vector;
using std::cout;

using json = nlohmann::json;

layer::layer() { } //Default constructor

void layer::initializeJSON(json& json_in, vessel *parent) { 

    parent_vessel = parent;
    int nts = parent->nts;
    a_h = json_in["inner_radius_homeostatic"];
    a_h = a_h * parent-> mm_to_m;
    h_h = json_in["thickness_homeostatic"];
    h_h = h_h * parent-> mm_to_m;
    rho_hat_h = json_in["whole_mass_density"];

    // Computed values
    a_mid_h = a_h + h_h / 2; //in vivo reference wall mid-point
    rhoR_h = rho_hat_h;
    a_e = 0.0;
    h_e = 0.0;
    lambda_th_curr = 1.0;
    J_di = 1.0;

    // Histories
    a = { 0.0 }, a_mid = { 0.0 }, h = { 0.0 }, A = { 0.0 }, A_mid = { 0.0 }, H = { 0.0 };
    a_pas = { 0.0 }, h_pas = {0.0}, a_mid_pas = {0.0};
    lambda_z_pre = { 0.0 }, rho = { 0.0 }, rhoR = { 0.0 }, ups_infl_p = { 0.0 }, ups_infl_d = {0.0};
    a.resize(nts);
    a_mid.resize(nts);
    h.resize(nts);
    a_pas.resize(nts);
    h_pas.resize(nts);
    a_mid_pas.resize(nts);
    A.resize(nts);
    A_mid.resize(nts);
    H.resize(nts);
    lambda_z_pre.resize(nts);
    rho.resize(nts);
    rhoR.resize(nts);
    ups_infl_p.resize(nts);
    ups_infl_d.resize(nts);
    a[0] = a_h;
    a_mid[0] = a_mid_h;
    h[0] = h_h;
    rho[0] = rho_hat_h;
    rhoR[0] = rho_hat_h;

    is_contacting = true;
    if (json_in.contains("initial_contact")) is_contacting = json_in["initial_contact"].get<bool>();
    
    // Stresses and stiffnesses
    sigma = {0.0}, sigma_h = {0.0}, Cbar = {0.0}, CC = {0.0};
    sigma.resize(3);
    sigma_h.resize(3);
    Cbar.resize(3);
    CC.resize(36);
    // update_sigma(this);
    sigma_h = sigma;
    sigma_inv_h = sigma[0] + sigma[1] + sigma[2];

    // Handshake stretch history
    F_s = {0.0}, Fi_curr = {0.0}, F_curr = {0.0};
    F_s.resize(9*nts);
    Fi_curr.resize(9);
    F_curr.resize(9);
    //Initialize to no stretch
    F_s[0] = 1.0;
    F_s[4] = 1.0;
    F_s[8] = 1.0;
    F_curr[0] = 1.0;
    F_curr[4] = 1.0;
    F_curr[8] = 1.0;
    Fi_curr[0] = 1.0;
    Fi_curr[4] = 1.0;
    Fi_curr[8] = 1.0;

    // Read constituents        
    json constituents_in = json_in["constituents"];
    constituents = {};
    int n_alpha = constituents_in.size();
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        constituent next_const;
        constituents.push_back(next_const);
        constituents[alpha].initializeJSON(constituents_in[alpha], this);
    } 
} //Default constructor

