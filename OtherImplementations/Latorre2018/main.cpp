//Reproduces Figure 1 from Lattore and Humphrey 2018, ZAMM
//NOTE: Stiffness calc does not include active stress

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

using std::string;
using std::vector;
using std::cout;

struct vessel
{
    string vessel_name = "default";

    //Time variables
    int nts = 0; //total number of time steps
    double dt = 0; //time increment
    int sn = 0; //current time step index
    double s = 0; //actual current time

    //Geometric quantities
    double A_h = 0, B_h = 0, H_h = 0; //Traction-free reference
    double A_mid_h = 0; //Tf Midpoint reference
    double a_h = 0, b_h = 0, h_h = 0; //In-vivo reference
    double a_mid_h = 0; //In vivo midpoint
    double lambda_z_h = 0; //Reference in vivo stretch

    vector<double> a = { 0 }, a_mid = { 0 }, b = { 0 }, h = { 0 }; //Evolving in vivo reference
    vector<double> A = { 0 }, A_mid = { 0 }, B = { 0 }, H = { 0 }, lambda_z_pre = { 0 }; //Evolving traction free

    //Number of constituents
    int n_alpha = 0;

    //Material properties
    vector<double> c_alpha_h = { 0 }, eta_alpha_h = { 0 }, g_alpha_h = { 0 }, G_alpha_h = { 0 };

    //Mass fractions, referential densities, kinetic quantities
    vector<double> phi_alpha_h = { 0 }, rhoR_alpha_h = { 0 }, mR_alpha_h = { 0 }, k_alpha_h = { 0 };
    vector<double> K_sigma_p_alpha_h = { 0 }, K_sigma_d_alpha_h = { 0 }, K_tauw_p_alpha_h = { 0 }, K_tauw_d_alpha_h = { 0 };
    double rhoR_h = 0;

    //Histories
    vector<double> rhoR = { 0 }, rhoR_alpha = { 0 }, mR_alpha = { 0 }, k_alpha = { 0 };

    //Reference loading quantities
    double P_h = 0, f_h = 0, bar_tauw_h = 0, Q_h = 0;
    vector<double> sigma_h = { 0 };

    //Current loading quantities
    double lambda_th_curr = 0, lambda_z_curr = 0;
    double P = 0, f = 0, bar_tauw = 0, Q = 0;
    vector<double> sigma = { 0 }, Cbar = { 0 }, lambda_alpha_tau = { 0 }, lambda_z_tau = { 0 };
    double mb_equil = 0; //Current mechanobiological equil. state

    //Active stress quantities
    vector<int> alpha_active = { 0 }; //boolean for active constituents
    vector<double> a_act = { 0 }; //active radius history
    double T_act = 0, T_act_h = 0; //Max active stress, homeostatic max active stress
    double k_act = 0; //Active remodeling parameters
    double lambda_0 = 0; //Min contractile stretch
    double lambda_m = 0; //Max contractile stretch
    double CB = 0; //Basal VC to VD ratio
    double CS = 0; //Scaling factor for VC to VD ratio

    //Mechanobiologically equilibrated quantities
    double a_e = 0; //equilibrated radius
    double h_e = 0; //equilibrated thickness
    double rho_c_e = 0; //equilibrated collagen density
    double rho_m_e = 0; //equilibrated smc density
    double f_z_e = 0; //equilibrated axial force
    double mb_equil_e = 0; //Current mechanobiological equil. state

    //Flags
    int num_exp_flag = 0; //indicates whether doing reg G&R step or a numerical experiment
    
};

void update_time_step(vessel &curr_vessel);
int find_tf_geom(void *curr_vessel);
int tf_obj_f(const gsl_vector * x, void* curr_vessel, gsl_vector* f);
int find_equil_geom(void* curr_vessel);
int equil_obj_f(const gsl_vector* x, void* curr_vessel, gsl_vector* f);
int print_state_mr(size_t iter, gsl_multiroot_fsolver* s);
int find_iv_geom(void* curr_vessel);
double iv_obj_f(double a_mid_guess, void* curr_vessel);
void update_kinetics(vessel &curr_vessel);
void update_sigma(void *curr_vessel);
vector<double> constitutive(void* curr_vessel, double lambda_alpha_s, int alpha, int ts, int dir);

int main()
{
    //All units meters, kg, days
    //Conversions
    double mm_to_m = pow(10, -3);
    double kPa_to_Pa = pow(10, 3);

    //Initilize the current vessel
    //Store quantities for the current vessel being simulated
    vessel curr_vessel;
    string vessel_name = "Cerebral"; //Type of vessel simulated
    curr_vessel.vessel_name = vessel_name;

    //Time parameters
    double n_days = 701.00; //days simulated
    curr_vessel.dt = 1.0; //time step size
    curr_vessel.nts = int(n_days / curr_vessel.dt); //number of G&R time steps
    curr_vessel.sn = 0; //Initialize current time index to zero;
    curr_vessel.s = 0; //Initialize the actual current time to zero

    //Geometric parameters
    //From Table 1, Latorre and Humphrey 2018, APL
    curr_vessel.a_h = 1.40 * mm_to_m; //in vivo reference inner radius
    curr_vessel.h_h = 0.12 * mm_to_m; //in vivo reference medial thickness

    //Derived geometric parameters
    //curr_vessel.b_h = curr_vessel.a_h + curr_vessel.h_h; //in vivo reference outer diameter
    curr_vessel.a_mid_h = curr_vessel.a_h + curr_vessel.h_h / 2; //in vivo reference wall mid-point

    //Axial stretch with in vivo reference
    curr_vessel.lambda_z_h = 1.0; //Reference in vivo stretch

    //Initialize the loaded geometry history
    curr_vessel.a.resize(curr_vessel.nts); //loaded inner radius history
    curr_vessel.a[0] = curr_vessel.a_h;
    //curr_vessel.b.resize(curr_vessel.nts); //loaded outer radius history
    //curr_vessel.b[0] = curr_vessel.b_h;
    curr_vessel.a_mid.resize(curr_vessel.nts); //loaded mid-radial history
    curr_vessel.a_mid[0] = curr_vessel.a_mid_h;
    curr_vessel.h.resize(curr_vessel.nts); //loaded thickness history
    curr_vessel.h[0] = curr_vessel.h_h;

    //Initialize the traction free geometry history
    //These are unkown with in vivo reference until solved for
    curr_vessel.A.resize(curr_vessel.nts); //TF inner radius history
    //curr_vessel.B.resize(curr_vessel.nts); //TF outer radius history
    curr_vessel.A_mid.resize(curr_vessel.nts); //TF mid-radial history
    curr_vessel.H.resize(curr_vessel.nts); //TF Thickness history
    curr_vessel.lambda_z_pre.resize(curr_vessel.nts); //Axial pre-stretch

    //Constituent material properties
    curr_vessel.n_alpha = 3; //number of constituents alpha

    //Medial and adventitial constituents have the same material behavior
    double c1_e = 70.6 * kPa_to_Pa;
    double c2_e = 0.0;
    double c1_m = 10.0 * kPa_to_Pa;
    double c2_m = 3.5;
    double c1_ct = 672.5 * kPa_to_Pa;
    double c2_ct = 22.0;

    curr_vessel.c_alpha_h = { c1_e, c2_e,
                             c1_m, c2_m,
                             c1_ct, c2_ct};

    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents
    double eta_e_h = -1.0;
    double eta_m_h = 90.0 * M_PI / 180.0;
    double eta_ct_h = 90.0 * M_PI / 180.0;

    curr_vessel.eta_alpha_h = { eta_e_h, eta_m_h, eta_ct_h };

    //Pre-stretch parameters
    double g_eh = 1.4;
    double g_mh = 1.30;
    double g_cth = 1.08;

    curr_vessel.g_alpha_h = { g_eh, g_mh, g_cth };

    double eta = 0;
    double g_alpha = 0;
    curr_vessel.G_alpha_h.resize(3 * curr_vessel.n_alpha);

    //Mass density parameters
    curr_vessel.rhoR_h = 1050.0; //density of tissue as a whole kg/m^3

    //Initializing time history storage
    curr_vessel.rhoR.resize(curr_vessel.nts); //Referential mass density history
    curr_vessel.rhoR[0] = 1050.0;
    curr_vessel.rhoR_alpha.resize(curr_vessel.n_alpha * curr_vessel.nts); //mass densities (time history)
    curr_vessel.mR_alpha.resize(curr_vessel.n_alpha * curr_vessel.nts); //mass production rate (time history)
    curr_vessel.k_alpha.resize(curr_vessel.n_alpha * curr_vessel.nts); //mass removal decay (time history)

    //Homeostatic mass fractions of constituents
    double phi_eh = 0.02; //elastin
    double phi_mh = 0.76; //SMCs
    double phi_cth = 0.22; //circ. col.

    curr_vessel.phi_alpha_h = { phi_eh, phi_mh, phi_cth };

    //Initial mass densities of constituents
    double rhoR_eh = phi_eh * curr_vessel.rhoR_h;
    double rhoR_mh = phi_mh * curr_vessel.rhoR_h;
    double rhoR_cth = phi_cth * curr_vessel.rhoR_h;

    curr_vessel.rhoR_alpha_h = { rhoR_eh, rhoR_mh, rhoR_cth };

    //Kinetic parameters
    //Degradation parameters
    double k_elh = 0.00; //days
    double k_mh = 1.0 / 14.0; //1.0 / 14.0; //days
    double k_cth = 1.0 / 10.0; //days
    curr_vessel.k_alpha_h = { k_elh, k_mh, k_cth };

    //Production parameters, intialized for each constituent in a loop below
    curr_vessel.mR_alpha_h.resize(curr_vessel.n_alpha);

    //Gain parameters
    //Stress mediated production
    double K_sigma_p_elh = 0.0;
    double K_sigma_p_mh = 2.0; // 2.0;
    double K_sigma_p_cth = 1.0; // 1.0
    curr_vessel.K_sigma_p_alpha_h = { K_sigma_p_elh, K_sigma_p_mh, K_sigma_p_cth };

    //Stress mediated degradation
    double K_sigma_d_elh = 0.0;
    double K_sigma_d_mh = 1.0;
    double K_sigma_d_cth = 1.0;
    curr_vessel.K_sigma_d_alpha_h = { K_sigma_d_elh, K_sigma_d_mh, K_sigma_d_cth };

    //Wall Shear Stress mediated production
    double K_tauw_p_elh = 0.0;
    double K_tauw_p_mh = 1.0; // 1.0;
    double K_tauw_p_cth = 0.5; // 0.5
    curr_vessel.K_tauw_p_alpha_h = { K_tauw_p_elh, K_tauw_p_mh, K_tauw_p_cth };

    //Wall Shear Stress mediated degradation
    double K_tauw_d_elh = 0.0;
    double K_tauw_d_mh = 0.0;
    double K_tauw_d_cth = 0.0;
    curr_vessel.K_tauw_d_alpha_h = { K_tauw_d_elh, K_tauw_d_mh, K_tauw_d_cth };

    //Initialize homeostatic loading variables
    //Pressure and Flow
    curr_vessel.P = curr_vessel.P_h = 14.18  * kPa_to_Pa; //Pressure in Pa
    curr_vessel.Q = curr_vessel.Q_h = 1.0;
    curr_vessel.bar_tauw = curr_vessel.bar_tauw_h = 1.0/pow(curr_vessel.a_h,3); //Wall shear stress in Pa

    //Initialize variables for load history
    curr_vessel.lambda_th_curr = 1.0;
    curr_vessel.lambda_z_curr = 1.0;
    curr_vessel.lambda_alpha_tau.resize(curr_vessel.nts * curr_vessel.n_alpha);
    curr_vessel.lambda_z_tau.resize(curr_vessel.nts); 
    curr_vessel.lambda_z_tau[0] = 1.0;

    //Initial stresses in each direction
    curr_vessel.sigma_h.resize(3);
    curr_vessel.sigma.resize(3);
    curr_vessel.Cbar.resize(3);

    //Active stress parameters
    curr_vessel.alpha_active = { 0, 1, 0 };
    curr_vessel.a_act.resize(curr_vessel.nts);
    curr_vessel.a_act[0] = curr_vessel.a_h;
    curr_vessel.T_act_h = 170.0 * kPa_to_Pa; //Max active stress
    curr_vessel.T_act = 170.0 * kPa_to_Pa; //Current active stress
    curr_vessel.k_act = 1.0 / 7.0; //Active remodeling parameters
    curr_vessel.lambda_0 = 0.4; //Min contractile stretch
    curr_vessel.lambda_m = 1.1; //Max contractile stretch
    curr_vessel.CB = 0.8326; //
    curr_vessel.CS = curr_vessel.CB / 2.0;

    //Initializion with looping through constituents
    for (int alpha = 0; alpha < curr_vessel.n_alpha; alpha++) {

        //Initialize the deposition tensor
        eta = curr_vessel.eta_alpha_h[alpha];
        g_alpha = curr_vessel.g_alpha_h[alpha];
        if (eta > 0) { //for anisotropic constituents
            curr_vessel.G_alpha_h[3 * alpha] = 0.0;
            curr_vessel.G_alpha_h[3 * alpha + 1] = g_alpha * sin(eta);
            curr_vessel.G_alpha_h[3 * alpha + 2] = g_alpha * cos(eta);
        }
        else { //for isotropic constituents
            curr_vessel.G_alpha_h[3 * alpha] = 1.0 / pow(g_alpha, 2);
            curr_vessel.G_alpha_h[3 * alpha + 1] = g_alpha;
            curr_vessel.G_alpha_h[3 * alpha + 2] = g_alpha;
        }

        //Initialize homeostatic mass productions
        curr_vessel.mR_alpha_h[alpha] = curr_vessel.k_alpha_h[alpha] * curr_vessel.rhoR_alpha_h[alpha];

        //Initialize time histories to their homeostatic values for native vessel
        curr_vessel.rhoR_alpha[curr_vessel.nts * alpha] = curr_vessel.rhoR_alpha_h[alpha];
        curr_vessel.mR_alpha[curr_vessel.nts * alpha] = curr_vessel.mR_alpha_h[alpha];
        curr_vessel.k_alpha[curr_vessel.nts * alpha] = curr_vessel.k_alpha_h[alpha];

        //Initilize stretch histories
        curr_vessel.lambda_alpha_tau[curr_vessel.nts * alpha] = 1.0;

    }

    //Find the true initial stress state
    int equil_check = find_iv_geom(&curr_vessel);
    printf("%s %f %s %f %s %f \n", "Inner radius: ", curr_vessel.a[0],
        "Thickness: ", curr_vessel.h[0],
        "Ref Density", curr_vessel.rhoR_alpha[curr_vessel.nts * 1]);

    //Store true homeostatic state
    update_sigma(&curr_vessel);
    
    //Store homeostatic stresses
    curr_vessel.sigma_h = curr_vessel.sigma;
    curr_vessel.bar_tauw_h = curr_vessel.Q_h / pow(curr_vessel.a[0], 3);
    curr_vessel.a_act[0] = curr_vessel.a[0];
    curr_vessel.f_h = M_PI * curr_vessel.h[0] * (2 * curr_vessel.a[0] + curr_vessel.h[0]) *
        curr_vessel.sigma[2];

    //---------------------------------------------------------------------------
    //Find the reference traction free geometry
    //curr_vessel.num_exp_flag = 1; //num exp flag on
    //int test_tf = find_tf_geom(&curr_vessel);
    curr_vessel.num_exp_flag = 0; //num exp flag off

    //curr_vessel.A_mid_h = curr_vessel.A_mid[0];
    //curr_vessel.H_h = curr_vessel.H[0];

    //Reset the initial stretch/stress conditions
    //curr_vessel.lambda_th_curr = 1.0;
    //curr_vessel.lambda_z_curr = 1.0;
    //curr_vessel.sigma = curr_vessel.sigma_h;
    //---------------------------------------------------------------------------

    //Setup file I/O for G&R output
    std::ofstream GnR_out("GnR_out.txt");

    //Find initial mechanobiological equilibrium state
    double eta_K = curr_vessel.K_sigma_p_alpha_h[1] /
        curr_vessel.K_tauw_p_alpha_h[1];

    curr_vessel.mb_equil = 1 + curr_vessel.K_sigma_p_alpha_h[2] * ((curr_vessel.sigma[1] + curr_vessel.sigma[2]) /
        (curr_vessel.sigma_h[1] + curr_vessel.sigma_h[2]) - 1) -
        curr_vessel.K_tauw_p_alpha_h[2] * ((curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1);

    //Write initial state to file
    GnR_out << curr_vessel.a[0] << "\t"<< curr_vessel.h[0] << "\t"
        << curr_vessel.rhoR_alpha[curr_vessel.nts * 1 + 0] << "\t" << curr_vessel.rhoR_alpha[curr_vessel.nts * 2 + 0] << "\t"
        << curr_vessel.mb_equil << "\t" << curr_vessel.P / curr_vessel.P_h << "\n";

    //Run the G&R time stespping
    double perturb_offset = 140.0; //when the mechanical perturbation occurs
    for (int sn = 1; sn < curr_vessel.nts; sn++) {

        curr_vessel.s = curr_vessel.dt * sn;

        //Mechanical perturbation to pressure, axial stretch, and flow
        if (curr_vessel.s > perturb_offset) {
            curr_vessel.P = curr_vessel.P_h * (1.0 + 0.5 * (1 - exp(-(curr_vessel.s - perturb_offset) / 10)));
            curr_vessel.lambda_z_curr = curr_vessel.lambda_z_h * (1.0 + 0.5 * (1 - exp(-(curr_vessel.s - perturb_offset) / 10)));
            curr_vessel.Q = curr_vessel.Q_h * (1.0 + 0.5 * (1 - exp(-(curr_vessel.s - perturb_offset) / 10)));
        }

        curr_vessel.sn = sn;
        update_time_step(curr_vessel);
        printf("%s \n", "---------------------------");

        //Write full model outputs
        GnR_out << curr_vessel.a[sn] << "\t" << curr_vessel.h[sn] << "\t"
            << curr_vessel.rhoR_alpha[curr_vessel.nts * 1 + sn] << "\t" << curr_vessel.rhoR_alpha[curr_vessel.nts * 2 + sn] << "\t"
            << curr_vessel.mb_equil << "\t" << curr_vessel.P / curr_vessel.P_h << "\n";

        //Store axial stretch history
        curr_vessel.lambda_z_tau[sn] = curr_vessel.lambda_z_curr;
    }

    //Find the equilibrated solution
    int equil_check2 = find_equil_geom(&curr_vessel);
    //Write equilibrated outputs
    GnR_out << curr_vessel.a_e << "\t" << curr_vessel.h_e << "\t"
        << curr_vessel.rho_m_e << "\t" << curr_vessel.rho_c_e << "\t"
        << curr_vessel.mb_equil_e << "\t" << curr_vessel.P / curr_vessel.P_h << "\n";

    //Print equilibrated state
    printf("%s %f %s %f \n", "a_e", curr_vessel.a_e, "h_e", curr_vessel.h_e);

    return 0;


}

void update_time_step(vessel &curr_vessel) {
    //Solves equilibrium equations at the current time point and updates kinetic variables
    //Find current time step
    int n_alpha = curr_vessel.n_alpha;
    int nts = curr_vessel.nts;
    int sn = curr_vessel.sn;
    double s = sn * curr_vessel.dt;
    double rhoR_s0 = 0, rhoR_s1 = 0; //Total mass at current step
    double tol = 1E-14; //Convergence tolerance
    int iter = 0, equil_check = 0;
    
    //Get initial mass density production prediction based on mass densities from previous time step
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        curr_vessel.rhoR_alpha[nts * alpha + sn] = curr_vessel.rhoR_alpha[nts * alpha + sn - 1];
    }
    //Get initial guess for radius
    curr_vessel.a_mid[sn] = curr_vessel.a_mid[sn - 1];
    curr_vessel.a_act[sn] = curr_vessel.a_act[sn - 1];
    update_kinetics(curr_vessel);
    equil_check = find_iv_geom(&curr_vessel);

    //Find real mass density production at the current time step iteratively
    double mass_check = 0.0;
    do {
        iter++;
        //Value from previous prediction
        rhoR_s0 = curr_vessel.rhoR[sn];

        //Update prediction
        update_kinetics(curr_vessel);
        equil_check = find_iv_geom(&curr_vessel);

        rhoR_s1 = curr_vessel.rhoR[sn];
        mass_check = abs((rhoR_s1 - rhoR_s0) / rhoR_s0);
    } while (mass_check > tol && iter < 100);

    //Find mechano-bio state
    double eta_K = curr_vessel.K_sigma_p_alpha_h[1] /
        curr_vessel.K_tauw_p_alpha_h[1];

    curr_vessel.mb_equil = 1 + curr_vessel.K_sigma_p_alpha_h[2] * ((curr_vessel.sigma[1] + curr_vessel.sigma[2]) /
        (curr_vessel.sigma_h[1] + curr_vessel.sigma_h[2]) - 1 ) -
        curr_vessel.K_tauw_p_alpha_h[2] * ((curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1 );

    //Print current state
    printf("%s %f %s %f %s %f %s %f\n", "Time:", s, "a: ", curr_vessel.a[sn], "h:", curr_vessel.h[sn],
        "Equil:", curr_vessel.mb_equil);

}

int find_tf_geom(void *curr_vessel) {
    //For geometries defined in the loaded configuration this code solves for the unloaded geometric 
    //variables at the current time point and stores them in the current v
    //vessel structure. Unloaded is also
    //referred to as traction free (but not stress free).

    //Get initial guesses from the loaded geometry
    double lambda_th_ul = 0.95;
    double lambda_z_ul = 0.95 * ((struct vessel*) curr_vessel)->lambda_z_curr;

    //Update vesseel loads to zero for traction free
    ((struct vessel*) curr_vessel)->P = 0.0;
    ((struct vessel*) curr_vessel)->f = 0.0;
    ((struct vessel*) curr_vessel)->T_act = 0.0;

    const gsl_multiroot_fsolver_type* T;
    gsl_multiroot_fsolver* s;

    int status;
    size_t iter = 0;

    const size_t n = 2;

    gsl_multiroot_function f = {&tf_obj_f, n, curr_vessel};
    double x_init[2] = { lambda_th_ul, lambda_z_ul };
    gsl_vector* x = gsl_vector_alloc(n);

    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, 2);

    gsl_multiroot_fsolver_set(s, &f, x);

    print_state_mr(iter, s);

    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        print_state_mr(iter, s);

        if (status)   // check if solver is stuck
            break;

        status = gsl_multiroot_test_residual(s->f, 1e-7);
    }
    while (status == GSL_CONTINUE && iter < 100);

    printf("status = %s\n", gsl_strerror(status));

    //Store the traction free results
    int sn = ((struct vessel*) curr_vessel)->sn;
    ((struct vessel*) curr_vessel)->A_mid[sn] = gsl_vector_get(s->x, 0) * ((struct vessel*) curr_vessel)->a_mid[sn];
    ((struct vessel*) curr_vessel)->lambda_z_pre[sn] = 1 / gsl_vector_get(s->x, 1);
    ((struct vessel*) curr_vessel)->H[sn] = 1.0 / (gsl_vector_get(s->x, 0) * gsl_vector_get(s->x, 0)) * ((struct vessel*) curr_vessel)->h[sn];

    //Return current vars to the loaded conditions
    ((struct vessel*) curr_vessel)->P = ((struct vessel*) curr_vessel)->P_h;
    ((struct vessel*) curr_vessel)->f = ((struct vessel*) curr_vessel)->f_h;
    ((struct vessel*) curr_vessel)->T_act = ((struct vessel*) curr_vessel)->T_act_h;
    ((struct vessel*) curr_vessel)->lambda_z_curr = ((struct vessel*) curr_vessel)->lambda_z_h;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;

}

int tf_obj_f(const gsl_vector* x, void *curr_vessel, gsl_vector* f) {

    //Seperate out inputs
    const double lambda_th_ul_guess = gsl_vector_get(x, 0);
    const double lambda_z_ul_guess = gsl_vector_get(x, 1);

    //Finds the difference in the theoretical stress from Laplace for deformed mixture
    //from the stress calculated from the mixture equations
    int sn = ((struct vessel*) curr_vessel) -> sn;

    //Reference config.
    double a_mid_0 = ((struct vessel*) curr_vessel)->a_mid[0];

    //Current loaded config
    double a_mid = ((struct vessel*) curr_vessel)->a_mid[sn];

    //Current stretches ref -> loaded
    double lambda_th_ref = a_mid / a_mid_0;
    double lambda_z_ref = ((struct vessel*) curr_vessel)->lambda_z_h;

    //double lambda_t_tf = A_mid_guess / a_mid;
    //double lambda_z_tf = lambda_ul_guess;

    //find det(F)
    double J_s = ((struct vessel*) curr_vessel)->rhoR[sn]  / ((struct vessel*) curr_vessel)->rhoR[0]; //Isochoric transient loading with volume changes from G&R

    //Update current total stretches
    ((struct vessel*) curr_vessel)->lambda_th_curr = lambda_th_ul_guess * lambda_th_ref;
    ((struct vessel*) curr_vessel)->lambda_z_curr = lambda_z_ul_guess * lambda_z_ref;

    update_sigma(curr_vessel);

    //Should be 0 for the traction-free configuration
    double J1 = ((struct vessel*) curr_vessel)->sigma[1];
    double J2 = ((struct vessel*) curr_vessel)->sigma[2];

    gsl_vector_set(f, 0, J1);
    gsl_vector_set(f, 1, J2);

    return GSL_SUCCESS;

}

int find_equil_geom(void* curr_vessel) {
    //Finds the mechanobiologically equilibrated geometry for a given set of loads inclduing
    //pressure, flow, and axial stretch with a set of G&R parameter values from the original
    //homeostatic state

    //Loading changes
    double gamma = ((struct vessel*) curr_vessel)->P 
        / ((struct vessel*) curr_vessel)->P_h; //Fold change in pressure from homeostatic
    double epsilon = ((struct vessel*) curr_vessel)->Q 
        / ((struct vessel*) curr_vessel)->Q_h; //Fold chance in flow from homeostatic
    double lambda = ((struct vessel*) curr_vessel)->lambda_z_curr 
        / ((struct vessel*) curr_vessel)->lambda_z_h; //Fold change in axial stretch from homeostatic

    //Homeostatic geometry
    double a_h = ((struct vessel*) curr_vessel)->a_h;
    double h_h = ((struct vessel*) curr_vessel)->h_h;

    //Initial guesses based on the loading changes
    double a_e_guess = pow(epsilon, 1.0 / 3.0) * a_h;
    double h_e_guess = gamma * pow(epsilon, 1.0 / 3.0) * h_h;
    double rho_c_e_guess = ((struct vessel*) curr_vessel)->rhoR_alpha_h[2];
    double f_z_e_guess = ((struct vessel*) curr_vessel)->f_h* (h_e_guess * (2 * a_e_guess + h_e_guess)) / (h_h * (2 * a_h + h_h));

    const gsl_multiroot_fsolver_type* T;
    gsl_multiroot_fsolver* s;

    int status;
    size_t iter = 0;

    const size_t n = 4;

    gsl_multiroot_function f = { &equil_obj_f, n, curr_vessel };
    double x_init[4] = {a_e_guess, h_e_guess, rho_c_e_guess, f_z_e_guess};
    gsl_vector* x = gsl_vector_alloc(n);

    for (int i = 0; i < n; i++) {
        gsl_vector_set(x, i, x_init[i]);
    }

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, n);

    gsl_multiroot_fsolver_set(s, &f, x);

    print_state_mr(iter, s);

    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        print_state_mr(iter, s);

        if (status)
            break;

        status = gsl_multiroot_test_residual(s->f, 1e-7);

    } while (status == GSL_CONTINUE && iter < 100);
    
    printf("status = %s \n", gsl_strerror(status));

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;

}

int equil_obj_f(const gsl_vector* x, void* curr_vessel, gsl_vector* f) {
    //Mechanobiologically equilibrated objective function
    //Unknown input variables
    const double a_e_guess = gsl_vector_get(x, 0);
    const double h_e_guess = gsl_vector_get(x, 1);
    const double rho_c_e_guess = gsl_vector_get(x, 2);
    const double f_z_e_guess = gsl_vector_get(x, 3);

    //Equations for J1
    //Stress values from equilibrium equations
    double sigma_e_th_lmb = ((struct vessel*) curr_vessel)->P * a_e_guess / h_e_guess;
    double sigma_e_z_lmb = f_z_e_guess / (M_PI * h_e_guess * (2 * a_e_guess + h_e_guess));

    //WSS from Pousielle flow, constant viscosity
    double bar_tauw_e = ((struct vessel*) curr_vessel)->Q / pow(a_e_guess, 3);

    //Ratio of stress:WSS mediated matrix production
    double eta_K = ((struct vessel*) curr_vessel)->K_sigma_p_alpha_h[1] /
        ((struct vessel*) curr_vessel)->K_tauw_p_alpha_h[1];

    //Stress and WSS deviations from HS state
    double delta_sigma = (sigma_e_th_lmb + sigma_e_z_lmb) /
        (((struct vessel*) curr_vessel)->sigma_h[1] + ((struct vessel*) curr_vessel)->sigma_h[2]) - 1;
    double delta_tauw = bar_tauw_e / ((struct vessel*) curr_vessel)->bar_tauw_h - 1;

    //Equations for J2
    //Homeostatic geometry
    double a_h = ((struct vessel*) curr_vessel)->a_h;
    double h_h = ((struct vessel*) curr_vessel)->h_h;

    //Equilibrated stretches
    double lambda_r_e = h_e_guess / h_h;
    double lambda_th_e = (a_e_guess + h_e_guess / 2) / (a_h + h_h / 2);
    double lambda_z_e = ((struct vessel*) curr_vessel)->lambda_z_curr;
    double F_e[3] = { lambda_r_e, lambda_th_e, lambda_z_e };

    //Equilibrated volume change
    double J_e = lambda_r_e * lambda_th_e * lambda_z_e; 

    //Equilibrated mass densities
    //Equilibrated elastin density
    double rho_el_e = ((struct vessel*) curr_vessel)->rhoR_alpha_h[0] / J_e; 

    //Ratio of degrdation rate for smc:col
    double eta_q = ((struct vessel*) curr_vessel)->k_alpha_h[1] / 
        ((struct vessel*) curr_vessel)->k_alpha_h[2];

    //Ratio of stress med prod for smc:col must be same for WSS
    double eta_ups = ((struct vessel*) curr_vessel)->K_sigma_p_alpha_h[1] /
        ((struct vessel*) curr_vessel)->K_sigma_p_alpha_h[2];

    //Equilibrated muscle density
    double rho_m_e = ((struct vessel*) curr_vessel)->rhoR_alpha_h[1] / J_e *
        pow(J_e * rho_c_e_guess / ((struct vessel*) curr_vessel)->rhoR_alpha_h[2], eta_q * eta_ups);

    //Array of equilibrated densities
    double rho_alpha[3] = { rho_el_e , rho_m_e , rho_c_e_guess };
    double rho_h = ((struct vessel*) curr_vessel)->rhoR_h;

    //Equations for J3 & J4
    //number of constituents
    int n_alpha = ((struct vessel*) curr_vessel)->n_alpha;

    //cauchy stress hat for each const. in each direction
    vector<double> hat_sigma_alpha_dir(3 * n_alpha, 0);
    vector<double> sigma_e_dir(3, 0);

    //equilibrated stresses
    vector<double> sigma_e(3, 0);

    //equilibrated active stress
    double C = ((struct vessel*) curr_vessel)->CB -
        ((struct vessel*) curr_vessel)->CS * delta_tauw;
    double lambda_act = 1.0;
    double parab_act = 1 - pow((((struct vessel*) curr_vessel)->lambda_m - lambda_act) /
        (((struct vessel*) curr_vessel)->lambda_m - ((struct vessel*) curr_vessel)->lambda_0), 2);
    double hat_sigma_act_e = ((struct vessel*) curr_vessel)->T_act * (1 - exp(-pow(C, 2))) * lambda_act * parab_act;

    //equilibrated constituent strech
    double lambda_alpha_ntau_s = 0.0;

    //equilbirated 2nd PK stress hat
    double hat_S_alpha = 0.0;

    for (int alpha = 0; alpha < n_alpha; alpha++) {
        for (int dir = 0; dir < 3; dir++) {

            //Check if vessel is anisotropic
            if (((struct vessel*) curr_vessel)->eta_alpha_h[alpha] > 0) {
                //Constituent stretch is equal to deposition stretch
                lambda_alpha_ntau_s = ((struct vessel*) curr_vessel)->g_alpha_h[alpha];

                //2nd PK stress hat at equilibrium
                hat_S_alpha = ((struct vessel*) curr_vessel)->c_alpha_h[2 * alpha] * (pow(lambda_alpha_ntau_s, 2) - 1) *
                    exp(((struct vessel*) curr_vessel)->c_alpha_h[2 * alpha + 1] * pow(pow(lambda_alpha_ntau_s, 2) - 1, 2));

                //Cauchy stress hat at equilibrium
                hat_sigma_alpha_dir[alpha * dir + dir] = ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + dir] *
                    hat_S_alpha * ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + dir];


            }
            else {
                //2nd PK stress hat at equilibrium
                hat_S_alpha = ((struct vessel*) curr_vessel)->c_alpha_h[2 * alpha];

                //Cauchy stress hat at equilibrium
                hat_sigma_alpha_dir[alpha * dir + dir] = ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + dir] * hat_S_alpha *
                    ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + dir];

                //Check if the consituent is present from the initial time point
                //Account for volume change and mixture deformation
                if (((struct vessel*) curr_vessel)->k_alpha_h[alpha] == 0) {
                    hat_sigma_alpha_dir[alpha * dir + dir] = F_e[dir] * hat_sigma_alpha_dir[alpha * dir + dir] * F_e[dir];
                }

            }

            sigma_e_dir[dir] += rho_alpha[alpha] / rho_h * hat_sigma_alpha_dir[alpha * dir + dir];

            if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1 && dir == 1) {
                sigma_e_dir[dir] += rho_alpha[alpha] / rho_h * hat_sigma_act_e;
            }
        }
    }

    //Four objective equations
    double J1 = eta_K*delta_sigma - delta_tauw; //mechano-mediated matrix production equation
    double J2 = rho_el_e + rho_m_e + rho_c_e_guess - rho_h; //mixture mass balance equation
    double J3 = sigma_e_dir[1] - sigma_e_dir[0] - sigma_e_th_lmb; //circumferential const - lmb
    double J4 = sigma_e_dir[2] - sigma_e_dir[0] - sigma_e_z_lmb; //axial const - lmb

    gsl_vector_set(f, 0, J1);
    gsl_vector_set(f, 1, J2);
    gsl_vector_set(f, 2, J3);
    gsl_vector_set(f, 3, J4);

    //Store equilibrated results        
    ((struct vessel*) curr_vessel)->a_e = a_e_guess;
    ((struct vessel*) curr_vessel)->h_e = h_e_guess;
    ((struct vessel*) curr_vessel)->rho_c_e = rho_c_e_guess * J_e;
    ((struct vessel*) curr_vessel)->rho_m_e = rho_m_e * J_e;
    ((struct vessel*) curr_vessel)->f_z_e = f_z_e_guess;
    ((struct vessel*) curr_vessel)->mb_equil_e = 1 + ((struct vessel*) curr_vessel)->K_sigma_p_alpha_h[2] * delta_sigma -
        ((struct vessel*) curr_vessel)->K_tauw_p_alpha_h[2] * delta_tauw;

    return GSL_SUCCESS;
}

int print_state_mr(size_t iter, gsl_multiroot_fsolver* s)
{

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

int find_iv_geom(void* curr_vessel) {
    //Finds the loaded (in vivo) conifiguration for a given pressure and axial stretch

    int status;
    int iter = 0;
    int max_iter = 100;

    int sn = ((struct vessel*) curr_vessel)->sn;

    const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    gsl_function f = { &iv_obj_f, curr_vessel };

    //Set search range for new mid radius
    double a_mid_act, a_mid_high, a_mid_low;
    a_mid_low = 0.95 * ((struct vessel*) curr_vessel)->a_mid[sn];
    a_mid_high = 1.05 * ((struct vessel*) curr_vessel)->a_mid[sn];

    gsl_root_fsolver_set( s, &f, a_mid_low, a_mid_high );

    //printf("Using %s method \n", gsl_root_fsolver_name(s));
    //printf("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err (est)");

    a_mid_act = 0.0;
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        a_mid_act = gsl_root_fsolver_root(s);
        a_mid_low = gsl_root_fsolver_x_lower(s);
        a_mid_high = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(a_mid_low, a_mid_high, 0, 0.00001);

        if (status == GSL_SUCCESS) {
            //printf("Loaded Config Converged:\n");
            //printf("%5d [%.7f, %.7f] %.7f %.7fs\n", iter, a_mid_low, a_mid_high, a_mid_act, a_mid_high - a_mid_low);
        }

        //printf("%5d [%.7f, %.7f] %.7f %.7fs\n", iter, a_mid_low, a_mid_high, a_mid_act, a_mid_high - a_mid_low);
    } while (status == GSL_CONTINUE && iter < max_iter);

    //int set_iv = iv_obj_f(a_mid_act, curr_vessel);
    //((struct vessel*) curr_vessel)->a_mid[sn] = a_mid_act;

    return status;
}

double iv_obj_f(double a_mid_guess, void *curr_vessel) {
    //Finds the difference in the theoretical stress from Laplace for deformed mixture
    //from the stress calculated from the mixture equations

    int sn = ((struct vessel*) curr_vessel)->sn;
    double a = 0.0, h = 0.0, lambda_t = 0.0, lambda_z = 0.0, J_s = 0.0;

    if (sn > 0) {
        lambda_t = a_mid_guess / ((struct vessel*) curr_vessel)->a_mid[0];
        lambda_z = ((struct vessel*) curr_vessel)->lambda_z_curr;
        J_s = ((struct vessel*) curr_vessel)->rhoR[sn] / ((struct vessel*) curr_vessel)->rhoR[0];

        //Update vessel geometry for calculation of next time step
        h = J_s / (lambda_t * lambda_z) * ((struct vessel*) curr_vessel)->h[0];
        a = a_mid_guess - h / 2;
        ((struct vessel*) curr_vessel)->a_mid[sn] = a_mid_guess;
        //((struct vessel*) curr_vessel)->a_act[sn] = ((struct vessel*) curr_vessel)->a_act[0];
        ((struct vessel*) curr_vessel)->a[sn] = a;
        ((struct vessel*) curr_vessel)->h[sn] = h;
        ((struct vessel*) curr_vessel)->lambda_th_curr = lambda_t;
        ((struct vessel*) curr_vessel)->lambda_z_curr = lambda_z;
        ((struct vessel*) curr_vessel)->bar_tauw = ((struct vessel*) curr_vessel)->Q / pow(a, 3);
        
    }
    else {
        //((struct vessel*) curr_vessel)->a_mid[0] = a_mid_guess;
        lambda_t = a_mid_guess / ((struct vessel*) curr_vessel)->a_mid_h;
        lambda_z = ((struct vessel*) curr_vessel)->lambda_z_curr;
        J_s = 1.0;

        //Update vessel geometry for calculation of next time step
        h = J_s / (lambda_t * lambda_z) * ((struct vessel*) curr_vessel)->h_h;
        a = a_mid_guess - h / 2;
        ((struct vessel*) curr_vessel)->a_mid[sn] = a_mid_guess;
        //((struct vessel*) curr_vessel)->a_act[sn] = ((struct vessel*) curr_vessel)->a_act[0];
        ((struct vessel*) curr_vessel)->a[sn] = a;
        ((struct vessel*) curr_vessel)->h[sn] = h;
        ((struct vessel*) curr_vessel)->lambda_th_curr = 1.0;
        ((struct vessel*) curr_vessel)->lambda_z_curr = lambda_z;
        ((struct vessel*) curr_vessel)->bar_tauw = ((struct vessel*) curr_vessel)->Q / pow(a, 3);
    }

    double sigma_t_th = ((struct vessel*) curr_vessel)->P * a / h;
    update_sigma(curr_vessel);

    ((struct vessel*) curr_vessel)->f = M_PI * h * (2 * a + h) *
        ((struct vessel*) curr_vessel)->sigma[2];

    double J = ((struct vessel*) curr_vessel)->sigma[1] - sigma_t_th;

    return J;
}

void update_kinetics(vessel &curr_vessel) {

    //This function updates the kinetics for G&R.
    int n_alpha = curr_vessel.n_alpha;
    int nts = curr_vessel.nts;
    double dt = curr_vessel.dt;
    double s = curr_vessel.s;
    int sn = curr_vessel.sn;
    int taun_min = 0;
    
    double tau_max = 10 * (1 / curr_vessel.k_alpha_h[1]); //max time of 10 half-lives

    //Differences in current mechanical state from the reference state
    double delta_sigma = (curr_vessel.sigma[1] + curr_vessel.sigma[2]) /
                        (curr_vessel.sigma_h[1] + curr_vessel.sigma_h[2]) - 1;

    double delta_tauw = (curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1;

    //Initialize pars for looping later
    double K_sigma_p = 0, K_tauw_p = 0, K_sigma_d = 0, K_tauw_d = 0;
    double upsilon_mech_p = 0, upsilon_mech_d = 0;

    double k_alpha_s = 0;
    double mR_alpha_s = 0;

    double mq_0 = 0, mq_1 = 0, mq_2;
    double q_0 = 0, q_1 = 0, q_2;
    double k_0 = 0, k_1 = 0, k_2;
    double rhoR_s = 0, rhoR_alpha_s = 0;

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
        deg_check = curr_vessel.k_alpha_h[alpha] > 0;

        if ( sn > 0 && deg_check) {

            //Get the gains for the current constituent
            K_sigma_p = curr_vessel.K_sigma_p_alpha_h[alpha];
            K_tauw_p = curr_vessel.K_tauw_p_alpha_h[alpha];

            K_sigma_d = curr_vessel.K_sigma_d_alpha_h[alpha];
            K_tauw_d = curr_vessel.K_tauw_d_alpha_h[alpha];

            //Update the stimulus functions for each constituent
            upsilon_mech_p = 1 + K_sigma_p * delta_sigma - K_tauw_p * delta_tauw;
            upsilon_mech_d = 1 + K_sigma_d * pow(delta_sigma, 2) + K_tauw_d * pow(delta_tauw, 2);

            //Reset these to zero for each constituent
            k_alpha_s = 0;
            mR_alpha_s = 0;
            rhoR_alpha_s = 0;

            //Update kinetic values for the current time
            k_alpha_s = curr_vessel.k_alpha_h[alpha] * upsilon_mech_d;
            mR_alpha_s = k_alpha_s * upsilon_mech_p * curr_vessel.rhoR_alpha[nts * alpha + sn];
            curr_vessel.k_alpha[nts * alpha + sn] = k_alpha_s;
            curr_vessel.mR_alpha[nts * alpha + sn] = mR_alpha_s;

            k_2 = curr_vessel.k_alpha[nts * alpha + sn];
            q_2 = 1.0;
            mq_2 = curr_vessel.mR_alpha[nts * alpha + sn] * q_2;

            //loop through and update constituent densities from previous time points
            //starting from the current time point and counting down is more efficient
            for (int taun = sn - 1; taun >= taun_min + 1; taun = taun - 2) {

               //Simpsons rule     
                k_1 = curr_vessel.k_alpha[nts * alpha + taun];
                q_1 = exp(-(k_2 + k_1) * dt / 2) * q_2;
                mq_1 = curr_vessel.mR_alpha[nts * alpha + taun] * q_1;

                k_0 = curr_vessel.k_alpha[nts * alpha + taun - 1];
                q_0 = exp(-(k_2 + 4.0 * k_1 + k_0) * dt / 3) * q_2;
                mq_0 = curr_vessel.mR_alpha[nts * alpha + taun - 1] * q_0;

                rhoR_alpha_s += (mq_2 + 4.0 * mq_1 + mq_0) * dt / 3;

                k_2 = k_0;
                q_2 = q_0;
                mq_2 = mq_0;

            }

            //At last time step, doing trapezoidal integration if even integration pts
            if (even_n) {
                k_0 = curr_vessel.k_alpha[nts * alpha + taun_min];
                q_0 = exp(-(k_2 + k_0) * dt / 2) * q_2;
                mq_0 = curr_vessel.mR_alpha[nts * alpha + taun_min] * q_0;

                rhoR_alpha_s += (mq_2 + mq_0) * dt / 2;
            }

            //Account for the cohort of material present initially
            if (taun_min == 0) {
                rhoR_alpha_s += curr_vessel.rhoR_alpha[nts * alpha + 0] * q_0;
            }

        }
        else {
            rhoR_alpha_s = curr_vessel.rhoR_alpha[nts * alpha + 0];
        }

        curr_vessel.rhoR_alpha[nts * alpha + sn] = rhoR_alpha_s;
        rhoR_s += rhoR_alpha_s;

    }
    
    curr_vessel.rhoR[sn] = rhoR_s;

}

void update_sigma(void* curr_vessel) {

    //Get current time index
    double s = ((struct vessel*)curr_vessel)->s;
    int sn = ((struct vessel*)curr_vessel)->sn;
    int nts = ((struct vessel*)curr_vessel)->nts;
    double dt = ((struct vessel*)curr_vessel)->dt;
    int taun_min = 0;

    double tau_max = 100 * (1 / ((struct vessel*)curr_vessel)->k_alpha_h[2]); //max time of 10 half-lives

    //Specify vessel geometry
    double a0 = ((struct vessel*)curr_vessel)->a[0];
    double h0 = ((struct vessel*)curr_vessel)->h[0];

    //Calculate vessel stretches
    double lambda_th_s = ((struct vessel*)curr_vessel)->lambda_th_curr;
    double lambda_z_s = ((struct vessel*)curr_vessel)->lambda_z_curr;

    //Calculate constituent specific stretches for evolving constituents at the current time
    int n_alpha = ((struct vessel*)curr_vessel)->n_alpha;
    vector<double> lambda_alpha_s(n_alpha, 0);
    double eta_alpha = 0;
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Check to see if constituent is isotropic
        eta_alpha = ((struct vessel*)curr_vessel)->eta_alpha_h[alpha];
        if (eta_alpha >= 0) {

            //Stretch is equal to the sqrt of I4
            lambda_alpha_s[alpha] = sqrt(pow(lambda_z_s * cos(eta_alpha), 2)
                + pow(lambda_th_s * sin(eta_alpha), 2));

            //Update stored current stretch if not numerical experiment
            if (((struct vessel*)curr_vessel)->num_exp_flag == 0) {
                ((struct vessel*)curr_vessel)->lambda_alpha_tau[nts * alpha + sn] = lambda_alpha_s[alpha];
            }
        }
    }

    //Find the current deformation gradient
    double J_s = ((struct vessel*)curr_vessel)->rhoR[sn] / ((struct vessel*)curr_vessel)->rhoR_h;
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
    double Cbar[3] = { 0 };
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
    double k_act = ((struct vessel*) curr_vessel)->k_act;
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

        //Trapz rule allows for fast heredity integral evaluation
        k_2 = ((struct vessel*)curr_vessel)->k_alpha[nts * alpha + sn];
        q_2 = 1.0;
        mq_2 = ((struct vessel*)curr_vessel)->mR_alpha[nts * alpha + sn];

        //Find active radius from current cohort
        if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
            a_2 = ((struct vessel*) curr_vessel)->a[sn];
            q_act_2 = 1.0;
        }

        //Kinematics
        F_tau[0] = F_s[0], F_tau[1] = F_s[1], F_tau[2] = F_s[2];
        J_tau = J_s;

        //Find stress from current cohort
        for (int dir = 0; dir < 3; dir++) {

            constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, sn, dir);
            hat_S_alpha = constitutive_return[0];
            hat_dSdC_alpha = constitutive_return[1];
            F_alpha_ntau_s = F_s[dir] / F_tau[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
            hat_sigma_2[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
            hat_Cbar_2[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;
        }

        //Boolean for whether the constituent increases ref mass density
        deg_check = ((struct vessel*)curr_vessel)->mR_alpha_h[alpha] > 0;

        //Check if during G&R or at initial time point
        if (sn > 0 && deg_check) {

            for (int taun = sn - 1; taun >= taun_min + 1; taun = taun - 2) {

                //Find the 1st intermediate deformation gradient
                a = ((struct vessel*)curr_vessel)->a[taun];
                h = ((struct vessel*)curr_vessel)->h[taun];
                lambda_th_tau = (a + h / 2) / (a0 + h0 / 2);
                lambda_z_tau = ((struct vessel*)curr_vessel)->lambda_z_tau[taun];
                J_tau = ((struct vessel*)curr_vessel)->rhoR[taun] / ((struct vessel*)curr_vessel)->rhoR_h;
                F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau);
                F_tau[1] = lambda_th_tau;
                F_tau[2] = lambda_z_tau;

                //Find 1st intermediate kinetics
                k_1 = ((struct vessel*)curr_vessel)->k_alpha[nts * alpha + taun];
                q_1 = exp(-(k_2 + k_1) * dt / 2) * q_2;
                mq_1 = ((struct vessel*)curr_vessel)->mR_alpha[nts * alpha + taun] * q_1;

                //Find intermediate active state
                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    a_1 = a;
                    q_act_1 = exp(-k_act * dt) * q_act_2;
                }

                //Find 1st intermeidate stress component in each direction
                for (int dir = 0; dir < 3; dir++) {

                    constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, taun, dir);
                    hat_S_alpha = constitutive_return[0];
                    hat_dSdC_alpha = constitutive_return[1];
                    F_alpha_ntau_s = F_s[dir] / F_tau[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
                    hat_sigma_1[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                    hat_Cbar_1[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;
                }

                //Find the 2nd intermediate deformation gradient
                a = ((struct vessel*)curr_vessel)->a[taun - 1];
                h = ((struct vessel*)curr_vessel)->h[taun - 1];
                lambda_th_tau = (a + h / 2) / (a0 + h0 / 2);
                lambda_z_tau = ((struct vessel*)curr_vessel)->lambda_z_tau[taun - 1];
                J_tau = ((struct vessel*)curr_vessel)->rhoR[taun - 1] / ((struct vessel*)curr_vessel)->rhoR_h;;
                F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau);
                F_tau[1] = lambda_th_tau;
                F_tau[2] = lambda_z_tau;

                //Find 2nd intermediate kinetics
                k_0 = ((struct vessel*)curr_vessel)->k_alpha[nts * alpha + taun - 1];
                q_0 = exp(-(k_2 + 4 * k_1 + k_0) * dt / 3) * q_2;
                mq_0 = ((struct vessel*)curr_vessel)->mR_alpha[nts * alpha + taun - 1] * q_0;

                //Find intermediate active state
                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    a_0 = a;
                    q_act_0 = exp(-k_act * dt) * q_act_1;
                }

                //Find component in each direction
                for (int dir = 0; dir < 3; dir++) {

                    constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, taun - 1, dir);
                    hat_S_alpha = constitutive_return[0];
                    hat_dSdC_alpha = constitutive_return[1];
                    F_alpha_ntau_s = F_s[dir] / F_tau[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
                    hat_sigma_0[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                    hat_Cbar_0[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;

                    //Add to the stress and stiffness contribution in the given direction
                    sigma[dir] += (mq_2 * hat_sigma_2[dir] + 4 * mq_1 * hat_sigma_1[dir] + mq_0 * hat_sigma_0[dir])
                        / ((struct vessel*)curr_vessel)->rhoR_h * dt / 3;
                    Cbar[dir] += (mq_2 * hat_Cbar_2[dir] + 4 * mq_1 * hat_Cbar_1[dir] + mq_0 * hat_Cbar_0[dir])
                        / ((struct vessel*)curr_vessel)->rhoR_h * dt / 3;
                }

                //Store active vars for next iteration
                //Find intermediate active state
                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    a_act += k_act * (q_act_2 * a_2 + 4 * q_act_1 * a_1 + q_act_0 * a_0) * dt / 3;
                    a_2 = a_0;
                    q_act_2 = q_act_0;
                }

                //Store intermediate kinetics for next iteration
                k_2 = k_0;
                q_2 = q_0;
                mq_2 = mq_0;

                //Store intermediate stress and stiffness for next iteration
                hat_sigma_2 = hat_sigma_0;
                hat_Cbar_2 = hat_Cbar_0;

            }

            if (even_n) {

                //Find the 2nd intermediate deformation gradient
                a = ((struct vessel*)curr_vessel)->a[taun_min];
                h = ((struct vessel*)curr_vessel)->h[taun_min];
                lambda_th_tau = (a + h / 2) / (a0 + h0 / 2);
                lambda_z_tau = ((struct vessel*)curr_vessel)->lambda_z_tau[taun_min];
                J_tau = ((struct vessel*)curr_vessel)->rhoR[taun_min] / ((struct vessel*)curr_vessel)->rhoR_h;
                F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau);
                F_tau[1] = lambda_th_tau;
                F_tau[2] = lambda_z_tau;

                //Find 2nd intermediate kinetics
                k_0 = ((struct vessel*)curr_vessel)->k_alpha[nts * alpha + taun_min];
                q_0 = exp(-(k_2 + k_0) * dt / 2) * q_2;
                mq_0 = ((struct vessel*)curr_vessel)->mR_alpha[nts * alpha + taun_min] * q_0;

                //Find intermediate active state
                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    a_0 = a;
                    q_act_0 = exp(-k_act * dt) * q_act_2;
                }

                //Find component in each direction
                for (int dir = 0; dir < 3; dir++) {

                    constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, taun_min, dir);
                    hat_S_alpha = constitutive_return[0];
                    hat_dSdC_alpha = constitutive_return[1];
                    F_alpha_ntau_s = F_s[dir] / F_tau[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
                    hat_sigma_0[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                    hat_Cbar_0[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;

                    //Add to the stress and stiffness contribution in the given direction
                    sigma[dir] += (mq_2 * hat_sigma_2[dir] + mq_0 * hat_sigma_0[dir])
                        / ((struct vessel*)curr_vessel)->rhoR_h * dt / 2;
                    Cbar[dir] += (mq_2 * hat_Cbar_2[dir] + mq_0 * hat_Cbar_0[dir])
                        / ((struct vessel*)curr_vessel)->rhoR_h * dt / 2;
                }

                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    a_act += k_act * (q_act_2 * a_2 + q_act_0 * a_0) * dt / 2;
                }
            }

            //Add in the stress and stiffness contributions of the initial material
            if (taun_min == 0) {
                for (int dir = 0; dir < 3; dir++) {
                    sigma[dir] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + 0]
                        / ((struct vessel*)curr_vessel)->rhoR_h * q_0 * hat_sigma_0[dir];
                    Cbar[dir] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + 0]
                        / ((struct vessel*)curr_vessel)->rhoR_h * q_0 * hat_Cbar_0[dir];
                }
            }

        }
        //Initial time point and constituents with prescribed degradation profiles
        else {
            //Find stress from initial cohort          
            for (int dir = 0; dir < 3; dir++) {

                constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, 0, dir);
                hat_S_alpha = constitutive_return[0];
                hat_dSdC_alpha = constitutive_return[1];
                F_alpha_ntau_s = F_s[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
                hat_sigma_2[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                hat_Cbar_2[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;

                sigma[dir] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + sn] /
                    ((struct vessel*)curr_vessel)->rhoR_h * hat_sigma_2[dir];
                Cbar[dir] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + sn] /
                    ((struct vessel*)curr_vessel)->rhoR_h * hat_Cbar_2[dir];
            }

        }

        if (taun_min == 0 && ((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
            a_act += ((struct vessel*) curr_vessel)->a_act[0] * q_act_0;
        }


    }

    //Find active stress contribtion
    //add in initial active stress radius contribution
    if (sn == 0) {
        a_act = ((struct vessel*) curr_vessel)->a_act[0];
    }

    C = ((struct vessel*) curr_vessel)->CB -
        ((struct vessel*) curr_vessel)->CS * (((struct vessel*) curr_vessel)->bar_tauw /
        ((struct vessel*) curr_vessel)->bar_tauw_h - 1);

    lambda_act = ((struct vessel*) curr_vessel)->a[sn] / a_act;

    parab_act = 1 - pow((((struct vessel*) curr_vessel)->lambda_m - lambda_act) /
        (((struct vessel*) curr_vessel)->lambda_m - ((struct vessel*) curr_vessel)->lambda_0), 2);

    hat_sigma_act = ((struct vessel*) curr_vessel)->T_act * (1 - exp(-pow(C, 2))) * lambda_act * parab_act;

    sigma_act = ((struct vessel*) curr_vessel)->rhoR_alpha[nts * 1 + sn] / J_s / ((struct vessel*)curr_vessel)->rhoR_h * hat_sigma_act;

    //The Lagrange multiplier is the radial stress component
    //subtract from each direction

    lagrange = sigma[0];
    for (int dir = 0; dir < 3; dir++) {

        Cbar[dir] = 2 * sigma[dir] + 2 * Cbar[dir];
        sigma[dir] = sigma[dir] - lagrange;
        //Accounting for active stress
        if (dir == 1) {
            sigma[dir] += sigma_act;
        }

        ((struct vessel*)curr_vessel)->sigma[dir] = sigma[dir];
        ((struct vessel*)curr_vessel)->Cbar[dir] = Cbar[dir];
    }

    //Save updated active radius
    ((struct vessel*) curr_vessel)->a_act[sn] = a_act;

}

vector<double> constitutive(void* curr_vessel, double lambda_alpha_s, int alpha, int ts, int dir) {

    double lambda_alpha_ntau_s = 0;
    double Q1 = 0;
    double Q2 = 0;
    double hat_S_alpha = 0;
    double hat_dSdC_alpha = 0;
    double pol_mod = 0;
    int nts = ((struct vessel*)curr_vessel)->nts;
    vector<double> return_constitutive = { 0, 0 };

    //Check if ansisotropic
    if (((struct vessel*)curr_vessel)->eta_alpha_h[alpha] >= 0) {

        lambda_alpha_ntau_s = ((struct vessel*)curr_vessel)->g_alpha_h[alpha] *
            lambda_alpha_s / ((struct vessel*)curr_vessel)->lambda_alpha_tau[nts * alpha + ts];

        if (lambda_alpha_ntau_s < 1) {
            lambda_alpha_ntau_s = 1;
        }

        Q1 = (pow(lambda_alpha_ntau_s, 2) - 1);
        Q2 = ((struct vessel*)curr_vessel)->c_alpha_h[2 * alpha + 1] * pow(Q1, 2);
        hat_S_alpha = ((struct vessel*)curr_vessel)->c_alpha_h[2 * alpha] * Q1 * exp(Q2);
        hat_dSdC_alpha = ((struct vessel*)curr_vessel)->c_alpha_h[2 * alpha] * exp(Q2) * (1 + 2 * Q2);

    }
    else {

        hat_S_alpha = ((struct vessel*)curr_vessel)->c_alpha_h[2 * alpha];
    }

    return_constitutive = { hat_S_alpha , hat_dSdC_alpha };

    return return_constitutive;

}

