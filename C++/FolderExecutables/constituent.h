// gnr_vessel.h
#ifndef CONSTITUENT
#define CONSTITUENT
#include "json_fwd.hpp"
using std::string;
using std::vector;
using std::cout;
using json = nlohmann::json;
class layer;
class constituent {
public:
    //Constituent inflammatory status
    int alpha_infl;
    // boolean that indicates if it's a polymer
    int is_pol; 

    //Material properties
    vector<double> c_alpha_h, G_alpha_h, g_alpha_h;
    double eta_alpha_h;
    double epsilon_pol_min;

    //Mass fractions, referential apparent densities, kinetic quantities
    double phi_alpha_h, rhoR_alpha_h, mR_alpha_h, k_alpha_h;
    double K_sigma_p_alpha_h, K_sigma_d_alpha_h, K_tauw_p_alpha_h, K_tauw_d_alpha_h;
    vector<double> K_infl_p_alpha;

    //True mass densities, initial volume fractions
    double rho_hat_alpha_h, epsilonR_alpha_0;
    
    //Histories
    vector<double> rhoR_alpha, mR_alpha, k_alpha;
    vector<double> epsilonR_alpha, epsilon_alpha;
    vector<double> degradation_profile; // for materials that do not grow (e.g. elastin, degradable sheath), specifies the proportion of remaining material at time t

    //Current loading quantities
    vector<double> lambda_alpha_tau;
    
    // Active Params
    vector<double> lambda_act; //active radius history

    //Active stress quantities
    int alpha_active; //boolean for active constituents
    vector<double> a_act; //active radius history
    double T_act, T_act_h; //Max active stress, homeostatic max active stress
    double k_act; //Active remodeling parameters
    double lambda_0; //Min contractile stretch
    double lambda_m; //Max contractile stretch
    double CB; //Basal VC to VD ratio
    double CS; //Scaling factor for VC to VD ratio

    //Mechanobiologically equilibrated quantities
    double rho_alpha_e; //equilibrated constituent density

    // Maintain reference to parent layer   
    layer *parent_layer;

    // Temp state vars
    double ups_saved;
 
    constituent();
    void initializeJSON(json& json_in, layer* parent_l); 
};



#endif /* GNR_VESSEL */
