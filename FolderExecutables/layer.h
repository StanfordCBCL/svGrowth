// gnr_vessel.h
#ifndef LAYER
#define LAYER

#include "json_fwd.hpp"
using std::string;
using std::vector;
using std::cout;
using json = nlohmann::json;
class constituent;
class vessel;
class layer {
public:

    //Geometric quantities
    double a_h, h_h; //In-vivo reference
    double a_mid_h; //In vivo midpoint

    vector<double> a, a_mid, h; //Evolving in vivo reference
    vector<double> a_pas, a_mid_pas, h_pas; //Evolving in vivo reference
    vector<double> A, A_mid, H, lambda_z_pre; //Evolving traction free

    //Reference true mass density
    double rhoR_h;

    //Histories
    vector<double> rhoR, rho; 
    vector<double> ups_infl_p, ups_infl_d;

    //Reference loading quantities
    vector<double> sigma_h;

    //Current loading quantities
    double lambda_th_curr;
    vector<double> sigma, Cbar, CC;

    double sigma_inv, sigma_inv_h;

    vector<double> F_s;
    double J_di;
    vector<double> F_curr;
    vector<double> Fi_curr;
    
    //Mechanobiologically equilibrated quantities
    double a_e; //equilibrated radius
    double h_e; //equilibrated thickness

    //Mass density parameters
    double rho_hat_h;

    bool is_contacting; // True for outer layers if contacting inner vessels, always true for inner layer. 

    // Maintain reference to parent vessel   
    vessel *parent_vessel;

    // Create array of constituents
    vector<constituent> constituents;

;

    
    layer();
    void initializeJSON(json& json_in, vessel *parent); 

};



#endif /* GNR_VESSEL */
