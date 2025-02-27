// gnr_vessel.h
#ifndef VESSEL
#define VESSEL

#include "json_fwd.hpp"
using std::string;
using std::vector;
using std::cout;
using json = nlohmann::json;
class constituent;
class layer;
class vessel {
public:
    string vessel_name;
    string file_name;
    string gnr_name;
    string equil_gnr_name;
    string exp_name;
    string hns_name;

    //Time variables
    int nts; //total number of time steps
    double dt; //time increment
    int sn; //current time step index
    double s; //actual current time

    //Geometric quantities
    double lambda_z_h; //Reference in vivo stretch

    //Reference loading quantities
    double P_h, f_h, bar_tauw_h, Q_h, P_prev;

    //Current loading quantities
    vector<double> P, Q, lambda_z_tau;
    double f, bar_tauw;
    double mb_equil;

    //Mechanobiologically equilibrated quantities
    double f_z_e; //equilibrated axial force
    double mb_equil_e; //Current mechanobiological equil. state

    //Flags
    int num_exp_flag; //indicates whether doing reg G&R step or a numerical experiment
    int pol_only_flag; //indicates whether other constituents are produced
    int wss_calc_flag; //indicates if GnR should update its own current WSS
    int app_visc_flag; //indicates whether to use the empirical correction for viscosity from Secomb 2017

    string variable_defs;
    string constant_defs; 


    //Conversions
    double um_to_m = pow(10, -6);
    double mm_to_m = pow(10, -3);
    double kPa_to_Pa = pow(10, 3);
        
    // Create vector of layers
    vector<layer> layers;

    std::ofstream GnR_out, Equil_GnR_out, Exp_out, Mat_out, HnS_out;

    vessel(); //Default constructor
    //Vessel(string file_name); //File name constructor ***ELS USE DELEGATING CONSTRUCTOR***
    //Vessel(string file_name, string vessel_type); //File name and type constructor ***ELS USE DELEGATING CONSTRUCTOR***
    void printNativeOutputs();
    void initializeJSON(string json_name, double n_days_inp = 10, double dt_inp = 1);
    void evaluate_expr(json& expression_in, vector<double>& result, double dt, int nts);
    //~Vessel() destructor

};



#endif /* GNR_VESSEL */
