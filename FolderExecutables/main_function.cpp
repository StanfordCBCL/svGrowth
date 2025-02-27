//Models the G&R of a TEVG
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

#include "vessel.h"
#include "functions.h"

using std::string;
using std::vector;
using std::cout;

//All units meters, kg, days
//Conversions
double um_to_m = pow(10, -6);
double mm_to_m = pow(10, -3);
double kPa_to_Pa = pow(10, 3);
#ifdef __cplusplus
extern "C"
#endif
int run(char* prefix_char, char* name_char, int restart_arg, int iter_arg, int gnr_arg, int num_days, double step_size,
        double sigma_arg, double tauw_arg, double anysm_arg, double tevg_arg, double F0, double F1, double F2, double F3, double F4,
        double F5, double F6, double F7, double F8) {

    //std::cout << string(prefix_char) << std::endl;
    //std::cout << string(name_char) << std::endl;
    //std::cout << restart_arg << std::endl;
    //std::cout << iter_arg << std::endl;
    //std::cout << gnr_arg << std::endl;
    //std::cout << num_days << std::endl;
    //std::cout << step_size << std::endl;
    //std::cout << sigma_arg << std::endl;

    string prefix_arg = string(prefix_char);
    string name_arg = string(name_char);
    vector<double> F_s = {F0,F1,F2,F3,F4,F5,F6,F7,F8};

    try {
        //std::cout << "Restarting simulation: " << restart_arg << "\n";
        //std::cout << "Filename suffix: " << name_arg << "\n";

        //Initialize the reference vessel for the simulation
        string native_file = "FolderVesselConfigurationFiles/Native_in_handshake_";// + name_arg;
        string immune_file = "FolderVesselConfigurationFiles/Immune_in_";// + name_arg;
        string scaffold_file = "FolderVesselConfigurationFiles/Scaffold_in_";// + name_arg;

        vessel simulation_vessel;

        //Initialize TEVG files if necessary
        if(tevg_arg > 0.0) {
            //Initialize TEVG
            simulation_vessel.initializeNativeHandshake(native_file,num_days,step_size);
            simulation_vessel.initializeTEVGHandshake(scaffold_file,immune_file,tevg_arg,num_days,step_size);
        } else {
            simulation_vessel.initializeNativeHandshake(native_file,num_days,step_size);
        }

        //Get all other input arguements and apply to TEVG
        simulation_vessel.gnr_name = prefix_arg + "/" + simulation_vessel.gnr_name + "_" + name_arg;
        simulation_vessel.exp_name = prefix_arg + "/" + simulation_vessel.exp_name + "_" + name_arg;
        simulation_vessel.file_name = prefix_arg + "/" + simulation_vessel.file_name + "_" + name_arg;
        simulation_vessel.hns_name = prefix_arg + "/" + simulation_vessel.hns_name + "_" + name_arg;

        //------------------------------------------------------------------------

        //For elastin degradation 
        double Q_p1;//, s;

        //double s_gnd_off = 30, k_gnd_h = 0.1;

        if(tevg_arg <= 0.0) {
            for (int sn = 1; sn < simulation_vessel.nts; sn++) {
                //Calculate elastin degradation
                //s = sn * simulation_vessel.dt;
                
                Q_p1 = 1.0;
                
                simulation_vessel.epsilonR_alpha[0 * simulation_vessel.nts + sn] = Q_p1 * simulation_vessel.epsilonR_alpha[0 * simulation_vessel.nts + 0];

                simulation_vessel.rhoR_alpha[0 * simulation_vessel.nts + sn] = simulation_vessel.epsilonR_alpha[0 * simulation_vessel.nts + sn] * 
                                                                    simulation_vessel.rho_hat_alpha_h[0];
            }
        }
        if(anysm_arg > 0.0) {
            //Change endothelial functioning to be proportional to damage
            for (int alpha = 0; alpha < simulation_vessel.n_alpha; alpha++) {
                simulation_vessel.K_tauw_p_alpha_h[alpha] = std::max(simulation_vessel.K_tauw_p_alpha_h[alpha] - anysm_arg, 0.0);
            }
            
            simulation_vessel.c_alpha_h[0] = simulation_vessel.c_alpha_h[0]*(1.0 - anysm_arg);
            
            /*

            for (int sn = 1; sn < simulation_vessel.nts; sn++) {
                //Calculate elastin degradation
                s = sn * simulation_vessel.dt;

                Q_p1 = (s > s_gnd_off) * (anysm_arg * exp(-k_gnd_h * (s - s_gnd_off)) + (1 - anysm_arg))
                + (s <= s_gnd_off) * 1;
                
                simulation_vessel.epsilonR_alpha[0 * simulation_vessel.nts + sn] = Q_p1 * simulation_vessel.epsilonR_alpha[0 * simulation_vessel.nts + 0];

                simulation_vessel.rhoR_alpha[0 * simulation_vessel.nts + sn] = simulation_vessel.epsilonR_alpha[0 * simulation_vessel.nts + sn] * 
                                                                    simulation_vessel.rho_hat_alpha_h[0];

                //std::cout << sn << ": " << simulation_vessel.epsilonR_alpha[0 * simulation_vessel.nts + sn] << std::endl;
            }
            */
        }

        //If restarting file, load previous vessel
        if(!restart_arg && !iter_arg) {
            //simulation_vessel.GnR_out.open(simulation_vessel.gnr_name);
            simulation_vessel.HnS_out.open(simulation_vessel.hns_name);
        } else {
            //simulation_vessel.GnR_out.open(simulation_vessel.gnr_name, std::ofstream::out | std::ofstream::app);
            simulation_vessel.HnS_out.open(simulation_vessel.hns_name, std::ofstream::out | std::ofstream::app);
            //Read vessel from file
            simulation_vessel.load();
            //std::cout << "Loading vessel" <<std::endl;
        }

        if (simulation_vessel.sn == 0){
            //Write initial state to file
            update_sigma_handshake(&simulation_vessel);
            update_kinetics(simulation_vessel);
            //printf("%s \n", "---------------------------");
            //fflush(stdout);
            simulation_vessel.printOutputsHandshake();
        }

        simulation_vessel.sigma_inv = sigma_arg;
        simulation_vessel.bar_tauw = tauw_arg;
        for (int i = 0; i < 9; i++) {
            simulation_vessel.F_curr[i] = F_s[i];
        }
            
        if(gnr_arg){
            update_time_step_handshake(simulation_vessel, iter_arg);
            //printf("%s \n", "---------------------------");
            //fflush(stdout);
        }
        //else{
            //Print current state
            //printf("%s %f", "Time:", simulation_vessel.s);
            //printf("%s \n", "xxxxxxxxxxxxxxxxxxxxxxxxxxx");
            //fflush(stdout); 
        //}

        //Write full model outputs
        simulation_vessel.printOutputsHandshake();
        //Print vessel to file
        simulation_vessel.save();

        //simulation_vessel.GnR_out.close();
        simulation_vessel.HnS_out.close();

    } catch(std::exception& e) {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;

}

int main(){}
