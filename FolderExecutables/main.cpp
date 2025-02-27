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
#include "layer.h"
#include "json.hpp"
#include "constituent.h"
#include "functions.h"

using std::string;
using std::vector;
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//All units meters, kg, days
//Conversions
double um_to_m = pow(10, -6);
double mm_to_m = pow(10, -3);
double kPa_to_Pa = pow(10, 3);

int main( int ac, char* av[] ) {

    try{

        //Input argument variables to store to
        //Set name, number steps, and new pressure and wss
        int restart_arg;
        int step_arg;
        double step_size;
        double P_arg;
        double tauw_arg;
        int gnr_arg;
        string name_arg;
        int opt;
        int portnum;
        int num_days;

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("restart,r", po::value<int>(&restart_arg)->default_value(0),
                "restarting simulation")
            ("step,s", po::value<int>(&step_arg), "number of timesteps to run simulation")
            ("time step size,d", po::value<double>(&step_size)->default_value(1.0), "size of each time step in days")
            ("name,n", po::value<string>(&name_arg)->default_value(""),
                "suffix of simulation files.")
            ("pressure,p", po::value<double>(&P_arg), "input pressure")
            ("wss,w", po::value<double>(&tauw_arg), "input wall shear stress")
            ("simulate,u", po::value<int>(&gnr_arg)->default_value(1), "execute simulation")
            ("max_days,m", po::value<int>(&num_days)->default_value(361), "maximum days to simulate")
        ;

        po::positional_options_description p;
        p.add("name", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).
                  options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << "Usage: options_description [options]\n";
            cout << desc;
            return 0;
        }

        std::cout << "Restarting simulation: " << restart_arg << "\n";
        std::cout << "Filename suffix: " << name_arg << "\n";

        //Initialize the reference vessel for the simulation
        string native_file = "../FolderVesselConfigurationFiles/JSON_in_" + name_arg + ".json";
        std::cout << "Filename native: " << native_file << "\n";

        vessel native_vessel;
        //native_vessel.initializeNative(native_file,num_days,step_size);
        native_vessel.initializeJSON(native_file,num_days,step_size);

        if (!vm.count("step")) {
            step_arg = int( num_days / step_size );
        }
        std::cout << "Steps to simulate: " << step_arg << "\n";
        // if (vm.count("pressure")){
        //     native_vessel.P = P_arg;
        //     std::cout << "Updated pressure: " << P_arg << std::endl;
        // }
        // //Apparent viscosity vars
        // double mu = 0.0;
        // if (vm.count("wss")){
        //     native_vessel.bar_tauw = tauw_arg;
        //     mu = get_app_visc(&native_vessel, 0);
        //     native_vessel.Q = native_vessel.bar_tauw / (4 * mu / (3.14159265 * pow(native_vessel.layers[0].a_h*100, 3)));
        //     std::cout << "Updated WSS: " << tauw_arg << std::endl;
        // }

        //Setup all other output files
        native_vessel.gnr_name = native_vessel.gnr_name + "_" + name_arg;
        native_vessel.equil_gnr_name = native_vessel.equil_gnr_name + "_" + name_arg;
        native_vessel.exp_name = native_vessel.exp_name + "_" + name_arg;
        native_vessel.file_name = native_vessel.file_name + "_" + name_arg;

        //------------------------------------------------------------------------

        //For elastin degradation 
        
        double Q_p1, s;
        int L_Z = 2;
        double s_gnd_off = 30, epsilonR_gnd_min = 0.00, k_gnd_h = 0.1;
        for (int layer = 0; layer < native_vessel.layers.size(); layer ++) {
            for (int alpha = 0; alpha < native_vessel.layers[layer].constituents.size(); alpha++ ) {
                if (native_vessel.layers[layer].constituents[alpha].k_alpha_h <= 0) {
                    for (int sn = 1; sn < native_vessel.nts; sn++) {
                        //Calculate elastin degradation
                        s = sn * native_vessel.dt;
                        Q_p1 = native_vessel.layers[layer].constituents[alpha].degradation_profile[sn]; // (1 + exp(-k_p1 * zeta_p1)) / (1 + exp(k_p1 * gamma_p_d1 * (s - zeta_p1 * gamma_p_d2)));
                        
                        //Uncomment for elastin damage by section
                        //if(L_Z == 2) {
                    //  Q_p1 = (s > s_gnd_off) * std::max(((1 - epsilonR_gnd_min) * exp(-k_gnd_h * (s - s_gnd_off)) + epsilonR_gnd_min),0.2)
                        //      + (s <= s_gnd_off) * 1;
                    // }

                        //Q_p1 = (s > s_gnd_off) * ((1 - epsilonR_gnd_min) * exp(-k_gnd_h * (s - s_gnd_off)) + epsilonR_gnd_min)
                        //+ (s <= s_gnd_off) * 1;
                        
                        native_vessel.GnR_out.flush();

                        native_vessel.layers[layer].constituents[alpha].epsilonR_alpha[sn] = Q_p1 * native_vessel.layers[layer].constituents[alpha].epsilonR_alpha[0];
                        //native_vessel.epsilonR_alpha[0 * native_vessel.nts + sn] = Q_p1 * native_vessel.epsilonR_alpha[0 * native_vessel.nts + 0];

                        native_vessel.layers[layer].constituents[alpha].rhoR_alpha[sn] = native_vessel.layers[layer].constituents[alpha].epsilonR_alpha[sn] * 
                                                                            native_vessel.layers[layer].constituents[alpha].rho_hat_alpha_h;

                    }
                }
            }
        }
        

         


        //PD Tests
        //double P_low = 0.0 * 133.33;
        //double P_high = 10.0 * 133.33;
        //double lambda_z_test = 1.0;
        //int pd_test_check = 0.0;
        //int pd_test_count = 0;
        //vector<int> pd_test_ind = { 1, 11, 91, 151 }; //, 91, 151

//                            printf("%s %f %s %f %s %f %s %f %s %f %s %f\n", "Time:", native_vessel.s, "a: ", native_vessel.a[native_vessel.sn], "h:", native_vessel.h[native_vessel.sn],
//                   "Cbar_r:", native_vessel.Cbar[0], "Cbar_t:", native_vessel.Cbar[1], "Cbar_z:", native_vessel.Cbar[2]);

        printf("%s %f\n", "Time:", native_vessel.s);
        for (int layer = 0; layer < native_vessel.layers.size(); layer++) {
        printf("%s %i %s %f %s %f %s %f %s %f %s %f","Layer #:", layer, "a: ", native_vessel.layers[layer].a[native_vessel.sn], "h:", native_vessel.layers[layer].h[native_vessel.sn], 
        "Cbar_r:", native_vessel.layers[layer].Cbar[0], "Cbar_t:", native_vessel.layers[layer].Cbar[1], "Cbar_z:", native_vessel.layers[layer].Cbar[2]);
        }
        printf("%s \n", "---------------------------");
        fflush(stdout);

        if(!restart_arg)
        {
            std::cout << "Initializing new simulation..." << std::endl;
            //Setup file I/O for G&R output
            native_vessel.GnR_out.open(native_vessel.gnr_name);
            native_vessel.Equil_GnR_out.open(native_vessel.equil_gnr_name);
            native_vessel.Exp_out.open(native_vessel.exp_name);

            //Write initial state to file
            int sn = 0;
//            native_vessel.printNativeOutputs();
            //Simulate number of timesteps
            //native_vessel.nts = step_arg; 
            //Set flag for WSS calculation
            if (!(vm.count("wss"))){
                native_vessel.wss_calc_flag = 1;
            }

            //Run the G&R time stepping
            for (int sn = 1; sn < std::min(step_arg,native_vessel.nts); sn++) {
                
                if(gnr_arg){

                    native_vessel.s = native_vessel.dt * sn;
                    native_vessel.sn = sn;
                    update_time_step(native_vessel);
                    // P/P0, a/a0, tv / tv0, 
                    // ts / ts0, sigmatt_v, sigmatt_s, sigmazz_v, sigmazz_s
                    // printf("%f, %f, %f, %f, %f, %f, %f, %f;\n", 
                    //     native_vessel.P[sn] / native_vessel.P[0], native_vessel.layers[0].a[sn] / native_vessel.layers[0].a[0], native_vessel.layers[0].h[sn] / native_vessel.layers[0].h[0], 
                    //     native_vessel.layers[1].h[sn] / native_vessel.layers[1].h[0], native_vessel.layers[0].sigma[1] / 1000, native_vessel.layers[1].sigma[1] / 1000, native_vessel.layers[0].sigma[2] / 1000, native_vessel.layers[1].sigma[2] / 1000);

                    printf("%s %f\n", "Time:", native_vessel.s);
                    for (int layer = 0; layer < native_vessel.layers.size(); layer++) {
                        printf("%s %i %s %f %s %f %s %f %s %f %s %f\n","Layer #:", layer, "a: ", native_vessel.layers[layer].a[native_vessel.sn], "h:", native_vessel.layers[layer].h[native_vessel.sn], 
                            "Cbar_r:", native_vessel.layers[layer].Cbar[0], "Cbar_t:", native_vessel.layers[layer].Cbar[1], "Cbar_z:", native_vessel.layers[layer].Cbar[2]);
                    }
                    printf("%s \n", "---------------------------");
                    native_vessel.printNativeOutputs();
                    fflush(stdout);

                }
                else{
                    //Print current state
                    printf("%s %f\n", "Time:", native_vessel.s);
                    for (int layer = 0; layer < native_vessel.layers.size(); layer++) {
                        printf("%s %i %s %f %s %f %s %f %s %f %s %f\n","Layer #:", layer, "a: ", native_vessel.layers[layer].a[native_vessel.sn], "h:", native_vessel.layers[layer].h[native_vessel.sn], 
                            "Cbar_r:", native_vessel.layers[layer].Cbar[0], "Cbar_t:", native_vessel.layers[layer].Cbar[1], "Cbar_z:", native_vessel.layers[layer].Cbar[2]);
                    }
                    printf("%s \n", "xxxxxxxxxxxxxxxxxxxxxxxxxxx");
                    fflush(stdout); 
                }

                //To check to run a pressure-diameter test
                //if (pd_test_count < pd_test_ind.size()) {
                    //if (sn == pd_test_ind[pd_test_count]) {
                        //pd_test_check = run_pd_test(TEVG, P_low, P_high, lambda_z_test);
                        //pd_test_count++;
                    //}
                //}

                //Write full model outputs
//                native_vessel.printNativeOutputs();
            }

            //Print vessel to file
  //          native_vessel.save();

            //Long-term equilibrated solution
		    //int sn = native_vessel.nts - 1;
            //double s = native_vessel.dt * sn;
		    //native_vessel.s = s;
		    //native_vessel.sn = sn;
				
            //int equil_solve = find_equil_geom(&native_vessel);
            //native_vessel.printNativeEquilibratedOutputs();
            //Print equilibrated solution
            //printf("%s %f %s %f %s %f %s %f\n", "Time:", s, "a_e: ", native_vessel.a_e, "h_e:", native_vessel.h_e, "mb_equil:", native_vessel.mb_equil_e);
            //fflush(stdout);

        }
        else
        {
            std::cout << "Continuing simulation from file..." << std::endl;
            //Setup file I/O for G&R output
            native_vessel.GnR_out.open(native_vessel.gnr_name, std::ofstream::out | std::ofstream::app);
            native_vessel.Exp_out.open(native_vessel.exp_name, std::ofstream::out | std::ofstream::app);

            //Read vessel from file
//            native_vessel.load();
            //Load pressure and WSS
            // if (vm.count("pressure")){
            //     native_vessel.P = P_arg;
            // }
            // if (vm.count("wss")){
            //     native_vessel.bar_tauw = tauw_arg;
            //     mu = get_app_visc(&native_vessel, native_vessel.sn);
            //     native_vessel.Q = native_vessel.bar_tauw / (4 * mu / (3.14159265 * pow(native_vessel.layers[0].a[native_vessel.sn]*100, 3)));
            // }

            if (native_vessel.sn == 0){
                //Write initial state to file
                int sn = 0;
//                native_vessel.printNativeOutputs();
            }

            if (!(vm.count("wss"))){
                native_vessel.wss_calc_flag = 1;
            }

            int csn = native_vessel.sn+1;

            //Run the G&R time stepping
            for (int sn = csn; sn < std::min(csn+step_arg,native_vessel.nts); sn++) {
                
                if(gnr_arg){
                    native_vessel.s = native_vessel.dt * sn;
                    native_vessel.sn = sn;
                    update_time_step(native_vessel);
                    printf("%s \n", "---------------------------");
                    fflush(stdout);
                }
                else{

                }
                
                //if (pd_test_count < pd_test_ind.size()) {
                    //if (sn == pd_test_ind[pd_test_count]) {
                        //pd_test_check = run_pd_test(TEVG, P_low, P_high, lambda_z_test);
                        //pd_test_count++;
                    //}
                //}

                //Write full model outputs
//                native_vessel.printNativeOutputs();

            }

            //Print vessel to file
//            native_vessel.save();

            //Long-term equilibrated solution
		    //int sn = native_vessel.nts - 1;
            //double s = native_vessel.dt * sn;
		    //native_vessel.s = s;
		    //native_vessel.sn = sn;
				
            //int equil_solve = find_equil_geom(&native_vessel);
            //native_vessel.printNativeEquilibratedOutputs();
            //Print equilibrated solution
            //printf("%s %f %s %f %s %f %s %f\n", "Time:", s, "a_e: ", native_vessel.a_e, "h_e:", native_vessel.h_e, "mb_equil:", native_vessel.mb_equil_e);
            //fflush(stdout);

        }

        native_vessel.GnR_out.close();
        native_vessel.Equil_GnR_out.close();
        native_vessel.Exp_out.close();

    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;

}
