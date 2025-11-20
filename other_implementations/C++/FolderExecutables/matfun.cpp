#include <cmath>
#include <vector>
#include <iostream>

#include "matfun.h"

using std::vector;
using std::cout;

template <typename T>
// C++ template to print vessel contents
std::ostream& operator<<(std::ostream& os, const vector<T>& v) {
    for (int i = 0; i < v.size(); ++i) { 
        os << v[i]; 
        if (i != v.size() - 1) 
            os << " "; 
    }
    return os; 
}

double mat_det(vector<double> F){
    double J = F[0]*(F[4]*F[8] - F[5]*F[7]) + F[1]*(F[5]*F[6] -
                    F[3]*F[8]) + F[2]*(F[3]*F[7] - F[4]*F[6]);
    return J;
}

double mat_ddot(vector<double> A, vector<double> B){
    double s = 0.0;


    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            s = s + A[(i*3)+j] * B[(i*3)+j];
        }
    }

    return s;

}

vector<double> mat_inv(vector<double> A){
    vector<double> Ainv(9,0);
    double d = mat_det(A);

    Ainv[0] = (A[4]*A[8]-A[5]*A[7]) / d;
    Ainv[1] = (A[2]*A[7]-A[1]*A[8]) / d;
    Ainv[2] = (A[1]*A[5]-A[2]*A[4]) / d;

    Ainv[3] = (A[5]*A[6]-A[3]*A[8]) / d;
    Ainv[4] = (A[0]*A[8]-A[2]*A[6]) / d;
    Ainv[5] = (A[2]*A[3]-A[0]*A[5]) / d;

    Ainv[6] = (A[3]*A[7]-A[4]*A[6]) / d;
    Ainv[7] = (A[1]*A[6]-A[0]*A[7]) / d;
    Ainv[8] = (A[0]*A[4]-A[1]*A[3]) / d;

    return Ainv;
}


vector<double> mat_trans(vector<double> A){
    vector<double> Atrans(9,0);

    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            Atrans[(j*3)+i]=A[(i*3)+j];
        }
    }

    return Atrans;
}


vector<double> mat_mul(vector<double> A, vector<double> B){
    vector<double> C(9,0);

    for (int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                C[(i*3)+j]+=A[(i*3)+k]*B[(k*3)+j];
            }
        }
    }

    return C;
}

vector<double> mat_pol_dec(vector<double> A){
    //This is literally the slowest way I could do this lol
    double tol = 1E-14; //Convergence tolerance
    double tol_check = 0.0; //Convergence tolerance
    int iter = 0;

    vector<double> Q1(9,0);
    vector<double> Q2(9,0);

    Q1 = A;

    do {
        tol_check = 0.0;
        iter++;

        for(int i = 0; i < 9; i++){

            Q2[i] = 0.5*Q1[i] + 0.5*mat_trans(mat_inv(Q1))[i];

            tol_check = tol_check + abs(Q1[i]-Q2[i]);
        }

        Q1 = Q2;

    } while ( iter < 1000);

    return Q2;
}

vector<double> ten_dyadprod(vector<double> A, vector<double> B){
    vector<double> C(27,0);

    for (int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                for(int l = 0; l < 3; l++){
                    C[(((i*3)+j)*3+k)*3+l] = A[(i*3)+j]*B[(k*3)+l];
                }
            }
        }
    }

    return C;
}

vector<double> ten_symmprod(vector<double> A, vector<double> B){
    vector<double> C(27,0);

    for (int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                for(int l = 0; l < 3; l++){
                    C[(((i*3)+j)*3+k)*3+l] = 0.5*(A[(i*3)+k]*B[(j*3)+l] + A[(i*3)+l]*B[(j*3)+k]);
                }
            }
        }
    }

    return C;
}

vector<double> scl_mul(vector<double> A, double b){
    vector<double> C(9,0);

    for (int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            A[(i*3)+j]=A[(i*3)+j]*b;
        }
    }

    return A;
}


vector<double> elem_mul(vector<double> A, vector<double> B){
    vector<double> C(9,0);

    for (int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
             C[(i*3)+j]=A[(i*3)+j]*B[(i*3)+j];
        }
    }

    return C;
}

void get_material_stiffness(vector<double> F_alpha_ntau_s, double hat_dS_dlambda2_alpha, double J_s, double hat_CC[3][3][3][3]){

    double hat_CC_temp[3][3][3][3] = {0};
    double hat_dS_dC_alpha_ABPQ = 0;
    double FiA = 0;
    double FjB = 0;
    double FkP = 0;
    double FlQ = 0;

    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){

                    for(int A = 0; A<3; A++){
                        for(int B = 0; B<3; B++){
                            for(int P = 0; P<3; P++){
                                for(int Q = 0; Q<3; Q++){
                                    //Current stretches (from reference configuration)

                                    FiA = F_alpha_ntau_s[i*3 + A];
                                    FjB = F_alpha_ntau_s[j*3 + B];
                                    FkP = F_alpha_ntau_s[k*3 + P];
                                    FlQ = F_alpha_ntau_s[l*3 + Q];

                                    hat_CC_temp[i][j][k][l]+= 2*FiA*FjB*FkP*FlQ*hat_dS_dlambda2_alpha/ J_s;


                                }
                            }
                        }
                    }
                    
                    hat_CC[i][j][k][l] = hat_CC_temp[i][j][k][l];

                }
            }
        }
    }

}

void get_active_material_stiffness(double S, double CC[3][3][3][3],vector<double> F_s, double J_s){

    double hat_CC_temp[3][3][3][3] = {0};
    double hat_dS_dC_alpha_ABPQ = 0;
    double FiA = 0;
    double FjB = 0;
    double FkP = 0;
    double FlQ = 0;

    hat_CC_temp[1][1][1][1] = S;

    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){

                    for(int A = 0; A<3; A++){
                        for(int B = 0; B<3; B++){
                            for(int P = 0; P<3; P++){
                                for(int Q = 0; Q<3; Q++){
                                    //Current stretches (from reference configuration)

                                    FiA = F_s[i*3 + A];
                                    FjB = F_s[j*3 + B];
                                    FkP = F_s[k*3 + P];
                                    FlQ = F_s[l*3 + Q];

                                    CC[i][j][k][l] += FiA*FjB*FkP*FlQ*hat_CC_temp[A][B][P][Q]/J_s;


                                }
                            }
                        }
                    }

                }
            }
        }
    }


}
