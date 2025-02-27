#ifndef MATFUN
#define MATFUN


using std::vector;

double mat_det(vector<double> F);
double mat_ddot(vector<double> A, vector<double> B);
vector<double> mat_inv(vector<double> A);
vector<double> mat_trans(vector<double> A);
vector<double> mat_mul(vector<double> A, vector<double> B);
vector<double> ten_dyadprod(vector<double> A, vector<double> B);
vector<double> ten_symmprod(vector<double> A, vector<double> B);
vector<double> scl_mul(vector<double> A, double b);
vector<double> elem_mul(vector<double> A, vector<double> B);
vector<double> mat_pol_dec(vector<double> A);

void get_material_stiffness(vector<double> F_alpha_ntau_s, double hat_dS_dlambda2_alpha, double J_s, double hat_CC[3][3][3][3]);
void get_active_material_stiffness(double S, double CC[3][3][3][3],vector<double> F_s, double J_s);

#endif /* MATFUN */