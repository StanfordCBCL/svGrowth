
class Vessel:
    def __init__(self):
        # Name/identifier
        self.vessel_name = "default"

        # Time variables
        self.nts = 0          # total number of time steps
        self.dt = 0.0         # time increment
        self.sn = 0           # current time step index
        self.s = 0.0          # actual current time

        # Geometric quantities (reference states)
        self.A_h = 0.0
        self.B_h = 0.0
        self.H_h = 0.0
        self.A_mid_h = 0.0
        self.a_h = 0.0
        self.b_h = 0.0
        self.h_h = 0.0
        self.a_mid_h = 0.0
        self.lambda_z_h = 0.0

        # Evolving in vivo references
        self.a = [0.0]
        self.a_mid = [0.0]
        self.b = [0.0]
        self.h = [0.0]

        # Evolving traction-free references
        self.A = [0.0]
        self.A_mid = [0.0]
        self.B = [0.0]
        self.H = [0.0]
        self.lambda_z_pre = [0.0]

        # Number of constituents
        self.n_alpha = 0

        # Material properties
        self.c_alpha_h = [0.0]
        self.eta_alpha_h = [0.0]
        self.g_alpha_h = [0.0]
        self.G_alpha_h = [0.0]

        # Mass fractions, referential densities, kinetic quantities
        self.phi_alpha_h = [0.0]
        self.rhoR_alpha_h = [0.0]
        self.mR_alpha_h = [0.0]
        self.k_alpha_h = [0.0]
        self.K_sigma_p_alpha_h = [0.0]
        self.K_sigma_d_alpha_h = [0.0]
        self.K_tauw_p_alpha_h = [0.0]
        self.K_tauw_d_alpha_h = [0.0]
        self.rhoR_h = 0.0

        # Histories
        self.rhoR = [0.0]
        self.rhoR_alpha = [0.0]
        self.mR_alpha = [0.0]
        self.k_alpha = [0.0]

        # Reference loading quantities
        self.P_h = 0.0
        self.f_h = 0.0
        self.bar_tauw_h = 0.0
        self.Q_h = 0.0
        self.sigma_h = [0.0]

        # Current loading quantities
        self.lambda_th_curr = 0.0
        self.lambda_z_curr = 0.0
        self.P = 0.0
        self.f = 0.0
        self.bar_tauw = 0.0
        self.Q = 0.0
        self.sigma = [0.0]
        self.Cbar = [0.0]
        self.lambda_alpha_tau = [0.0]
        self.lambda_z_tau = [0.0]
        self.mb_equil = 0.0

        # Active stress quantities
        self.alpha_active = [0]   # boolean flags for active constituents
        self.a_act = [0.0]        # active radius history
        self.T_act = 0.0
        self.T_act_h = 0.0
        self.k_act = 0.0
        self.lambda_0 = 0.0
        self.lambda_m = 0.0
        self.CB = 0.0
        self.CS = 0.0

        # Mechanobiologically equilibrated quantities
        self.a_e = 0.0
        self.h_e = 0.0
        self.rho_c_e = 0.0
        self.rho_m_e = 0.0
        self.f_z_e = 0.0
        self.mb_equil_e = 0.0

        # Flags
        self.num_exp_flag = 0

