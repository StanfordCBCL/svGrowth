# This file is the main run file for the Equilibrated G&R implementation from Lattore and Humphrey 2018
import math
from vessel import Vessel
from geometry import find_iv_geom, find_equil_geom
from stress import update_sigma
from timestep import update_time_step_py


def main():
    """
    Sets up the Vessel object, runs the G&R simulation, and writes
    results to 'GnR_out.txt'.
    """

    # ------------------------------------------------------------------
    # 1) Unit conversions and initial setup
    # ------------------------------------------------------------------
    mm_to_m   = 1.0e-3  # mllimeters to meters
    kPa_to_Pa = 1.0e3   # kiloPascals to Pascals

    # Create a new Vessel object (assuming you have a Vessel class)
    curr_vessel = Vessel()
    curr_vessel.vessel_name = "Cerebral"  # name label for this vessel

    # Time parameters
    n_days            = 701.0  # total time (in days) for simulation
    curr_vessel.dt    = 1.0  # time-step size (days)
    curr_vessel.nts   = int(n_days / curr_vessel.dt)  # number of time steps
    curr_vessel.sn    = 0  # current time-step index
    curr_vessel.s     = 0.0  # current physical time (days)

    # ------------------------------------------------------------------
    # 2) Geometric parameters (in vivo reference)
    # ------------------------------------------------------------------
    curr_vessel.a_h      = 1.40 * mm_to_m   # in-vivo reference inner radius (m)
    curr_vessel.h_h      = 0.12 * mm_to_m   # in-vivo reference medial thickness (m)
    curr_vessel.a_mid_h  = curr_vessel.a_h + 0.5 * curr_vessel.h_h  # mid-wall radius (m)
    curr_vessel.lambda_z_h = 1.0  # homeostatic axial stretch (dimensionless)

    # Initialize arrays/lists for loaded geometry history
    curr_vessel.a      = [0.0] * curr_vessel.nts  # loaded inner radius over time
    curr_vessel.a_mid  = [0.0] * curr_vessel.nts  # loaded mid-radius over time
    curr_vessel.h      = [0.0] * curr_vessel.nts  # loaded wall thickness over time


    # Set initial loaded geometry to the in-vivo reference
    curr_vessel.a[0]     = curr_vessel.a_h
    curr_vessel.a_mid[0] = curr_vessel.a_mid_h
    curr_vessel.h[0]     = curr_vessel.h_h

    # Initialize traction-free geometry arrays
    curr_vessel.A      = [0.0] * curr_vessel.nts  # traction-free inner radius
    curr_vessel.A_mid  = [0.0] * curr_vessel.nts  # traction-free mid-radius
    curr_vessel.H      = [0.0] * curr_vessel.nts  # traction-free wall thickness
    curr_vessel.lambda_z_pre = [0.0] * curr_vessel.nts  # axial pre-stretch array

    # ------------------------------------------------------------------
    # 3) Constituent material properties
    # ------------------------------------------------------------------
    curr_vessel.n_alpha = 3  # number of constituent families (e.g., elastin, muscle, collagen)

    # c_alpha_h: material parameters for each constituent
    c1_e  =  70.6 * kPa_to_Pa  # elastin stiffness scale
    c2_e  =   0.0              # elastin exponent
    c1_m  =  10.0 * kPa_to_Pa  # muscle stiffness scale
    c2_m  =   3.5              # muscle exponent
    c1_ct = 672.5 * kPa_to_Pa  # collagen stiffness scale
    c2_ct =  22.0              # collagen exponent
    curr_vessel.c_alpha_h = [
        c1_e,  c2_e,   # elastin  c1, c2
        c1_m,  c2_m,   # muscle   c1, c2
        c1_ct, c2_ct   # collagen c1, c2
    ]

    # Orientation angles (in radians); negative or -1 => isotropic
    eta_e_h  = -1.0                    # elastin isotropic marker
    eta_m_h  = 90.0 * math.pi / 180.0  # SMCs at 90 deg
    eta_ct_h = 90.0 * math.pi / 180.0  # collagen at 90 deg
    curr_vessel.eta_alpha_h = [eta_e_h, eta_m_h, eta_ct_h]

    # Deposition stretches
    g_eh  = 1.40  # elastin deposition stretch
    g_mh  = 1.30  # muscle deposition stretch
    g_cth = 1.08  # collagen deposition stretch
    curr_vessel.g_alpha_h = [g_eh, g_mh, g_cth]

    # We'll build G_alpha_h (3 per alpha) below after we define the vessel
    curr_vessel.G_alpha_h = [0.0]*(3*curr_vessel.n_alpha)

    # Tissue density
    curr_vessel.rhoR_h = 1050.0  # kg/m^3

    # Initialize referential density storage
    curr_vessel.rhoR       = [0.0]*curr_vessel.nts  # total referential density over time
    curr_vessel.rhoR[0]    = 1050.0                 # set initial density
    curr_vessel.rhoR_alpha = [0.0]*(curr_vessel.n_alpha * curr_vessel.nts)  # alpha-specific densities
    curr_vessel.mR_alpha   = [0.0]*(curr_vessel.n_alpha * curr_vessel.nts)  # alpha-specific production rates
    curr_vessel.k_alpha    = [0.0]*(curr_vessel.n_alpha * curr_vessel.nts)  # alpha-specific degradation rates

    # Homeostatic mass fractions
    phi_eh  = 0.02  # elastin fraction
    phi_mh  = 0.76  # SMC fraction
    phi_cth = 0.22  # collagen fraction
    curr_vessel.phi_alpha_h = [phi_eh, phi_mh, phi_cth]

    # Initial mass densities
    rhoR_eh  = phi_eh  * curr_vessel.rhoR_h  # elastin mass density
    rhoR_mh  = phi_mh  * curr_vessel.rhoR_h  # SMC mass density
    rhoR_cth = phi_cth * curr_vessel.rhoR_h  # collagen mass density
    curr_vessel.rhoR_alpha_h = [rhoR_eh, rhoR_mh, rhoR_cth]

    # Kinetic parameters
    # Degradation rates (1/day) for elastin, muscle, collagen
    k_elh = 0.0
    k_mh  = 1.0/14.0
    k_cth = 1.0/10.0
    curr_vessel.k_alpha_h = [k_elh, k_mh, k_cth]

    # Production parameters (initialized below in a loop)
    curr_vessel.mR_alpha_h = [0.0]*curr_vessel.n_alpha  # set below in loop

    # Gains (coefficients in stress-/WSS-mediated production and degradation)
    # Stress-mediated production
    K_sigma_p_elh = 0.0
    K_sigma_p_mh  = 2.0
    K_sigma_p_cth = 1.0
    curr_vessel.K_sigma_p_alpha_h = [K_sigma_p_elh, K_sigma_p_mh, K_sigma_p_cth]

    # Stress-mediated degradation
    K_sigma_d_elh = 0.0
    K_sigma_d_mh  = 1.0
    K_sigma_d_cth = 1.0
    curr_vessel.K_sigma_d_alpha_h = [K_sigma_d_elh, K_sigma_d_mh, K_sigma_d_cth]

    # WSS-mediated production
    K_tauw_p_elh = 0.0
    K_tauw_p_mh  = 1.0
    K_tauw_p_cth = 0.5
    curr_vessel.K_tauw_p_alpha_h = [K_tauw_p_elh, K_tauw_p_mh, K_tauw_p_cth]

    # WSS-mediated degradation
    K_tauw_d_elh = 0.0
    K_tauw_d_mh  = 0.0
    K_tauw_d_cth = 0.0
    curr_vessel.K_tauw_d_alpha_h = [K_tauw_d_elh, K_tauw_d_mh, K_tauw_d_cth]

    # ------------------------------------------------------------------
    # 4) Loading variables (pressure, flow, etc.)
    # ------------------------------------------------------------------
    curr_vessel.P     = 14.18 * kPa_to_Pa  # current luminal pressure (Pa)
    curr_vessel.P_h   = 14.18 * kPa_to_Pa  # homeostatic pressure (Pa)
    curr_vessel.Q     = 1.0                # current volumetric flow (m^3/day, or dimensionless scale)
    curr_vessel.Q_h   = 1.0                # homeostatic flow
    # WSS in Pa ( bar_tauw = Q / a^3 )
    curr_vessel.bar_tauw   = 1.0 / (curr_vessel.a_h**3)  # current wall shear stress (Pa)
    curr_vessel.bar_tauw_h = 1.0 / (curr_vessel.a_h**3)  # homeostatic wall shear stress (Pa)

    curr_vessel.lambda_th_curr = 1.0  # current circumferential stretch
    curr_vessel.lambda_z_curr  = 1.0  # current axial stretch

    curr_vessel.lambda_alpha_tau = [0.0]*(curr_vessel.nts * curr_vessel.n_alpha)  # alpha-specific stretch histories
    curr_vessel.lambda_z_tau     = [0.0]*curr_vessel.nts  # axial stretch history
    curr_vessel.lambda_z_tau[0]  = 1.0

    # ------------------------------------------------------------------
    # 5) Stresses, active stress, etc.
    # ------------------------------------------------------------------
    curr_vessel.sigma_h = [0.0, 0.0, 0.0]  # homeostatic stress in radial/circ/axial
    curr_vessel.sigma    = [0.0, 0.0, 0.0]  # current mixture stress
    curr_vessel.Cbar     = [0.0, 0.0, 0.0]  # a stiffness-like measure for each direction

    # Active stress parameters
    curr_vessel.alpha_active = [0, 1, 0]  # elastin=0, muscle=1, collagen=0
    curr_vessel.a_act = [0.0]*curr_vessel.nts
    curr_vessel.a_act[0] = curr_vessel.a_h  # reference radius for active muscle tone

    curr_vessel.T_act_h = 170.0 * kPa_to_Pa  # homeostatic max active stress
    curr_vessel.T_act   = 170.0 * kPa_to_Pa  # current active stress
    curr_vessel.k_act   = 1.0 / 7.0          # active remodeling rate (1/day)
    curr_vessel.lambda_0 = 0.4               # minimum contractile stretch
    curr_vessel.lambda_m = 1.1               # maximum contractile stretch
    curr_vessel.CB       = 0.8326            # basal tone coefficient
    curr_vessel.CS       = curr_vessel.CB / 2.0  # shear-sensitivity coefficient

    # ------------------------------------------------------------------
    # 6) Initialize each constituent
    # ------------------------------------------------------------------
    for alpha in range(curr_vessel.n_alpha):
        eta_alpha   = curr_vessel.eta_alpha_h[alpha]
        g_alpha_val = curr_vessel.g_alpha_h[alpha]

        # Deposition tensor G_alpha_h
        # 3 directions: radial (0), circumferential (1), axial (2)
        if eta_alpha > 0.0:
            # anisotropic
            curr_vessel.G_alpha_h[3*alpha + 0] = 0.0
            curr_vessel.G_alpha_h[3*alpha + 1] = g_alpha_val * math.sin(eta_alpha)
            curr_vessel.G_alpha_h[3*alpha + 2] = g_alpha_val * math.cos(eta_alpha)
        else:
            # isotropic
            curr_vessel.G_alpha_h[3*alpha + 0] = 1.0 / (g_alpha_val**2)
            curr_vessel.G_alpha_h[3*alpha + 1] = g_alpha_val
            curr_vessel.G_alpha_h[3*alpha + 2] = g_alpha_val

        # Homeostatic mass production
        curr_vessel.mR_alpha_h[alpha] = (
            curr_vessel.k_alpha_h[alpha] * curr_vessel.rhoR_alpha_h[alpha]
        )

        # Initialize time histories
        curr_vessel.rhoR_alpha[curr_vessel.nts*alpha + 0] = curr_vessel.rhoR_alpha_h[alpha]
        curr_vessel.mR_alpha[curr_vessel.nts*alpha + 0]   = curr_vessel.mR_alpha_h[alpha]
        curr_vessel.k_alpha[curr_vessel.nts*alpha + 0]    = curr_vessel.k_alpha_h[alpha]

        # Stretch histories
        curr_vessel.lambda_alpha_tau[curr_vessel.nts*alpha + 0] = 1.0

    # ------------------------------------------------------------------
    # 7) Find the true initial stress state in vivo
    # ------------------------------------------------------------------
    print("Inner radius:", curr_vessel.a[0],
          "Thickness:", curr_vessel.h[0],
          "Ref Density:", curr_vessel.rhoR_alpha[curr_vessel.nts * 1 + 0])

    find_iv_geom(curr_vessel)

    print("Inner radius:", curr_vessel.a[0],
          "Thickness:", curr_vessel.h[0],
          "Ref Density:", curr_vessel.rhoR_alpha[curr_vessel.nts * 1 + 0])

    # Then update mixture-based stress
    update_sigma(curr_vessel)

    # Store homeostatic stresses
    curr_vessel.sigma_h    = curr_vessel.sigma[:]  # copy
    curr_vessel.bar_tauw_h = curr_vessel.Q_h / (curr_vessel.a[0]**3)
    curr_vessel.a_act[0]   = curr_vessel.a[0]
    curr_vessel.f_h        = (math.pi * curr_vessel.h[0] *
                              (2.0 * curr_vessel.a[0] + curr_vessel.h[0]) *
                              curr_vessel.sigma[2])

    # Optionally find the traction-free geometry, if needed
    # curr_vessel.num_exp_flag = 1
    # find_tf_geom_py(curr_vessel)
    curr_vessel.num_exp_flag = 0

    # ------------------------------------------------------------------
    # 8) Setup file I/O for G&R output
    # ------------------------------------------------------------------
    GnR_out = open("GnR_out.txt", "w")

    # Initial mechanobiological equilibrium
    sigma_rel = ((curr_vessel.sigma[1] + curr_vessel.sigma[2]) /
                 (curr_vessel.sigma_h[1] + curr_vessel.sigma_h[2]))
    tauw_rel  = curr_vessel.bar_tauw / curr_vessel.bar_tauw_h
    curr_vessel.mb_equil = 1.0 \
        + curr_vessel.K_sigma_p_alpha_h[2] * (sigma_rel - 1.0) \
        - curr_vessel.K_tauw_p_alpha_h[2]  * (tauw_rel - 1.0)

    # Write initial state to file
    GnR_out.write(f"{curr_vessel.a[0]}\t{curr_vessel.h[0]}\t"
                  f"{curr_vessel.rhoR_alpha[curr_vessel.nts*1 + 0]}\t"
                  f"{curr_vessel.rhoR_alpha[curr_vessel.nts*2 + 0]}\t"
                  f"{curr_vessel.mb_equil}\t"
                  f"{curr_vessel.P / curr_vessel.P_h}\n")

    # ------------------------------------------------------------------
    # 9) Run the G&R time stepping
    # ------------------------------------------------------------------
    perturb_offset = 140.0
    for sn in range(1, curr_vessel.nts):
        curr_vessel.s = curr_vessel.dt * sn

        # Mechanical perturbation after 'perturb_offset' days
        if curr_vessel.s > perturb_offset:
            factor = (1.0 - math.exp(-(curr_vessel.s - perturb_offset)/10.0))
            curr_vessel.P = curr_vessel.P_h * (1.0 + 0.5 * factor)
            curr_vessel.lambda_z_curr = curr_vessel.lambda_z_h * (1.0 + 0.5 * factor)
            curr_vessel.Q = curr_vessel.Q_h * (1.0 + 0.5 * factor)

        # Update vessel index and solve
        curr_vessel.sn = sn
        update_time_step_py(curr_vessel)
        print("---------------------------")

        # Write model outputs
        GnR_out.write(f"{curr_vessel.a[sn]}\t{curr_vessel.h[sn]}\t"
                      f"{curr_vessel.rhoR_alpha[curr_vessel.nts*1 + sn]}\t"
                      f"{curr_vessel.rhoR_alpha[curr_vessel.nts*2 + sn]}\t"
                      f"{curr_vessel.mb_equil}\t"
                      f"{curr_vessel.P / curr_vessel.P_h}\n")

        # Store axial stretch history
        curr_vessel.lambda_z_tau[sn] = curr_vessel.lambda_z_curr

    # ------------------------------------------------------------------
    # 10) Find the final equilibrated solution
    # ------------------------------------------------------------------
    find_equil_geom(curr_vessel)

    # Write equilibrium geometry
    GnR_out.write(f"{curr_vessel.a_e}\t{curr_vessel.h_e}\t"
                  f"{curr_vessel.rho_m_e}\t{curr_vessel.rho_c_e}\t"
                  f"{curr_vessel.mb_equil_e}\t"
                  f"{curr_vessel.P / curr_vessel.P_h}\n")
    GnR_out.close()

    # Print the final equilibrium state
    print("a_e", curr_vessel.a_e, "h_e", curr_vessel.h_e)
