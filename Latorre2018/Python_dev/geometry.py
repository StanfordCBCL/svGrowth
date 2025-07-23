import math
import numpy as np
from scipy.optimize import root, brentq, toms748
from stress import update_sigma

def equil_obj_f_py(vars, curr_vessel):
    """
    This function:
      - Reads in the unknowns (a_e_guess, h_e_guess, rho_c_e_guess, f_z_e_guess).
      - Computes 4 residual equations (J1, J2, J3, J4).
      - Updates 'curr_vessel' fields for the current guess.
      - Returns these 4 residuals as a NumPy array of shape (4,).

    Parameters
    ----------
    vars : list or array-like of length 4
        [a_e_guess, h_e_guess, rho_c_e_guess, f_z_e_guess].
    curr_vessel : Vessel
        Python Vessel object containing fields with current vessel parameters and histories.

    Returns
    -------
    residuals : np.ndarray of length 4
        [J1, J2, J3, J4] as described in the original code.
    """

    # ----------------------------------------------------------------
    # 1) Extract the unknown variables from 'vars'
    # ----------------------------------------------------------------
    a_e_guess, h_e_guess, rho_c_e_guess, f_z_e_guess = vars

    # ----------------------------------------------------------------
    # 2) Compute intermediate mechanical quantities
    # ----------------------------------------------------------------
    # (a) "Sigma-e-th-lmb" and "Sigma-e-z-lmb" from equilibrium relationships
    sigma_e_th_lmb = curr_vessel.P * a_e_guess / h_e_guess
    sigma_e_z_lmb  = f_z_e_guess / (math.pi * h_e_guess * (2.0 * a_e_guess + h_e_guess))

    # (b) Wall shear stress from Poiseuille flow assumption
    bar_tauw_e = curr_vessel.Q / (a_e_guess**3)  # WSS ~ Q / a^3 for constant viscosity

    # (c) Ratio of stress to WSS-mediated matrix production
    #     i.e., eta_K = K_sigma_p_alpha_h[1] / K_tauw_p_alpha_h[1]
    eta_K = (curr_vessel.K_sigma_p_alpha_h[1] /
             curr_vessel.K_tauw_p_alpha_h[1])

    # (d) Deviations from homeostatic stress and homeostatic WSS
    delta_sigma = ((sigma_e_th_lmb + sigma_e_z_lmb) /
                   (curr_vessel.sigma_h[1] + curr_vessel.sigma_h[2])
                   - 1.0)
    delta_tauw = (bar_tauw_e / curr_vessel.bar_tauw_h) - 1.0

    # ----------------------------------------------------------------
    # 3) Compute geometry-based parameters and volumetric stretch (J_e)
    # ----------------------------------------------------------------
    a_h = curr_vessel.a_h
    h_h = curr_vessel.h_h

    # "Radial" stretch is h_e / h_h, "circumferential" is (a_e + h_e/2)/(a_h + h_h/2)
    lambda_r_e  = h_e_guess / h_h
    lambda_th_e = (a_e_guess + h_e_guess / 2.0) / (a_h + h_h / 2.0)
    lambda_z_e  = curr_vessel.lambda_z_curr  # from the vessel object
    F_e = [lambda_r_e, lambda_th_e, lambda_z_e]

    # Total volumetric ratio
    J_e = lambda_r_e * lambda_th_e * lambda_z_e

    # ----------------------------------------------------------------
    # 4) Compute equilibrated mass densities
    # ----------------------------------------------------------------
    # (a) Elastin is assumed to have no turnover, so it's scaled by 1/J_e
    rho_el_e = curr_vessel.rhoR_alpha_h[0] / J_e

    # (b) Ratios for muscle vs collagen turnover
    eta_q   = (curr_vessel.k_alpha_h[1] /
               curr_vessel.k_alpha_h[2])  # ratio of SMC:collagen degradation rates
    eta_ups = (curr_vessel.K_sigma_p_alpha_h[1] /
               curr_vessel.K_sigma_p_alpha_h[2])  # ratio of stress-mediated production gains

    # (c) Muscle density at equilibrium depends on ratio of J_e*rho_c : rho_c_h, raised to (eta_q*eta_ups)
    rho_m_e = (curr_vessel.rhoR_alpha_h[1] / J_e *
               (J_e * rho_c_e_guess / curr_vessel.rhoR_alpha_h[2])**(eta_q * eta_ups))

    # Collect partial densities: [elastin, muscle, collagen].
    # Collagen is directly the guess (rho_c_e_guess), i.e. already divided by J_e in the code logic.
    rho_alpha = [rho_el_e, rho_m_e, rho_c_e_guess]

    # Homeostatic total density (sum of all constituents in homeostatic state)
    rho_h = curr_vessel.rhoR_h

    # ----------------------------------------------------------------
    # 5) Compute the 3-directional Cauchy stress from each constituent
    # ----------------------------------------------------------------
    n_alpha = curr_vessel.n_alpha

    # We'll store partial Cauchy stress for each direction in sigma_e_dir
    sigma_e_dir = [0.0, 0.0, 0.0]

    # Active stress contribution
    #  C = CB - CS*(delta_tauw)
    C = (curr_vessel.CB -
         curr_vessel.CS * delta_tauw)

    # For now, the code sets lambda_act=1 for equilibrium,
    # and uses a "parab_act" expression with (lambda_m - lambda_0).
    lambda_act = 1.0
    parab_act  = 1.0 - ((curr_vessel.lambda_m - lambda_act) /
                        (curr_vessel.lambda_m - curr_vessel.lambda_0))**2

    # Then the "hat_sigma_act_e"
    hat_sigma_act_e = (curr_vessel.T_act *
                       (1.0 - math.exp(-C**2)) *
                       lambda_act * parab_act)

    # We'll build up sigma_e_dir for each direction by summing the constituent contributions.
    # (In the original code, hat_sigma_alpha_dir is dimension 3 * n_alpha, but we only need partial sums.)
    for alpha in range(n_alpha):
        for dir_ in range(3):
            # If anisotropic and eta_alpha_h[alpha] > 0, we use "deposition stretch" g_alpha_h[alpha].
            if curr_vessel.eta_alpha_h[alpha] > 0.0:
                lambda_alpha_ntau_s = curr_vessel.g_alpha_h[alpha]
                # 2nd PK stress "hat" at equilibrium
                L2m1 = (lambda_alpha_ntau_s**2) - 1.0
                exponent_part = (curr_vessel.c_alpha_h[2 * alpha + 1] *
                                 (L2m1**2))

                hat_S_alpha = (curr_vessel.c_alpha_h[2 * alpha] *
                               L2m1 *
                               math.exp(exponent_part))

                # Convert 2nd PK to Cauchy with G_alpha_h
                G_dir = curr_vessel.G_alpha_h[3 * alpha + dir_]
                hat_sigma_alpha_dir = G_dir * hat_S_alpha * G_dir

            else:
                # If "eta_alpha_h[alpha] <= 0", treat as isotropic with constant c_alpha_h
                hat_S_alpha = curr_vessel.c_alpha_h[2 * alpha]
                # Again, G_dir factor:
                G_dir = curr_vessel.G_alpha_h[3 * alpha + dir_]
                hat_sigma_alpha_dir = G_dir * hat_S_alpha * G_dir

                # If k_alpha_h[alpha] == 0 => no turnover => scale by F_e[dir]^2
                # (the code multiplies by F_e[dir]*... to account for mixture deformation).
                if curr_vessel.k_alpha_h[alpha] == 0:
                    hat_sigma_alpha_dir *= (F_e[dir_]**2)

            # Add in this constituent's partial stress
            sigma_e_dir[dir_] += (rho_alpha[alpha] / rho_h) * hat_sigma_alpha_dir

            # If this alpha is "active" and direction is circumferential (dir_=1),
            # add the active stress
            if curr_vessel.alpha_active[alpha] == 1 and dir_ == 1:
                sigma_e_dir[dir_] += (rho_alpha[alpha] / rho_h) * hat_sigma_act_e

    # ----------------------------------------------------------------
    # 6) Form the 4 residual equations: J1, J2, J3, J4
    # ----------------------------------------------------------------
    # J1: mechano-mediated matrix production => eta_K * delta_sigma - delta_tauw
    J1 = eta_K * delta_sigma - delta_tauw

    # J2: mixture mass balance => sum of eq. densities minus total homeostatic density
    #     elastin + muscle + guessed collagen = rho_h
    J2 = (rho_el_e + rho_m_e + rho_c_e_guess) - rho_h

    # J3: difference in circumferential stress minus radial (Lagrange) stress
    #     minus the 'lmb' quantity from pressure
    #     => sigma_e_circ - sigma_e_rad - sigma_e_th_lmb
    # in code: sigma_e_dir[1] - sigma_e_dir[0] - sigma_e_th_lmb
    J3 = sigma_e_dir[1] - sigma_e_dir[0] - sigma_e_th_lmb

    # J4: difference in axial stress minus radial stress minus the 'lmb' from f_z
    #     => sigma_e_ax - sigma_e_rad - sigma_e_z_lmb
    J4 = sigma_e_dir[2] - sigma_e_dir[0] - sigma_e_z_lmb

    # ----------------------------------------------------------------
    # 7) Store the "equilibrated" results back into 'curr_vessel'
    # ----------------------------------------------------------------
    curr_vessel.a_e     = a_e_guess
    curr_vessel.h_e     = h_e_guess
    curr_vessel.rho_c_e = rho_c_e_guess * J_e  # re-scale by J_e
    curr_vessel.rho_m_e = rho_m_e * J_e
    curr_vessel.f_z_e   = f_z_e_guess

    # Mechanobiological equilibrium measure
    #   mb_equil_e = 1 + K_sigma_p_alpha_h[2]*delta_sigma - K_tauw_p_alpha_h[2]*delta_tauw
    curr_vessel.mb_equil_e = (
        1.0
        + curr_vessel.K_sigma_p_alpha_h[2] * delta_sigma
        - curr_vessel.K_tauw_p_alpha_h[2] * delta_tauw
    )

    # ----------------------------------------------------------------
    # 8) Return the 4 residuals in a NumPy array
    # ----------------------------------------------------------------
    return np.array([J1, J2, J3, J4], dtype=float)


def find_equil_geom(curr_vessel):
    """
    Uses scipy.optimize.root to solve for
      a_e, h_e, rho_c_e, and f_z_e.
    """

    # -----------------------------------------------------------------
    # 1. Compute the fold changes (gamma, epsilon, lambda)
    #    from the Vessel object's current vs. homeostatic values
    # -----------------------------------------------------------------
    gamma = curr_vessel.P / curr_vessel.P_h  # fold-change in pressure
    epsilon = curr_vessel.Q / curr_vessel.Q_h  # fold-change in flow
    lam = curr_vessel.lambda_z_curr / curr_vessel.lambda_z_h  # fold-change in axial stretch

    # -----------------------------------------------------------------
    # 2. Gather homeostatic geometry
    # -----------------------------------------------------------------
    a_h = curr_vessel.a_h
    h_h = curr_vessel.h_h

    # -----------------------------------------------------------------
    # 3. Create initial guesses
    # -----------------------------------------------------------------
    a_e_guess = (epsilon ** (1.0 / 3.0)) * a_h
    h_e_guess = gamma * (epsilon ** (1.0 / 3.0)) * h_h
    rho_c_e_guess = curr_vessel.rhoR_alpha_h[2]  # example from code
    f_z_e_guess = (curr_vessel.f_h * (h_e_guess * (2 * a_e_guess + h_e_guess)) /
                   (h_h * (2 * a_h + h_h)))

    x0 = [a_e_guess, h_e_guess, rho_c_e_guess, f_z_e_guess]

    # -----------------------------------------------------------------
    # 4. Define the wrapper objective for scipy's solver
    # -----------------------------------------------------------------
    def objective(vars):
        return equil_obj_f_py(vars, curr_vessel)

    # -----------------------------------------------------------------
    # 5. Call scipy's root-finding method
    # -----------------------------------------------------------------
    # You can choose method='hybr' to replicate GSL's "hybrids" approach,
    # or try method='lm' (Levenberg–Marquardt), etc.
    sol = root(objective, x0, method='hybr')

    # -----------------------------------------------------------------
    # 6. Check convergence
    # -----------------------------------------------------------------
    if sol.success:
        print("Solver converged successfully. Message:", sol.message)
        # Extract solution
        a_e_sol, h_e_sol, rho_c_e_sol, f_z_e_sol = sol.x
        print(f"Solution: a_e={a_e_sol}, h_e={h_e_sol}, rho_c_e={rho_c_e_sol}, f_z_e={f_z_e_sol}")
    else:
        print("Solver failed to converge. Message:", sol.message)
        a_e_sol, h_e_sol, rho_c_e_sol, f_z_e_sol = [None] * 4

    return 0  # or return sol

def iv_obj_f_py(a_mid_guess, curr_vessel):
    """
    This function calculates the difference between:
      - The theoretical Laplace-based circumferential stress (sigma_t_th = P * a / h)
      - The mixture-based circumferential stress (curr_vessel.sigma[1])
    and returns that difference.

    Parameters
    ----------
    a_mid_guess : float
        Proposed mid-radius for the vessel at the current time step.
    curr_vessel : Vessel
        Python object with attributes corresponding the current vessel state (including histories)

    Returns
    -------
    float
        J = (mixture-based circumferential stress) - (Laplace-based stress)
    """
    sn = curr_vessel.sn

    # Initialize local variables
    a = 0.0
    h = 0.0
    lambda_t = 0.0
    lambda_z = 0.0
    J_s = 0.0

    # Update geometry according to the current time step 'sn'
    if sn > 0:
        # Deposition stretch in circumferential direction
        lambda_t = a_mid_guess / curr_vessel.a_mid[0]
        # Axial stretch remains the same
        lambda_z = curr_vessel.lambda_z_curr
        # Ratio of referential densities at sn vs. 0
        J_s = curr_vessel.rhoR[sn] / curr_vessel.rhoR[0]

        # New thickness, radius, mid-radius
        h = (J_s / (lambda_t * lambda_z)) * curr_vessel.h[0]
        a = a_mid_guess - (h / 2.0)
        curr_vessel.a_mid[sn] = a_mid_guess
        curr_vessel.a[sn] = a
        curr_vessel.h[sn] = h

        # Current in-vivo stretches
        curr_vessel.lambda_th_curr = lambda_t
        curr_vessel.lambda_z_curr  = lambda_z

        # Update WSS (bar_tauw)
        curr_vessel.bar_tauw = curr_vessel.Q / (a**3)

    else:
        # For sn=0, use the homeostatic mid radius 'a_mid_h'
        lambda_t = a_mid_guess / curr_vessel.a_mid_h
        lambda_z = curr_vessel.lambda_z_curr
        J_s      = 1.0  # Because it's the reference (no density change from 0 to 0)

        h = (J_s / (lambda_t * lambda_z)) * curr_vessel.h_h
        a = a_mid_guess - (h / 2.0)
        curr_vessel.a_mid[sn] = a_mid_guess
        curr_vessel.a[sn] = a
        curr_vessel.h[sn] = h

        # Here the code sets the "lambda_th_curr" to 1.0 for the initial guess
        curr_vessel.lambda_th_curr = 1.0
        curr_vessel.lambda_z_curr  = lambda_z

        curr_vessel.bar_tauw = curr_vessel.Q / (a**3)

    # Theoretical Laplace-based circumferential stress
    sigma_t_th = curr_vessel.P * (a / h)

    # Update mixture-based stress using your Python 'update_sigma'
    update_sigma(curr_vessel)

    # Update the axial force f = π * h * (2a + h) * sigma_axial
    # In your convention: sigma[2] = axial direction
    curr_vessel.f = math.pi * h * (2.0 * a + h) * curr_vessel.sigma[2]

    # The returned value: difference between mixture-based (sigma[1]) and Laplace-based (sigma_t_th)
    J = curr_vessel.sigma[1] - sigma_t_th

    return J

def find_iv_geom(curr_vessel):
    """
    This finds the in vivo (loaded) configuration by solving
    for a_mid (the mid radius) such that iv_obj_f == 0.
    """

    # 1) Extract the current time index
    sn = curr_vessel.sn

    # 2) Define the search interval around the current guess (±5%)
    a_mid_guess = curr_vessel.a_mid[sn]
    a_mid_low   = 0.95 * a_mid_guess
    a_mid_high  = 1.05 * a_mid_guess

    # 3) Define a Python wrapper for the objective function
    def objective(a_mid):
        return iv_obj_f_py(a_mid, curr_vessel)

    # 4) Use 'brentq' to find a root within [a_mid_low, a_mid_high]
    try:
        root = toms748(objective, a_mid_low, a_mid_high, xtol=1e-5)

        # On success, store the solution
        # curr_vessel.a_mid[sn] = root

        # We can mimic returning a status code: 0 -> success
        return 0

    except ValueError:
        # brentq throws an exception if it fails (e.g. sign of f(low)*f(high) > 0)
        # For GSL compatibility, we might return a non-zero status.
        return 1

def tf_obj_f_py(vars, curr_vessel):
    """
    Parameters
    ----------
    vars : array-like of length 2
        [lambda_th_ul_guess, lambda_z_ul_guess]
    curr_vessel : Vessel
        Python object that contains current vessel state information and histories

    Returns
    -------
    residuals : np.ndarray of length 2
        [J1, J2]
    """

    # Unpack the unknowns
    lambda_th_ul_guess, lambda_z_ul_guess = vars

    # Current time index
    sn = curr_vessel.sn

    # Homeostatic mid-radius vs. current mid-radius
    a_mid_0 = curr_vessel.a_mid[0]
    a_mid   = curr_vessel.a_mid[sn]

    # Reference stretches from the "reference" to "loaded" state
    lambda_th_ref = a_mid / a_mid_0
    lambda_z_ref  = curr_vessel.lambda_z_h

    # The code also computes a volumetric ratio J_s,
    # but it doesn't use it directly in the final residuals:
    # J_s = curr_vessel.rhoR[sn] / curr_vessel.rhoR[0]

    # Update the current total stretches in the vessel
    curr_vessel.lambda_th_curr = lambda_th_ul_guess * lambda_th_ref
    curr_vessel.lambda_z_curr  = lambda_z_ul_guess  * lambda_z_ref

    # Recompute the mixture-based stresses with the new stretches
    update_sigma(curr_vessel)

    # Residuals:
    # J1 = sigma[1], J2 = sigma[2]
    J1 = curr_vessel.sigma[1]
    J2 = curr_vessel.sigma[2]

    # Return them as a NumPy array
    return np.array([J1, J2], dtype=float)


def find_tf_geom(curr_vessel):
    """
    1) Sets the vessel's loads (P, f, T_act) to zero to find the traction-free
       geometry for the current time point 'sn'.
    2) Solves for [lambda_th_ul, lambda_z_ul] using the tf_obj_f_py objective.
    3) Stores the traction-free geometry (A_mid, lambda_z_pre, H) in
       curr_vessel for the current time step 'sn'.
    4) Resets the loads to their homeostatic values.

    Returns
    -------
    int
        (0 on success).
    """

    # -----------------------------------------------------------
    # 1) Gather initial guesses for the unknowns
    # -----------------------------------------------------------
    sn = curr_vessel.sn
    lambda_th_ul_guess = 0.95
    lambda_z_ul_guess = 0.95 * curr_vessel.lambda_z_curr

    # -----------------------------------------------------------
    # 2) Temporarily set vessel loads to zero for traction-free
    # -----------------------------------------------------------
    old_P = curr_vessel.P
    old_f = curr_vessel.f
    old_T_act = curr_vessel.T_act
    old_lz_curr = curr_vessel.lambda_z_curr

    curr_vessel.P = 0.0
    curr_vessel.f = 0.0
    curr_vessel.T_act = 0.0

    # (The original sets 'curr_vessel->lambda_z_curr' to ???
    #  but the code does not override it here explicitly,
    #  so we leave it as is, just use the guesses for the solver.)

    # -----------------------------------------------------------
    # 3) Define the system of equations using tf_obj_f_py
    #    (already translated in previous steps)
    # -----------------------------------------------------------
    def system(vars_):
        # tf_obj_f_py expects [lambda_th_ul_guess, lambda_z_ul_guess]
        # and returns [J1, J2]
        return tf_obj_f_py(vars_, curr_vessel)

    # -----------------------------------------------------------
    # 4) Solve the 2x2 system with SciPy's root
    # -----------------------------------------------------------
    x0 = [lambda_th_ul_guess, lambda_z_ul_guess]
    sol = root(system, x0, method='hybr', tol=1e-7)

    # (Optional) check iteration or debug info:
    # print(f"Iteration count: {sol.nfev}, status: {sol.message}")

    # -----------------------------------------------------------
    # 5) On success, store traction-free geometry in 'curr_vessel'
    # -----------------------------------------------------------
    if sol.success:
        lambda_th_ul_sol, lambda_z_ul_sol = sol.x
        #   A_mid[sn]      = x[0] * a_mid[sn]
        #   lambda_z_pre[sn] = 1 / x[1]
        #   H[sn]          = 1.0 / (x[0]^2) * h[sn]
        # Note: x[0] = lambda_th_ul, x[1] = lambda_z_ul
        curr_vessel.A_mid[sn] = lambda_th_ul_sol * curr_vessel.a_mid[sn]
        curr_vessel.lambda_z_pre[sn] = 1.0 / lambda_z_ul_sol
        curr_vessel.H[sn] = (1.0 / (lambda_th_ul_sol ** 2)) * curr_vessel.h[sn]

    else:
        print("Solver failed to converge:", sol.message)
        # You could handle error or return a non-zero code here

    # -----------------------------------------------------------
    # 6) Return the vessel loads to their homeostatic (original) values
    # -----------------------------------------------------------
    curr_vessel.P = curr_vessel.P_h
    curr_vessel.f = curr_vessel.f_h
    curr_vessel.T_act = curr_vessel.T_act_h
    curr_vessel.lambda_z_curr = curr_vessel.lambda_z_h

    return 0
