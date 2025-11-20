import math


def constitutive(curr_vessel, lambda_alpha_s, alpha, ts, dir_):
    """
    Parameters
    ----------
    curr_vessel : Vessel
        The vessel object containing fields like eta_alpha_h, g_alpha_h, etc.
    lambda_alpha_s : float
        The current in-plane stretch for constituent alpha.
    alpha : int
        Index of the constituent.
    ts : int
        Current time step index.
    dir_ : int
        Direction index (radial=0, circumferential=1, axial=2).
        Note that in this snippet, dir_ is not actually used.

    Returns
    -------
    list of two floats [hat_S_alpha, hat_dSdC_alpha]
    """

    # Local variables
    lambda_alpha_ntau_s = 0.0
    Q1 = 0.0
    Q2 = 0.0
    hat_S_alpha = 0.0
    hat_dSdC_alpha = 0.0
    pol_mod = 0.0  # not used in this snippet, but preserved for completeness

    nts = curr_vessel.nts

    # Check if anisotropic (eta_alpha_h[alpha] >= 0 means we have an orientation)
    if curr_vessel.eta_alpha_h[alpha] >= 0.0:
        # Compute ratio of current stretch to reference
        lambda_alpha_ntau_s = (
            curr_vessel.g_alpha_h[alpha] * lambda_alpha_s
            / curr_vessel.lambda_alpha_tau[nts * alpha + ts]
        )

        # Clamp to 1 if below 1
        if lambda_alpha_ntau_s < 1.0:
            lambda_alpha_ntau_s = 1.0

        # Compute the exponents
        Q1 = (lambda_alpha_ntau_s**2 - 1.0)   # (lam^2 - 1)
        Q2 = (curr_vessel.c_alpha_h[2 * alpha + 1] * (Q1**2))
        # Compute partial stress (hat_S_alpha) and derivative (hat_dSdC_alpha)
        hat_S_alpha    = curr_vessel.c_alpha_h[2 * alpha] * Q1 * math.exp(Q2)
        hat_dSdC_alpha = curr_vessel.c_alpha_h[2 * alpha] * math.exp(Q2) * (1.0 + 2.0 * Q2)

    else:
        # If eta_alpha < 0, we treat it as isotropic with a constant c_alpha_h[2*alpha].
        # No stretch-based exponential term is used in this case.
        hat_S_alpha = curr_vessel.c_alpha_h[2 * alpha]
        hat_dSdC_alpha = 0.0  # No dependence on stretch in this model branch

    return [hat_S_alpha, hat_dSdC_alpha]


def update_sigma(curr_vessel):
    """
    Parameters
    ----------
    curr_vessel : Vessel
        The Vessel object (in Python) containing all of the
        fields that describe the current state and histories of the vessel
    """

    # ----------------------------------------------------
    # 1. Extract some attributes for convenience
    # ----------------------------------------------------
    s = curr_vessel.s
    sn = curr_vessel.sn
    nts = curr_vessel.nts
    dt = curr_vessel.dt

    # We set the default for taun_min to 0
    taun_min = 0

    # tau_max = 100 * (1 / k_alpha_h[2])
    tau_max = 100.0 * (1.0 / curr_vessel.k_alpha_h[2])

    # Vessel geometry at the initial reference
    a0 = curr_vessel.a[0]
    h0 = curr_vessel.h[0]

    # Current "global" stretches
    lambda_th_s = curr_vessel.lambda_th_curr
    lambda_z_s = curr_vessel.lambda_z_curr

    # Number of constituents
    n_alpha = curr_vessel.n_alpha

    # ----------------------------------------------------
    # 2. Compute the constituent-specific current stretches
    # ----------------------------------------------------
    lambda_alpha_s = [0.0] * n_alpha
    for alpha in range(n_alpha):
        eta_alpha = curr_vessel.eta_alpha_h[alpha]
        # If eta_alpha >= 0, it means the constituent is anisotropic
        # but defined by an angle 'eta_alpha' (or possibly isotropic if 0).
        if eta_alpha >= 0.0:
            # In-plane stretch = sqrt( (lambda_z * cos(eta))^2 + (lambda_th * sin(eta))^2 )
            lam_alpha = math.sqrt((lambda_z_s * math.cos(eta_alpha)) ** 2 +
                                  (lambda_th_s * math.sin(eta_alpha)) ** 2)
            lambda_alpha_s[alpha] = lam_alpha

            # Update the stored current stretch if not doing a numerical experiment
            if curr_vessel.num_exp_flag == 0:
                curr_vessel.lambda_alpha_tau[nts * alpha + sn] = lam_alpha

    # ----------------------------------------------------
    # 3. Current deformation gradient (F_s) and Jacobian (J_s)
    #    J_s = ratio of total referential density at current time
    #          to the "homeostatic" density
    # ----------------------------------------------------
    J_s = curr_vessel.rhoR[sn] / curr_vessel.rhoR_h
    F_s0 = J_s / (lambda_th_s * lambda_z_s)  # radial stretch component
    F_s1 = lambda_th_s  # circumferential stretch
    F_s2 = lambda_z_s  # axial stretch
    F_s = [F_s0, F_s1, F_s2]

    # ----------------------------------------------------
    # 4. Initialize local accumulators
    # ----------------------------------------------------
    sigma = [0.0, 0.0, 0.0]  # total Cauchy stress in radial/circ/axial directions
    Cbar = [0.0, 0.0, 0.0]  # "stiffness-like" measure

    # Variables for active stress
    a_act = 0.0
    q_act_0 = 0.0
    q_act_1 = 0.0
    q_act_2 = 0.0
    k_act = curr_vessel.k_act

    # Temporary arrays for the constituent-based partial stresses/stiffness
    hat_sigma_0 = [0.0, 0.0, 0.0]
    hat_sigma_1 = [0.0, 0.0, 0.0]
    hat_sigma_2 = [0.0, 0.0, 0.0]
    hat_Cbar_0 = [0.0, 0.0, 0.0]
    hat_Cbar_1 = [0.0, 0.0, 0.0]
    hat_Cbar_2 = [0.0, 0.0, 0.0]

    # ----------------------------------------------------
    # 5. Determine if we are beyond the initial time history
    #    and set up integration bounds
    # ----------------------------------------------------
    if s > tau_max:
        taun_min = sn - int(tau_max / dt)
    else:
        taun_min = 0

    n = (sn - taun_min) + 1  # number of integration points
    even_n = (n % 2 == 0)  # check if even for trapezoid step

    # ----------------------------------------------------
    # 6. Loop over each constituent and integrate history
    # ----------------------------------------------------
    for alpha in range(n_alpha):

        # Current kinetics from the newest "cohort"
        k_2 = curr_vessel.k_alpha[nts * alpha + sn]
        q_2 = 1.0
        mq_2 = curr_vessel.mR_alpha[nts * alpha + sn]

        # If this constituent is "active"
        if curr_vessel.alpha_active[alpha] == 1:
            a_2 = curr_vessel.a[sn]  # radius at current time
            q_act_2 = 1.0

        # We'll use the current stretch quantities as the "cohort" reference:
        F_tau = [F_s0, F_s1, F_s2]
        J_tau = J_s

        # ------------------------------------------------
        # 6.1. Evaluate stress from the newly formed cohort
        # ------------------------------------------------
        for dir_ in range(3):
            # The constitutive() function is assumed to be available.
            # It should return [hat_S_alpha, hat_dSdC_alpha].
            #   - hat_S_alpha    = partial stress measure
            #   - hat_dSdC_alpha = partial derivative measure
            [hat_S_alpha, hat_dSdC_alpha] = constitutive(curr_vessel,
                                                         lambda_alpha_s[alpha],
                                                         alpha,
                                                         sn,  # time index
                                                         dir_)  # direction

            # F_alpha_ntau_s = (F_s[dir] / F_tau[dir]) * G_alpha_h[3*alpha + dir]
            F_alpha_ntau_s = (F_s[dir_] / F_tau[dir_]) * curr_vessel.G_alpha_h[3 * alpha + dir_]

            hat_sigma_2[dir_] = (F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s) / J_s
            hat_Cbar_2[dir_] = (F_alpha_ntau_s ** 2) * hat_dSdC_alpha * (F_alpha_ntau_s ** 2) / J_s

        # Check if this constituent has a positive mass formation rate
        deg_check = (curr_vessel.mR_alpha_h[alpha] > 0.0)

        # If we have multiple time steps (sn>0) and the constituent actually remodels:
        if sn > 0 and deg_check:

            # ----------------------------------------------------
            # 6.2. Simpson (and trapezoid) integration backward in time
            # ----------------------------------------------------
            for taun in range(sn - 1, taun_min, -2):
                #  (a) 1st intermediate deformation gradient (time = taun)
                a_t = curr_vessel.a[taun]
                h_t = curr_vessel.h[taun]
                lambda_th_tau = (a_t + h_t / 2.0) / (a0 + h0 / 2.0)
                lambda_z_tau = curr_vessel.lambda_z_tau[taun]
                J_tau = curr_vessel.rhoR[taun] / curr_vessel.rhoR_h

                F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau)
                F_tau[1] = lambda_th_tau
                F_tau[2] = lambda_z_tau

                #  (b) 1st intermediate kinetics
                k_1 = curr_vessel.k_alpha[nts * alpha + taun]
                q_1 = math.exp(-(k_2 + k_1) * dt / 2.0) * q_2
                mq_1 = curr_vessel.mR_alpha[nts * alpha + taun] * q_1

                #  (c) If alpha is active, gather radius for active stress
                if curr_vessel.alpha_active[alpha] == 1:
                    a_1 = a_t
                    q_act_1 = math.exp(-k_act * dt) * q_act_2


                #  (d) 1st intermediate stress and stiffness
                for dir_ in range(3):
                    [hat_S_alpha, hat_dSdC_alpha] = constitutive(curr_vessel,
                                                                 lambda_alpha_s[alpha],
                                                                 alpha,
                                                                 taun,
                                                                 dir_)
                    F_alpha_ntau_s = (F_s[dir_] / F_tau[dir_]) * curr_vessel.G_alpha_h[3 * alpha + dir_]
                    hat_sigma_1[dir_] = (F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s) / J_s
                    hat_Cbar_1[dir_] = (F_alpha_ntau_s ** 2) * hat_dSdC_alpha * (F_alpha_ntau_s ** 2) / J_s

                #  (e) 2nd intermediate deformation gradient (time = taun-1)
                a_tm1 = curr_vessel.a[taun - 1]
                h_tm1 = curr_vessel.h[taun - 1]
                lambda_th_tau = (a_tm1 + h_tm1 / 2.0) / (a0 + h0 / 2.0)
                lambda_z_tau = curr_vessel.lambda_z_tau[taun - 1]
                J_tau = curr_vessel.rhoR[taun - 1] / curr_vessel.rhoR_h

                F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau)
                F_tau[1] = lambda_th_tau
                F_tau[2] = lambda_z_tau

                #  (f) 2nd intermediate kinetics
                k_0 = curr_vessel.k_alpha[nts * alpha + (taun - 1)]
                q_0 = math.exp(-(k_2 + 4.0 * k_1 + k_0) * dt / 3.0) * q_2
                mq_0 = curr_vessel.mR_alpha[nts * alpha + (taun - 1)] * q_0

                #  (g) Active radius
                if curr_vessel.alpha_active[alpha] == 1:
                    a_0 = a_tm1
                    q_act_0 = math.exp(-k_act * dt) * q_act_1


                #  (h) 2nd intermediate stress and stiffness
                for dir_ in range(3):
                    [hat_S_alpha, hat_dSdC_alpha] = constitutive(curr_vessel,
                                                                 lambda_alpha_s[alpha],
                                                                 alpha,
                                                                 taun - 1,
                                                                 dir_)
                    F_alpha_ntau_s = (F_s[dir_] / F_tau[dir_]) * curr_vessel.G_alpha_h[3 * alpha + dir_]
                    hat_sigma_0[dir_] = (F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s) / J_s
                    hat_Cbar_0[dir_] = (F_alpha_ntau_s ** 2) * hat_dSdC_alpha * (F_alpha_ntau_s ** 2) / J_s

                    # Simpson's rule integration of stress & stiffness
                    sigma[dir_] += ((mq_2 * hat_sigma_2[dir_] +
                                     4.0 * mq_1 * hat_sigma_1[dir_] +
                                     mq_0 * hat_sigma_0[dir_]) /
                                    curr_vessel.rhoR_h) * (dt / 3.0)

                    Cbar[dir_] += ((mq_2 * hat_Cbar_2[dir_] +
                                    4.0 * mq_1 * hat_Cbar_1[dir_] +
                                    mq_0 * hat_Cbar_0[dir_]) /
                                   curr_vessel.rhoR_h) * (dt / 3.0)

                # (i) Store active variables for next iteration
                if curr_vessel.alpha_active[alpha] == 1:
                    a_act += k_act * (q_act_2 * a_2 + 4.0 * q_act_1 * a_1 + q_act_0 * a_0) * (dt / 3.0)
                    a_2 = a_0
                    q_act_2 = q_act_0

                # (j) Store intermediate kinetics for next iteration
                k_2 = k_0
                q_2 = q_0
                mq_2 = mq_0

                # (k) Copy the final “hat_sigma_0” into “hat_sigma_2” for next iteration
                for dir_ in range(3):
                    hat_sigma_2[dir_] = hat_sigma_0[dir_]
                    hat_Cbar_2[dir_] = hat_Cbar_0[dir_]

            # -----------------------------------------------
            # 6.3. If we have an even number of points, do a trapezoid step
            # -----------------------------------------------
            if even_n:
                a_t = curr_vessel.a[taun_min]
                h_t = curr_vessel.h[taun_min]
                lambda_th_tau = (a_t + h_t / 2.0) / (a0 + h0 / 2.0)
                lambda_z_tau = curr_vessel.lambda_z_tau[taun_min]
                J_tau = curr_vessel.rhoR[taun_min] / curr_vessel.rhoR_h

                F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau)
                F_tau[1] = lambda_th_tau
                F_tau[2] = lambda_z_tau

                k_0 = curr_vessel.k_alpha[nts * alpha + taun_min]
                q_0 = math.exp(-(k_2 + k_0) * dt / 2.0) * q_2
                mq_0 = curr_vessel.mR_alpha[nts * alpha + taun_min] * q_0

                if curr_vessel.alpha_active[alpha] == 1:
                    a_0 = a_t
                    q_act_0 = math.exp(-k_act * dt) * q_act_2


                # Compute stress/stiffness from this "taun_min" state
                for dir_ in range(3):
                    [hat_S_alpha, hat_dSdC_alpha] = constitutive(curr_vessel,
                                                                 lambda_alpha_s[alpha],
                                                                 alpha,
                                                                 taun_min,
                                                                 dir_)
                    F_alpha_ntau_s = (F_s[dir_] / F_tau[dir_]) * curr_vessel.G_alpha_h[3 * alpha + dir_]

                    hat_sigma_0[dir_] = (F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s) / J_s
                    hat_Cbar_0[dir_] = (F_alpha_ntau_s ** 2) * hat_dSdC_alpha * (F_alpha_ntau_s ** 2) / J_s

                    # Trapezoid integration
                    sigma[dir_] += ((mq_2 * hat_sigma_2[dir_] +
                                     mq_0 * hat_sigma_0[dir_]) /
                                    curr_vessel.rhoR_h) * (dt / 2.0)
                    Cbar[dir_] += ((mq_2 * hat_Cbar_2[dir_] +
                                    mq_0 * hat_Cbar_0[dir_]) /
                                   curr_vessel.rhoR_h) * (dt / 2.0)

                if curr_vessel.alpha_active[alpha] == 1:
                    a_act += k_act * (q_act_2 * a_2 + q_act_0 * a_0) * (dt / 2.0)

            # -----------------------------------------------
            # 6.4. Include the initial material if taun_min==0
            # -----------------------------------------------
            if taun_min == 0:
                for dir_ in range(3):
                    sigma[dir_] += (curr_vessel.rhoR_alpha[nts * alpha + 0] /
                                    curr_vessel.rhoR_h) * q_0 * hat_sigma_0[dir_]
                    Cbar[dir_] += (curr_vessel.rhoR_alpha[nts * alpha + 0] /
                                   curr_vessel.rhoR_h) * q_0 * hat_Cbar_0[dir_]

        else:
            # ---------------------------------------------
            # 6.5. If sn==0 or if the constituent doesn't remodel
            # ---------------------------------------------
            for dir_ in range(3):
                [hat_S_alpha, hat_dSdC_alpha] = constitutive(curr_vessel,
                                                             lambda_alpha_s[alpha],
                                                             alpha,
                                                             0,  # reference time
                                                             dir_)
                F_alpha_ntau_s = F_s[dir_] * curr_vessel.G_alpha_h[3 * alpha + dir_]
                hat_sigma_2[dir_] = (F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s) / J_s
                hat_Cbar_2[dir_] = (F_alpha_ntau_s ** 2) * hat_dSdC_alpha * (F_alpha_ntau_s ** 2) / J_s

                sigma[dir_] += (curr_vessel.rhoR_alpha[nts * alpha + sn] /
                                curr_vessel.rhoR_h) * hat_sigma_2[dir_]

                Cbar[dir_] += (curr_vessel.rhoR_alpha[nts * alpha + sn] /
                               curr_vessel.rhoR_h) * hat_Cbar_2[dir_]

        # (Optional) if taun_min==0 and the constituent is active,
        # add the initial active radius
        if taun_min == 0 and curr_vessel.alpha_active[alpha] == 1:
            # q_act_0 was never changed in this branch, so if you strictly need it,
            # you might define it as 1.0 by default. We'll assume 1.0 if not set.
            a_act += curr_vessel.a_act[0] * q_act_0  # times q_act_0 if needed

    # ----------------------------------------------------
    # 7. Compute the final active stress contribution
    # ----------------------------------------------------
    # If sn==0, we just set a_act to the initial a_act[0]
    if sn == 0:
        a_act = curr_vessel.a_act[0]

    # C is the "vasomotor" control variable:
    C = (curr_vessel.CB
         - curr_vessel.CS * (curr_vessel.bar_tauw / curr_vessel.bar_tauw_h - 1.0))

    lambda_act = curr_vessel.a[sn] / a_act if a_act != 0.0 else 0.0

    # Quadratic "parab_act" factor:  1 - [ (lambda_m - lambda_act)/(lambda_m-lambda_0) ]^2
    # clipped if necessary
    lm = curr_vessel.lambda_m
    l0 = curr_vessel.lambda_0
    top = (lm - lambda_act)
    bot = (lm - l0) if (lm != l0) else 1.0
    parab_act = 1.0 - (top / bot) ** 2

    # Hat-sigma-active
    T_act = curr_vessel.T_act
    hat_sigma_act = T_act * (1.0 - math.exp(-(C ** 2))) * lambda_act * parab_act

    # We assume the 'active' alpha is at index=1 based on the original code
    sigma_act = (curr_vessel.rhoR_alpha[nts * 1 + sn] / (J_s * curr_vessel.rhoR_h)
                 * hat_sigma_act)

    # ----------------------------------------------------
    # 8. The Lagrange multiplier is the radial stress (sigma[0])
    #    Subtract it from each direction. Then add active stress
    #    in the circumferential direction (dir == 1).
    # ----------------------------------------------------
    lagrange = sigma[0]
    for dir_ in range(3):
        # The original code: Cbar[dir] = 2 * sigma[dir] + 2 * Cbar[dir];
        # We'll replicate that:
        Cbar[dir_] = 2.0 * sigma[dir_] + 2.0 * Cbar[dir_]

        # Subtract radial stress from each direction
        sigma[dir_] = sigma[dir_] - lagrange

        # Add active stress in the circumferential direction (dir==1)
        if dir_ == 1:
            sigma[dir_] += sigma_act

        # Update the vessel fields
        curr_vessel.sigma[dir_] = sigma[dir_]
        curr_vessel.Cbar[dir_] = Cbar[dir_]

    # ----------------------------------------------------
    # 9. Finally, save the updated active radius
    # ----------------------------------------------------
    curr_vessel.a_act[sn] = a_act

