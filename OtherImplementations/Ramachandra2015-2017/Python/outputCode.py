
import numpy as np
import math

def GnRadaptationvein(percent_pressure, percent_flow, filename):
    """
    Main program for step change simulations
    -------------------------------------------------------------
    Units: [cm, g, sec]
    -------------------------------------------------------------
    Function takes as input - 
    percent change in pressure - float/double/integer
    percent change inflow -float/double/integer
    name of output file -string
    It writes to disk a file with the time(days), pressure (dyn/cm^2), flow (ml/s),
    radius(cm), thicness(cm)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    example input - for a 10% increase in flow and 50% increase in pressure
    >> GnRadaptationvein(0.5,0.1,'foo_simulation_output')
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    <^>
    -----
    (~ ~)
    /(` 
      )
    Code to simulate basilar artery adaptations for a step change in load using following constituents:
      --4 fiber collagen family (axial, circumferential, +/- 45deg),
      --circumferential smooth muscle
      -- elastin
    ref paper - Valentin, A., et al. Journal of The Royal Society Interface 6.32 (2009): 293-306.
    questions? - email - amarsden@stanford.edu
    for class purposes only
    """
    
    comp_t = 1500
    delta_P = 1 + percent_pressure  # convert percentage to fold increase
    delta_Q = 1 + percent_flow      # convert percentage to fold increase

    # Additional flags
    is_collagen_angle_const = False
    loading = "step"  # 'step' or 'gradual'

    K_c1 = 0.010000        # rate parameter - collagen, intramural
    K_c2 = 8.006000        # rate parameter - collagen, shear
    K_m1 = 2.009000        # rate parameter - smooth muscle, intramural
    K_m2 = 6.007000        # rate parameter - smooth muscle, shear
    C_ratio = 20
    G_ch = 1.039000
    G_mh = 1.310000
    G_et = 1.75
    G_ez = 1.75
    kq_m = 0.901000
    kq_c = 0.010000
    T_M = 109000.000000
    
    # appended = strcat(filename, str_a, str_b, str_c, str_d, str_e, str_f);
    appended = filename
    
    # Open the main output file
    fid = open(appended, 'w')

    # Days to computer pressure-diameter curves
    # format short e (Python doesn't need this)

    # Transmural pressure
    P_h = 5.0 * 133.32237 * 10.0  # (dyn/cm s^2)

    # Vessel geometry
    a_M = 0.2           # (cm)
    # L_M = 10              # (cm)

    # Viscosity and density of blood
    mu = 0.037             # (g/cm.s)
    # rho_f = 1.050         # (g/cm^3)

    # K_f = 0.5

    # Homeostatic wall shear stress
    # tau_wh = 50.6;%7.426;         # (g/cm.s^2) or (dynes/cm^2) or (0.1 Pa)
    tau_wh = 6
    
    # Calculate Homeostatic flow rate and blood velocity
    Q_Mh = tau_wh * math.pi * a_M**3 / (4.0 * mu)    # (mL/s)
    # v_Mh = Q_Mh / (pi* a_M^2);              #(cm/s)

    # -----------------------------------------------
    # Solid material parameters
    # -----------------------------------------------
    nfiberfly = 4
    
    # Mass density of solid constituents
    rho_s = 1.050  # (g/cm^3)

    # Fraction of collagen fibers (axial, cir, helical, helical)
    c_frac = np.array([0.1, 0.1, 0.4, 0.4])

    # Alignment of collagen fibers
    alpha_ckh = np.array([0.0, 90.0, 45.0, 135.0]) * (math.pi / 180.0)

    # Homeostatic circumferential stress
    # sigma_h = 100.0*1000.0*10.0;      #N/cm^2
    sigma_h = 18.0 * 1000.0 * 10.0      # N/cm^2

    # Calculate homeostatic thickness
    h_h = P_h * a_M / sigma_h

    # Mass fractions
    phi_c = 0.42  # was 0.22
    phi_ck = phi_c * c_frac
    phi_e = 0.1   # was 0.02
    phi_m = 1.0 - phi_c - phi_e

    # Homeostatic masses
    M_h = h_h * rho_s
    M_ckh = M_h * phi_ck
    # M_ch = sum(M_ckh)
    M_mh = phi_m * M_h
    M_eh = phi_e * M_h

    # sigma_ch = 100.0*1000.0*10.0;
    # sigma_mh = 100.0*1000.0*10.0;  # pas + act

    sigma_ch = 18.0 * 1000.0 * 10.0
    sigma_mh = 18.0 * 1000.0 * 10.0  # pas + act
    sigma_eh = (1.0 / phi_e) * (sigma_h - phi_c * sigma_ch - phi_m * sigma_mh)

    # Muscle activation parameters
    Lambda_M = 1.35
    Lambda_0 = 0.50
    C_basal = 0.68
    T_S0 = C_ratio * C_basal
    #      T_M = 58.0*1000.0*10.0;
    # Passive response parameters and free parameter calculation

    # Calculate remaining material parameters (g/cm.s^2)
    c_m3 = 0.05025
    c_c2 = np.zeros(4)
    c_c2[0] = 24022.65 / phi_ck[0]
    c_c2[1] = 1000 / phi_ck[1]
    c_c2[2] = 1.0
    c_c2[3] = 1.0
    
    c_c3 = np.zeros(4)
    c_c3[0] = 0.1
    c_c3[1] = 0.05025
    c_c3[2] = 1.035
    c_c3[3] = 1.035

    c_e = sigma_eh / ((G_et**2.0) - (((G_et**2.0) * (G_ez**2.0))**(-1)))
    
    #      sigma_m_pas_h=c_m2*((G_mh**2.0)*(G_mh**(2.0)-1.0))*
    #     2 EXP(c_m3*((G_mh**(2.0)-1.0)**2.0))  
    #      sigma_act_h=sigma_mh-sigma_m_pas_h
    #      temp=1.0-((Lambda_M-1.0)/(Lambda_M-Lambda_0))**2.0
    #      T_M=sigma_act_h/((1.0-EXP(-(C_basal**2.0)))*temp)
    
    sigma_act_h = T_M * (1.0 - math.exp(-C_basal**2.0)) * (1.0 - ((Lambda_M - 1.0) / (Lambda_M - Lambda_0))**2.0)
    print(sigma_act_h)
    
    if sigma_act_h < 0.0:
        sigma_act_h = 0.0  # only tension is possible

    temp = c_m3 * ((G_mh**2.0 - 1.0)**2.0)
    c_m2 = (sigma_mh - sigma_act_h) / ((G_mh**2.0) * (G_mh**2.0 - 1.0) * math.exp(temp))
    print(c_m2)
    
    if c_m2 < 0.0:
        print('c_m2 is negative')
        exit()

    str_ch = np.zeros(4)
    for i in range(nfiberfly):
        str_ch[i] = c_frac[i] * c_c2[i] * (G_ch**2.0) * ((G_ch**2.0) - 1.0) * math.exp(c_c3[i] * ((G_ch**2.0) - 1.0)**2.0) * (math.sin(alpha_ckh[i]))**2.0

    c_c2[3] = (sigma_ch - str_ch[0] - str_ch[1]) / (2.0 * str_ch[3])
    c_c2[2] = c_c2[3]
    c_c2_c = c_m2
    c_c3_c = c_m3
    c_m2_c = c_m2
    c_m3_c = c_m3

    # checking whether the equlibrium equation is satisfied with material
    # parameters

    dwdLt_c = 0.0
    for i in range(4):
        dwdLt_c = dwdLt_c + phi_ck[i] * h_h * dWkdLn(G_ch, c_c2[i], c_c3[i], c_c2_c, c_c3_c) * G_ch * (math.sin(alpha_ckh[i]))**2

    dwdLt_m = h_h * phi_m * dWmdLn(G_mh, c_m2, c_m3, c_m2_c, c_m3_c) * G_mh
    dwdLt_e = h_h * phi_e * dWedLn(G_et, G_ez, 1, c_e) * G_et
    dwdLt = dwdLt_c + dwdLt_m + dwdLt_e

    T_act = sigma_act_h * phi_m * h_h
    F = dwdLt + T_act - P_h * a_M  # this should be zero at equilibrium

    # sprintf('%0.5f, %0.5f, %0.5f, %0.5f\n', ...
    #         F, c_c2(1), c_m2, c_e)

    # --------------------------------------------------
    #   Find a stretch to satisfy:
    #   dWkdLn(Ln)/dWkdLn(G_ck) = yielding criteria (about 10)
    #   dWmdLn(Ln)/dWndLn(G_mh) = yielding criteria
    #    (Newton-Raptson method)
    # --------------------------------------------------

    y_Lkn = 3.0
    y_Lmn = 3.0

    # --------------------------------------------------
    #  Other kinetic parameters (stress-mediation, etc)
    # --------------------------------------------------
    # Degredation rates (day^-1)

    # kq_m = 1.0/80.0;  # for muscle
    # kq_c = 1.0/80.0;  # for collagen

    # Rate parameter for changing in the ref config of active VSM tone
    k_act = 1.0 / 20.0  # (day^-1)

    # Algorithm parameters
    Max_error = 2.0 * 10**(-4)
    Max_error_2 = 1.0 * 10**(-6)
    Max_it = 200

    # Number of time step per day
    n_t_step = 10.0
    dt = 1 / n_t_step  # (day)

    # Maximum life span of collagen and smooth muscle (days)
    age_max = 1000.0

    # number of data per life span
    num_DL = int(n_t_step * age_max + 1)
    num_t = int(comp_t * n_t_step + 1)

    DQ_m = np.zeros(num_DL)
    DQ_c = np.zeros(num_DL)
    Da = np.zeros(num_t)
    Dm_c = np.zeros((num_t, 4))     # rate of mass production of collagen
    Dm_m = np.zeros(num_t)
    Dalpha = np.zeros(num_t)        # fiber angle of newly produced collagen fibers

    DQ2_c = np.zeros((num_t, 4))    # tension dependent degradation
    DQ2_m = np.zeros(num_t)

    # Survival fractions Q(t)
    DQ_m[0] = 0.0
    DQ_c[0] = 0.0

    for i in range(1, num_DL):
        t = dt * i
        DQ_m[i] = DQ_m[i - 1] + 0.5 * dt * (q_i(t - dt, kq_m) + q_i(t, kq_m))
        DQ_c[i] = DQ_c[i - 1] + 0.5 * dt * (q_i(t - dt, kq_c) + q_i(t, kq_c))

    mean_age_m = DQ_m[num_DL - 1]
    mean_age_c = DQ_c[num_DL - 1]

    for i in range(num_DL):
        DQ_m[i] = (mean_age_m - DQ_m[i]) / mean_age_m
        DQ_c[i] = (mean_age_c - DQ_c[i]) / mean_age_c

    m_basal_ck = M_ckh / mean_age_c
    m_basal_m = M_mh / mean_age_m

    # Data for nn=1 at t=0
    Da[0] = a_M
    # a_act = a_M    #reference radius for active SMT
    # da_act_p = 0 #d a_act/dt
    a_act_p = a_M
    da_act_p = 0.0

    # Rate of mass production of collagen (4 families) and muscle
    Dm_c[0, :] = m_basal_ck
    Dm_m[0] = m_basal_m

    # Fiber angle of newly produced collagen fibers
    Dalpha[0] = alpha_ckh[2]
    # alpha_ck = alpha_ckh

    for i in range(num_t):
        DQ2_c[i, :] = np.array([1.0, 1.0, 1.0, 1.0])
        DQ2_m[i] = 1.0
    
    sigma_ck = np.array([1.0, 1.0, 1.0, 1.0]) * sigma_ch

    # Initialize output arrays for checking
    Sig_ = np.zeros(num_t)
    C_ = np.zeros(num_t)
    Lam_z = np.zeros(num_t)
    Lam_t = np.zeros(num_t)
    Masses_ = np.zeros(num_t)

    # --------------------------------------------------
    # Begin marching through time
    # --------------------------------------------------
    for n_t in range(1, num_t):
        t = n_t * dt
        # initial guess (predictors)
        beta = 0.3  # an adjusting parameter for a Newton-Rapson method

        if loading == "step":
            if t > 2:
                P = P_h * delta_P
                Q_M = Q_Mh * delta_Q
            else:
                P = P_h
                Q_M = Q_Mh
        elif loading == "gradual":
            if t <= 2:
                P = P_h
                Q_M = Q_Mh
            elif t > 2 and t < 32:
                p_slope = (P_h * delta_P - P_h) / (32 - 2)
                # Define the linear gradual loading equation
                P = P_h + p_slope * (t - 2)
                Q_M = Q_Mh * delta_Q
            else:
                P = P_h * delta_P
                Q_M = Q_Mh * delta_Q

        if n_t > 2:
            if abs(Da[n_t - 1] - Da[n_t - 2]) / Da[n_t - 1] < 0.1:
                Da[n_t] = 2.0 * Da[n_t - 1] - Da[n_t - 2]
                Dm_c[n_t, :] = 2.0 * Dm_c[n_t - 1, :] - Dm_c[n_t - 2, :]
                Dm_m[n_t] = 2.0 * Dm_m[n_t - 1] - Dm_m[n_t - 2]
            else:
                Da[n_t] = Da[n_t - 1]
                Dm_c[n_t, :] = Dm_c[n_t - 1, :]
                Dm_m[n_t] = Dm_m[n_t - 1]
        else:
            Da[n_t] = Da[n_t - 1]
            Dm_c[n_t, :] = Dm_c[n_t - 1, :]
            Dm_m[n_t] = Dm_m[n_t - 1]

        a_t = Da[n_t]

        # predictor for a_act
        a_act = a_act_p + 0.5 * dt * (da_act_p + k_act * (a_t - (a_act_p + dt * da_act_p)))

        tol_m = 100.0
        tau_w = tau_wh

        num_it1 = 0
        # num_it2 = 0

        while (tol_m > Max_error_2) and (num_it1 < Max_it):
            num_it1 = num_it1 + 1
            tol_a = 100.0
            num_it2 = 0

            while (tol_a > Max_error) and (num_it2 < Max_it):
                num_it2 = num_it2 + 1

                # Caclulate the current wall shear stress
                tau_w = 4.0 * mu * Q_M / (math.pi * a_t**3)

                T_S = T_S0
                C_t = C_basal - T_S * (tau_w / tau_wh - 1.0)
                dn_tau = (tau_w / tau_wh - 1.0)

                if C_t < 0.0:
                    C_t = 0.0

                # 2D stretches
                L_z = 1.0          # Constant axial stretch
                L_t = a_t / a_M    # Changing circumferential stretch # Inner radius instead of mid-point radius!

                # Alignment of a new collagen
                # Note that is defined by stretches NOT stresses
                if is_collagen_angle_const:
                    Dalpha[n_t] = alpha_ckh[2]
                    # disp(Dalpha[n_t])
                else:
                    Dalpha[n_t] = math.atan(L_t / L_z * math.tan(alpha_ckh[2]))
                    # disp(Dalpha[n_t])

                dwdLt_c = 0.0
                dwdLz_c = 0.0
                dwdLt_m = 0.0
                # dwdLz_m = 0.0
                ddwddLt_c = 0.0
                # ddwddLz_c = 0.0
                # ddwdLtdLz_c = 0.0
                ddwddLt_m = 0.0
                # ddwddLz_m = 0.0
                # ddwdLtdLz_m = 0.0

                M_ck = np.zeros(4)
                dWck_dLn = np.zeros(4)
                M_m = 0.0

                if n_t < num_DL:
                    M_ck = M_ckh * DQ_c[n_t] * DQ2_c[0, :]
                    sQ_m = DQ_m[n_t]  # *DQ2_m[0]
                    M_m = M_mh * sQ_m
                    Lc_k = np.sqrt(L_t**2 * np.sin(alpha_ckh)**2 + L_z**2 * np.cos(alpha_ckh)**2)

                    # Due to initial collagen
                    for i in range(4):
                        sQ_ck = DQ_c[n_t] * DQ2_c[0, i]
                        # Lc_kn = Lc_k[i]*G_ch

                        dLndLt = G_ch * L_t * math.sin(alpha_ckh[i])**2 / Lc_k[i]
                        dLndLz = G_ch * L_z * math.cos(alpha_ckh[i])**2 / Lc_k[i]
                        ddLnddLt = G_ch * (L_z * math.sin(alpha_ckh[i]) * math.cos(alpha_ckh[i]))**2 / Lc_k[i]**3
                        # ddLnddLz = G_ch*(L_t*sin(alpha_ckh[i])*cos(alpha_ckh[i]))**2/Lc_k[i]**3
                        # ddLndLtdLz = -G_ch*L_t*L_z*(sin(alpha_ckh[i])*cos(alpha_ckh[i]))**2/Lc_k[i]**3

                        dWck_dLn[i] = dWkdLn(Lc_k[i] * G_ch, c_c2[i], c_c3[i], c_c2_c, c_c3_c)

                        dwdLt_c = dwdLt_c + (M_ckh[i] / rho_s) * sQ_ck * dWck_dLn[i] * dLndLt
                        dwdLz_c = dwdLz_c + (M_ckh[i] / rho_s) * sQ_ck * dWck_dLn[i] * dLndLz
                        ddwddLt_c = ddwddLt_c + (M_ckh[i] / rho_s) * sQ_ck * \
                            (ddWkddLn(Lc_k[i] * G_ch, c_c2[i], c_c3[i], c_c2_c, c_c3_c) * (dLndLt)**2 + \
                            dWck_dLn[i] * ddLnddLt)

                    # Smooth Muscle presented at time zero
                    Lm_n = L_t * G_mh

                    dWm_dLn = dWmdLn(Lm_n, c_m2, c_m3, c_m2_c, c_m3_c)

                    dwdLt_m = dwdLt_m + (M_mh / rho_s) * sQ_m * dWm_dLn * G_mh
                    ddwddLt_m = ddwddLt_m + (M_mh / rho_s) * sQ_m * ddWmddLn(Lm_n, c_m2, c_m3, c_m2_c, c_m3_c) * G_mh**2

                if n_t < num_DL:
                    tn0 = 0
                else:  # calculate only for maximum life span of constituent
                    tn0 = n_t - num_DL

                for n_tau in range(tn0, n_t):
                    if n_tau == tn0 or n_tau == n_t - 1:
                        wt = 0.5 * dt  # for integration
                    else:
                        wt = dt

                    alpha_tau = np.array([alpha_ckh[0], alpha_ckh[1], Dalpha[n_tau], 2.0 * math.pi - Dalpha[n_tau]])
                    Lt_tau = Da[n_tau] / a_M
                    Lz_tau = 1.0
                    Lc_k_tau = np.sqrt(Lt_tau**2 * np.sin(alpha_tau)**2 + Lz_tau**2 * np.cos(alpha_tau)**2)
                    Lc_k = np.sqrt(L_t**2 * np.sin(alpha_tau)**2 + L_z**2 * np.cos(alpha_tau)**2)
                    Lc_kn = G_ch * Lc_k / Lc_k_tau

                    for i in range(4):
                        mc_tau = Dm_c[n_tau, i]
                        # if the stretch is larger than yiedling value,
                        # then it degradated. so it will not included in mass

                        if Lc_kn[i] <= y_Lkn:
                            sq_ck = q_i((n_t - n_tau) * dt, kq_c) * DQ2_c[n_tau, i]
                            M_ck[i] = M_ck[i] + mc_tau * sq_ck * wt
                            dLndLt = (G_ch / Lc_k_tau[i]) * L_t * math.sin(alpha_tau[i])**2 / Lc_k[i]
                            dLndLz = (G_ch / Lc_k_tau[i]) * L_z * math.cos(alpha_tau[i])**2 / Lc_k[i]
                            ddLnddLt = (G_ch / Lc_k_tau[i]) * (L_z * math.sin(alpha_tau[i]) * \
                                math.cos(alpha_tau[i]))**2 / Lc_k[i]**3

                            dWck_dLn[i] = dWkdLn(Lc_kn[i], c_c2[i], c_c3[i], c_c2_c, c_c3_c)

                            dwdLt_c = dwdLt_c + (mc_tau / rho_s) * sq_ck * dWck_dLn[i] * dLndLt * wt
                            dwdLz_c = dwdLz_c + (mc_tau / rho_s) * sq_ck * dWck_dLn[i] * dLndLz * wt
                            ddwddLt_c = ddwddLt_c + (mc_tau / rho_s) * sq_ck * \
                                (ddWkddLn(Lc_kn[i], c_c2[i], c_c3[i], c_c2_c, c_c3_c) * (dLndLt)**2 + \
                                dWck_dLn[i] * ddLnddLt) * wt

                    mm_tau = Dm_m[n_tau]
                    Lm_n = G_mh * L_t / Lt_tau

                    if Lm_n <= y_Lmn:
                        sq_m = q_i((n_t - n_tau) * dt, kq_m) * DQ2_m[n_tau]
                        M_m = M_m + mm_tau * sq_m * wt

                        dWm_dLn = dWmdLn(Lm_n, c_m2, c_m3, c_m2_c, c_m3_c)

                        dwdLt_m = dwdLt_m + (mm_tau / rho_s) * sq_m * dWm_dLn * (G_mh / Lt_tau) * wt
                        ddwddLt_m = ddwddLt_m + (mm_tau / rho_s) * sq_m * \
                            ddWmddLn(Lm_n, c_m2, c_m3, c_m2_c, c_m3_c) * (G_mh / Lt_tau)**2 * wt

                M_c = np.sum(M_ck)

                # strain energy due to elastin
                M_e = M_eh
                dwdLt_e = (M_e / rho_s) * dWedLn(L_t * G_et, L_z * G_ez, 1, c_e) * G_et
                # dwdLz_e = (M_e/rho_s)*dWedLn(L_t*G_et, L_z*G_ez, 2, c_e)*G_ez
                ddwddLt_e = (M_e / rho_s) * ddWeddLn(L_t * G_et, L_z * G_ez, 1, c_e) * G_et**2
                # ddwdLtdLz_e = (M_e/rho_s)*ddWeddLn(L_t*G_et, L_z*G_ez, 3, c_e)*G_et*G_ez

                # strain energy due to active SMT
                # T_B should be given as a function of time
                L_m_act = a_t / a_act
                J = L_t * L_z

                T_act = T_M * M_m / (rho_s * J) * (1.0 - math.exp(-C_t**2)) * L_m_act * \
                    (1.0 - ((Lambda_M - L_m_act) / (Lambda_M - Lambda_0))**2)

                # Active muscle nan only generate tensions
                if T_act < 0.0:
                    T_act = 0.0

                # Recall C_t = C_basal - T_S * (tau_w / tau_wh - 1);
                dC_tda = T_S * 12.0 * mu * Q_M / (tau_wh * math.pi * a_t**4)

                dT_actda = -T_M * M_m * L_z / (rho_s * J**2 * a_t) * (1.0 - math.exp(-C_t**2)) * L_m_act * (1.0 - ((Lambda_M - L_m_act) / (Lambda_M - Lambda_0))**2) + \
                    T_M * M_m / (rho_s * J) * (1.0 - math.exp(-C_t**2)) / a_act * (1.0 - ((Lambda_M - L_m_act) / (Lambda_M - Lambda_0))**2) + \
                    T_M * M_m / (rho_s * J) * (1.0 - math.exp(-C_t**2)) * L_m_act * (Lambda_M - L_m_act) / ((Lambda_M - Lambda_0)**2 * a_act) + \
                    T_M * M_m / (rho_s * J) * 2.0 * C_t * dC_tda * math.exp(-C_t**2) * L_m_act * (Lambda_M - L_m_act) / ((Lambda_M - Lambda_0)**2 * a_act)

                total_M = M_c + M_e + M_m
                dwdLt = dwdLt_c + dwdLt_m + dwdLt_e
                ddwddLt = ddwddLt_c + ddwddLt_m + ddwddLt_e

                F_a = (1.0 / L_z) * dwdLt + T_act - P * a_t
                dF_ada = (1.0 / L_z) * ddwddLt * (1.0 / a_M) + dT_actda - P  # - dPda*a_t

                a_t = a_t - beta * F_a / dF_ada
                tol_a = math.sqrt((a_t - Da[n_t])**2 / Da[n_t]**2)

                Da[n_t] = a_t
                a_act = 1.0 / (1.0 + 0.5 * dt * k_act) * (a_act_p + 0.5 * dt * (da_act_p + k_act * a_t))
                L_t = a_t / a_M
                
                if is_collagen_angle_const:
                    Dalpha[n_t] = alpha_ckh[2]
                else:
                    Dalpha[n_t] = math.atan(L_t / L_z * math.tan(alpha_ckh[2]))

            # adjust beta
            # if NR method don't converge
            if num_it2 == Max_it:
                beta = 0.1

            T_t_c = (rho_s * L_t * L_z / M_c) * (dwdLt_c) / L_z
            T_z_c = (rho_s * L_t * L_z / M_c) * (dwdLz_c) / L_t
            T_t_m = (rho_s * L_t * L_z / M_m) * ((dwdLt_m) / L_z + T_act)

            alpha_ck = np.array([0.0, math.pi / 2.0, Dalpha[n_t], 2.0 * math.pi - Dalpha[n_t]])

            for i in range(4):
                sigma_ck[i] = math.sqrt(T_z_c**2 * math.cos(alpha_ck[i])**2 + T_t_c**2 * math.sin(alpha_ck[i])**2)

            ### outputs to check sethu
            Sig_[n_t] = T_t_m
            C_[n_t] = C_t
            Lam_z[n_t] = L_z
            Lam_t[n_t] = L_t
            Masses_[n_t] = M_m
            
            # n_c = M_c/M_ch
            # n_SMC = M_m/M_mh
            nc_stress = (sigma_ck / sigma_ch - 1.0)
            nm_stress = (T_t_m / sigma_mh - 1.0)
            ncauchy_stress = ((sigma_ck[0] + sigma_ck[1] + sigma_ck[2] + sigma_ck[3] + T_t_m) / (sigma_ch + sigma_mh)) - 1.0

# =============================================================================
#         % Mass production functions % 1's used to be n_c's and n_SMC's
# %          dn_C      = (C_t/C_basal - 1.0);
# %          m_c = [rm_c(1, m_basal_ck(1), nc_stress(1), dn_C, K_c1, K_c2),...
# %              rm_c(1, m_basal_ck(2), nc_stress(2), dn_C, K_c1, K_c2),...
# %              rm_c(1, m_basal_ck(3), nc_stress(3), dn_C, K_c1, K_c2), ...
# %              rm_c(1, m_basal_ck(4), nc_stress(4), dn_C, K_c1, K_c2)];
# %          m_m = rm_m(1, m_basal_m, nm_stress, dn_C, K_m1, K_m2);
#         
# %        m_c = [M_ck(1)/M_ckh(1)*rm_c_dtauversion(1, m_basal_ck(1), nc_stress(1), dn_tau, K_c1, K_c2),...
#  %          M_ck(2)/M_ckh(2)*rm_c_dtauversion(1, m_basal_ck(2), nc_stress(2), dn_tau, K_c1, K_c2),...
#   %         M_ck(3)/M_ckh(3)*rm_c_dtauversion(1, m_basal_ck(3), nc_stress(3), dn_tau, K_c1, K_c2), ...
#    %        M_ck(4)/M_ckh(4)*rm_c_dtauversion(1, m_basal_ck(4), nc_stress(4), dn_tau, K_c1, K_c2)];
#     %    m_m = M_m/M_mh*rm_m_dtauversion(1, m_basal_m, nm_stress, dn_tau, K_m1, K_m2);
#     
# =============================================================================
            # Mass production functions # 1's used to be n_c's and n_SMC's
            m_c = np.array([
                M_ck[0] / M_ckh[0] * rm_c_dtauversion(1, m_basal_ck[0], nc_stress[0], dn_tau, K_c1, K_c2),
                M_ck[1] / M_ckh[1] * rm_c_dtauversion(1, m_basal_ck[1], nc_stress[1], dn_tau, K_c1, K_c2),
                M_ck[2] / M_ckh[2] * rm_c_dtauversion(1, m_basal_ck[2], nc_stress[2], dn_tau, K_c1, K_c2),
                M_ck[3] / M_ckh[3] * rm_c_dtauversion(1, m_basal_ck[3], nc_stress[3], dn_tau, K_c1, K_c2)
            ])
            m_m = M_m / M_mh * rm_m_dtauversion(1, m_basal_m, nm_stress, dn_tau, K_m1, K_m2)
            
            # Make sure we don't have negative production rates
            for i in range(4):
                if m_c[i] < 0.0:
                    m_c[i] = 0.0

            if m_m < 0.0:
                m_m = 0.0

            error = np.sum((m_c - Dm_c[n_t, :])**2) + (m_m - Dm_m[n_t])**2
            tol_m = math.sqrt(error / (np.sum(m_c**2) + m_m**2))
            Dm_c[n_t, :] = m_c
            Dm_m[n_t] = m_m

        a_act_p = a_act
        da_act_p = k_act * (a_t - a_act)

        h = total_M / (rho_s * L_t * L_z)

        # dWcdLn_h = dWkdLn(G_ch, c_c2, c_c3, c_c2_c, c_c3_c)
        # dWmdLn_h = dWmdLn(G_mh, c_m2, c_m3, c_m2_c, c_m3_c)

        # Output some data to console to monitor progress
        # sprintf('%0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f\n', ...
        #        t,     a_t/a_M, nc_stress[0],  nc_stress[1],   nc_stress[2],   nc_stress[3])

        # sprintf('%0.5f, %0.5f, %0.5f, %0.5f, %0.5f\n', ...
        #         C_t, C_basal, T_S, tau_w, tau_wh)

        # Write data to file
        # 1-5
        fid.write(f'{t:10.8f}\t   {P:10.8f}\t    {Q_M:10.8f}\t                {a_t:10.8f}\t                {h:10.8f}\n')

        # update DQ2_c and DQ2_m
        # tension dependent degradation
        L_t = a_t / a_M

        if n_t < num_DL:
            tn0 = 0
        else:
            tn0 = n_t - num_DL

        dWcdLn_h = np.zeros(4)
        for i in range(4):
            dWcdLn_h[i] = dWkdLn(G_ch, c_c2[i], c_c3[i], c_c2_c, c_c3_c)
        
        dWmdLn_h = dWmdLn(G_mh, c_m2, c_m3, c_m2_c, c_m3_c)

        for n_tau in range(tn0, n_t):
            alpha_tau = np.array([alpha_ckh[0], alpha_ckh[1], Dalpha[n_tau], 2.0 * math.pi - Dalpha[n_tau]])
            Lt_tau = Da[n_tau] / a_M

            Lz_tau = 1.0

            Lc_k_tau = np.sqrt(Lt_tau**2 * np.sin(alpha_tau)**2 + Lz_tau**2 * np.cos(alpha_tau)**2)
            Lc_k = np.sqrt(L_t**2 * np.sin(alpha_tau)**2 + L_z**2 * np.cos(alpha_tau)**2)
            Lc_kn = G_ch * (Lc_k / Lc_k_tau)

            for i in range(4):
                # Tension dependent degradation rate
                beta_c = dWkdLn(Lc_kn[i], c_c2[i], c_c3[i], c_c2_c, c_c3_c) / dWcdLn_h[i]
                add_remv = math.exp(-1.0 * f_beta(beta_c, kq_c) * dt)
                DQ2_c[n_tau, i] = add_remv  # *DQ2_c[n_tau, i]

            Lm_n = G_mh * (L_t / Lt_tau)
            beta_m = dWmdLn(Lm_n, c_m2, c_m3, c_m2_c, c_m3_c) / dWmdLn_h
            add_remv = math.exp(-1.0 * f_beta(beta_m, kq_m) * dt)
            DQ2_m[n_tau] = add_remv  # *DQ2_m[n_tau]

        if math.isnan(h):
            break

    fid.close()
    rad = a_t
    thick = h
    stress_mat = np.array([tau_w / tau_wh - 1, nc_stress[0], nc_stress[1], nc_stress[2], nc_stress[3], nm_stress])

    return rad, thick

# Helper functions (implementations based on context and typical usage)

def q_i(t, kq):
    """Survival function"""
    return math.exp(-kq * t)

def dWkdLn(Ln, c_c2, c_c3, c_c2_c, c_c3_c):
    if Ln >= 1:
        y = c_c2 * Ln * (Ln**2 - 1) * np.exp(c_c3 * (Ln**2 - 1)**2)
    else:
        y = c_c2_c * Ln * (Ln**2 - 1) * np.exp(c_c3_c * (Ln**2 - 1)**2)
    return y

def dWmdLn(Ln, c_m2, c_m3, c_m2_c, c_m3_c):
    if Ln >= 1:
        y = c_m2 * Ln * (Ln**2 - 1) * np.exp(c_m3 * (Ln**2 - 1)**2)
    else:
        y = c_m2_c * Ln * (Ln**2 - 1.0) * np.exp(c_m3_c * (Ln**2 - 1.0)**2)
    return y

def dWedLn(Ln_t, Ln_z, axis, c_e):
    if axis == 1:  # circumferential
        y = c_e * (Ln_t - 1 / (Ln_t**3 * Ln_z**2))
    elif axis == 2:  # axial
        y = c_e * (Ln_z - 1 / (Ln_z**3 * Ln_t**2))
    else:
        raise ValueError('wrong parameter in dWedLn')
    return y

def ddWmddLn(Ln, c_m2, c_m3, c_m2_c, c_m3_c):
    if Ln >= 1:
        y = c_m2 * (3 * Ln**2 - 1 + 4 * c_m3 * (Ln * (Ln**2 - 1))**2) * np.exp(c_m3 * (Ln**2 - 1)**2)
    else:
        exp_Q = np.exp(c_m3_c * (Ln**2 - 1.0)**2)
        y = c_m2_c * (3 * Ln**2 - 1.0 + 4 * c_m3_c * (Ln * (Ln**2 - 1.0))**2) * exp_Q
    return y

def ddWkddLn(Ln, c_c2, c_c3, c_c2_c, c_c3_c):
    if Ln >= 1:
        y = c_c2 * (3 * Ln**2 - 1 + 4 * c_c3 * (Ln * (Ln**2 - 1))**2) * np.exp(c_c3 * (Ln**2 - 1)**2)
    else:
        exp_Q = np.exp(c_c3_c * (Ln**2 - 1.0)**2)
        y = c_c2_c * (3 * Ln**2 - 1.0 + 4 * c_c3_c * (Ln * (Ln**2 - 1.0))**2) * exp_Q
    return y

def ddWeddLn(Ln_t, Ln_z, f, c_e):
    if f == 1:
        y = c_e * (1 + 3 / (Ln_t**4 * Ln_z**2))
    elif f == 2:
        y = c_e * (1 + 3 / (Ln_z**4 * Ln_t**2))
    elif f == 3:
        y = c_e * 2 / (Ln_t * Ln_z)**3
    else:
        raise ValueError('Wrong parameter for ddWeddLn')
    return y

def f_beta(beta, kq):
# =============================================================================
#     gamma = 0.1
#     y_f = 5.0  # was 5.0
#     
#     if beta < y_f:
#         y = 0.0
#     else:
#         y = 9 * kq * (1 - math.exp(-gamma * (beta - y_f) ** 2))
# =============================================================================
    
    # The last line in the natural language prompt seems to be an alternative calculation
    # Uncomment the following line if you want to use that instead
    y = abs(beta - 1) * kq
    
    return y

def rm_m_dtauversion(n_SMC, m_basal, dn_stress, dn_tau, Zeta_m, Km2):
    y = n_SMC * m_basal * (Zeta_m * dn_stress - Km2 * dn_tau + 1)
    return y

def rm_c_dtauversion(n_c, m_basal, dn_stress, dn_tau, Zeta_c, Kc2):
    y = n_c * m_basal * (Zeta_c * dn_stress - Kc2 * dn_tau + 1)
    return y

