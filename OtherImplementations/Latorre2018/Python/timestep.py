from kinetics import update_kinetics
from geometry import find_iv_geom


def update_time_step_py(curr_vessel):
    """
    Python equivalent of:
        void update_time_step(vessel &curr_vessel)

    1. Sets initial guesses for radius and mass densities based on previous time step.
    2. Calls update_kinetics and find_iv_geom to solve for the new in vivo (loaded) geometry.
    3. Iterates until the total referential mass density converges within tolerance or hits
       a maximum iteration count.
    4. Computes the mechano-biological equilibrium state (mb_equil).
    5. Prints the current time, radius, thickness, and equilibrium state.
    """

    # 1) Basic indexing and parameters
    n_alpha = curr_vessel.n_alpha
    nts     = curr_vessel.nts
    sn      = curr_vessel.sn
    s       = sn * curr_vessel.dt   # current time in "real" units
    tol     = 1.0e-14              # convergence tolerance
    iter_   = 0
    max_iter= 100

    # 2) Initialize mass density for each constituent at this time step
    #    (copy from previous time step)
    for alpha in range(n_alpha):
        curr_vessel.rhoR_alpha[nts * alpha + sn] = curr_vessel.rhoR_alpha[nts * alpha + sn - 1]

    # 3) Initialize guess for geometry from previous time step
    curr_vessel.a_mid[sn] = curr_vessel.a_mid[sn - 1]
    curr_vessel.a_act[sn] = curr_vessel.a_act[sn - 1]

    # 4) Update kinetics with these guesses and solve for the in vivo geometry
    update_kinetics(curr_vessel)
    find_iv_geom(curr_vessel)  # previously translated function

    # 5) Now iteratively refine the mass production until total mass density converges
    while True:
        iter_ += 1

        # (a) Store the old total referential mass density
        rhoR_s0 = curr_vessel.rhoR[sn]

        # (b) Update the kinetic variables again
        update_kinetics(curr_vessel)

        # (c) Re-solve for in vivo geometry
        find_iv_geom(curr_vessel)

        # (d) Compute new total referential mass density
        rhoR_s1 = curr_vessel.rhoR[sn]

        # (e) Check relative convergence in total mass density
        mass_check = abs((rhoR_s1 - rhoR_s0) / rhoR_s0 if rhoR_s0 != 0 else 0)

        if (mass_check <= tol) or (iter_ >= max_iter):
            break

    # 6) Compute the current mechano-biological equilibrium state (mb_equil)
    #    mb_equil = 1 + K_sigma_p_alpha_h[2] * ( sigma_rel - 1 ) - K_tauw_p_alpha_h[2]* ( tauw_rel - 1 )
    sigma_rel = (curr_vessel.sigma[1] + curr_vessel.sigma[2]) / (curr_vessel.sigma_h[1] + curr_vessel.sigma_h[2])
    tauw_rel  = curr_vessel.bar_tauw / curr_vessel.bar_tauw_h

    curr_vessel.mb_equil = (
        1.0
        + curr_vessel.K_sigma_p_alpha_h[2] * (sigma_rel - 1.0)
        - curr_vessel.K_tauw_p_alpha_h[2]  * (tauw_rel - 1.0)
    )

    # 7) Print the current state
    #    e.g.: "Time: {s} a: {curr_vessel.a[sn]} h: {curr_vessel.h[sn]} Equil: {curr_vessel.mb_equil}"
    print(f"Time: {s:.5g}  a: {curr_vessel.a[sn]:.5g}  h: {curr_vessel.h[sn]:.5g}  Equil: {curr_vessel.mb_equil:.5g}")

