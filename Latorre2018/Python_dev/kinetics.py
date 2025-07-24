import math
import numpy as np
from vessel import Vessel
import matplotlib.pyplot as plt


def apply_perturbation(current_time, perturbation_time, homeostatic_value, perturbation_type="step", perturbation_percentage=None, duration=None):
    """
    Applies a perturbation to a given quantity based on the specified parameters and perturbation type.

    Parameters:
    - current_time (float): The current time in the simulation.
    - perturbation_time (float): The time at which the perturbation is applied.
    - homeostatic_value (float): The homeostatic value of the quantity.
    - perturbation_type (str): The type of perturbation ("step", "latorre2018", or "linear").
    - perturbation_percentage (float, optional): The magnitude of the perturbation as a percentage of the homeostatic value.
                                       Positive for increase, negative for decrease.
    - duration (float, optional): The duration over which the linear perturbation is applied.

    Returns:
    - float: The perturbed value of the quantity.
    """
    
    if current_time < perturbation_time:
        return homeostatic_value
    else:
        if perturbation_type == "step":
            if perturbation_percentage is None:
                raise ValueError("Perturbation percentage must be specified for linear perturbation.")
            perturbation_magnitude = homeostatic_value * (perturbation_percentage / 100)
            return homeostatic_value + perturbation_magnitude
   
        elif perturbation_type == "latorre2018":
                factor = (1.0 - math.exp(-(current_time - perturbation_time) / 10.0))
                return homeostatic_value * (1.0 + 0.5 * factor)
  
        elif perturbation_type == "linear":
            if perturbation_percentage is None:
                raise ValueError("Perturbation percentage must be specified for linear perturbation.")
            if duration is None:
                raise ValueError("Duration must be specified for linear perturbation.")
            elapsed_time = current_time - perturbation_time
            perturbation_magnitude = homeostatic_value * (perturbation_percentage / 100)
            perturbed_value = homeostatic_value + perturbation_magnitude
            p_slope = (perturbed_value - homeostatic_value) / duration
            if elapsed_time < duration:
                return homeostatic_value + p_slope * (elapsed_time)
            else:
                return perturbed_value
            
        else:
            raise ValueError("Invalid perturbation type. Choose 'step', 'latorre2018', or 'linear'.")


def update_kinetics(curr_vessel):
    """
    Updates the G&R kinetics for the given vessel instance.

    Parameters
    ----------
    curr_vessel : Vessel
        An instance of the Vessel class that has all the required
        data arrays (rhoR_alpha, k_alpha, mR_alpha, etc.) and attributes
        (n_alpha, nts, dt, etc.).
    """

    # Unpack needed variables
    n_alpha = curr_vessel.n_alpha
    nts = curr_vessel.nts
    dt = curr_vessel.dt
    s = curr_vessel.s
    sn = curr_vessel.sn

    # Determine the earliest time index (taun_min) for integration
    # (using a maximum 'time of interest' tau_max = 10 half-lives)
    tau_max = 10.0 * (1.0 / curr_vessel.k_alpha_h[1])
    if s > tau_max:
        taun_min = sn - int(tau_max / dt)
    else:
        taun_min = 0

    # Compute number of integration points for Simpson's rule
    n = (sn - taun_min) + 1
    even_n = (n % 2 == 0)

    # Calculate mechanical-state differences from reference
    delta_sigma = (
        (curr_vessel.sigma[1] + curr_vessel.sigma[2]) /
        (curr_vessel.sigma_h[1] + curr_vessel.sigma_h[2])
        - 1.0
    )
    delta_tauw = (curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1.0

    # Initialize total current referential density
    rhoR_s = 0.0

    # Loop over each constituent
    for alpha in range(n_alpha):
        deg_check = (curr_vessel.k_alpha_h[alpha] > 0)

        if sn > 0 and deg_check:
            # Unpack gains for current constituent
            K_sigma_p = curr_vessel.K_sigma_p_alpha_h[alpha]
            K_tauw_p = curr_vessel.K_tauw_p_alpha_h[alpha]
            K_sigma_d = curr_vessel.K_sigma_d_alpha_h[alpha]
            K_tauw_d = curr_vessel.K_tauw_d_alpha_h[alpha]

            # Update stimulus functions
            upsilon_mech_p = 1.0 + K_sigma_p * delta_sigma - K_tauw_p * delta_tauw
            upsilon_mech_d = 1.0 + K_sigma_d * (delta_sigma ** 2) + K_tauw_d * (delta_tauw ** 2) 

            k_alpha_s = 0
            mR_alpha_s = 0
            rhoR_alpha_s = 0

            # Current degradation rate and mass-production rate for alpha
            k_alpha_s = curr_vessel.k_alpha_h[alpha] * upsilon_mech_d # TODO: looks inconsistent with Latorre 2018, eq.42
            mR_alpha_s = (k_alpha_s *
                          upsilon_mech_p *
                          curr_vessel.rhoR_alpha[nts * alpha + sn])

            # Store these updated kinetic values
            curr_vessel.k_alpha[nts * alpha + sn] = k_alpha_s
            curr_vessel.mR_alpha[nts * alpha + sn] = mR_alpha_s

            # Setup for Simpson integration
            k_2 = k_alpha_s
            q_2 = 1.0
            mq_2 = mR_alpha_s * q_2
            rhoR_alpha_s = 0.0

            # Simpson's rule stepping (going backwards from sn-1 down to taun_min+1 by steps of 2)
            for taun in range(sn - 1, taun_min, -2):
                k_1 = curr_vessel.k_alpha[nts * alpha + taun]
                q_1 = math.exp(-(k_2 + k_1) * dt / 2.0) * q_2
                mq_1 = curr_vessel.mR_alpha[nts * alpha + taun] * q_1

                k_0 = curr_vessel.k_alpha[nts * alpha + (taun - 1)]
                q_0 = math.exp(-(k_2 + 4.0 * k_1 + k_0) * dt / 3.0) * q_2
                mq_0 = curr_vessel.mR_alpha[nts * alpha + (taun - 1)] * q_0

                # Simpson's rule accumulator
                rhoR_alpha_s += (mq_2 + 4.0 * mq_1 + mq_0) * dt / 3.0

                # Shift variables for next iteration
                k_2 = k_0
                q_2 = q_0
                mq_2 = mq_0

            # If we have an even number of integration points,
            # do a trapezoid step at the end
            if even_n:
                k_0 = curr_vessel.k_alpha[nts * alpha + taun_min]
                q_0 = math.exp(-(k_2 + k_0) * dt / 2.0) * q_2
                mq_0 = curr_vessel.mR_alpha[nts * alpha + taun_min] * q_0

                rhoR_alpha_s += (mq_2 + mq_0) * dt / 2.0

            # Add the initially present material if taun_min == 0
            if taun_min == 0:
                rhoR_alpha_s += (curr_vessel.rhoR_alpha[nts * alpha] * q_0)

        else:
            # If no degradation or this is the first time step,
            # density just remains what's initially in the system
            rhoR_alpha_s = curr_vessel.rhoR_alpha[nts * alpha]

        # Update the referential density for alpha at current time sn
        curr_vessel.rhoR_alpha[nts * alpha + sn] = rhoR_alpha_s
        # Sum to total referential density
        rhoR_s += rhoR_alpha_s

    # Update total referential density
    curr_vessel.rhoR[sn] = rhoR_s
    


def plot_perturbation_evolution(perturbation_time, homeostatic_value, perturbation_percentage, perturbation_type="step", duration=None, total_time=20):
    """
    Plots the evolution of a quantity over time with a specified perturbation.

    Parameters:
    - perturbation_time (float): The time at which the perturbation is applied.
    - homeostatic_value (float): The homeostatic value of the quantity.
    - perturbation_percentage (float): The magnitude of the perturbation as a percentage of the homeostatic value.
    - perturbation_type (str): The type of perturbation ("step", "latorre2018", or "linear").
    - duration (float, optional): The duration over which the linear perturbation is applied.
    - total_time (float): The total time to simulate and plot.
    """
    times = np.arange(0, total_time + 0.1, 0.1)
    values = [apply_perturbation(t, perturbation_time, homeostatic_value, perturbation_percentage, perturbation_type, duration) for t in times]

    plt.figure(figsize=(10, 6))
    plt.plot(times, values, label=f'{perturbation_type.capitalize()} Perturbation')
    plt.axvline(x=perturbation_time, color='r', linestyle='--', label='Perturbation Time')
    plt.xlabel('Time')
    plt.ylabel('Quantity Value')
    plt.title(f'Evolution of Quantity with {perturbation_type.capitalize()} Perturbation')
    plt.legend()
    plt.grid(True)
    plt.show()

# =============================================================================
# # Example usage for plotting perturbations:
# perturbation_time = 5.0
# homeostatic_pressure = 100.0  # Example homeostatic value for pressure
# perturbation_percentage = 10.0  # Increase by 10%
# duration = 8.0  # Duration for linear perturbation
# 
# # Plot step perturbation
# plot_perturbation_evolution(perturbation_time, homeostatic_pressure, perturbation_percentage, "step")
# 
# # Plot latorre2018 perturbation
# plot_perturbation_evolution(perturbation_time, homeostatic_pressure, perturbation_percentage, "latorre2018")
# 
# # Plot linear perturbation
# plot_perturbation_evolution(perturbation_time, homeostatic_pressure, perturbation_percentage, "linear", duration)
# 
# =============================================================================
