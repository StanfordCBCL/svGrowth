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