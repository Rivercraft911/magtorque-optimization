import numpy as np
from scipy.optimize import minimize

# Constants and Constraints
rho_copper = 1.68e-8        # Copper resistivity (Ω·m)
oz_to_m = 0.035e-3         # Copper weight to meters conversion (m)
J_max = 2.5e8                # Maximum current density (A/m²)
min_trace_width = 0.0001    # Minimum trace width (m)
max_trace_width = 0.0100    # Maximum trace width (m)
min_trace_spacing = 0.0001  # Minimum trace spacing (m)
voltage_V = float(input("Enter voltage supplied (V): "))
max_power = float(input("Enter max power (W): "))

# User Inputs
inner_length_mm = float(input("Enter inner length (mm): "))
inner_width_mm = float(input("Enter inner width (mm): "))
outer_length_mm = float(input("Enter outer length (mm): "))
outer_width_mm = float(input("Enter outer width (mm): "))
copper_weight_oz = float(input("Enter copper weight (oz): "))
number_of_layers = int(input("Enter number of layers: "))

# Convert dimensions to meters
inner_length = inner_length_mm / 1000
inner_width = inner_width_mm / 1000
outer_length = outer_length_mm / 1000
outer_width = outer_width_mm / 1000

copper_thickness = copper_weight_oz * oz_to_m
coil_layers = number_of_layers - 1  # One layer connections

# Define functions
def calculate_max_turns(trace_width):
    if trace_width <= 0:
        return 0

    min_inner_clearance = trace_width + 2 * min_trace_spacing
    effective_inner_length = inner_length + 2 * min_inner_clearance
    effective_inner_width = inner_width + 2 * min_inner_clearance
    available_length = (outer_length - effective_inner_length) / 2
    available_width = (outer_width - effective_inner_width) / 2

    if available_length <= 0 or available_width <= 0:
        return 0

    turn_width_total = trace_width + min_trace_spacing
    max_turns_length = int(available_length / turn_width_total)
    max_turns_width = int(available_width / turn_width_total)

    return max(1, min(max_turns_length, max_turns_width))

def calculate_turn_length(turn_number, trace_width, num_turns):
    offset = turn_number * (trace_width + min_trace_spacing)
    current_length = outer_length - 2 * offset
    current_width = outer_width - 2 * offset
    perimeter = 2 * (current_length + current_width)

    if turn_number < num_turns - 1:
        connection_length = trace_width + min_trace_spacing
    else:
        connection_length = 0  # No connection for the last turn

    return perimeter + connection_length

def calculate_total_length(trace_width, num_turns):
    total_length = sum(calculate_turn_length(n, trace_width, num_turns) for n in range(num_turns))
    total_length *= coil_layers
    return total_length

def calculate_area(turn_number, trace_width):
    offset = turn_number * (trace_width + min_trace_spacing)
    length = outer_length - 2 * offset
    width = outer_width - 2 * offset
    return length * width

def calculate_total_area(trace_width, num_turns):
    total_area = sum(calculate_area(n, trace_width) for n in range(num_turns))
    total_area *= coil_layers
    return total_area

def calculate_resistance(trace_width, total_length):
    A_cross = trace_width * copper_thickness
    if A_cross == 0:
        return np.inf
    R = (rho_copper * total_length) / A_cross
    return R

def calculate_current(R, trace_width):
    if R == 0:
        return 0

    I_P = max_power / voltage_V
    A_cross = trace_width * copper_thickness
    I_J = J_max * A_cross
    I_V = voltage_V / R
    I = min(I_P, I_J, I_V)
    return I

def calculate_magnetic_moment(total_area, current):
    m = total_area * current
    return m

def check_constraints(trace_width):
    if trace_width < min_trace_width or trace_width > max_trace_width:
        return False

    num_turns = calculate_max_turns(trace_width)
    if num_turns == 0:
        return False

    total_length = calculate_total_length(trace_width, num_turns)
    R = calculate_resistance(trace_width, total_length)
    current = calculate_current(R, trace_width)
    A_cross = trace_width * copper_thickness
    current_density = current / A_cross
    if current_density > J_max:
        return False

    power = current * voltage_V
    if power > max_power:
        return False

    return True

def objective_function(x):
    trace_width = x[0]
    if not check_constraints(trace_width):
        return 0

    num_turns = calculate_max_turns(trace_width)
    total_length = calculate_total_length(trace_width, num_turns)
    R = calculate_resistance(trace_width, total_length)
    current = calculate_current(R, trace_width)
    total_area = calculate_total_area(trace_width, num_turns)
    m = calculate_magnetic_moment(total_area, current)

    return -m  # Negative because we are minimizing

# Optimization
result = minimize(
    objective_function,
    x0=[min_trace_width],
    bounds=[(min_trace_width, max_trace_width)],
    method='L-BFGS-B'
)

# Output Results
optimal_trace_width = result.x[0]
if not check_constraints(optimal_trace_width):
    print("No feasible design found")
else:
    num_turns = calculate_max_turns(optimal_trace_width)
    total_length = calculate_total_length(optimal_trace_width, num_turns)
    R = calculate_resistance(optimal_trace_width, total_length)
    current = calculate_current(R, optimal_trace_width)
    total_area = calculate_total_area(optimal_trace_width, num_turns)
    m = calculate_magnetic_moment(total_area, current)
    A_cross = optimal_trace_width * copper_thickness
    current_density = current / A_cross

    print("\nOptimal Design Parameters:")
    print(f"Trace Width: {optimal_trace_width * 1000:.3f} mm")
    print(f"Magnetic Moment: {m:.6f} A·m²")
    print(f"Power Consumption: {current * voltage_V:.6f} W")
    print(f"Current: {current:.6f} A")
    print(f"Resistance: {R:.6f} Ω")
    print(f"Current Density: {current_density / 1e6:.2f} A/mm²")
    print(f"Total Turns: {num_turns * coil_layers}")
    print(f"Turns per Layer: {num_turns}")
    print(f"Total Length of Coil: {total_length:.3f} m")
    print(f"Outer Dimensions: {outer_length_mm} mm x {outer_width_mm} mm")
