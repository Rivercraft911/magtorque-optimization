import numpy as np
from scipy.optimize import minimize
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Constants and Constraints
rho_copper = 1.68e-8        # Copper resistivity (Ω·m)
oz_to_m = 0.035e-3         # Copper weight to meters conversion (m)
J_max = 2.5e8                # Maximum current density (A/m²)
min_trace_width = 0.0001    # Minimum trace width (m)
max_trace_width = 0.0100    # Maximum trace width (m)
min_trace_spacing = 0.0001  # Minimum trace spacing (m)

# Define functions (calculations, as in your script)
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

def calculate_optimal_trace_width():
    result = minimize(
        objective_function,
        x0=[min_trace_width],
        bounds=[(min_trace_width, max_trace_width)],
        method='SLSQP'
    )
    return result.x[0]

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


def update_graphs():
    global voltage_V, max_power, inner_length, inner_width, outer_length, outer_width, copper_thickness, coil_layers

    # Fetch user input values from the GUI variables
    voltage_V = voltage_var.get()
    max_power = power_var.get()
    inner_length = inner_length_var.get() / 1000
    inner_width = inner_width_var.get() / 1000
    outer_length = outer_length_var.get() / 1000
    outer_width = outer_width_var.get() / 1000
    copper_weight_oz = copper_weight_var.get()
    copper_thickness = copper_weight_oz * oz_to_m
    coil_layers = layers_var.get() - 1  # Adjust for coil layers

    # Calculate optimal trace width using the updated parameters
    optimal_trace_width = calculate_optimal_trace_width()
    num_turns = calculate_max_turns(optimal_trace_width)
    total_length = calculate_total_length(optimal_trace_width, num_turns)
    R = calculate_resistance(optimal_trace_width, total_length)
    current = calculate_current(R, optimal_trace_width)
    total_area = calculate_total_area(optimal_trace_width, num_turns)
    magnetic_moment = calculate_magnetic_moment(total_area, current)
    power_to_moment_ratio = current * voltage_V / magnetic_moment if magnetic_moment != 0 else np.inf

    # Update Graphs
    trace_widths = np.linspace(min_trace_width, max_trace_width, 100)
    magnetic_moments = [
        calculate_magnetic_moment(calculate_total_area(w, calculate_max_turns(w)),
                                  calculate_current(
                                      calculate_resistance(w, calculate_total_length(w, calculate_max_turns(w))), w))
        for w in trace_widths
    ]
    power_to_moment_ratios = [
        (calculate_current(calculate_resistance(w, calculate_total_length(w, calculate_max_turns(w))),
                           w) * voltage_V / m if m != 0 else np.inf)
        for w, m in zip(trace_widths, magnetic_moments)
    ]
    ax1.clear()
    ax1.plot(trace_widths, magnetic_moments, label="Magnetic Moment")
    ax1.axvline(optimal_trace_width, color="r", linestyle="--", label="Optimal Width")
    ax1.set_title("Magnetic Moment vs. Trace Width")
    ax1.set_xlabel("Trace Width (m)")
    ax1.set_ylabel("Magnetic Moment (A·m²)")
    ax1.legend()

    ax2.clear()
    ax2.plot(trace_widths, power_to_moment_ratios, label="Power-to-Magnetic Moment Ratio")
    ax2.axvline(optimal_trace_width, color="r", linestyle="--", label="Optimal Width")
    ax2.set_title("Power/Magnetic Moment Ratio vs. Trace Width")
    ax2.set_xlabel("Trace Width (m)")
    ax2.set_ylabel("Efficiency (W/A·m²)")
    ax2.legend()

    canvas.draw()

# Tkinter GUI setup
root = tk.Tk()
root.title("Magnetic Coil Optimization Tool")

# Input fields for user parameters
voltage_var = tk.DoubleVar(value=50.0)
power_var = tk.DoubleVar(value=10.0)
inner_length_var = tk.DoubleVar(value=30.0)
inner_width_var = tk.DoubleVar(value=30.0)
outer_length_var = tk.DoubleVar(value=100.0)
outer_width_var = tk.DoubleVar(value=100.0)
copper_weight_var = tk.DoubleVar(value=1.0)
layers_var = tk.IntVar(value=2)

# Layout inputs with labels
tk.Label(root, text="Voltage (V):").pack()
tk.Entry(root, textvariable=voltage_var).pack()
tk.Label(root, text="Max Power (W):").pack()
tk.Entry(root, textvariable=power_var).pack()
tk.Label(root, text="Inner Length (mm):").pack()
tk.Entry(root, textvariable=inner_length_var).pack()
tk.Label(root, text="Inner Width (mm):").pack()
tk.Entry(root, textvariable=inner_width_var).pack()
tk.Label(root, text="Outer Length (mm):").pack()
tk.Entry(root, textvariable=outer_length_var).pack()
tk.Label(root, text="Outer Width (mm):").pack()
tk.Entry(root, textvariable=outer_width_var).pack()
tk.Label(root, text="Copper Weight (oz):").pack()
tk.Entry(root, textvariable=copper_weight_var).pack()
tk.Label(root, text="Number of Layers:").pack()
tk.Entry(root, textvariable=layers_var).pack()

# Update button
update_button = ttk.Button(root, text="Update Graphs", command=update_graphs)
update_button.pack()

# Set up Matplotlib figures with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack()

root.mainloop()
