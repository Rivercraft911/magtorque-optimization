import numpy as np
from scipy.optimize import minimize
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Constants and Constraints
rho_copper = 1.68e-8        # Copper resistivity (Ω·m)
oz_to_mm = 0.035            # Copper weight to millimeters conversion (mm)
J_max = 2.5e8                # Maximum current density (A/m²)
min_trace_width = 0.1       # Minimum trace width (mm)
max_trace_width = 10.0      # Maximum trace width (mm)
min_trace_spacing = 0.1     # Minimum trace spacing (mm)

def calculate_max_turns(trace_width, outer_length, inner_length):
    if trace_width <= 0:
        return 0

    min_inner_clearance = trace_width + 2 * min_trace_spacing
    effective_inner_length = inner_length + 2 * min_inner_clearance
    available_length = (outer_length - effective_inner_length) / 2

    if available_length <= 0:
        return 0

    turn_width_total = trace_width + min_trace_spacing
    max_turns_length = int(available_length / turn_width_total)

    return max(1, max_turns_length)

def calculate_turn_length(turn_number, trace_width, num_turns, outer_length):
    offset = turn_number * (trace_width + min_trace_spacing)
    current_length = outer_length - 2 * offset
    perimeter = 2 * current_length

    if turn_number < num_turns - 1:
        connection_length = trace_width + min_trace_spacing
    else:
        connection_length = 0  # No connection for the last turn

    return perimeter + connection_length

def calculate_total_length(trace_width, num_turns, outer_length):
    total_length = sum(calculate_turn_length(n, trace_width, num_turns, outer_length) for n in range(num_turns))
    total_length *= coil_layers
    return total_length / 1000  # Convert mm to meters

def calculate_area(turn_number, trace_width, outer_length):
    offset = turn_number * (trace_width + min_trace_spacing)
    length = outer_length - 2 * offset
    return (length ** 2) / 1e6  # Convert mm² to m²

def calculate_total_area(trace_width, num_turns, outer_length):
    total_area = sum(calculate_area(n, trace_width, outer_length) for n in range(num_turns))
    total_area *= coil_layers
    return total_area

def calculate_resistance(trace_width, total_length):
    A_cross = (trace_width * copper_thickness) / 1e6  # Convert mm² to m²
    if A_cross == 0:
        return np.inf
    R = (rho_copper * total_length) / A_cross
    return R

def calculate_current(R, trace_width):
    if R == 0:
        return 0

    I_P = max_power / voltage_V
    A_cross = (trace_width * copper_thickness) / 1e6  # Convert mm² to m²
    I_J = J_max * A_cross
    I_V = voltage_V / R
    I = min(I_P, I_J, I_V)
    return I

def calculate_magnetic_moment(total_area, current):
    m = total_area * current
    return m

def check_constraints(trace_width, outer_length, inner_length):
    if trace_width < min_trace_width or trace_width > max_trace_width:
        return False

    num_turns = calculate_max_turns(trace_width, outer_length, inner_length)
    if num_turns == 0:
        return False

    total_length = calculate_total_length(trace_width, num_turns, outer_length)
    R = calculate_resistance(trace_width, total_length)
    current = calculate_current(R, trace_width)
    A_cross = (trace_width * copper_thickness) / 1e6  # Convert mm² to m²
    current_density = current / A_cross
    if current_density > J_max:
        return False

    power = current * voltage_V
    if power > max_power:
        return False

    return True

def objective_function(x, outer_length, inner_length):
    trace_width = x[0]
    if not check_constraints(trace_width, outer_length, inner_length):
        return 0

    num_turns = calculate_max_turns(trace_width, outer_length, inner_length)
    total_length = calculate_total_length(trace_width, num_turns, outer_length)
    R = calculate_resistance(trace_width, total_length)
    current = calculate_current(R, trace_width)
    total_area = calculate_total_area(trace_width, num_turns, outer_length)
    m = calculate_magnetic_moment(total_area, current)

    return -m  # Negative because we are minimizing

def graph_optimal_magnetic_moments(outer_length):
    inner_lengths = np.linspace(1, outer_length - 1, 50)  # Inner lengths in mm
    optimal_magnetic_moments = []
    optimal_trace_widths = []

    for inner_length in inner_lengths:
        result = minimize(
            objective_function,
            x0=[min_trace_width],
            bounds=[(min_trace_width, max_trace_width)],
            method='L-BFGS-B',
            args=(outer_length, inner_length)
        )
        optimal_trace_width = result.x[0]
        num_turns = calculate_max_turns(optimal_trace_width, outer_length, inner_length)
        total_length = calculate_total_length(optimal_trace_width, num_turns, outer_length)
        R = calculate_resistance(optimal_trace_width, total_length)
        current = calculate_current(R, optimal_trace_width)
        total_area = calculate_total_area(optimal_trace_width, num_turns, outer_length)
        magnetic_moment = calculate_magnetic_moment(total_area, current)
        optimal_magnetic_moments.append(magnetic_moment)
        optimal_trace_widths.append(optimal_trace_width)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    ax1.plot(inner_lengths, optimal_magnetic_moments, label="Optimal Magnetic Moment")
    ax1.set_xlabel("Inner Length (mm)")
    ax1.set_ylabel("Magnetic Moment (A·mm²)")
    ax1.set_title("Optimal Magnetic Moment vs. Inner Length")
    ax1.legend()
    ax1.grid()

    ax2.plot(inner_lengths, optimal_trace_widths, label="Optimal Trace Width", color='orange')
    ax2.set_xlabel("Inner Length (mm)")
    ax2.set_ylabel("Trace Width (mm)")
    ax2.set_title("Optimal Trace Width vs. Inner Length")
    ax2.legend()
    ax2.grid()

    plt.tight_layout()
    plt.show()

# GUI for user input
root = tk.Tk()
root.title("Magnetic Coil Optimization Tool")

# Input fields for user parameters
voltage_var = tk.DoubleVar(value=8.0)
power_var = tk.DoubleVar(value=0.5)
outer_length_var = tk.DoubleVar(value=80.0)
copper_weight_var = tk.DoubleVar(value=2.0)
layers_var = tk.IntVar(value=6)

input_frame = tk.Frame(root)
input_frame.pack(side=tk.LEFT, padx=10, pady=10)

tk.Label(input_frame, text="Voltage (V):").pack()
tk.Entry(input_frame, textvariable=voltage_var).pack()
tk.Label(input_frame, text="Max Power (W):").pack()
tk.Entry(input_frame, textvariable=power_var).pack()
tk.Label(input_frame, text="Outer Length (mm):").pack()
tk.Entry(input_frame, textvariable=outer_length_var).pack()
tk.Label(input_frame, text="Copper Weight (oz):").pack()
tk.Entry(input_frame, textvariable=copper_weight_var).pack()
tk.Label(input_frame, text="Number of Layers:").pack()
tk.Entry(input_frame, textvariable=layers_var).pack()

def update_graphs():
    global voltage_V, max_power, outer_length, copper_thickness, coil_layers

    # Fetch user input values from the GUI variables
    voltage_V = voltage_var.get()
    max_power = power_var.get()
    outer_length = outer_length_var.get()  # Outer length in mm
    copper_weight_oz = copper_weight_var.get()
    copper_thickness = copper_weight_oz * oz_to_mm  # Copper thickness in mm
    coil_layers = layers_var.get()

    # Plot optimal magnetic moments for different inner lengths
    graph_optimal_magnetic_moments(outer_length)

# Update button
update_button = ttk.Button(input_frame, text="Update Graphs", command=update_graphs)
update_button.pack(pady=10)

root.mainloop()
