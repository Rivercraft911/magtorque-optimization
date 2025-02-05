import numpy as np
from scipy.optimize import minimize
import tkinter as tk
from tkinter import ttk, font
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# Constants and Constraints
rho_copper = 1.68e-8        # Copper resistivity (Ω·m)
oz_to_m = 0.035e-3         # Copper weight to meters conversion (m)
J_max = 2.5e8                # Maximum current density (A/m²)
min_trace_width = 0.00018    # Minimum trace width (m)
max_trace_width = 0.0100    # Maximum trace width (m)
min_trace_spacing = 0.0002  # Minimum trace spacing (m)

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
        method='L-BFGS-B'
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

class MagnetorquerOptimizer(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Magnetorquer Coil Optimization Tool")
        self.configure(bg='#f0f0f0')
        
        # Initialize variables first
        self.initialize_variables()
        
        # Create matplotlib figure before creating frames
        self.create_figure()
        
        # Set up the main container with grid
        self.main_container = ttk.Frame(self)
        self.main_container.grid(row=0, column=0, sticky="nsew", padx=20, pady=20)
        
        # Configure style
        style = ttk.Style()
        style.configure('Title.TLabel', font=('Helvetica', 14, 'bold'))
        style.configure('Results.TLabel', font=('Courier', 10))
        
        # Create frames in the correct order
        self.create_input_frame()
        self.create_results_frame()
        self.create_graphs_frame()
        
        # Initial update
        self.update_graphs()

    def initialize_variables(self):
        self.voltage_var = tk.DoubleVar(value=8.0)
        self.power_var = tk.DoubleVar(value=0.5)
        self.inner_length_var = tk.DoubleVar(value=40.0)
        self.inner_width_var = tk.DoubleVar(value=40.0)
        self.outer_length_var = tk.DoubleVar(value=80.0)
        self.outer_width_var = tk.DoubleVar(value=80.0)
        self.copper_weight_var = tk.DoubleVar(value=2.0)
        self.layers_var = tk.IntVar(value=6)

    def create_input_frame(self):
        input_frame = ttk.LabelFrame(self.main_container, text="Input Parameters", padding="10")
        input_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

        parameters = [
            ("Voltage (V)", self.voltage_var),
            ("Max Power (W)", self.power_var),
            ("Inner Length (mm)", self.inner_length_var),
            ("Inner Width (mm)", self.inner_width_var),
            ("Outer Length (mm)", self.outer_length_var),
            ("Outer Width (mm)", self.outer_width_var),
            ("Copper Weight (oz)", self.copper_weight_var),
            ("Number of Layers", self.layers_var)
        ]

        for i, (label_text, var) in enumerate(parameters):
            ttk.Label(input_frame, text=label_text).grid(row=i, column=0, sticky="w", pady=2)
            entry = ttk.Entry(input_frame, width=10, textvariable=var)
            entry.grid(row=i, column=1, padx=5, pady=2)

        ttk.Button(input_frame, text="Update", command=self.update_graphs).grid(
            row=len(parameters), column=0, columnspan=2, pady=10)

    def create_results_frame(self):
        results_frame = ttk.LabelFrame(self.main_container, text="Results", padding="10")
        results_frame.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
        
        self.result_label = ttk.Label(results_frame, style='Results.TLabel', justify="left")
        self.result_label.grid(row=0, column=0, sticky="w")

    def create_figure(self):
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        plt.style.use('classic')
        self.fig.patch.set_facecolor('#f0f0f0')

    def create_graphs_frame(self):
        graphs_frame = ttk.LabelFrame(self.main_container, text="Visualization", padding="10")
        graphs_frame.grid(row=0, column=1, rowspan=2, sticky="nsew", padx=5, pady=5)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=graphs_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def update_graphs(self):
        # Get values from GUI
        global voltage_V, max_power, inner_length, inner_width, outer_length, outer_width, copper_thickness, coil_layers
        
        voltage_V = self.voltage_var.get()
        max_power = self.power_var.get()
        inner_length = self.inner_length_var.get() / 1000
        inner_width = self.inner_width_var.get() / 1000
        outer_length = self.outer_length_var.get() / 1000
        outer_width = self.outer_width_var.get() / 1000
        copper_weight_oz = self.copper_weight_var.get()
        copper_thickness = copper_weight_oz * oz_to_m
        coil_layers = self.layers_var.get()

        try:
            # Calculate optimal values
            optimal_trace_width = calculate_optimal_trace_width()
            num_turns = calculate_max_turns(optimal_trace_width)
            total_length = calculate_total_length(optimal_trace_width, num_turns)
            R = calculate_resistance(optimal_trace_width, total_length)
            current = calculate_current(R, optimal_trace_width)
            total_area = calculate_total_area(optimal_trace_width, num_turns)
            magnetic_moment = calculate_magnetic_moment(total_area, current)

            # Update results display
            result_text = (
                f"Optimal Design Parameters:\n"
                f"━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                f"Trace Width:      {optimal_trace_width * 1000:.3f} mm\n"
                f"Magnetic Moment:  {magnetic_moment:.6f} A·m²\n"
                f"Power Usage:      {current * voltage_V:.3f} W\n"
                f"Current:          {current:.3f} A\n"
                f"Resistance:       {R:.3f} Ω\n"
                f"Current Density:  {current / (optimal_trace_width * copper_thickness) / 1e6:.2f} A/mm²\n"
                f"Total Turns:      {num_turns * coil_layers}\n"
                f"Turns per Layer:  {num_turns}\n"
                f"Total Length:     {total_length:.3f} m"
            )
            self.result_label.config(text=result_text)

            # Clear previous plots
            self.ax1.clear()
            self.ax2.clear()

            # Generate data for plots
            trace_widths = np.linspace(min_trace_width, max_trace_width, 200)
            magnetic_moments = []
            power_to_moment_ratios = []

            for w in trace_widths:
                n_turns = calculate_max_turns(w)
                t_length = calculate_total_length(w, n_turns)
                r = calculate_resistance(w, t_length)
                i = calculate_current(r, w)
                t_area = calculate_total_area(w, n_turns)
                m = calculate_magnetic_moment(t_area, i)
                magnetic_moments.append(m)
                power_to_moment_ratios.append((i * voltage_V / m) if m != 0 else np.inf)

            # Plot magnetic moment
            self.ax1.plot(trace_widths * 1000, magnetic_moments, 'b-', label="Magnetic Moment")
            self.ax1.axvline(optimal_trace_width * 1000, color='r', linestyle='--', label="Optimal Width")
            self.ax1.set_title("Magnetic Moment vs. Trace Width", pad=10)
            self.ax1.set_xlabel("Trace Width (mm)")
            self.ax1.set_ylabel("Magnetic Moment (A·m²)")
            self.ax1.grid(True, linestyle='--', alpha=0.7)
            self.ax1.legend()

            # Plot power to moment ratio
            self.ax2.plot(trace_widths * 1000, power_to_moment_ratios, 'g-', label="Power/Moment Ratio")
            self.ax2.axvline(optimal_trace_width * 1000, color='r', linestyle='--', label="Optimal Width")
            self.ax2.set_title("Power/Magnetic Moment Ratio vs. Trace Width", pad=10)
            self.ax2.set_xlabel("Trace Width (mm)")
            self.ax2.set_ylabel("Power/Moment Ratio (W/A·m²)")
            self.ax2.grid(True, linestyle='--', alpha=0.7)
            self.ax2.legend()
            self.fig.tight_layout()
            self.canvas.draw()

        except Exception as e:
            self.result_label.config(text=f"Error in calculations: {str(e)}")

if __name__ == "__main__":
    app = MagnetorquerOptimizer()
    app.mainloop()