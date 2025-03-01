import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from typing import Tuple, Dict, List


class MagnetorquerOptimizer:
    """
    A tool for optimizing PCB magnetorquer coils by calculating the optimal trace width
    to maximize magnetic moment while considering various constraints.
    """

    # Constants
    COPPER_RESISTIVITY = 1.68e-8  # Ohm-m
    TRACE_SPACING_MIN = 0.15  # mm (minimum trace spacing)
    TRACE_WIDTH_MIN = 0.15  # mm (minimum trace width)
    MM_TO_M = 0.001  # Conversion from mm to m
    OZ_TO_M = 0.0000347  # Conversion from oz/ft² to m (thickness)

    def __init__(self,
                 voltage: float,
                 max_power: float,
                 inner_length: float,
                 inner_width: float,
                 outer_length: float,
                 outer_width: float,
                 copper_weight: float,
                 num_layers: int):
        """
        Initialize the magnetorquer optimizer with design parameters.

        Args:
            voltage: Supply voltage (V)
            max_power: Maximum power (W)
            inner_length: Inner rectangle length (mm)
            inner_width: Inner rectangle width (mm)
            outer_length: Outer rectangle length (mm)
            outer_width: Outer rectangle width (mm)
            copper_weight: Copper weight (oz/ft²)
            num_layers: Number of PCB layers
        """
        self.voltage = voltage
        self.max_power = max_power
        self.inner_length = inner_length * self.MM_TO_M  # Convert to m
        self.inner_width = inner_width * self.MM_TO_M  # Convert to m
        self.outer_length = outer_length * self.MM_TO_M  # Convert to m
        self.outer_width = outer_width * self.MM_TO_M  # Convert to m
        self.copper_thickness = copper_weight * self.OZ_TO_M  # Convert to m
        self.num_layers = num_layers

        # Calculate available winding area (in m²)
        self.winding_area = (self.outer_length * self.outer_width) - (self.inner_length * self.inner_width)

        # Display parameters
        self._display_input_parameters()

    def _display_input_parameters(self):
        """Display the input parameters in a formatted way."""
        print("=== Magnetorquer Optimization Parameters ===")
        print(f"Voltage: {self.voltage} V")
        print(f"Max Power: {self.max_power} W")
        print(f"Inner Dimensions: {self.inner_length / self.MM_TO_M:.2f} × {self.inner_width / self.MM_TO_M:.2f} mm")
        print(f"Outer Dimensions: {self.outer_length / self.MM_TO_M:.2f} × {self.outer_width / self.MM_TO_M:.2f} mm")
        print(f"Copper Weight: {self.copper_thickness / self.OZ_TO_M:.2f} oz")
        print(f"Number of Layers: {self.num_layers}")
        print(f"Winding Area: {self.winding_area * 1e6:.2f} mm²")
        print("============================================")

    def calculate_trace_parameters(self, trace_width_mm: float) -> Dict:
        """
        Calculate coil parameters for a given trace width.

        Args:
            trace_width_mm: Trace width in mm

        Returns:
            Dictionary containing calculated parameters
        """
        # Convert trace width to meters
        trace_width = trace_width_mm * self.MM_TO_M

        # Minimum spacing in meters
        trace_spacing = self.TRACE_SPACING_MIN * self.MM_TO_M

        # Calculate perimeter of one turn (average of inner and outer)
        avg_length = (self.inner_length + self.outer_length) / 2
        avg_width = (self.inner_width + self.outer_width) / 2
        perimeter = 2 * (avg_length + avg_width)

        # Calculate number of turns per layer
        available_width = (self.outer_length - self.inner_length) / 2
        turns_per_side = int(available_width / (trace_width + trace_spacing))

        # Total number of turns accounting for all layers
        total_turns = turns_per_side * self.num_layers

        # Calculate trace cross-sectional area
        trace_area = trace_width * self.copper_thickness

        # Calculate resistance
        trace_length = total_turns * perimeter
        resistance = self.COPPER_RESISTIVITY * trace_length / trace_area

        # Calculate current based on voltage and resistance
        current = self.voltage / resistance if resistance > 0 else 0

        # Calculate power
        power = self.voltage * current

        # Check if power exceeds maximum
        if power > self.max_power:
            # Recalculate current based on max power
            current = self.max_power / self.voltage
            # Recalculate effective resistance
            resistance = self.voltage / current if current > 0 else float('inf')
            # Update power
            power = self.max_power

        # Calculate current density
        current_density = current / trace_area if trace_area > 0 else 0

        # Calculate magnetic moment (A·m²)
        # For a rectangular coil, magnetic moment = current * area * number of turns
        coil_area = self.inner_length * self.inner_width
        magnetic_moment = current * coil_area * total_turns

        return {
            'trace_width_mm': trace_width_mm,
            'trace_spacing_mm': self.TRACE_SPACING_MIN,
            'turns_per_layer': turns_per_side,
            'total_turns': total_turns,
            'resistance': resistance,
            'current': current,
            'power': power,
            'current_density': current_density,
            'magnetic_moment': magnetic_moment,
            'trace_length': trace_length
        }

    def find_optimal_trace_width(self,
                                 step_size: float = 0.01,
                                 min_width: float = None,
                                 max_width: float = None) -> Tuple[float, Dict]:
        """
        Find the optimal trace width that maximizes the magnetic moment.

        Args:
            step_size: Step size for trace width sweep (mm)
            min_width: Minimum trace width to consider (mm), defaults to minimum allowed
            max_width: Maximum trace width to consider (mm), defaults to half of available width

        Returns:
            Tuple of (optimal_trace_width, parameters)
        """
        if min_width is None:
            min_width = self.TRACE_WIDTH_MIN

        if max_width is None:
            # Default to half the available width
            available_width = (self.outer_length - self.inner_length) / 2
            max_width = (available_width / self.MM_TO_M) / 2

        trace_widths = np.arange(min_width, max_width, step_size)
        results = []

        for width in trace_widths:
            params = self.calculate_trace_parameters(width)
            results.append(params)

        # Find the width that maximizes magnetic moment
        magnetic_moments = [r['magnetic_moment'] for r in results]
        if not magnetic_moments:
            return 0, {}

        optimal_index = np.argmax(magnetic_moments)
        optimal_params = results[optimal_index]
        optimal_width = optimal_params['trace_width_mm']

        return optimal_width, optimal_params

    def plot_optimization_results(self, min_width: float = None, max_width: float = None, step_size: float = 0.01):
        """
        Plot optimization results showing how parameters vary with trace width.

        Args:
            min_width: Minimum trace width to consider (mm)
            max_width: Maximum trace width to consider (mm)
            step_size: Step size for trace width sweep (mm)
        """
        if min_width is None:
            min_width = self.TRACE_WIDTH_MIN

        if max_width is None:
            # Default to half the available width
            available_width = (self.outer_length - self.inner_length) / 2
            max_width = (available_width / self.MM_TO_M) / 2

        trace_widths = np.arange(min_width, max_width, step_size)

        # Calculate parameters for each trace width
        magnetic_moments = []
        currents = []
        resistances = []
        powers = []
        total_turns = []

        for width in trace_widths:
            params = self.calculate_trace_parameters(width)
            magnetic_moments.append(params['magnetic_moment'])
            currents.append(params['current'])
            resistances.append(params['resistance'])
            powers.append(params['power'])
            total_turns.append(params['total_turns'])

        # Find optimal width
        optimal_index = np.argmax(magnetic_moments)
        optimal_width = trace_widths[optimal_index]
        optimal_params = self.calculate_trace_parameters(optimal_width)

        # Create figure with subplots
        fig = plt.figure(figsize=(15, 12))
        gs = gridspec.GridSpec(3, 2, height_ratios=[2, 1, 1])

        # Plot 1: Magnetic Moment vs Trace Width
        ax1 = plt.subplot(gs[0, 0])
        ax1.plot(trace_widths, magnetic_moments, 'b-')
        ax1.scatter([optimal_width], [magnetic_moments[optimal_index]], color='red', s=100, zorder=5)
        ax1.set_xlabel('Trace Width (mm)')
        ax1.set_ylabel('Magnetic Moment (A·m²)')
        ax1.set_title('Magnetic Moment vs Trace Width')
        ax1.grid(True)
        ax1.annotate(f'Optimal Width: {optimal_width:.2f} mm\nMax Moment: {magnetic_moments[optimal_index]:.2e} A·m²',
                     xy=(optimal_width, magnetic_moments[optimal_index]),
                     xytext=(optimal_width + 0.1, magnetic_moments[optimal_index] * 0.9),
                     arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                     fontsize=10)

        # Plot 2: Current, Power vs Trace Width
        ax2 = plt.subplot(gs[0, 1])
        line1 = ax2.plot(trace_widths, currents, 'g-', label='Current (A)')
        ax2.set_xlabel('Trace Width (mm)')
        ax2.set_ylabel('Current (A)', color='g')
        ax2.tick_params(axis='y', labelcolor='g')

        ax2_twin = ax2.twinx()
        line2 = ax2_twin.plot(trace_widths, powers, 'r-', label='Power (W)')
        ax2_twin.set_ylabel('Power (W)', color='r')
        ax2_twin.tick_params(axis='y', labelcolor='r')

        # Add horizontal line for max power
        ax2_twin.axhline(y=self.max_power, color='r', linestyle='--', alpha=0.5, label=f'Max Power: {self.max_power}W')

        # Combine legends
        lines = line1 + line2
        labels = [l.get_label() for l in lines]
        ax2.legend(lines, labels, loc='upper left')

        ax2.set_title('Current and Power vs Trace Width')
        ax2.grid(True)

        # Plot 3: Resistance vs Trace Width
        ax3 = plt.subplot(gs[1, 0])
        ax3.plot(trace_widths, resistances, 'purple')
        ax3.set_xlabel('Trace Width (mm)')
        ax3.set_ylabel('Resistance (Ω)')
        ax3.set_title('Resistance vs Trace Width')
        ax3.grid(True)

        # Plot 4: Number of Turns vs Trace Width
        ax4 = plt.subplot(gs[1, 1])
        ax4.plot(trace_widths, total_turns, 'orange')
        ax4.set_xlabel('Trace Width (mm)')
        ax4.set_ylabel('Total Number of Turns')
        ax4.set_title('Total Turns vs Trace Width')
        ax4.grid(True)

        # Plot 5: Visualization of coil
        ax5 = plt.subplot(gs[2, :])
        self._draw_coil_visualization(ax5, optimal_params)

        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        fig.suptitle(
            f'Magnetorquer Optimization Results (Copper: {self.copper_thickness / self.OZ_TO_M:.1f}oz, Layers: {self.num_layers})',
            fontsize=16)

        # Display optimal parameters as text
        optimal_info = (
            f"Optimal Trace Width: {optimal_width:.2f} mm\n"
            f"Trace Spacing: {optimal_params['trace_spacing_mm']:.2f} mm\n"
            f"Total Turns: {optimal_params['total_turns']}\n"
            f"Resistance: {optimal_params['resistance']:.2f} Ω\n"
            f"Current: {optimal_params['current']:.3f} A\n"
            f"Power: {optimal_params['power']:.2f} W\n"
            f"Magnetic Moment: {optimal_params['magnetic_moment']:.2e} A·m²\n"
            f"Total Trace Length: {optimal_params['trace_length']:.2f} m"
        )

        fig.text(0.5, 0.02, optimal_info, ha='center',
                 bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))

        plt.show()

        return optimal_width, optimal_params

    def _draw_coil_visualization(self, ax, params):
        """Draw a visualization of the coil with the given parameters."""
        # Convert back to mm for visualization
        inner_length_mm = self.inner_length / self.MM_TO_M
        inner_width_mm = self.inner_width / self.MM_TO_M
        outer_length_mm = self.outer_length / self.MM_TO_M
        outer_width_mm = self.outer_width / self.MM_TO_M
        trace_width_mm = params['trace_width_mm']
        trace_spacing_mm = params['trace_spacing_mm']

        # Set up the plot
        ax.set_xlim(-outer_length_mm / 2 - 5, outer_length_mm / 2 + 5)
        ax.set_ylim(-outer_width_mm / 2 - 5, outer_width_mm / 2 + 5)
        ax.set_aspect('equal')
        ax.set_title('Magnetorquer Coil Visualization (Top Layer)')

        # Draw outer boundary
        outer_rect = Rectangle((-outer_length_mm / 2, -outer_width_mm / 2),
                               outer_length_mm, outer_width_mm,
                               linewidth=1, edgecolor='black', facecolor='none')
        ax.add_patch(outer_rect)

        # Draw inner boundary
        inner_rect = Rectangle((-inner_length_mm / 2, -inner_width_mm / 2),
                               inner_length_mm, inner_width_mm,
                               linewidth=1, edgecolor='black', facecolor='none')
        ax.add_patch(inner_rect)

        # Simulate trace with spacing (simplified for visualization)
        # Only show a few turns to avoid cluttering
        combined_width = trace_width_mm + trace_spacing_mm
        total_turns_per_layer = min(params['turns_per_layer'], 5)  # Limit to 5 turns for clarity

        # Draw traces from outside to inside
        for i in range(total_turns_per_layer):
            offset = i * combined_width

            # Current turn dimensions
            current_length = outer_length_mm - 2 * offset
            current_width = outer_width_mm - 2 * offset

            # Skip if we've reached the inner boundary
            if (current_length < inner_length_mm or current_width < inner_width_mm):
                break

            rect = Rectangle((-current_length / 2, -current_width / 2),
                             current_length, current_width,
                             linewidth=trace_width_mm, edgecolor='red', facecolor='none')
            ax.add_patch(rect)

        # Add annotations
        ax.text(0, -outer_width_mm / 2 - 3, f"Trace Width: {trace_width_mm:.2f} mm",
                ha='center', va='center')
        ax.text(0, outer_width_mm / 2 + 3, f"Trace Spacing: {trace_spacing_mm:.2f} mm",
                ha='center', va='center')

        # Show total turns info
        if params['turns_per_layer'] > total_turns_per_layer:
            ax.text(0, 0, f"(Showing {total_turns_per_layer} of {params['turns_per_layer']} turns per layer)",
                    ha='center', va='center', fontsize=9, color='blue')

        ax.set_xlabel('Length (mm)')
        ax.set_ylabel('Width (mm)')
        ax.grid(True, linestyle='--', alpha=0.7)


def main():
    """Main function to run the magnetorquer optimizer."""
    # Default parameters
    voltage = 5.0  # V
    max_power = 1.0  # W
    inner_length = 40  # mm
    inner_width = 20  # mm
    outer_length = 60  # mm
    outer_width = 40  # mm
    copper_weight = 2.0  # oz
    num_layers = 6  # layers

    # Create the optimizer
    optimizer = MagnetorquerOptimizer(
        voltage=voltage,
        max_power=max_power,
        inner_length=inner_length,
        inner_width=inner_width,
        outer_length=outer_length,
        outer_width=outer_width,
        copper_weight=copper_weight,
        num_layers=num_layers
    )

    # Find and plot optimal results
    optimal_width, optimal_params = optimizer.plot_optimization_results()

    print("\n=== Optimization Results ===")
    print(f"Optimal Trace Width: {optimal_width:.2f} mm")
    print(f"Magnetic Moment: {optimal_params['magnetic_moment']:.2e} A·m²")
    print(f"Current: {optimal_params['current']:.3f} A")
    print(f"Power: {optimal_params['power']:.2f} W")
    print(f"Resistance: {optimal_params['resistance']:.2f} Ω")
    print(f"Total Turns: {optimal_params['total_turns']}")
    print(f"Current Density: {optimal_params['current_density']:.2e} A/m²")
    print(f"Total Trace Length: {optimal_params['trace_length']:.2f} m")

    # Interactive mode for exploring different parameters
    explore_more = input("\nWould you like to explore different parameters? (y/n): ")
    if explore_more.lower() == 'y':
        while True:
            print("\n=== Parameter Exploration ===")
            print("Enter new parameters (or press Enter to keep current value):")

            try:
                new_voltage = input(f"Voltage [{voltage} V]: ")
                voltage = float(new_voltage) if new_voltage else voltage

                new_max_power = input(f"Max Power [{max_power} W]: ")
                max_power = float(new_max_power) if new_max_power else max_power

                new_inner_length = input(f"Inner Length [{inner_length} mm]: ")
                inner_length = float(new_inner_length) if new_inner_length else inner_length

                new_inner_width = input(f"Inner Width [{inner_width} mm]: ")
                inner_width = float(new_inner_width) if new_inner_width else inner_width

                new_outer_length = input(f"Outer Length [{outer_length} mm]: ")
                outer_length = float(new_outer_length) if new_outer_length else outer_length

                new_outer_width = input(f"Outer Width [{outer_width} mm]: ")
                outer_width = float(new_outer_width) if new_outer_width else outer_width

                new_copper_weight = input(f"Copper Weight [{copper_weight} oz]: ")
                copper_weight = float(new_copper_weight) if new_copper_weight else copper_weight

                new_num_layers = input(f"Number of Layers [{num_layers}]: ")
                num_layers = int(new_num_layers) if new_num_layers else num_layers

                # Create new optimizer with updated parameters
                optimizer = MagnetorquerOptimizer(
                    voltage=voltage,
                    max_power=max_power,
                    inner_length=inner_length,
                    inner_width=inner_width,
                    outer_length=outer_length,
                    outer_width=outer_width,
                    copper_weight=copper_weight,
                    num_layers=num_layers
                )

                # Find and plot optimal results
                optimal_width, optimal_params = optimizer.plot_optimization_results()

                continue_exploring = input("\nContinue exploring? (y/n): ")
                if continue_exploring.lower() != 'y':
                    break

            except ValueError as e:
                print(f"Error: {e}. Please enter valid numbers.")
                continue

    print("\nOptimization complete. Thank you for using the Magnetorquer Optimization Tool!")


if __name__ == "__main__":
    main()