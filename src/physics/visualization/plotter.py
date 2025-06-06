"""
Plotting and visualization tools for SEMITIP simulation results.

This module provides functions to visualize potential distributions,
band diagrams, current-voltage curves, and other simulation outputs.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from typing import List, Optional, Tuple, Dict
import pickle

from ...simulation.multint import SimulationResults
from ..core.potential import PotentialProfile
from ..solvers.grid import Grid3D


class SEMITIPPlotter:
    """Main plotting class for SEMITIP visualization."""
    
    def __init__(self, style: str = 'default'):
        """
        Initialize plotter with matplotlib style.
        
        Args:
            style: Matplotlib style name or 'default'
        """
        if style != 'default':
            plt.style.use(style)
        
        # Set up default colors
        self.colors = {
            'potential': 'blue',
            'vb': 'red',
            'cb': 'green',
            'current': 'black',
            'surface': 'orange'
        }
        
        # Figure settings
        self.fig_dpi = 300
        self.fig_size = (8, 6)
    
    def plot_potential_profile(self, profile: PotentialProfile,
                             show_band_edges: bool = True,
                             band_gap: float = 1.42,
                             fermi_level: float = 0.0) -> plt.Figure:
        """
        Plot 1D potential profile along z-axis.
        
        Args:
            profile: Potential profile from simulation
            show_band_edges: Whether to show VB/CB edges
            band_gap: Semiconductor band gap (eV)
            fermi_level: Fermi level position (eV)
            
        Returns:
            Matplotlib figure
        """
        fig, ax = plt.subplots(figsize=self.fig_size)
        
        # Get combined profile
        z, pot = profile.get_combined_profile()
        
        # Plot potential
        ax.plot(z, pot, 'b-', linewidth=2, label='Potential')
        
        # Add band edges if requested
        if show_band_edges:
            vb_edge = pot - band_gap
            cb_edge = pot
            ax.plot(z, vb_edge, 'r--', linewidth=1.5, label='VB edge')
            ax.plot(z, cb_edge, 'g--', linewidth=1.5, label='CB edge')
        
        # Add Fermi level
        ax.axhline(y=fermi_level, color='k', linestyle=':', 
                  linewidth=1, label='Fermi level')
        
        # Mark interface
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax.text(0.02, 0.95, 'Surface', transform=ax.transAxes,
               verticalalignment='top')
        
        # Labels and formatting
        ax.set_xlabel('Position z (nm)')
        ax.set_ylabel('Energy (eV)')
        ax.set_title(f'Potential Profile at r={profile.r_position:.1f} nm')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Set reasonable y-limits
        ax.set_ylim(min(pot) - 0.5, max(pot) + 0.5)
        
        plt.tight_layout()
        return fig
    
    def plot_potential_2d(self, grid: Grid3D, potential_3d: np.ndarray,
                        phi_slice: int = 0,
                        contour_levels: int = 20) -> plt.Figure:
        """
        Plot 2D cross-section of potential distribution.
        
        Args:
            grid: Computational grid
            potential_3d: 3D potential array
            phi_slice: Angular slice index
            contour_levels: Number of contour levels
            
        Returns:
            Matplotlib figure
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Extract 2D slice
        potential_2d = potential_3d[:, :, phi_slice]
        
        # Create meshgrid
        R, Z = np.meshgrid(grid.r, grid.z, indexing='ij')
        
        # Plot filled contours
        levels = np.linspace(potential_2d.min(), potential_2d.max(), contour_levels)
        cs = ax.contourf(R, Z, potential_2d, levels=levels, cmap='RdBu_r')
        
        # Add contour lines
        ax.contour(R, Z, potential_2d, levels=levels, colors='black', 
                  linewidths=0.5, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(cs, ax=ax, label='Potential (V)')
        
        # Mark surface
        ax.axhline(y=0, color='black', linewidth=2)
        ax.text(grid.r[-1]*0.9, 2, 'Vacuum', ha='right')
        ax.text(grid.r[-1]*0.9, -2, 'Semiconductor', ha='right')
        
        # Labels and formatting
        ax.set_xlabel('Radial position r (nm)')
        ax.set_ylabel('Vertical position z (nm)')
        ax.set_title(f'Potential Distribution (φ = {grid.phi[phi_slice]*180/np.pi:.0f}°)')
        
        # Set aspect ratio
        ax.set_aspect('equal')
        
        plt.tight_layout()
        return fig
    
    def plot_current_voltage(self, results: List[SimulationResults],
                           plot_components: bool = True) -> plt.Figure:
        """
        Plot current-voltage characteristics.
        
        Args:
            results: List of simulation results at different biases
            plot_components: Whether to show VB/CB components
            
        Returns:
            Matplotlib figure
        """
        # Extract data
        biases = [r.bias_voltage for r in results]
        currents = [r.current for r in results]
        
        # Create figure
        fig, ax = plt.subplots(figsize=self.fig_size)
        
        # Plot total current
        ax.plot(biases, currents, 'ko-', linewidth=2, markersize=6,
               label='Total current')
        
        # Plot components if available and requested
        if plot_components and hasattr(results[0], 'current_components'):
            vb_currents = [r.current_components.get('vb', 0) for r in results]
            cb_currents = [r.current_components.get('cb', 0) for r in results]
            
            ax.plot(biases, vb_currents, 'r^--', linewidth=1.5,
                   markersize=5, label='VB current')
            ax.plot(biases, cb_currents, 'gs--', linewidth=1.5,
                   markersize=5, label='CB current')
        
        # Formatting
        ax.set_xlabel('Bias Voltage (V)')
        ax.set_ylabel('Current (A)')
        ax.set_title('Tunneling Current vs Bias Voltage')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Use log scale if current spans many orders of magnitude
        if len(currents) > 1 and max(np.abs(currents)) / min(np.abs(currents)) > 100:
            ax.set_yscale('log')
            ax.set_ylabel('Current (A) - log scale')
        
        plt.tight_layout()
        return fig
    
    def plot_band_bending(self, results: List[SimulationResults]) -> plt.Figure:
        """
        Plot band bending vs bias voltage.
        
        Args:
            results: List of simulation results
            
        Returns:
            Matplotlib figure
        """
        biases = [r.bias_voltage for r in results]
        band_bendings = [r.band_bending for r in results]
        
        fig, ax = plt.subplots(figsize=self.fig_size)
        
        ax.plot(biases, band_bendings, 'bo-', linewidth=2, markersize=6)
        
        ax.set_xlabel('Bias Voltage (V)')
        ax.set_ylabel('Band Bending (eV)')
        ax.set_title('Surface Band Bending vs Bias')
        ax.grid(True, alpha=0.3)
        
        # Add zero line
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        return fig
    
    def plot_convergence(self, convergence_info: Dict) -> plt.Figure:
        """
        Plot convergence history of Poisson solver.
        
        Args:
            convergence_info: Convergence information from solver
            
        Returns:
            Matplotlib figure
        """
        fig, ax = plt.subplots(figsize=self.fig_size)
        
        history = convergence_info['convergence_history']
        iterations = np.arange(1, len(history) + 1)
        
        ax.semilogy(iterations, history, 'b-', linewidth=2)
        
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Maximum Change (V)')
        ax.set_title('Poisson Solver Convergence')
        ax.grid(True, alpha=0.3, which='both')
        
        # Mark convergence point
        if convergence_info['converged']:
            ax.axvline(x=convergence_info['iterations'], 
                      color='green', linestyle='--',
                      label=f"Converged at {convergence_info['iterations']}")
            ax.legend()
        
        plt.tight_layout()
        return fig
    
    def create_summary_figure(self, result: SimulationResults,
                            band_gap: float = 1.42) -> plt.Figure:
        """
        Create a summary figure with multiple panels.
        
        Args:
            result: Single simulation result
            band_gap: Semiconductor band gap
            
        Returns:
            Matplotlib figure with multiple subplots
        """
        fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(2, 2, figure=fig)
        
        # Panel 1: Potential profile
        ax1 = fig.add_subplot(gs[0, 0])
        if hasattr(result, 'potential_profile') and result.potential_profile:
            try:
                z, pot = result.potential_profile.get_combined_profile()
                # Ensure arrays have same length
                min_len = min(len(z), len(pot))
                z = z[:min_len]
                pot = pot[:min_len]
                ax1.plot(z, pot, 'b-', linewidth=2)
                ax1.set_xlabel('Position z (nm)')
                ax1.set_ylabel('Potential (eV)')
                ax1.set_title('Potential Profile')
                ax1.grid(True, alpha=0.3)
                ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
            except Exception as e:
                # Skip plotting if there's an error
                ax1.text(0.5, 0.5, f'Error plotting profile:\n{str(e)}', 
                        transform=ax1.transAxes, ha='center', va='center')
        
        # Panel 2: Band diagram
        ax2 = fig.add_subplot(gs[0, 1])
        if hasattr(result, 'potential_profile') and result.potential_profile:
            try:
                z, pot = result.potential_profile.get_combined_profile()
                # Ensure arrays have same length
                min_len = min(len(z), len(pot))
                z = z[:min_len]
                pot = pot[:min_len]
                vb = pot - band_gap
                cb = pot
                ax2.plot(z, vb, 'r-', linewidth=2, label='VB')
                ax2.plot(z, cb, 'g-', linewidth=2, label='CB')
                ax2.axhline(y=0, color='k', linestyle=':', linewidth=1, label='EF')
                ax2.set_xlabel('Position z (nm)')
                ax2.set_ylabel('Energy (eV)')
                ax2.set_title('Band Diagram')
                ax2.legend()
                ax2.grid(True, alpha=0.3)
                ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
            except Exception as e:
                ax2.text(0.5, 0.5, f'Error plotting bands:\n{str(e)}', 
                        transform=ax2.transAxes, ha='center', va='center')
        
        # Panel 3: Convergence
        ax3 = fig.add_subplot(gs[1, 0])
        if 'convergence_history' in result.convergence_info:
            history = result.convergence_info['convergence_history']
            ax3.semilogy(history, 'b-', linewidth=2)
            ax3.set_xlabel('Iteration')
            ax3.set_ylabel('Error (V)')
            ax3.set_title('Convergence History')
            ax3.grid(True, alpha=0.3, which='both')
        
        # Panel 4: Information text
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.axis('off')
        info_text = f"""Simulation Results:
        
Bias Voltage: {result.bias_voltage:.3f} V
Band Bending: {result.band_bending:.3f} eV
Depletion Width: {result.depletion_width:.1f} nm
Current: {result.current:.3e} A

Convergence:
Iterations: {result.convergence_info['iterations']}
Final Error: {result.convergence_info['final_error']:.2e} V
Time: {result.convergence_info['time']:.2f} s"""
        
        ax4.text(0.1, 0.9, info_text, transform=ax4.transAxes,
                verticalalignment='top', fontfamily='monospace')
        
        # Overall title
        fig.suptitle(f'SEMITIP Simulation Summary - Bias = {result.bias_voltage:.2f} V',
                    fontsize=14)
        
        plt.tight_layout()
        return fig
    
    def save_all_figures(self, results: List[SimulationResults],
                        output_dir: str = 'figures'):
        """
        Generate and save all standard figures.
        
        Args:
            results: List of simulation results
            output_dir: Directory to save figures
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        # I-V curve
        fig_iv = self.plot_current_voltage(results)
        fig_iv.savefig(f'{output_dir}/current_voltage.png', dpi=self.fig_dpi)
        plt.close(fig_iv)
        
        # Band bending
        fig_bb = self.plot_band_bending(results)
        fig_bb.savefig(f'{output_dir}/band_bending.png', dpi=self.fig_dpi)
        plt.close(fig_bb)
        
        # Summary for first and last bias
        for i, idx in enumerate([0, -1]):
            if idx < len(results):
                fig_summary = self.create_summary_figure(results[idx])
                fig_summary.savefig(f'{output_dir}/summary_bias_{i+1}.png', 
                                  dpi=self.fig_dpi)
                plt.close(fig_summary)
        
        print(f"Figures saved to {output_dir}/")


def plot_simulation_results(results_file: str, output_dir: Optional[str] = None):
    """
    Load and plot simulation results from file.
    
    Args:
        results_file: Path to pickled results file
        output_dir: Directory to save figures (optional)
    """
    # Load results
    with open(results_file, 'rb') as f:
        data = pickle.load(f)
    
    results = data['results']
    config = data.get('config', None)
    
    # Create plotter
    plotter = SEMITIPPlotter()
    
    # Generate plots
    if output_dir:
        plotter.save_all_figures(results, output_dir)
    else:
        # Interactive plots
        fig_iv = plotter.plot_current_voltage(results)
        fig_bb = plotter.plot_band_bending(results)
        
        if results:
            fig_summary = plotter.create_summary_figure(results[0])
        
        plt.show()


if __name__ == "__main__":
    # Example usage
    import sys
    if len(sys.argv) > 1:
        plot_simulation_results(sys.argv[1])