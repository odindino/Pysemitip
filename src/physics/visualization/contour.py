"""
Contour plotting for SEMITIP simulations.

This module implements contour plot generation similar to CONTR3 in the
Fortran code, creating publication-quality potential contour plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.path import Path
import matplotlib.patheffects as path_effects
from typing import Tuple, Optional, List

from ..solvers.grid import Grid3D
from ..materials.tip import TipModel


class ContourPlotter:
    """
    Generate contour plots of potential distributions.
    
    Implements functionality similar to CONTR3 routine.
    """
    
    def __init__(self, grid: Grid3D, tip: TipModel):
        """
        Initialize contour plotter.
        
        Args:
            grid: Computational grid
            tip: Tip model for geometry
        """
        self.grid = grid
        self.tip = tip
    
    def plot_potential_contours(self, potential_3d: np.ndarray,
                              r_max: Optional[float] = None,
                              z_range: Optional[Tuple[float, float]] = None,
                              n_contours: int = 20,
                              potential_range: Optional[Tuple[float, float]] = None,
                              phi_slice: int = 0,
                              show_tip: bool = True,
                              save_path: Optional[str] = None) -> plt.Figure:
        """
        Create contour plot of potential distribution.
        
        Args:
            potential_3d: 3D potential array
            r_max: Maximum radius to plot (nm)
            z_range: (z_min, z_max) range to plot (nm)
            n_contours: Number of contour levels
            potential_range: (min, max) potential values for contours
            phi_slice: Angular slice to plot
            show_tip: Whether to show tip outline
            save_path: Path to save figure
            
        Returns:
            Matplotlib figure
        """
        # Set plot limits
        if r_max is None:
            r_max = self.grid.params.rmax / 2  # Default to half grid
        if z_range is None:
            z_range = (-self.grid.params.smax / 2, self.grid.params.vmax / 2)
        
        # Extract 2D slice
        potential_2d = potential_3d[:, :, phi_slice]
        
        # Create coordinate arrays
        r = self.grid.r
        z = self.grid.z
        
        # Limit to plot range
        r_mask = r <= r_max
        z_mask = (z >= z_range[0]) & (z <= z_range[1])
        
        r_plot = r[r_mask]
        z_plot = z[z_mask]
        potential_plot = potential_2d[np.ix_(r_mask, z_mask)]
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Determine contour levels
        if potential_range is None:
            v_min, v_max = potential_plot.min(), potential_plot.max()
        else:
            v_min, v_max = potential_range
        
        # Ensure v_min < v_max and levels are valid
        if v_max <= v_min:
            v_range = max(abs(v_min), abs(v_max), 1e-6)  # Fallback to small range
            v_min = -v_range
            v_max = v_range
        
        levels = np.linspace(v_min, v_max, n_contours)
        
        # Ensure levels are strictly increasing (remove duplicates)
        levels = np.unique(levels)
        if len(levels) < 3:  # Need at least 3 levels for contour
            v_center = (v_min + v_max) / 2
            v_range = max(abs(v_max - v_min), 1e-6)
            levels = np.linspace(v_center - v_range/2, v_center + v_range/2, n_contours)
        
        # Create meshgrid
        R, Z = np.meshgrid(r_plot, z_plot, indexing='ij')
        
        # Plot filled contours
        cs = ax.contourf(R, Z, potential_plot, levels=levels, 
                        cmap='RdBu_r', extend='both')
        
        # Add contour lines
        contour_lines = ax.contour(R, Z, potential_plot, levels=levels,
                                  colors='black', linewidths=0.5, alpha=0.5)
        
        # Label contours
        ax.clabel(contour_lines, inline=True, fontsize=8, fmt='%.3f')
        
        # Add colorbar
        cbar = plt.colorbar(cs, ax=ax, label='Potential (V)', 
                          shrink=0.9, pad=0.02)
        
        # Draw surface
        ax.axhline(y=0, color='black', linewidth=2)
        ax.fill_between(r_plot, -100, 0, color='lightgray', alpha=0.3)
        
        # Draw tip if requested
        if show_tip:
            self._draw_tip_outline(ax, r_plot, z_plot)
        
        # Labels and formatting
        ax.set_xlabel('Radial position r (nm)', fontsize=12)
        ax.set_ylabel('Vertical position z (nm)', fontsize=12)
        ax.set_title('Potential Distribution Contours', fontsize=14)
        
        # Set aspect ratio
        ax.set_aspect('equal')
        
        # Grid
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Annotations
        ax.text(0.02, 0.98, 'Vacuum', transform=ax.transAxes,
               verticalalignment='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax.text(0.02, 0.02, 'Semiconductor', transform=ax.transAxes,
               verticalalignment='bottom', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.tight_layout()
        
        # Save if requested
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Contour plot saved to {save_path}")
        
        return fig
    
    def _draw_tip_outline(self, ax, r_range: np.ndarray, z_range: np.ndarray):
        """Draw tip outline on the plot."""
        # Get tip parameters
        eta, a, z0, c = self.tip.hyperboloid_parameters()
        
        # Tip apex position
        z_apex = self.tip.separation + self.tip.protrusion_radius
        
        # Draw hyperbolic tip shape
        z_tip = np.linspace(z_apex, z_range[-1], 100)
        r_tip = np.zeros_like(z_tip)
        
        for i, z in enumerate(z_tip):
            if z > z_apex:
                # Hyperboloid equation
                r_tip[i] = a * np.sqrt((z - z0)**2 / a**2 - 1)
        
        # Account for tip position
        r_tip_shifted = r_tip + abs(self.tip.position[0])
        
        # Draw tip outline
        ax.plot(r_tip_shifted, z_tip, 'k-', linewidth=2, label='Tip')
        ax.plot(-r_tip_shifted, z_tip, 'k-', linewidth=2)
        
        # Fill tip region
        ax.fill_betweenx(z_tip, -r_tip_shifted, r_tip_shifted,
                        color='gray', alpha=0.5)
        
        # Draw tip apex
        ax.plot(0, z_apex, 'ko', markersize=6)
    
    def plot_field_magnitude(self, potential_3d: np.ndarray,
                           phi_slice: int = 0,
                           log_scale: bool = True) -> plt.Figure:
        """
        Plot electric field magnitude.
        
        Args:
            potential_3d: 3D potential array
            phi_slice: Angular slice to plot
            log_scale: Use logarithmic color scale
            
        Returns:
            Matplotlib figure
        """
        from ..core.potential import PotentialProcessor
        
        # Calculate electric field
        processor = PotentialProcessor(self.grid)
        E_magnitude = processor.calculate_field_magnitude(potential_3d)
        
        # Extract 2D slice
        E_2d = E_magnitude[:, :, phi_slice]
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create meshgrid
        R, Z = np.meshgrid(self.grid.r, self.grid.z, indexing='ij')
        
        # Plot field magnitude
        if log_scale:
            # Use logarithmic scale
            E_plot = np.log10(E_2d + 1e-10)  # Avoid log(0)
            im = ax.pcolormesh(R, Z, E_plot, cmap='hot', shading='auto')
            cbar = plt.colorbar(im, ax=ax, label='log₁₀(|E| / V/nm)')
        else:
            im = ax.pcolormesh(R, Z, E_2d, cmap='hot', shading='auto')
            cbar = plt.colorbar(im, ax=ax, label='|E| (V/nm)')
        
        # Draw surface
        ax.axhline(y=0, color='white', linewidth=2)
        
        # Labels
        ax.set_xlabel('Radial position r (nm)')
        ax.set_ylabel('Vertical position z (nm)')
        ax.set_title('Electric Field Magnitude')
        
        plt.tight_layout()
        return fig
    
    def create_publication_figure(self, potential_3d: np.ndarray,
                                band_bending: float,
                                current: float,
                                bias: float) -> plt.Figure:
        """
        Create publication-quality figure with contours and annotations.
        
        Args:
            potential_3d: 3D potential array
            band_bending: Band bending value (eV)
            current: Tunneling current (A)
            bias: Bias voltage (V)
            
        Returns:
            Matplotlib figure
        """
        # Create figure with specific size for publication
        fig = plt.figure(figsize=(8, 10))
        
        # Main contour plot
        ax_main = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
        
        # Extract potential slice
        potential_2d = potential_3d[:, :, 0]
        
        # Plot contours
        r_max = 30  # nm
        z_range = (-30, 20)  # nm
        
        r_mask = self.grid.r <= r_max
        z_mask = (self.grid.z >= z_range[0]) & (self.grid.z <= z_range[1])
        
        r_plot = self.grid.r[r_mask]
        z_plot = self.grid.z[z_mask]
        potential_plot = potential_2d[np.ix_(r_mask, z_mask)]
        
        R, Z = np.meshgrid(r_plot, z_plot, indexing='ij')
        
        # Contour levels
        v_min, v_max = potential_plot.min(), potential_plot.max()
        
        # Ensure valid range
        if v_max <= v_min:
            v_range = max(abs(v_min), abs(v_max), 1e-6)
            v_min = -v_range
            v_max = v_range
        
        levels = np.linspace(v_min, v_max, 20)
        levels = np.unique(levels)  # Remove duplicates
        
        if len(levels) < 3:
            v_center = (v_min + v_max) / 2
            v_range = max(abs(v_max - v_min), 1e-6)
            levels = np.linspace(v_center - v_range/2, v_center + v_range/2, 20)
        
        cs = ax_main.contourf(R, Z, potential_plot, levels=levels,
                            cmap='RdBu_r', extend='both')
        contours = ax_main.contour(R, Z, potential_plot, levels=levels,
                                  colors='black', linewidths=0.5, alpha=0.5)
        
        # Surface and tip
        ax_main.axhline(y=0, color='black', linewidth=2)
        self._draw_tip_outline(ax_main, r_plot, z_plot)
        
        # Labels
        ax_main.set_xlabel('r (nm)', fontsize=12)
        ax_main.set_ylabel('z (nm)', fontsize=12)
        ax_main.set_title(f'V = {bias:.2f} V, BB = {band_bending:.3f} eV',
                         fontsize=14)
        
        # Colorbar
        cbar = plt.colorbar(cs, ax=ax_main, label='Potential (V)')
        
        # Bottom panel - 1D profile
        ax_profile = plt.subplot2grid((4, 1), (3, 0))
        
        # Extract profile at r=0
        r0_idx = 0
        z_profile = self.grid.z[z_mask]
        pot_profile = potential_2d[r0_idx, z_mask]
        
        ax_profile.plot(z_profile, pot_profile, 'b-', linewidth=2)
        ax_profile.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax_profile.set_xlabel('z (nm)')
        ax_profile.set_ylabel('Potential (V)')
        ax_profile.set_xlim(z_range)
        ax_profile.grid(True, alpha=0.3)
        
        # Add current annotation
        current_text = f'I = {current:.2e} A'
        ax_main.text(0.95, 0.95, current_text, transform=ax_main.transAxes,
                    verticalalignment='top', horizontalalignment='right',
                    fontsize=12, bbox=dict(boxstyle='round', 
                                         facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        return fig


def create_contour_plots(results_file: str, output_dir: str = 'contours'):
    """
    Create contour plots from simulation results.
    
    Args:
        results_file: Path to simulation results
        output_dir: Directory to save contour plots
    """
    import pickle
    import os
    
    # Load results
    with open(results_file, 'rb') as f:
        data = pickle.load(f)
    
    results = data['results']
    config = data['config']
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create grid and tip from config
    from ..materials import create_tip_from_config
    from ..solvers import create_grid_from_config
    
    tip = create_tip_from_config(config)
    grid = create_grid_from_config(config, tip)
    
    # Create contour plotter
    plotter = ContourPlotter(grid, tip)
    
    # Plot contours for each bias
    for i, result in enumerate(results):
        # Standard contour plot
        fig = plotter.plot_potential_contours(
            result.potential_3d,
            n_contours=20,
            save_path=f'{output_dir}/contour_bias_{i+1}.png'
        )
        plt.close(fig)
        
        # Publication figure
        fig_pub = plotter.create_publication_figure(
            result.potential_3d,
            result.band_bending,
            result.current,
            result.bias_voltage
        )
        fig_pub.savefig(f'{output_dir}/publication_bias_{i+1}.png',
                       dpi=300, bbox_inches='tight')
        plt.close(fig_pub)
    
    print(f"Contour plots saved to {output_dir}/")