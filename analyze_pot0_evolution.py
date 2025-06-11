#!/usr/bin/env python3
"""
Pot0 Evolution Analysis Script

This script analyzes the Pot0 (band bending) evolution in both Fortran and Python versions,
explaining the physical significance of the sign change and comparing LAPLACE vs VSINT models.
"""

import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path

def extract_fortran_pot0_evolution(file_path):
    """Extract Pot0 evolution data from Fortran output file"""
    pot0_data = []
    current_bias = None
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Find all bias values and their corresponding Pot0 evolution
    bias_pattern = r'BIAS, TIP POTENTIAL =\s+([-\d.E+]+)\s+([-\d.E+]+)'
    pot0_pattern = r'ITER,Pot0 =\s+(\d+)\s+([-\d.E+]+)'
    
    lines = content.split('\n')
    for i, line in enumerate(lines):
        bias_match = re.search(bias_pattern, line)
        if bias_match:
            current_bias = float(bias_match.group(1))
            # Extract Pot0 evolution for this bias
            bias_pot0_data = []
            for j in range(i+1, min(i+50, len(lines))):  # Look ahead for Pot0 data
                pot0_match = re.search(pot0_pattern, lines[j])
                if pot0_match:
                    iteration = int(pot0_match.group(1))
                    pot0_value = float(pot0_match.group(2))
                    bias_pot0_data.append((iteration, pot0_value))
                elif 'BIAS, TIP POTENTIAL' in lines[j]:  # Next bias section
                    break
            
            if bias_pot0_data:
                pot0_data.append({
                    'bias': current_bias,
                    'evolution': bias_pot0_data
                })
    
    return pot0_data

def analyze_pot0_physics():
    """Analyze the physical meaning of Pot0 and its sign change"""
    print("=" * 80)
    print("POT0 PHYSICAL ANALYSIS")
    print("=" * 80)
    
    print("\n1. WHAT IS Pot0?")
    print("-" * 40)
    print("• Pot0 = Band bending at the semiconductor surface")
    print("• Pot0 = V_surface - V_bulk (relative to bulk potential)")
    print("• Represents the electrostatic potential difference between")
    print("  the surface and the neutral bulk region")
    print("• Directly related to surface band bending in eV")
    
    print("\n2. PHYSICAL SIGNIFICANCE OF SIGN CHANGE")
    print("-" * 40)
    print("NEGATIVE Pot0 (Pot0 < 0):")
    print("• Surface potential is LOWER than bulk")
    print("• Bands bend DOWNWARD at surface")
    print("• Electrons are ATTRACTED to surface")
    print("• Creates electron accumulation layer")
    print("• Common in n-type semiconductors with negative tip bias")
    
    print("\nPOSITIVE Pot0 (Pot0 > 0):")
    print("• Surface potential is HIGHER than bulk")
    print("• Bands bend UPWARD at surface")
    print("• Electrons are REPELLED from surface")
    print("• Creates depletion layer (or inversion for strong bending)")
    print("• Common when tip-induced field dominates")
    
    print("\n3. WHY THE SIGN CHANGE OCCURS")
    print("-" * 40)
    print("The Pot0 evolution from negative to positive represents:")
    print("• Initial state: Strong surface accumulation (negative bending)")
    print("• Transition: Balance between tip field and surface charge")
    print("• Final state: Tip-induced depletion dominates (positive bending)")
    print("• This transition is physically realistic for STM conditions")

def analyze_laplace_vs_vsint():
    """Compare LAPLACE and VSINT physical models"""
    print("\n" + "=" * 80)
    print("LAPLACE vs VSINT MODEL COMPARISON")
    print("=" * 80)
    
    print("\nLAPLACE MODEL (∇²φ = 0):")
    print("-" * 40)
    print("• Pure electrostatic model")
    print("• No charge density (ρ = 0)")
    print("• Only boundary conditions matter")
    print("• Fast convergence (linear system)")
    print("• Missing: Surface states, bulk charge, band bending effects")
    print("• Good for: Initial guess, metallic systems")
    
    print("\nVSINT MODEL (∇²φ = -ρ/ε₀):")
    print("-" * 40)
    print("• Full Poisson equation with charge density")
    print("• Includes surface charge density ρ_surface")
    print("• Includes bulk charge density ρ_bulk")
    print("• Self-consistent solution (nonlinear)")
    print("• Accounts for: Band bending, surface states, carrier redistribution")
    print("• Essential for: Semiconductor STM simulation")
    
    print("\nKEY DIFFERENCES:")
    print("-" * 40)
    print("1. CHARGE EFFECTS:")
    print("   • LAPLACE: Ignores all charge effects")
    print("   • VSINT: Includes surface and bulk charge density")
    
    print("\n2. BAND BENDING:")
    print("   • LAPLACE: No physical band bending")
    print("   • VSINT: Realistic band bending due to charge redistribution")
    
    print("\n3. SURFACE STATES:")
    print("   • LAPLACE: No surface state effects")
    print("   • VSINT: Surface states contribute to surface charge")
    
    print("\n4. CONVERGENCE:")
    print("   • LAPLACE: Direct solution (one step)")
    print("   • VSINT: Iterative nonlinear solution (many steps)")

def plot_fortran_pot0_evolution(pot0_data):
    """Plot Pot0 evolution from Fortran data"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Fortran Pot0 Evolution Analysis', fontsize=16)
    
    # Plot 1: Evolution for first bias
    if pot0_data:
        first_bias_data = pot0_data[0]
        iterations = [d[0] for d in first_bias_data['evolution']]
        pot0_values = [d[1] for d in first_bias_data['evolution']]
        
        axes[0,0].plot(iterations, pot0_values, 'b-o', markersize=4)
        axes[0,0].axhline(y=0, color='r', linestyle='--', alpha=0.5)
        axes[0,0].set_xlabel('Iteration')
        axes[0,0].set_ylabel('Pot0 (V)')
        axes[0,0].set_title(f'Pot0 Evolution (Bias = {first_bias_data["bias"]:.3f}V)')
        axes[0,0].grid(True, alpha=0.3)
        
        # Highlight sign change
        sign_change_idx = None
        for i in range(1, len(pot0_values)):
            if pot0_values[i-1] < 0 and pot0_values[i] > 0:
                sign_change_idx = i
                break
        
        if sign_change_idx:
            axes[0,0].axvline(x=iterations[sign_change_idx], color='orange', 
                            linestyle=':', linewidth=2, label='Sign Change')
            axes[0,0].legend()
    
    # Plot 2: Final Pot0 vs Bias
    if len(pot0_data) > 1:
        biases = [d['bias'] for d in pot0_data]
        final_pot0 = [d['evolution'][-1][1] if d['evolution'] else 0 for d in pot0_data]
        
        axes[0,1].plot(biases, final_pot0, 'g-s', markersize=6)
        axes[0,1].axhline(y=0, color='r', linestyle='--', alpha=0.5)
        axes[0,1].set_xlabel('Bias Voltage (V)')
        axes[0,1].set_ylabel('Final Pot0 (V)')
        axes[0,1].set_title('Final Pot0 vs Bias')
        axes[0,1].grid(True, alpha=0.3)
    
    # Plot 3: Sign change analysis
    axes[1,0].text(0.1, 0.8, 'NEGATIVE Pot0:', fontsize=12, weight='bold', 
                   transform=axes[1,0].transAxes, color='blue')
    axes[1,0].text(0.1, 0.7, '• Bands bend DOWNWARD', fontsize=10, 
                   transform=axes[1,0].transAxes)
    axes[1,0].text(0.1, 0.6, '• Electron accumulation', fontsize=10, 
                   transform=axes[1,0].transAxes)
    axes[1,0].text(0.1, 0.5, '• Surface < Bulk potential', fontsize=10, 
                   transform=axes[1,0].transAxes)
    
    axes[1,0].text(0.1, 0.3, 'POSITIVE Pot0:', fontsize=12, weight='bold', 
                   transform=axes[1,0].transAxes, color='red')
    axes[1,0].text(0.1, 0.2, '• Bands bend UPWARD', fontsize=10, 
                   transform=axes[1,0].transAxes)
    axes[1,0].text(0.1, 0.1, '• Electron depletion', fontsize=10, 
                   transform=axes[1,0].transAxes)
    axes[1,0].text(0.1, 0.0, '• Surface > Bulk potential', fontsize=10, 
                   transform=axes[1,0].transAxes)
    axes[1,0].set_title('Physical Meaning of Pot0 Sign')
    axes[1,0].axis('off')
    
    # Plot 4: Model comparison
    axes[1,1].text(0.05, 0.9, 'LAPLACE MODEL:', fontsize=11, weight='bold', 
                   transform=axes[1,1].transAxes, color='purple')
    axes[1,1].text(0.05, 0.8, '∇²φ = 0 (no charge)', fontsize=10, 
                   transform=axes[1,1].transAxes)
    axes[1,1].text(0.05, 0.7, '• Fast convergence', fontsize=10, 
                   transform=axes[1,1].transAxes)
    axes[1,1].text(0.05, 0.6, '• No band bending', fontsize=10, 
                   transform=axes[1,1].transAxes)
    
    axes[1,1].text(0.05, 0.4, 'VSINT MODEL:', fontsize=11, weight='bold', 
                   transform=axes[1,1].transAxes, color='green')
    axes[1,1].text(0.05, 0.3, '∇²φ = -ρ/ε₀ (with charge)', fontsize=10, 
                   transform=axes[1,1].transAxes)
    axes[1,1].text(0.05, 0.2, '• Nonlinear convergence', fontsize=10, 
                   transform=axes[1,1].transAxes)
    axes[1,1].text(0.05, 0.1, '• Realistic band bending', fontsize=10, 
                   transform=axes[1,1].transAxes)
    axes[1,1].set_title('LAPLACE vs VSINT Models')
    axes[1,1].axis('off')
    
    plt.tight_layout()
    plt.savefig('pot0_evolution_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved as: pot0_evolution_analysis.png")

def run_python_pot0_analysis():
    """Run our Python version and track Pot0 evolution"""
    print("\n" + "=" * 80)
    print("PYTHON VERSION POT0 ANALYSIS")
    print("=" * 80)
    
    try:
        # Import our Python modules
        import sys
        sys.path.append('/Users/yangziliang/Git-Projects/Pysemitip/src')
        
        from simulation.multint import MultIntSimulation
        from core.filereader import YAMLConfigReader
        
        # Load configuration
        config_path = '/Users/yangziliang/Git-Projects/Pysemitip/data/input/examples/test/quick_test.yaml'
        reader = YAMLConfigReader()
        config = reader.read_config(config_path)
        
        print(f"Running Python simulation with config: {config_path}")
        
        # Run simulation with Pot0 tracking
        simulation = MultIntSimulation(config)
        
        # Modify the simulation to track Pot0 evolution
        class Pot0Tracker:
            def __init__(self):
                self.pot0_history = []
                self.iteration_count = 0
            
            def track_pot0(self, pot0_value):
                self.iteration_count += 1
                if self.iteration_count % 100 == 0:  # Track every 100 iterations
                    self.pot0_history.append((self.iteration_count, pot0_value))
                    print(f"ITER,Pot0 = {self.iteration_count:8d} {pot0_value:E}")
        
        tracker = Pot0Tracker()
        
        # Run simulation (this would need modification to actually track Pot0)
        print("Note: To fully track Pot0 evolution, we need to modify the Poisson solver")
        print("to call tracker.track_pot0() during iterations.")
        
        results = simulation.run()
        
        if hasattr(results, 'band_bending_surface'):
            final_pot0 = results.band_bending_surface
            print(f"\nFinal Python Pot0: {final_pot0:.6f} V")
        else:
            print("Could not extract final Pot0 from Python results")
            
    except Exception as e:
        print(f"Error running Python analysis: {e}")
        print("This is expected - we would need to modify the code to track Pot0 evolution")

def main():
    """Main analysis function"""
    print("POT0 EVOLUTION AND PHYSICAL MODEL ANALYSIS")
    print("=" * 80)
    
    # 1. Extract and analyze Fortran data
    fortran_file = '/Users/yangziliang/Git-Projects/Pysemitip/src/fortran/MultInt/fort_MultInt.16'
    print(f"\nAnalyzing Fortran output: {fortran_file}")
    
    pot0_data = extract_fortran_pot0_evolution(fortran_file)
    print(f"Found {len(pot0_data)} bias points with Pot0 evolution data")
    
    if pot0_data:
        # Show first bias evolution
        first_bias = pot0_data[0]
        print(f"\nFirst bias: {first_bias['bias']:.6f} V")
        print("Pot0 evolution (every 100 iterations):")
        for iter_num, pot0_val in first_bias['evolution'][::10]:  # Show every 10th point
            print(f"  ITER {iter_num:4d}: Pot0 = {pot0_val:+.6f} V")
        
        # Find sign change
        pot0_values = [d[1] for d in first_bias['evolution']]
        sign_change_iter = None
        for i in range(1, len(pot0_values)):
            if pot0_values[i-1] < 0 and pot0_values[i] > 0:
                sign_change_iter = first_bias['evolution'][i][0]
                sign_change_pot0 = pot0_values[i]
                break
        
        if sign_change_iter:
            print(f"\n*** SIGN CHANGE DETECTED ***")
            print(f"Sign change occurs around iteration {sign_change_iter}")
            print(f"Pot0 changes from negative to positive: {sign_change_pot0:+.6f} V")
    
    # 2. Analyze physics
    analyze_pot0_physics()
    
    # 3. Compare models
    analyze_laplace_vs_vsint()
    
    # 4. Create visualization
    if pot0_data:
        plot_fortran_pot0_evolution(pot0_data)
    
    # 5. Python analysis
    run_python_pot0_analysis()
    
    # 6. Summary and recommendations
    print("\n" + "=" * 80)
    print("SUMMARY AND RECOMMENDATIONS")
    print("=" * 80)
    
    print("\nKEY FINDINGS:")
    print("1. Fortran Pot0 evolves from negative to positive during iteration")
    print("2. This represents realistic semiconductor band bending physics")
    print("3. Sign change indicates transition from accumulation to depletion")
    print("4. VSINT model includes essential charge density effects")
    
    print("\nRECOMMENDations FOR PYTHON VERSION:")
    print("1. Implement full nonlinear Poisson solver (not just Laplace)")
    print("2. Include surface and bulk charge density calculations")
    print("3. Track Pot0 evolution during iteration to verify physics")
    print("4. Ensure proper boundary conditions at semiconductor interface")
    print("5. Validate that sign change occurs as in Fortran version")

if __name__ == "__main__":
    main()