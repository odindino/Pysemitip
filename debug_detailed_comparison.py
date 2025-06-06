#!/usr/bin/env python3
"""
詳細對比Python和Fortran的差異
"""

import numpy as np
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

from src.core.filereader import load_yaml_config
from src.simulation.multint import MultIntSimulation

def debug_detailed_comparison():
    """詳細對比Python和Fortran的參數計算"""
    
    print("="*60)
    print("詳細對比Python和Fortran的參數計算")
    print("="*60)
    
    # Load configuration
    config_path = "data/input/examples/quick_test.yaml"
    config = load_yaml_config(config_path)
    
    # Create simulation
    sim = MultIntSimulation(config)
    
    # 對比1: Hyperboloid參數
    print("\n1. Hyperboloid參數對比:")
    eta, a, z0, c = sim.tip.hyperboloid_parameters()
    print(f"   Python: eta={eta:.8f}, a={a:.8f}, z0={z0:.8e}, c={c:.8e}")
    print(f"   Fortran: eta=0.70710677, a=1.4142135, z0=5.96046448E-08, c=5.96046519E-08")
    
    # 對比2: 基本參數
    print(f"\n2. 基本參數對比:")
    print(f"   Tip radius: {sim.tip.radius}")
    print(f"   Tip slope: {sim.tip.slope}")
    print(f"   Separation: {sim.tip.separation}")
    
    # 對比3: Bias voltage計算
    print(f"\n3. Bias voltage 分析:")
    test_bias = -2.0
    sim.tip.bias_voltage = test_bias
    tip_potential = sim.tip.tip_potential
    print(f"   輸入bias: {test_bias}")
    print(f"   Python tip potential: {tip_potential}")
    print(f"   Fortran tip potential: -2.0707107")
    print(f"   差異: {abs(tip_potential - (-2.0707107)):.8f}")
    
    # 檢查是否有modulation_voltage相關計算
    modulation_v = config.voltage_scan.modulation_voltage
    print(f"   Modulation voltage: {modulation_v}")
    print(f"   測試: {test_bias + modulation_v}")
    print(f"   測試: {test_bias - modulation_v}")
    
    # 對比4: Depletion width計算
    print(f"\n4. Depletion width 分析:")
    depl_width = sim._estimate_depletion_width_1d(test_bias)
    print(f"   Python: {depl_width:.6f} nm")
    print(f"   Fortran: 54.337425 nm")
    print(f"   比例: {depl_width / 54.337425:.3f}")
    
    # 詳細分析depletion width計算步驟
    region = sim.semiconductor_regions[0]
    vbi = sim.fermi_level - region.valence_band_edge()
    total_potential = abs(vbi - test_bias)
    
    print(f"   計算步驟:")
    print(f"     Fermi level: {sim.fermi_level:.6f} eV")
    print(f"     VB edge: {region.valence_band_edge():.6f} eV")
    print(f"     Built-in potential: {vbi:.6f} V")
    print(f"     Total potential: {total_potential:.6f} V")
    print(f"     Net doping: {region.net_doping:.2e} cm^-3")
    print(f"     Permittivity: {region.permittivity}")
    
    # 對比5: Energy range計算
    print(f"\n5. Energy range 分析:")
    energy_range = sim._calculate_energy_range(test_bias)
    print(f"   Python: {energy_range[0]:.6f} to {energy_range[1]:.6f} eV")
    print(f"   Fortran: -6.6036396 to 10.218132 eV")
    print(f"   Range Python: {energy_range[1] - energy_range[0]:.6f}")
    print(f"   Range Fortran: {10.218132 - (-6.6036396):.6f}")
    
    # 對比6: Grid參數
    print(f"\n6. Grid 參數對比:")
    grid = sim.grid
    print(f"   Python: nr={grid.params.nr}, nv={grid.params.nv}, ns={grid.params.ns}, np={grid.params.np}")
    print(f"   Fortran: nr=16, nv=4, ns=16, np=8")
    print(f"   Python delr={grid.params.delr:.6f}, delv={grid.params.delv:.6f}")
    print(f"   Python dels={grid.params.dels:.6f}, delp={grid.params.delp:.6f}")
    print(f"   Fortran delr=0.50000, delv=0.25000, dels=0.50000, delp=0.39270")
    
    # 計算√2因子相關
    print(f"\n7. √2 因子分析:")
    sqrt2 = np.sqrt(2)
    print(f"   √2 = {sqrt2:.8f}")
    print(f"   1/√2 = {1/sqrt2:.8f}")
    print(f"   -2 * √2 / 2 = {-2 * sqrt2 / 2:.8f}")
    print(f"   差異 vs Fortran: {abs(-2 * sqrt2 / 2 - (-2.0707107)):.8f}")
    
    # 檢查是否有調制相關的計算
    print(f"\n8. 可能的調制計算:")
    bias_mod = test_bias * sqrt2 / 2
    print(f"   bias * √2 / 2 = {bias_mod:.8f}")
    print(f"   bias - modulation = {test_bias - modulation_v:.8f}")
    print(f"   bias + modulation = {test_bias + modulation_v:.8f}")

if __name__ == "__main__":
    debug_detailed_comparison()