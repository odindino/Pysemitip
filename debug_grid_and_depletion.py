#!/usr/bin/env python3
"""
調試網格間距和depletion width計算問題
"""

import numpy as np
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

from src.core.filereader import load_yaml_config
from src.simulation.multint import MultIntSimulation
from src.utils.constants import PhysicalConstants as PC

def debug_grid_and_depletion():
    """調試網格和depletion width問題"""
    
    print("="*60)
    print("調試網格間距和depletion width計算問題")
    print("="*60)
    
    # Load configuration
    config_path = "data/input/examples/quick_test.yaml"
    config = load_yaml_config(config_path)
    
    # Create simulation
    sim = MultIntSimulation(config)
    
    # 檢查網格生成
    print("\n1. 網格參數詳細分析:")
    grid = sim.grid
    print(f"   Grid dimensions: nr={grid.params.nr}, nv={grid.params.nv}, ns={grid.params.ns}, np={grid.params.np}")
    print(f"   Grid extents: rmax={grid.params.rmax}, vmax={grid.params.vmax}, smax={grid.params.smax}")
    print(f"   Grid spacing: delr={grid.params.delr:.6f}, delv={grid.params.delv:.6f}")
    print(f"                 dels={grid.params.dels:.6f}, delp={grid.params.delp:.6f}")
    
    # 比較與Fortran預期值
    print(f"\n   與Fortran比較:")
    print(f"   delr: {grid.params.delr:.6f} vs 0.50000 (比例: {grid.params.delr/0.5:.3f})")
    print(f"   delv: {grid.params.delv:.6f} vs 0.25000 (比例: {grid.params.delv/0.25:.3f})")
    print(f"   dels: {grid.params.dels:.6f} vs 0.50000 (比例: {grid.params.dels/0.5:.3f})")
    print(f"   delp: {grid.params.delp:.6f} vs 0.39270 (比例: {grid.params.delp/0.393:.3f})")
    
    # 檢查網格點實際值
    print(f"\n2. 網格點分佈:")
    print(f"   R values (first 5): {grid.r[:5]}")
    print(f"   R values (last 5): {grid.r[-5:]}")
    if hasattr(grid, 'zv'):
        print(f"   Z vacuum (first 5): {grid.zv[:5]}")
    if hasattr(grid, 'zs'):
        print(f"   Z semiconductor (first 5): {grid.zs[:5]}")
    
    # 詳細分析depletion width計算
    print(f"\n3. Depletion width 詳細計算:")
    test_bias = -2.0
    
    # 手動計算步驟
    region = sim.semiconductor_regions[0]
    fermi_level = sim.fermi_level
    
    print(f"   輸入參數:")
    print(f"     Bias voltage: {test_bias} V")
    print(f"     Fermi level: {fermi_level:.6f} eV")
    print(f"     Semiconductor parameters:")
    print(f"       Band gap: {region.band_gap} eV")
    print(f"       VB offset: {region.valence_band_offset} eV")
    print(f"       CB offset: derived from band gap")
    print(f"       Permittivity: {region.permittivity}")
    print(f"       Donor concentration: {region.donor_concentration:.2e} cm^-3")
    print(f"       Acceptor concentration: {region.acceptor_concentration:.2e} cm^-3")
    print(f"       Net doping: {region.net_doping:.2e} cm^-3")
    
    # 計算band edges
    vb_edge = region.valence_band_edge()
    cb_edge = region.conduction_band_edge()
    print(f"     Band edges:")
    print(f"       VB edge: {vb_edge:.6f} eV")
    print(f"       CB edge: {cb_edge:.6f} eV")
    
    # Built-in potential計算
    vbi = fermi_level - vb_edge
    print(f"     Built-in potential: {vbi:.6f} V")
    
    # Total potential
    total_potential = abs(vbi - test_bias)
    print(f"     Total potential: {total_potential:.6f} V")
    
    # 物理常數
    print(f"   物理常數:")
    print(f"     Electron charge: {PC.E:.6e} C")
    print(f"     Vacuum permittivity: {PC.EPSILON0:.6e} F/m")
    
    # Depletion width公式逐步計算
    eps = region.permittivity * PC.EPSILON0
    doping = abs(region.net_doping) * 1e6  # Convert to m^-3
    
    print(f"   Depletion width計算:")
    print(f"     Permittivity: {region.permittivity}")
    print(f"     ε = ε_r × ε_0 = {eps:.6e} F/m")
    print(f"     Doping: {abs(region.net_doping):.2e} cm^-3 = {doping:.2e} m^-3")
    
    # 公式: w = sqrt(2 * ε * V / (e * N))
    numerator = 2 * eps * total_potential
    denominator = PC.E * doping
    w_squared = numerator / denominator
    w_meters = np.sqrt(w_squared)
    w_nm = w_meters * 1e9
    
    print(f"     分子: 2εV = {numerator:.6e}")
    print(f"     分母: eN = {denominator:.6e}")
    print(f"     w² = {w_squared:.6e} m²")
    print(f"     w = {w_meters:.6e} m = {w_nm:.6f} nm")
    
    # 與Fortran比較
    fortran_w = 54.337425
    print(f"   結果比較:")
    print(f"     Python: {w_nm:.6f} nm")
    print(f"     Fortran: {fortran_w:.6f} nm")
    print(f"     比例: {w_nm/fortran_w:.3f}")
    print(f"     差異: {abs(w_nm - fortran_w):.3f} nm")
    
    # 檢查可能的單位或公式問題
    print(f"\n4. 可能的修正:")
    
    # 試驗不同的公式變形
    test_corrections = [
        ("Original", w_nm),
        ("除以√2", w_nm / np.sqrt(2)),
        ("除以2", w_nm / 2),
        ("乘以√2", w_nm * np.sqrt(2)),
        ("使用不同VBI", np.sqrt(2 * eps * abs(fermi_level - test_bias) / (PC.E * doping)) * 1e9),
        ("只用bias", np.sqrt(2 * eps * abs(test_bias) / (PC.E * doping)) * 1e9),
    ]
    
    for name, value in test_corrections:
        ratio = value / fortran_w
        print(f"     {name}: {value:.3f} nm (比例: {ratio:.3f})")

if __name__ == "__main__":
    debug_grid_and_depletion()