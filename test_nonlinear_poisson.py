#!/usr/bin/env python3
"""
快速測試非線性 Poisson 求解器是否正確實現
"""
import sys
import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid
from src.physics.core.charge_density import ChargeDensityCalculator

# 設置日誌
logging.basicConfig(level=logging.INFO, format='%(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_nonlinear_poisson():
    """測試非線性 Poisson 求解器"""
    print("Testing nonlinear Poisson solver with GSECT method...")
    
    try:
        # 創建小網格以加快測試
        grid = HyperbolicGrid(N_eta=8, N_nu=6, R=1.0, Z_TS=1.0, shank_slope=1.0)
        
        # 創建物理參數
        class MinimalProps:
            def __init__(self):
                class SemiconductorProps:
                    epsilon_r = 12.9
                    Ev_offset_eV = -5.17
                self.semiconductor_props = SemiconductorProps()
        
        props = MinimalProps()
        
        # 創建求解器
        solver = PoissonSOREquation(grid, props)
        
        # 創建模擬電荷密度計算器
        class MockChargeDensityCalculator:
            def get_charge_density_C_m3(self, ef_rel_vb_eV):
                # 模擬非線性電荷密度響應
                # 在費米能級附近有強烈的非線性
                if ef_rel_vb_eV > 0.5:
                    return 1e16 * np.tanh(ef_rel_vb_eV - 0.5)  # 強非線性
                elif ef_rel_vb_eV < -0.5:
                    return -1e16 * np.tanh(-ef_rel_vb_eV - 0.5)  # 強非線性
                else:
                    return 1e12 * ef_rel_vb_eV  # 弱線性
        
        charge_calc = MockChargeDensityCalculator()
        
        # 測試參數
        V_tip = -2.0
        V_sample = 0.0
        system_fermi = 1.42  # eV
        
        print(f"Grid size: {grid.N_eta} x {grid.N_nu}")
        print(f"Test conditions: V_tip={V_tip}V, V_sample={V_sample}V, EF={system_fermi}eV")
        
        # 運行非線性求解 (測試收斂性)
        potential, iterations, converged = solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=charge_calc,
            system_fermi_level_E_F_main_eV=system_fermi,
            max_iterations=800,  # 增加迭代次數測試收斂
            tolerance_Volts=5e-4,  # 更嚴格的容差
            omega=0.7  # 更保守的 omega
        )
        
        print(f"✅ Test completed:")
        print(f"  Iterations: {iterations}")
        print(f"  Converged: {converged}")
        
        # 計算 Pot0
        pot0 = solver._calculate_pot0_fortran_style(potential)
        print(f"  Pot0 (band bending): {pot0:.6f} V")
        
        # 檢查數值穩定性
        if np.any(np.isnan(potential)) or np.any(np.isinf(potential)):
            print("❌ FAILED: NaN or Inf values detected")
            return False
        
        if np.max(np.abs(potential)) > 50:
            print("❌ FAILED: Potential values too large (>50V)")
            return False
        
        # 檢查電位變化
        potential_range = np.max(potential) - np.min(potential)
        print(f"  Potential range: [{np.min(potential):.3f}, {np.max(potential):.3f}] V (span: {potential_range:.3f}V)")
        
        if potential_range < 0.1:
            print("⚠️  WARNING: Very small potential variation")
        elif potential_range > 20:
            print("⚠️  WARNING: Very large potential variation")
        else:
            print("✅ PASSED: Reasonable potential variation")
        
        # 檢查 Pot0 合理性
        if abs(pot0) < 0.01:
            print("⚠️  WARNING: Pot0 is very small, may indicate no band bending")
        elif abs(pot0) > 10:
            print("⚠️  WARNING: Pot0 magnitude is very large")
        else:
            print("✅ PASSED: Reasonable Pot0 magnitude")
            
        return True
        
    except Exception as e:
        print(f"❌ FAILED: Exception occurred: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_nonlinear_poisson()
    print(f"\n{'='*50}")
    if success:
        print("✅ Nonlinear Poisson solver test PASSED")
    else:
        print("❌ Nonlinear Poisson solver test FAILED")
    print(f"{'='*50}")
    sys.exit(0 if success else 1)