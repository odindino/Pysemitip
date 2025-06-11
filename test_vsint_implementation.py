#!/usr/bin/env python3
"""
測試新的VSINT實現對Pot0計算的改善
"""
import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid
from src.utils.constants import PhysicalConstants as PC

logging.basicConfig(level=logging.INFO, format='%(message)s')

def test_vsint_implementation():
    """測試VSINT實現對Pot0的改善效果"""
    print("測試新的VSINT實現")
    print("="*80)
    
    # 創建測試環境
    grid = HyperbolicGrid(N_eta=16, N_nu=8, R=1.0, Z_TS=1.0, shank_slope=1.0)
    
    class Props:
        def __init__(self):
            class SemiconductorProps:
                epsilon_r = 12.9
                Ev_offset_eV = -5.17
            self.semiconductor_props = SemiconductorProps()
    
    props = Props()
    solver = PoissonSOREquation(grid, props)
    
    # 測試參數
    V_tip = -2.07
    V_sample = 0.0
    system_fermi_level = 1.4187
    
    print(f"測試條件:")
    print(f"  V_tip = {V_tip} V")
    print(f"  V_sample = {V_sample} V")
    print(f"  System Fermi level = {system_fermi_level} eV")
    print(f"  Grid: N_eta={grid.N_eta}, N_nu={grid.N_nu}")
    
    # 創建簡化的電荷密度計算器 (用於測試)
    class MockChargeDensityCalculator:
        def __init__(self):
            self.system_fermi_level = system_fermi_level
            
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            """簡化的電荷密度計算 - 返回合理的值用於測試"""
            # 使用簡單的指數函數模擬費米狄拉克分佈
            kT = 0.0259  # 300K時的熱能 (eV)
            
            # 電子密度 (負電荷)
            n_electrons = 1e17 / (1 + np.exp(-ef_rel_vb_eV / kT))
            
            # 離子化施體密度 (正電荷)  
            Nd = 1e18  # cm^-3
            n_donors = Nd * 1e6  # 轉換為 m^-3
            
            # 電荷中性條件
            charge_density = PC.E * (n_donors - n_electrons)  # C/m^3
            
            return charge_density
    
    charge_calculator = MockChargeDensityCalculator()
    
    print(f"\n" + "="*60)
    print(f"測試1: 基本Laplace求解 (作為基準)")
    print(f"="*60)
    
    potential_laplace, iterations_l, error_l = solver.solve_laplace(
        V_tip, V_sample, max_iterations=200, tolerance=1e-4
    )
    
    pot0_laplace = solver._calculate_pot0_fortran_style(potential_laplace)
    print(f"Laplace 結果:")
    print(f"  迭代次數: {iterations_l}")
    print(f"  Pot0 = {pot0_laplace:.6f} V")
    print(f"  與Fortran (-0.08V) 差異: {abs(pot0_laplace - (-0.08)):.6f} V")
    
    print(f"\n" + "="*60)
    print(f"測試2: 使用VSINT的非線性求解")
    print(f"="*60)
    
    try:
        potential_vsint, iterations_vs, converged_vs = solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=charge_calculator,
            system_fermi_level_E_F_main_eV=system_fermi_level,
            max_iterations=300,  # 短一些的測試
            tolerance_Volts=1e-4
        )
        
        # 計算VSINT方法的Pot0
        vsint_array = solver._initialize_vsint_array()
        vsint_array = solver._update_vsint_with_surface_charge(
            vsint_array, potential_vsint, charge_calculator,
            system_fermi_level, V_tip)
        
        pot0_vsint = solver._calculate_pot0_fortran_style(
            potential_vsint, use_vsint=True, vsint_array=vsint_array)
        pot0_regular = solver._calculate_pot0_fortran_style(
            potential_vsint, use_vsint=False)
        
        print(f"VSINT 非線性求解結果:")
        print(f"  迭代次數: {iterations_vs}")
        print(f"  收斂狀態: {'是' if converged_vs else '否'}")
        print(f"  VSINT Pot0 = {pot0_vsint:.6f} V")
        print(f"  Regular Pot0 = {pot0_regular:.6f} V")
        print(f"  VSINT與Fortran差異: {abs(pot0_vsint - (-0.08)):.6f} V")
        print(f"  Regular與Fortran差異: {abs(pot0_regular - (-0.08)):.6f} V")
        
        # 分析VSINT陣列
        print(f"\nVSINT 陣列分析:")
        print(f"  VSINT[0,0] = {vsint_array[0,0]:.6f} V")
        print(f"  VSINT[1,0] = {vsint_array[1,0]:.6f} V")
        print(f"  對應的Regular potential[0,nu_max] = {potential_vsint[0, grid.N_nu-1]:.6f} V")
        print(f"  對應的Regular potential[1,nu_max] = {potential_vsint[1, grid.N_nu-1]:.6f} V")
        
        # 手動驗證PCENT公式
        manual_vsint = (9.0 * vsint_array[0,0] - vsint_array[1,0]) / 8.0
        manual_regular = (9.0 * potential_vsint[0, grid.N_nu-1] - potential_vsint[1, grid.N_nu-1]) / 8.0
        
        print(f"\n手動PCENT計算驗證:")
        print(f"  VSINT手動: (9*{vsint_array[0,0]:.4f} - {vsint_array[1,0]:.4f})/8 = {manual_vsint:.6f} V")
        print(f"  Regular手動: (9*{potential_vsint[0, grid.N_nu-1]:.4f} - {potential_vsint[1, grid.N_nu-1]:.4f})/8 = {manual_regular:.6f} V")
        
        print(f"\n" + "="*60)
        print(f"比較分析")
        print(f"="*60)
        
        improvement_vsint = abs(pot0_laplace - (-0.08)) - abs(pot0_vsint - (-0.08))
        improvement_regular = abs(pot0_laplace - (-0.08)) - abs(pot0_regular - (-0.08))
        
        print(f"結果比較:")
        print(f"  Laplace Pot0:     {pot0_laplace:.6f} V (差異: {abs(pot0_laplace - (-0.08)):.6f} V)")
        print(f"  VSINT Pot0:       {pot0_vsint:.6f} V (差異: {abs(pot0_vsint - (-0.08)):.6f} V)")
        print(f"  Regular Pot0:     {pot0_regular:.6f} V (差異: {abs(pot0_regular - (-0.08)):.6f} V)")
        print(f"  Fortran 目標:     -0.08 V")
        print(f"  VSINT改善:        {improvement_vsint:.6f} V ({'改善' if improvement_vsint > 0 else '惡化'})")
        print(f"  Regular改善:      {improvement_regular:.6f} V ({'改善' if improvement_regular > 0 else '惡化'})")
        
        if abs(pot0_vsint - (-0.08)) < 0.1:
            print(f"🎉 VSINT Pot0 已接近Fortran結果 (差異 < 0.1V)")
        elif abs(pot0_vsint - (-0.08)) < abs(pot0_regular - (-0.08)):
            print(f"✅ VSINT方法比Regular方法更好")
        else:
            print(f"⚠️  VSINT方法沒有明顯改善")
            
        return {
            'pot0_laplace': pot0_laplace,
            'pot0_vsint': pot0_vsint,
            'pot0_regular': pot0_regular,
            'vsint_array': vsint_array,
            'potential': potential_vsint
        }
        
    except Exception as e:
        print(f"❌ VSINT非線性求解失敗: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    results = test_vsint_implementation()
    
    if results:
        print(f"\n" + "="*80)
        print(f"總結:")
        print(f"="*80)
        
        fortran_target = -0.08
        print(f"Pot0 結果與Fortran目標({fortran_target}V)比較:")
        print(f"  Laplace:     {results['pot0_laplace']:.6f} V (差異: {abs(results['pot0_laplace'] - fortran_target):.6f} V)")
        print(f"  VSINT:       {results['pot0_vsint']:.6f} V (差異: {abs(results['pot0_vsint'] - fortran_target):.6f} V)")
        print(f"  Regular:     {results['pot0_regular']:.6f} V (差異: {abs(results['pot0_regular'] - fortran_target):.6f} V)")
        
        best_method = min([
            (abs(results['pot0_laplace'] - fortran_target), 'Laplace'),
            (abs(results['pot0_vsint'] - fortran_target), 'VSINT'),
            (abs(results['pot0_regular'] - fortran_target), 'Regular')
        ])
        
        print(f"\n最佳方法: {best_method[1]} (差異: {best_method[0]:.6f} V)")
        
        if best_method[1] == 'VSINT':
            print(f"🎉 VSINT實現成功改善了Pot0計算！")
        else:
            print(f"⚠️  需要進一步優化VSINT實現")
    else:
        print(f"❌ 測試失敗，需要修復VSINT實現")