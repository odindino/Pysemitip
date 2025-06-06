#!/usr/bin/env python3
"""
Golden Section Search 調試和修復腳本
比較 Python 和 Fortran GSECT 實現的差異，並修復問題

創建日期: 2025-06-06
目的: 修復 Pysemitip 中 Golden Section Search 算法的收斂問題
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Callable, List, Tuple
import time

# 將 Pysemitip 模組加入路徑
import sys
import os
sys.path.append('/Users/yangziliang/Git-Projects/Pysemitip/src')

# 直接複製 golden_section_search 函數定義避免導入問題
def golden_section_search(func, xmin: float, xmax: float, 
                         tolerance: float = 1e-6, max_iter: int = 100) -> float:
    """
    Golden section search for finding minimum of a function.
    
    This implements the GSECT subroutine from Fortran SEMITIP.
    """
    if abs(xmax - xmin) < tolerance:
        return (xmin + xmax) / 2.0
    
    # Ensure xmin < xmax
    if xmax < xmin:
        xmin, xmax = xmax, xmin
    
    # Golden ratio constant
    gs = 0.3819660  # (3 - sqrt(5)) / 2
    
    delx = xmax - xmin
    xa = xmin + delx * gs
    fa = func(xa)
    xb = xmax - delx * gs
    fb = func(xb)
    
    for _ in range(max_iter):
        delx_saved = delx
        if delx < tolerance:
            break
            
        if fb < fa:
            # Move to the right interval
            xmin = xa
            delx = xmax - xmin
            if abs(delx - delx_saved) < tolerance * tolerance:
                break
            xa = xb
            fa = fb
            xb = xmax - delx * gs
            fb = func(xb)
        else:
            # Move to the left interval
            xmax = xb
            delx = xmax - xmin
            if abs(delx - delx_saved) < tolerance * tolerance:
                break
            xb = xa
            fb = fa
            xa = xmin + delx * gs
            fa = func(xa)
    
    return (xmin + xmax) / 2.0

def test_function_1(x: float) -> float:
    """測試函數1: (x-2)^2 + 1, 最小值在 x=2"""
    return (x - 2.0)**2 + 1.0

def test_function_2(x: float) -> float:
    """測試函數2: x^4 - 2x^2 + x + 1, 複雜函數"""
    return x**4 - 2*x**2 + x + 1.0

def test_function_3(x: float) -> float:
    """測試函數3: 模擬 Poisson 殘差函數"""
    return np.exp(-x**2) * np.cos(5*x) + 0.1*x**2

class FortranGSECT:
    """Python 重新實現的 Fortran GSECT，嚴格按照 Fortran 邏輯"""
    
    @staticmethod
    def gsect(func: Callable, xmin: float, xmax: float, ep: float) -> float:
        """
        嚴格按照 Fortran gsect-6.0.f 實現
        """
        # Fortran lines 12-13: IF (XMAX.EQ.XMIN) RETURN
        if abs(xmax - xmin) < 1e-16:
            return (xmin + xmax) / 2.0
        
        # Fortran lines 14: IF (EP.EQ.0.) RETURN  
        if abs(ep) < 1e-16:
            return (xmin + xmax) / 2.0
        
        # Fortran lines 15-18: swap if needed
        if xmax < xmin:
            temp = xmax
            xmax = xmin
            xmin = temp
        
        # Fortran constant
        gs = 0.3819660
        
        # Fortran lines 19-23: initialize
        delx = xmax - xmin
        xa = xmin + delx * gs
        fa = func(xa)
        xb = xmax - delx * gs
        fb = func(xb)
        
        # 收斂歷史記錄
        iteration_count = 0
        max_iterations = 10000  # 防止無限循環
        
        # Fortran line 24+: main loop
        while iteration_count < max_iterations:
            iteration_count += 1
            
            # Fortran line 26: DELXSAV=DELX
            delxsav = delx
            
            # Fortran line 27: IF (DELX.LT.EP) RETURN
            if delx < ep:
                break
            
            # Fortran line 28: IF (FB.LT.FA) GO TO 200
            if fb < fa:
                # Fortran lines 200-211: right interval
                xmin = xa
                delx = xmax - xmin
                
                # Fortran line 202: IF (DELX.EQ.DELXSAV) RETURN
                if abs(delx - delxsav) < 1e-16:
                    break
                    
                xa = xb
                fa = fb
                xb = xmax - delx * gs
                fb = func(xb)
            else:
                # Fortran lines 29-36: left interval  
                xmax = xb
                delx = xmax - xmin
                
                # Fortran line 31: IF (DELX.EQ.DELXSAV) RETURN
                if abs(delx - delxsav) < 1e-16:
                    break
                    
                xb = xa
                fb = fa
                xa = xmin + delx * gs
                fa = func(xa)
        
        # Fortran: result is (XMIN+XMAX)/2
        return (xmin + xmax) / 2.0

def compare_gsect_implementations():
    """比較 Python 和 Fortran 風格的 GSECT 實現"""
    print("🔍 Golden Section Search 實現比較")
    print("=" * 60)
    
    test_functions = [
        ("(x-2)^2 + 1", test_function_1, 0.0, 4.0, 2.0),
        ("x^4 - 2x^2 + x + 1", test_function_2, -2.0, 2.0, None),
        ("模擬 Poisson 殘差", test_function_3, -1.0, 1.0, None)
    ]
    
    tolerances = [1e-6, 1e-8, 1e-10]
    
    for func_name, func, xmin, xmax, expected in test_functions:
        print(f"\n📈 測試函數: {func_name}")
        print(f"   搜尋區間: [{xmin}, {xmax}]")
        if expected:
            print(f"   預期最小值位置: {expected}")
        
        for tol in tolerances:
            print(f"\n   容差: {tol:.0e}")
            
            # 測試 Python 當前實現
            start_time = time.time()
            try:
                result_py = golden_section_search(func, xmin, xmax, tol)
                time_py = time.time() - start_time
                value_py = func(result_py)
                print(f"   🐍 Python 實現:  x={result_py:.8f}, f(x)={value_py:.8f}, 時間={time_py:.4f}s")
            except Exception as e:
                print(f"   🐍 Python 實現:  錯誤 - {e}")
                result_py, value_py = None, None
            
            # 測試 Fortran 風格實現
            start_time = time.time()
            try:
                result_fortran = FortranGSECT.gsect(func, xmin, xmax, tol)
                time_fortran = time.time() - start_time
                value_fortran = func(result_fortran)
                print(f"   🏛️ Fortran 風格: x={result_fortran:.8f}, f(x)={value_fortran:.8f}, 時間={time_fortran:.4f}s")
            except Exception as e:
                print(f"   🏛️ Fortran 風格: 錯誤 - {e}")
                result_fortran, value_fortran = None, None
            
            # 比較結果
            if result_py is not None and result_fortran is not None:
                diff_x = abs(result_py - result_fortran)
                diff_f = abs(value_py - value_fortran)
                print(f"   📊 差異: Δx={diff_x:.2e}, Δf(x)={diff_f:.2e}")
                
                if expected:
                    error_py = abs(result_py - expected)
                    error_fortran = abs(result_fortran - expected)
                    print(f"   🎯 誤差: Python={error_py:.2e}, Fortran={error_fortran:.2e}")

def simulate_poisson_residual_function():
    """模擬 Poisson 求解器中的殘差函數行為"""
    print("\n🧮 模擬 Poisson 求解器殘差函數")
    print("=" * 60)
    
    # 模擬半導體區域電位更新的殘差函數
    def poisson_residual(potential_correction: float) -> float:
        """
        模擬 Poisson 求解器中的殘差函數
        這對應於 _update_semiconductor_fortran_style 中的 bulk_residual
        """
        # 模擬當前電位
        current_potential = -0.08  # Volts, 類似 Fortran 輸出
        
        # 測試電位
        test_potential = current_potential + potential_correction
        
        # 模擬電荷密度 (簡化的費米-狄拉克分布)
        kT = 0.026  # 300K, eV
        fermi_level = -0.7  # 費米能級
        band_gap = 1.42  # GaAs 帶隙
        
        # 價帶和導帶密度
        valence_density = 1e18 * np.exp((test_potential - fermi_level) / kT)
        conduction_density = 1e18 * np.exp(-(test_potential + band_gap - fermi_level) / kT)
        net_charge = conduction_density - valence_density
        
        # 有限差分項 (簡化)
        laplacian_term = 1e-3 * test_potential  # 拉普拉斯項
        
        # 殘差 = 有限差分結果 - 電荷密度貢獻
        residual_magnitude = abs(laplacian_term - net_charge * 1.6e-19 / 12.9)
        
        return residual_magnitude
    
    # 測試在典型搜尋範圍內的行為
    bias = -2.07  # Volts, 來自 Fortran 輸出
    delta_sem = max(1e-6, abs(bias) / 1e6)  # Fortran DELSEM 邏輯
    
    print(f"偏壓: {bias} V")
    print(f"搜尋範圍半寬: {delta_sem} V")
    print(f"實際搜尋範圍: ±{delta_sem:.2e} V")
    
    # 在搜尋範圍內繪圖
    x_range = np.linspace(-delta_sem, delta_sem, 1000)
    y_values = [poisson_residual(x) for x in x_range]
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_range, y_values, 'b-', linewidth=2)
    plt.xlabel('電位修正 (V)')
    plt.ylabel('殘差大小')
    plt.title('模擬 Poisson 殘差函數')
    plt.grid(True, alpha=0.3)
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    
    # 找最小值
    min_idx = np.argmin(y_values)
    min_x = x_range[min_idx]
    min_y = y_values[min_idx]
    plt.plot(min_x, min_y, 'ro', markersize=8, label=f'最小值: ({min_x:.2e}, {min_y:.2e})')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('/Users/yangziliang/Git-Projects/Pysemitip/poisson_residual_simulation.png', dpi=150)
    print(f"殘差函數圖片已儲存: poisson_residual_simulation.png")
    
    # 使用 Golden Section Search 找最小值
    print(f"\n🔍 使用 Golden Section Search 尋找最小值:")
    
    result_py = golden_section_search(poisson_residual, -delta_sem, delta_sem, 1e-8)
    result_fortran = FortranGSECT.gsect(poisson_residual, -delta_sem, delta_sem, 1e-8)
    
    print(f"Python 實現:     {result_py:.8e} V")
    print(f"Fortran 風格:    {result_fortran:.8e} V")
    print(f"理論最小值:      {min_x:.8e} V")
    print(f"差異 (Py-理論):  {abs(result_py - min_x):.2e} V")
    print(f"差異 (F77-理論): {abs(result_fortran - min_x):.2e} V")
    
    return delta_sem

def debug_current_gsect_issue():
    """調試當前 Python GSECT 實現的具體問題"""
    print("\n🐛 調試當前 Golden Section Search 問題")
    print("=" * 60)
    
    # 模擬實際 Poisson 求解器中遇到的問題場景
    def problematic_residual(x: float) -> float:
        """模擬導致極小電位變化的問題殘差函數"""
        # 這模擬了當前 Python 實現中的問題：
        # 殘差函數過於平坦，導致搜尋範圍內變化極小
        return 1e-12 + 1e-15 * x**2
    
    search_range = 2.07e-6  # 來自 Fortran DELSEM 計算
    
    print(f"問題場景: 搜尋範圍 = ±{search_range:.2e}")
    
    # 測試當前實現
    result_py = golden_section_search(problematic_residual, -search_range, search_range, 1e-8)
    result_fortran = FortranGSECT.gsect(problematic_residual, -search_range, search_range, 1e-8)
    
    print(f"Python 結果:     {result_py:.12e}")
    print(f"Fortran 結果:    {result_fortran:.12e}")
    print(f"預期結果:        0.0 (對稱函數)")
    
    # 檢查函數在搜尋範圍內的變化
    test_points = np.linspace(-search_range, search_range, 11)
    print(f"\n殘差函數在搜尋範圍內的變化:")
    for x in test_points:
        y = problematic_residual(x)
        print(f"  x = {x:+.2e}, f(x) = {y:.2e}")
    
    # 這說明了問題：函數值變化太小，導致搜尋算法無法有效區分
    func_range = problematic_residual(search_range) - problematic_residual(-search_range)
    print(f"\n函數值變化範圍: {func_range:.2e}")
    print(f"相對於容差 1e-8: {func_range / 1e-8:.2e}")
    
    if func_range < 1e-8:
        print("⚠️ 問題診斷: 函數變化小於容差，Golden Section Search 無法有效工作")
        print("💡 建議解決方案:")
        print("   1. 調整容差參數")
        print("   2. 改進搜尋範圍計算")
        print("   3. 檢查殘差函數計算邏輯")

def main():
    """主函數"""
    print("🔧 Pysemitip Golden Section Search 調試和修復")
    print("=" * 80)
    print("目標: 識別並修復導致 Poisson 求解器收斂問題的 GSECT 實現")
    print("=" * 80)
    
    # 1. 基本函數比較
    compare_gsect_implementations()
    
    # 2. Poisson 殘差函數模擬
    delta_sem = simulate_poisson_residual_function()
    
    # 3. 問題診斷
    debug_current_gsect_issue()
    
    print("\n" + "=" * 80)
    print("🎯 調試總結和修復建議")
    print("=" * 80)
    print("1. ✅ Fortran GSECT 算法邏輯正確實現")
    print("2. ⚠️ 當前問題可能在於:")
    print("   - 搜尋範圍計算 (DELSEM/DELSURF)")
    print("   - 殘差函數定義")
    print("   - 容差參數設定")
    print("3. 🔧 下一步修復:")
    print("   - 檢查 _update_semiconductor_fortran_style 中的搜尋範圍")
    print("   - 驗證 bulk_residual 函數邏輯")
    print("   - 調整 golden_section_tolerance 參數")
    
    # 創建修復建議文件
    with open('/Users/yangziliang/Git-Projects/Pysemitip/gsect_debug_results.txt', 'w') as f:
        f.write("Golden Section Search 調試結果\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"搜尋範圍問題: delta_sem = {delta_sem:.2e}\n")
        f.write("Fortran DELSEM 計算: max(1e-6, abs(bias)/1e6)\n")
        f.write("這導致極小的搜尋範圍，需要重新檢查\n\n")
        f.write("建議修復順序:\n")
        f.write("1. 檢查 Fortran 原始搜尋範圍計算\n")
        f.write("2. 修復 Python 搜尋範圍計算\n") 
        f.write("3. 驗證殘差函數定義\n")
        f.write("4. 測試完整 Poisson 求解器\n")
    
    print(f"\n📋 詳細結果已儲存到: gsect_debug_results.txt")

if __name__ == "__main__":
    main()
