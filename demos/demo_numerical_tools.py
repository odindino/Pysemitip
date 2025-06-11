"""
Pysemitip數值工具演示

此腳本展示第一階段實現的基礎數值工具的主要功能和應用場景。
"""

import numpy as np
import sys
import os

# 添加src路徑
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from utils.numerical import (
    golden_section_search,
    fermi_dirac_integral,
    numerical_derivative,
    adaptive_quadrature
)

from utils.interpolation import (
    linear_interpolation,
    cubic_spline_interpolation,
    bilinear_interpolation
)

# 嘗試導入matplotlib，如果沒有就跳過繪圖
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def demo_golden_section():
    """演示黃金分割優化算法"""
    print("1. 黃金分割優化算法演示")
    print("-" * 40)
    
    # 模擬SEMITIP中的電子親和力優化問題
    def electron_affinity_objective(chi):
        """
        模擬SEMMIN中的電子親和力優化目標函數
        chi: 電子親和力 (eV)
        """
        # 簡化的目標函數：最小化實驗與理論tunneling current的差異
        optimal_chi = 4.2  # 假設的最優值
        return (chi - optimal_chi)**2 + 0.1 * np.sin(5 * chi)
    
    # 搜尋最優電子親和力
    chi_opt, f_opt, iters = golden_section_search(
        electron_affinity_objective, 3.0, 5.0, precision=1e-6
    )
    
    print(f"  電子親和力優化結果:")
    print(f"  最優值: χ = {chi_opt:.6f} eV")
    print(f"  目標函數值: {f_opt:.8f}")
    print(f"  迭代次數: {iters}")
    
    # 可視化優化過程
    chi_range = np.linspace(3.0, 5.0, 1000)
    f_values = [electron_affinity_objective(chi) for chi in chi_range]
    
    if HAS_MATPLOTLIB:
        try:
            plt.figure(figsize=(10, 6))
            plt.plot(chi_range, f_values, 'b-', linewidth=2, label='Objective Function')
            plt.axvline(chi_opt, color='r', linestyle='--', linewidth=2, 
                       label=f'Optimal Solution: χ = {chi_opt:.6f} eV')
            plt.xlabel('Electron Affinity χ (eV)')
            plt.ylabel('Objective Function Value')
            plt.title('SEMITIP Electron Affinity Optimization Demo')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig('demos/gsect_demo.png', dpi=150, bbox_inches='tight')
            print(f"  Optimization plot saved as demos/gsect_demo.png")
        except Exception as e:
            print(f"  (圖表生成失敗: {e})")
    else:
        print(f"  (無法生成圖表，請安裝matplotlib)")


def demo_fermi_dirac():
    """演示費米-狄拉克積分"""
    print("\n2. 費米-狄拉克積分演示")
    print("-" * 40)
    
    # 計算不同溫度和費米能級下的載流子濃度
    eta_values = np.linspace(-5, 5, 11)
    
    print(f"  半導體載流子濃度計算 (使用F_1/2積分):")
    print(f"  約化化學勢η   F_1/2(η)      載流子濃度比例")
    print(f"  " + "-" * 50)
    
    for eta in eta_values:
        f_val = fermi_dirac_integral(0.5, eta)
        # 載流子濃度與F_1/2成正比
        carrier_ratio = f_val / fermi_dirac_integral(0.5, 0)  # 相對於η=0的比例
        
        print(f"  {eta:8.1f}      {f_val:8.6f}     {carrier_ratio:8.4f}")
    
    # 溫度依賴性演示
    print(f"\n  溫度對費米-狄拉克分佈的影響:")
    temperatures = [77, 300, 500]  # K
    k_B = 8.617e-5  # eV/K
    E_F = 0.5  # 費米能級相對於導帶底 (eV)
    
    for T in temperatures:
        eta = E_F / (k_B * T)
        f_val = fermi_dirac_integral(0.5, eta)
        print(f"  T = {T:3d} K: η = {eta:6.2f}, F_1/2(η) = {f_val:8.6f}")


def demo_interpolation():
    """演示插值算法"""
    print("\n3. 插值算法演示")
    print("-" * 40)
    
    # 模擬SEMITIP中的電位分佈插值
    print(f"  電位場空間插值演示:")
    
    # 一維線性插值：沿z軸的電位分佈
    z_data = np.array([0, 5, 10, 15, 20])  # nm
    V_data = np.array([0, -0.5, -1.2, -1.8, -2.0])  # V
    
    z_interp = np.linspace(0, 20, 41)
    V_linear = linear_interpolation(z_data, V_data, z_interp)
    V_spline = cubic_spline_interpolation(z_data, V_data, z_interp)
    
    print(f"  z軸位置(nm)  已知電位(V)  線性插值(V)  三次樣條(V)")
    print(f"  " + "-" * 55)
    
    for i in range(0, len(z_interp), 8):
        z = z_interp[i]
        v_lin = V_linear[i]
        v_spl = V_spline[i]
        
        # 檢查是否是已知數據點
        if z in z_data:
            idx = np.where(z_data == z)[0][0]
            v_known = V_data[idx]
            print(f"  {z:8.1f}      {v_known:8.3f}     {v_lin:8.3f}     {v_spl:8.3f}")
        else:
            print(f"  {z:8.1f}         --       {v_lin:8.3f}     {v_spl:8.3f}")
    
    # 二維雙線性插值：表面電位分佈
    print(f"\n  表面電位分佈2D插值:")
    x_coords = np.array([0, 10, 20])  # nm
    y_coords = np.array([0, 10])      # nm
    V_surface = np.array([
        [0.0, -0.5, -1.0],    # y=0
        [-0.2, -0.7, -1.2]    # y=10
    ])
    
    # 插值到新的點
    test_points = [(5, 5), (15, 3), (8, 7)]
    print(f"  位置(nm)      插值電位(V)")
    print(f"  " + "-" * 25)
    
    for x, y in test_points:
        V_interp = bilinear_interpolation(x_coords, y_coords, V_surface, x, y)
        print(f"  ({x:2d}, {y:2d})        {V_interp:8.3f}")


def demo_numerical_calculus():
    """演示數值微積分"""
    print("\n4. 數值微積分演示")
    print("-" * 40)
    
    # 電場計算：E = -∇V
    def potential_function(z):
        """模擬一維電位函數 V(z)"""
        return -2.0 * np.exp(-z/10) + 0.1 * z
    
    positions = [0, 5, 10, 15, 20]
    print(f"  電場計算 E_z = -dV/dz:")
    print(f"  位置z(nm)   電位V(V)    電場E_z(V/nm)")
    print(f"  " + "-" * 40)
    
    for z in positions:
        V = potential_function(z)
        E_z = -numerical_derivative(potential_function, z, h=1e-3)
        print(f"  {z:6.1f}     {V:8.3f}     {E_z:8.4f}")
    
    # 數值積分：計算總電荷
    print(f"\n  電荷密度積分計算:")
    
    def charge_density(x):
        """模擬電荷密度分佈 ρ(x)"""
        return 1.6e-19 * np.exp(-(x-10)**2/25)  # C/nm³
    
    # 計算總電荷
    total_charge = adaptive_quadrature(charge_density, 0, 20, tolerance=1e-12)
    print(f"  積分區間: [0, 20] nm")
    print(f"  總電荷: Q = {total_charge:.6e} C")
    print(f"  (每nm³的電荷密度在x=10nm處最大)")


def main():
    """主演示函數"""
    print("Pysemitip數值工具功能演示")
    print("=" * 60)
    print("展示第一階段基礎工具層的核心功能")
    print("適用於掃描穿隧顯微鏡(STM)的3D有限差分求解器")
    print("=" * 60)
    
    demo_golden_section()
    demo_fermi_dirac()
    demo_interpolation()
    demo_numerical_calculus()
    
    print("\n" + "=" * 60)
    print("🎯 演示總結:")
    print("✅ 黃金分割優化 - 已實現並驗證，可用於SEMMIN/SURFMIN")
    print("✅ 費米-狄拉克積分 - 已實現載流子濃度計算")
    print("✅ 插值算法 - 支援1D/2D電位場插值")
    print("✅ 數值微積分 - 支援電場計算和電荷積分")
    print("=" * 60)
    print("\n🚀 第一階段基礎工具層實現完成！")
    print("📝 已建立:")
    print("   - utils/numerical.py (核心數值計算)")
    print("   - utils/interpolation.py (插值工具)")
    print("   - tests/test_numerical_tools.py (驗證測試)")
    print("   - tests/test_gsect_fortran_compatibility.py (Fortran兼容性)")
    print("\n🎯 準備進入第二階段: 物理模型層實現")


if __name__ == "__main__":
    main()
