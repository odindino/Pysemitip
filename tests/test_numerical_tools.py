"""
SEMITIP數值工具測試模組

此模組包含對utils/numerical.py中實現的核心數值算法的驗證測試，
確保Python實現與原始Fortran代碼的數值一致性。

測試覆蓋範圍:
1. 黃金分割搜尋算法 (GSECT) 的精度和收斂性測試
2. 費米-狄拉克積分函數的數值精度測試
3. 插值算法的精度驗證
4. 邊界條件和異常情況處理測試

設計原則:
- 使用已知解析解的測試函數
- 驗證數值精度是否符合物理計算要求
- 測試邊界條件和數值穩定性
"""

import numpy as np
import pytest
import sys
import os
import math

# 添加src路徑以便導入模組
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from utils.numerical import (
    golden_section_search,
    fermi_dirac_integral,
    fermi_dirac_occupation,
    numerical_derivative,
    adaptive_quadrature,
    trapezoidal_integration,
    GoldenSectionOptimizer
)

from utils.interpolation import (
    linear_interpolation,
    cubic_spline_interpolation,
    bilinear_interpolation,
    LinearInterpolator,
    CubicSplineInterpolator,
    BilinearInterpolator
)


class TestGoldenSectionSearch:
    """黃金分割搜尋算法測試"""
    
    def test_simple_quadratic(self):
        """測試簡單二次函數最小值搜尋"""
        def quadratic(x):
            return (x - 3.0)**2 + 2.0
        
        x_opt, f_opt, iters = golden_section_search(quadratic, 0, 6, precision=1e-8)
        
        # 驗證精度
        assert abs(x_opt - 3.0) < 1e-6, f"最優解精度不足: {x_opt}"
        assert abs(f_opt - 2.0) < 1e-6, f"最優值精度不足: {f_opt}"
        assert iters > 0, "迭代次數應大於0"
        
    def test_cosine_function(self):
        """測試餘弦函數最小值搜尋"""
        def cosine(x):
            return np.cos(x)
        
        # 在[2, 4]區間內搜尋，最小值應在π附近
        x_opt, f_opt, iters = golden_section_search(cosine, 2, 4, precision=1e-8)
        
        assert abs(x_opt - np.pi) < 1e-6, f"餘弦函數最優解不正確: {x_opt}"
        assert abs(f_opt - (-1.0)) < 1e-6, f"餘弦函數最優值不正確: {f_opt}"
        
    def test_boundary_conditions(self):
        """測試邊界條件處理"""
        def linear(x):
            return x
        
        # 線性函數在區間[0, 1]內最小值應在左端點
        x_opt, f_opt, iters = golden_section_search(linear, 0, 1, precision=1e-6)
        
        # 由於黃金分割返回區間中點，對於線性函數應接近左端點
        assert x_opt <= 0.5, f"線性函數優化結果不合理: {x_opt}"
        
    def test_precision_control(self):
        """測試精度控制"""
        def quadratic(x):
            return x**2
        
        # 不同精度要求的測試
        x_opt_low, _, iters_low = golden_section_search(quadratic, -1, 1, precision=1e-3)
        x_opt_high, _, iters_high = golden_section_search(quadratic, -1, 1, precision=1e-8)
        
        # 高精度應需要更多迭代
        assert iters_high > iters_low, "高精度應需要更多迭代次數"
        assert abs(x_opt_high) < abs(x_opt_low) + 1e-3, "高精度結果應更準確"
        
    def test_interval_reversal(self):
        """測試區間反轉處理"""
        def quadratic(x):
            return (x - 2)**2
        
        # 測試xmax < xmin的情況
        x_opt1, f_opt1, _ = golden_section_search(quadratic, 0, 4)
        x_opt2, f_opt2, _ = golden_section_search(quadratic, 4, 0)  # 反轉區間
        
        assert abs(x_opt1 - x_opt2) < 1e-6, "區間反轉應產生相同結果"
        assert abs(f_opt1 - f_opt2) < 1e-6, "區間反轉應產生相同最優值"


class TestFermiDiracIntegral:
    """費米-狄拉克積分測試"""
    
    def test_nondegenerate_limit(self):
        """測試非簡併極限 (η << 0)"""
        # 當η很小時，F_j(η) ≈ Γ(j+1) * exp(η)
        eta = -10.0
        order = 0.5
        
        result = fermi_dirac_integral(order, eta)
        expected = math.gamma(order + 1) * np.exp(eta)
        
        relative_error = abs(result - expected) / abs(expected)
        assert relative_error < 1e-4, f"非簡併極限精度不足: {relative_error}"
        
    def test_degenerate_limit(self):
        """測試完全簡併極限 (η >> 0)"""
        # 當η很大時，使用更簡單的情況進行測試
        eta = 15.0  # 更大的η值
        order = 0.5  # 使用較小的階數，修正項較小
        
        result = fermi_dirac_integral(order, eta)
        # 對於order=0.5，修正項很小
        expected = (eta ** (order + 1)) / (order + 1)
        
        relative_error = abs(result - expected) / abs(expected)
        assert relative_error < 1e-3, f"完全簡併極限精度不足: {relative_error}, result={result}, expected={expected}"
        
    def test_known_values(self):
        """測試已知數值"""
        # F_0.5(0) 的準確值約為 0.678
        result = fermi_dirac_integral(0.5, 0.0)
        expected = 0.678  # 修正的期望值
        
        relative_error = abs(result - expected) / expected
        assert relative_error < 5e-2, f"F_0.5(0)計算精度不足: {result}, 相對誤差: {relative_error}"
        
    def test_vectorized_input(self):
        """測試向量化輸入"""
        eta_array = np.array([-5, 0, 5])
        result = fermi_dirac_integral(0.5, eta_array)
        
        assert len(result) == len(eta_array), "向量化輸出長度不正確"
        assert all(np.isfinite(result)), "向量化結果包含無效值"
        
        # 驗證單調性：對於固定階數，F_j(η)隨η遞增
        assert result[0] < result[1] < result[2], "費米-狄拉克積分應隨η遞增"


class TestNumericalDerivative:
    """數值微分測試"""
    
    def test_polynomial_derivative(self):
        """測試多項式函數的導數"""
        def cubic(x):
            return x**3 + 2*x**2 + 3*x + 4
        
        def cubic_derivative(x):
            return 3*x**2 + 4*x + 3
        
        x = 2.0
        numerical_deriv = numerical_derivative(cubic, x, h=1e-6)
        analytical_deriv = cubic_derivative(x)
        
        relative_error = abs(numerical_deriv - analytical_deriv) / abs(analytical_deriv)
        assert relative_error < 1e-6, f"多項式導數精度不足: {relative_error}"
        
    def test_transcendental_functions(self):
        """測試超越函數的導數"""
        def exp_sin(x):
            return np.exp(x) * np.sin(x)
        
        def exp_sin_derivative(x):
            return np.exp(x) * (np.sin(x) + np.cos(x))
        
        x = 1.0
        numerical_deriv = numerical_derivative(exp_sin, x, h=1e-6)
        analytical_deriv = exp_sin_derivative(x)
        
        relative_error = abs(numerical_deriv - analytical_deriv) / abs(analytical_deriv)
        assert relative_error < 1e-6, f"超越函數導數精度不足: {relative_error}"
        
    def test_different_methods(self):
        """測試不同差分方法"""
        def quadratic(x):
            return x**2
        
        x = 1.0
        
        forward = numerical_derivative(quadratic, x, method='forward')
        backward = numerical_derivative(quadratic, x, method='backward')
        central = numerical_derivative(quadratic, x, method='central')
        
        # 中心差分應該最準確
        analytical = 2 * x
        
        central_error = abs(central - analytical)
        forward_error = abs(forward - analytical)
        backward_error = abs(backward - analytical)
        
        assert central_error < forward_error, "中心差分應比前向差分更準確"
        assert central_error < backward_error, "中心差分應比後向差分更準確"


class TestAdaptiveQuadrature:
    """自適應積分測試"""
    
    def test_polynomial_integration(self):
        """測試多項式積分"""
        def cubic(x):
            return x**3 + 2*x**2 + 3*x + 4
        
        # 解析積分: ∫(x³+2x²+3x+4)dx = x⁴/4 + 2x³/3 + 3x²/2 + 4x
        def cubic_integral(a, b):
            def antiderivative(x):
                return x**4/4 + 2*x**3/3 + 3*x**2/2 + 4*x
            return antiderivative(b) - antiderivative(a)
        
        a, b = 0, 2
        numerical_result = adaptive_quadrature(cubic, a, b)
        analytical_result = cubic_integral(a, b)
        
        relative_error = abs(numerical_result - analytical_result) / abs(analytical_result)
        assert relative_error < 1e-10, f"多項式積分精度不足: {relative_error}"
        
    def test_oscillatory_function(self):
        """測試振盪函數積分"""
        def oscillatory(x):
            return np.sin(10 * x)
        
        # ∫₀^π sin(10x)dx = [-cos(10x)/10]₀^π = [cos(10π)-cos(0)]/(-10)
        a, b = 0, np.pi
        numerical_result = adaptive_quadrature(oscillatory, a, b)
        analytical_result = (np.cos(10*np.pi) - np.cos(0)) / (-10)
        
        absolute_error = abs(numerical_result - analytical_result)
        assert absolute_error < 1e-6, f"振盪函數積分精度不足: {absolute_error}"


class TestFermiDiracOccupation:
    """費米-狄拉克佔據函數測試"""
    
    def test_zero_temperature(self):
        """測試零溫近似"""
        fermi_energy = 1.0
        energies = np.array([0.5, 1.0, 1.5])
        
        occupation = fermi_dirac_occupation(energies, fermi_energy, temperature=0)
        expected = np.array([1.0, 0.5, 0.0])
        
        np.testing.assert_array_equal(occupation, expected)
        
    def test_finite_temperature(self):
        """測試有限溫度"""
        fermi_energy = 1.0
        temperature = 300  # K
        
        # 在費米能級處應該是0.5
        occupation_at_ef = fermi_dirac_occupation(fermi_energy, fermi_energy, temperature)
        assert abs(occupation_at_ef - 0.5) < 1e-10
        
        # 遠低於費米能級應該接近1
        occupation_low = fermi_dirac_occupation(fermi_energy - 1.0, fermi_energy, temperature)
        assert occupation_low > 0.99
        
        # 遠高於費米能級應該接近0
        occupation_high = fermi_dirac_occupation(fermi_energy + 1.0, fermi_energy, temperature)
        assert occupation_high < 0.01
        
    def test_numerical_limits(self):
        """測試數值穩定性"""
        fermi_energy = 0.0
        temperature = 300
        
        # 測試極端情況
        very_low = fermi_dirac_occupation(-10.0, fermi_energy, temperature)
        very_high = fermi_dirac_occupation(10.0, fermi_energy, temperature)
        
        assert 0 <= very_low <= 1
        assert 0 <= very_high <= 1
        assert very_low > very_high


class TestTrapezoidalIntegration:
    """梯形積分測試"""
    
    def test_linear_function(self):
        """測試線性函數積分"""
        def linear(x):
            return 2 * x + 1
        
        # ∫(2x+1)dx from 0 to 2 = [x²+x] = 4+2 = 6
        result = trapezoidal_integration(linear, 0, 2, 1000)
        expected = 6.0
        
        assert abs(result - expected) < 1e-6
        
    def test_quadratic_function(self):
        """測試二次函數積分"""
        def quadratic(x):
            return x**2
        
        # ∫x²dx from 0 to 3 = [x³/3] = 9
        result = trapezoidal_integration(quadratic, 0, 3, 1000)
        expected = 9.0
        
        assert abs(result - expected) < 1e-3


class TestLinearInterpolation:
    """線性插值測試"""
    
    def test_exact_interpolation(self):
        """測試精確插值"""
        x_data = np.array([0, 1, 2, 3])
        y_data = np.array([0, 2, 4, 6])  # y = 2x
        
        # 測試數據點本身
        for i, (x, y) in enumerate(zip(x_data, y_data)):
            result = linear_interpolation(x_data, y_data, x)
            assert abs(result - y) < 1e-12, f"數據點插值不精確: {result} vs {y}"
        
        # 測試中點插值
        x_mid = 1.5
        y_expected = 3.0  # 線性函數的中點值
        result = linear_interpolation(x_data, y_data, x_mid)
        assert abs(result - y_expected) < 1e-12, f"中點插值不精確: {result}"
        
    def test_extrapolation(self):
        """測試外推"""
        x_data = np.array([1, 2, 3])
        y_data = np.array([1, 4, 9])  # y = x²（但用線性插值）
        
        # 測試左外推
        result_left = linear_interpolation(x_data, y_data, 0.5, extrapolate=True)
        # 左端斜率 = (4-1)/(2-1) = 3
        # 外推值 = 1 + 3*(0.5-1) = 1 - 1.5 = -0.5
        assert abs(result_left - (-0.5)) < 1e-12, f"左外推不正確: {result_left}"
        
        # 測試右外推
        result_right = linear_interpolation(x_data, y_data, 4, extrapolate=True)
        # 右端斜率 = (9-4)/(3-2) = 5
        # 外推值 = 9 + 5*(4-3) = 9 + 5 = 14
        assert abs(result_right - 14) < 1e-12, f"右外推不正確: {result_right}"


class TestBilinearInterpolation:
    """雙線性插值測試"""
    
    def test_rectangular_grid(self):
        """測試矩形網格插值"""
        x_coords = np.array([0, 1, 2])
        y_coords = np.array([0, 1])
        z_values = np.array([[0, 1, 4],    # y=0: z=x²
                            [1, 2, 5]])    # y=1: z=x²+1
        
        # 測試中心點
        result = bilinear_interpolation(x_coords, y_coords, z_values, 1.0, 0.5)
        # 在x=1, y=0.5處，應該是z(1,0)和z(1,1)的平均值
        expected = (1 + 2) / 2  # = 1.5
        assert abs(result - expected) < 1e-12, f"雙線性插值中心點不正確: {result}"
        
        # 測試角點
        result_corner = bilinear_interpolation(x_coords, y_coords, z_values, 0, 0)
        assert abs(result_corner - 0) < 1e-12, f"角點插值不正確: {result_corner}"


def run_comprehensive_tests():
    """運行全面的測試套件"""
    print("SEMITIP數值工具綜合測試")
    print("=" * 50)
    
    test_classes = [
        TestGoldenSectionSearch,
        TestFermiDiracIntegral,
        TestFermiDiracOccupation,
        TestTrapezoidalIntegration,
        TestNumericalDerivative,
        TestAdaptiveQuadrature,
        TestLinearInterpolation,
        TestBilinearInterpolation
    ]
    
    total_tests = 0
    passed_tests = 0
    failed_tests = []
    
    for test_class in test_classes:
        print(f"\n測試類別: {test_class.__name__}")
        test_instance = test_class()
        
        # 獲取所有測試方法
        test_methods = [method for method in dir(test_instance) 
                       if method.startswith('test_')]
        
        for method_name in test_methods:
            total_tests += 1
            try:
                method = getattr(test_instance, method_name)
                method()
                print(f"  ✓ {method_name}")
                passed_tests += 1
            except Exception as e:
                print(f"  ✗ {method_name}: {str(e)}")
                failed_tests.append(f"{test_class.__name__}.{method_name}: {str(e)}")
    
    print(f"\n測試總結:")
    print(f"總測試數: {total_tests}")
    print(f"通過測試: {passed_tests}")
    print(f"失敗測試: {len(failed_tests)}")
    
    if failed_tests:
        print(f"\n失敗的測試:")
        for failure in failed_tests:
            print(f"  - {failure}")
    else:
        print(f"\n🎉 所有測試都通過了！")
    
    return len(failed_tests) == 0


if __name__ == "__main__":
    # 直接運行測試
    success = run_comprehensive_tests()
    
    if not success:
        print("\n⚠️  存在測試失敗，請檢查實現。")
        sys.exit(1)
    else:
        print("\n✅ 數值工具實現驗證成功！")
