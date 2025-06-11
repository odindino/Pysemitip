"""
核心數值計算模組

此模組包含SEMITIP專案的核心數值計算函數，特別是從Fortran移植的關鍵算法。

主要功能:
1. 黃金分割優化算法 (GSECT) - 一維函數最小值搜尋
2. 費米-狄拉克積分函數 (FJINT) - 統計力學積分
3. 數值微分工具
4. 自適應積分算法

設計原則:
- 保持與原始Fortran代碼的數值一致性
- 提供現代Python介面
- 支援向量化操作
- 包含詳細的物理意義註解
"""

import numpy as np
from typing import Callable, Tuple, Optional, Union
import warnings
import math
from functools import wraps
from scipy import special


def validate_numerical_inputs(func):
    """數值輸入驗證裝飾器"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        # 基本的數值穩定性檢查
        for arg in args:
            if isinstance(arg, (int, float)):
                if not np.isfinite(arg):
                    raise ValueError(f"輸入包含非有限值: {arg}")
        return func(*args, **kwargs)
    return wrapper


class GoldenSectionOptimizer:
    """
    黃金分割優化器 - 從Fortran GSECT-6.0移植
    
    這是一個一維函數最小值搜尋算法，在SEMITIP中用於：
    - 電子親和力優化 (SEMMIN)
    - 表面位勢優化 (SURFMIN)
    
    算法原理:
    黃金分割搜尋使用黃金比例 φ = (√5-1)/2 ≈ 0.618 來逐步縮小搜尋區間，
    確保在每次迭代中保持最佳的區間分割。
    """
    
    def __init__(self):
        # 黃金分割比例 (與Fortran GSECT中的GS=0.3819660對應)
        self.golden_ratio = 0.3819660  # (3-√5)/2, 黃金分割的補角
        
    @validate_numerical_inputs
    def optimize(self, 
                func: Callable[[float], float],
                xmin: float, 
                xmax: float, 
                precision: float = 1e-6,
                max_iterations: int = 1000) -> Tuple[float, float, int]:
        """
        黃金分割搜尋最小值
        
        參數:
            func: 目標函數 f(x)
            xmin: 搜尋區間下界
            xmax: 搜尋區間上界  
            precision: 收斂精度 (對應Fortran中的EP)
            max_iterations: 最大迭代次數
            
        返回:
            tuple: (最優解x, 最優值f(x), 迭代次數)
            
        注意: 此實現直接對應Fortran GSECT-6.0的邏輯
        """
        # 處理邊界情況
        if abs(xmax - xmin) < precision:
            x_opt = (xmin + xmax) / 2
            return x_opt, func(x_opt), 0
            
        if precision <= 0:
            raise ValueError("精度必須大於0")
            
        # 確保xmin < xmax
        if xmax < xmin:
            xmin, xmax = xmax, xmin
            
        # 初始化
        delx = xmax - xmin
        xa = xmin + delx * self.golden_ratio
        fa = func(xa)
        xb = xmax - delx * self.golden_ratio  
        fb = func(xb)
        
        iterations = 0
        
        # 主搜尋循環 (對應Fortran中的100標籤)
        while iterations < max_iterations:
            delx_save = delx
            
            # 收斂檢查
            if delx < precision:
                break
                
            if fb < fa:
                # 移動到右半區間 (對應Fortran中的200標籤)
                xmin = xa
                delx = xmax - xmin
                
                # 防止無限循環
                if abs(delx - delx_save) < np.finfo(float).eps:
                    break
                    
                xa = xb
                fa = fb
                xb = xmax - delx * self.golden_ratio
                fb = func(xb)
            else:
                # 移動到左半區間
                xmax = xb
                delx = xmax - xmin
                
                # 防止無限循環  
                if abs(delx - delx_save) < np.finfo(float).eps:
                    break
                    
                xb = xa
                fb = fa
                xa = xmin + delx * self.golden_ratio
                fa = func(xa)
                
            iterations += 1
            
        # 返回區間中點作為最優解
        x_opt = (xmin + xmax) / 2
        f_opt = func(x_opt)
        
        if iterations >= max_iterations:
            warnings.warn(f"黃金分割搜尋達到最大迭代次數 {max_iterations}")
            
        return x_opt, f_opt, iterations


# 全域優化器實例
_golden_optimizer = GoldenSectionOptimizer()

def golden_section_search(func: Callable[[float], float],
                         xmin: float,
                         xmax: float, 
                         precision: float = 1e-6,
                         max_iterations: int = 1000) -> Tuple[float, float, int]:
    """
    黃金分割搜尋函數介面
    
    這是對GoldenSectionOptimizer的便捷封裝，直接對應原始Fortran GSECT子程序。
    
    參數:
        func: 目標函數
        xmin: 搜尋下界
        xmax: 搜尋上界
        precision: 收斂精度
        max_iterations: 最大迭代次數
        
    返回:
        tuple: (最優解, 最優值, 迭代次數)
        
    範例:
        >>> def parabola(x):
        ...     return (x - 2)**2 + 1
        >>> x_opt, f_opt, iters = golden_section_search(parabola, 0, 4)
        >>> print(f"最優解: x = {x_opt:.6f}, f(x) = {f_opt:.6f}")
    """
    return _golden_optimizer.optimize(func, xmin, xmax, precision, max_iterations)


def _fermi_dirac_numerical_integration(order: float, eta: float, precision: float) -> float:
    """
    費米-狄拉克積分的高精度數值計算
    
    使用專門針對費米-狄拉克積分優化的數值方法
    """
    def integrand(x):
        # 避免數值溢出
        if x - eta > 700:  # exp(700) 接近浮點數上限
            return 0.0
        elif x - eta < -700:
            return x ** order
        else:
            return (x ** order) / (np.exp(x - eta) + 1)
    
    # 分段積分以提高精度
    # 第一部分：[0, eta + 50]
    integral_1 = adaptive_quadrature(integrand, 0, max(eta + 50, 50), precision/2)
    
    # 第二部分：[eta + 50, ∞] 使用漸近展開
    if eta > -50:
        # 對於大x的漸近行為
        def asymptotic_integrand(x):
            return (x ** order) * np.exp(eta - x)
        
        integral_2 = adaptive_quadrature(asymptotic_integrand, max(eta + 50, 50), 
                                       max(eta + 100, 100), precision/2)
    else:
        integral_2 = 0.0
    
    return integral_1 + integral_2


class FermiDiracIntegrator:
    """
    費米-狄拉克積分器
    
    計算費米-狄拉克分佈相關的積分，這在半導體物理和統計力學中至關重要。
    
    費米-狄拉克積分定義為:
    F_j(η) = ∫[0,∞] (x^j)/(exp(x-η)+1) dx
    
    其中:
    - j 是積分階數
    - η 是約化化學勢 (μ-E_c)/k_B T
    """
    
    @staticmethod
    @validate_numerical_inputs
    def fermi_dirac_integral(order: float, 
                           eta: Union[float, np.ndarray],
                           precision: float = 1e-10) -> Union[float, np.ndarray]:
        """
        計算費米-狄拉克積分 F_j(η)
        
        參數:
            order: 積分階數 j
            eta: 約化化學勢 η
            precision: 積分精度
            
        返回:
            費米-狄拉克積分值
            
        注意: 
        此函數使用數值積分方法，對於特殊情況會使用解析近似。
        在SEMITIP中主要用於計算載流子濃度。
        """
        eta = np.asarray(eta)
        scalar_input = eta.ndim == 0
        eta = np.atleast_1d(eta)
        
        result = np.zeros_like(eta, dtype=float)
        
        for i, eta_val in enumerate(eta):
            if eta_val < -10:
                # 非簡併情況: F_j(η) ≈ Γ(j+1) * exp(η)
                result[i] = math.gamma(order + 1) * np.exp(eta_val)
            elif eta_val > 10:
                # 簡併情況: F_j(η) ≈ η^(j+1)/(j+1) + π²(j)(j-1)η^(j-1)/6 + ...
                # 主項
                main_term = (eta_val ** (order + 1)) / (order + 1)
                # 修正項（當order > 1時）
                if order > 1:
                    correction = (math.pi**2 / 6) * order * (order - 1) * (eta_val ** (order - 1))
                    result[i] = main_term + correction
                else:
                    result[i] = main_term
            else:
                # 一般情況: 使用更精確的數值積分
                result[i] = _fermi_dirac_numerical_integration(order, eta_val, precision)
        
        return result[0] if scalar_input else result


def fermi_dirac_integral(order: float, 
                        eta: Union[float, np.ndarray],
                        precision: float = 1e-10) -> Union[float, np.ndarray]:
    """費米-狄拉克積分的便捷介面"""
    return FermiDiracIntegrator.fermi_dirac_integral(order, eta, precision)


@validate_numerical_inputs
def numerical_derivative(func: Callable[[float], float],
                        x: float,
                        h: float = 1e-5,
                        method: str = 'central') -> float:
    """
    數值微分計算
    
    參數:
        func: 目標函數
        x: 求導點
        h: 步長
        method: 差分方法 ('forward', 'backward', 'central')
        
    返回:
        數值導數值
    """
    if method == 'forward':
        return (func(x + h) - func(x)) / h
    elif method == 'backward':
        return (func(x) - func(x - h)) / h
    elif method == 'central':
        return (func(x + h) - func(x - h)) / (2 * h)
    else:
        raise ValueError(f"未知的差分方法: {method}")


@validate_numerical_inputs  
def adaptive_quadrature(func: Callable[[float], float],
                       a: float,
                       b: float, 
                       tolerance: float = 1e-10,
                       max_depth: int = 15) -> float:
    """
    自適應積分算法
    
    使用Simpson規則的自適應版本進行數值積分。
    
    參數:
        func: 被積函數
        a: 積分下界
        b: 積分上界
        tolerance: 容許誤差
        max_depth: 最大遞歸深度
        
    返回:
        積分值
    """
    def simpson_rule(f, x0, x2):
        """Simpson 1/3規則"""
        x1 = (x0 + x2) / 2
        h = (x2 - x0) / 6
        return h * (f(x0) + 4*f(x1) + f(x2))
    
    def adaptive_simpson(f, x0, x2, tol, depth):
        """遞歸自適應Simpson積分"""
        if depth > max_depth:
            return simpson_rule(f, x0, x2)
            
        x1 = (x0 + x2) / 2
        
        # 計算整個區間和兩個子區間的Simpson積分
        s_whole = simpson_rule(f, x0, x2)
        s_left = simpson_rule(f, x0, x1) 
        s_right = simpson_rule(f, x1, x2)
        s_split = s_left + s_right
        
        # 檢查誤差
        error = abs(s_split - s_whole)
        if error < 15 * tol:  # Richardson外推誤差估計
            return s_split + (s_split - s_whole) / 15
        else:
            # 遞歸細分
            return (adaptive_simpson(f, x0, x1, tol/2, depth+1) + 
                   adaptive_simpson(f, x1, x2, tol/2, depth+1))
    
    return adaptive_simpson(func, a, b, tolerance, 0)


def compute_derivative_matrix(grid_points: np.ndarray,
                            order: int = 2,
                            boundary_conditions: str = 'periodic') -> np.ndarray:
    """
    構建有限差分導數矩陣
    
    這在SEMITIP的泊松求解器中用於構建線性系統。
    
    參數:
        grid_points: 網格點
        order: 差分階數
        boundary_conditions: 邊界條件類型
        
    返回:
        導數矩陣
    """
    n = len(grid_points)
    matrix = np.zeros((n, n))
    
    if order == 2:
        # 二階中心差分
        for i in range(1, n-1):
            h = grid_points[i+1] - grid_points[i]
            matrix[i, i-1] = 1/h**2
            matrix[i, i] = -2/h**2  
            matrix[i, i+1] = 1/h**2
            
        # 處理邊界條件
        if boundary_conditions == 'periodic':
            matrix[0, -1] = matrix[1, 0]
            matrix[0, 0] = matrix[1, 1]
            matrix[0, 1] = matrix[1, 2]
            matrix[-1, -2] = matrix[-2, -3]
            matrix[-1, -1] = matrix[-2, -2]
            matrix[-1, 0] = matrix[-2, -1]
    
    return matrix


def fermi_dirac_occupation(energy: Union[float, np.ndarray], 
                          fermi_energy: float,
                          temperature: float) -> Union[float, np.ndarray]:
    """
    費米-狄拉克佔據函數
    
    對應Fortran中的FD函數：f(E) = 1/(1 + exp((E-EF)/(kB*T)))
    
    參數:
        energy: 電子能量 (eV)
        fermi_energy: 費米能級 (eV)
        temperature: 溫度 (K)，如果為0則使用T=0近似
        
    返回:
        佔據機率 [0,1]
        
    注意: 此函數直接對應Fortran FD函數的邏輯
    """
    energy = np.asarray(energy)
    scalar_input = energy.ndim == 0
    energy = np.atleast_1d(energy)
    
    result = np.zeros_like(energy, dtype=float)
    
    if temperature == 0:
        # T=0 近似
        for i, E in enumerate(energy):
            if E == fermi_energy:
                result[i] = 0.5
            elif E < fermi_energy:
                result[i] = 1.0
            else:
                result[i] = 0.0
    else:
        # 有限溫度
        k_B = 8.617e-5  # eV/K
        for i, E in enumerate(energy):
            x = (E - fermi_energy) / (k_B * temperature)
            if x > 40:
                result[i] = 0.0
            elif x < -40:
                result[i] = 1.0
            else:
                result[i] = 1.0 / (1.0 + np.exp(x))
    
    return result[0] if scalar_input else result


def trapezoidal_integration(func: Callable[[float], float],
                          x_min: float,
                          x_max: float,
                          n_steps: int = 1000) -> float:
    """
    梯形積分法
    
    對應Fortran中的TRAP函數，提供簡單穩定的數值積分。
    
    參數:
        func: 被積函數
        x_min: 積分下界
        x_max: 積分上界
        n_steps: 積分步數
        
    返回:
        積分值
    """
    if n_steps < 1:
        n_steps = 1
        
    if n_steps == 1:
        return (func(x_min) + func(x_max)) * (x_max - x_min) / 2.0
    
    h = (x_max - x_min) / n_steps
    total = (func(x_min) + func(x_max)) / 2.0
    
    for i in range(1, n_steps):
        x = x_min + i * h
        total += func(x)
    
    return total * h


# 便捷別名，對應Fortran函數名
fd = fermi_dirac_occupation  # 對應Fortran FD函數
trap = trapezoidal_integration  # 對應Fortran TRAP函數


if __name__ == "__main__":
    # 基本測試
    print("Pysemitip數值工具模組測試")
    print("=" * 40)
    
    # 測試黃金分割搜尋
    print("1. 黃金分割搜尋測試:")
    def test_func(x):
        return (x - 2.5)**2 + 1.5
    
    x_opt, f_opt, iters = golden_section_search(test_func, 0, 5)
    print(f"   目標函數: (x-2.5)² + 1.5")
    print(f"   最優解: x = {x_opt:.8f}")
    print(f"   最優值: f(x) = {f_opt:.8f}")
    print(f"   迭代次數: {iters}")
    print(f"   理論最優解: x = 2.5, f(x) = 1.5")
    
    # 測試費米-狄拉克積分
    print("\n2. 費米-狄拉克積分測試:")
    eta_values = [-5, 0, 5]
    for eta in eta_values:
        f_val = fermi_dirac_integral(0.5, eta)
        print(f"   F_0.5({eta}) = {f_val:.8f}")
    
    print("\n測試完成！")
