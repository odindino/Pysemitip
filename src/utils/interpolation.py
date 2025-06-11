"""
插值工具模組

此模組提供SEMITIP專案所需的各種插值算法，用於處理離散數據點之間的數值插值。

主要功能:
1. 線性插值 - 一維線性插值
2. 三次樣條插值 - 平滑的三次多項式插值  
3. 雙線性插值 - 二維矩形網格插值
4. 三維三線性插值 - 三維網格插值

這些插值方法在SEMITIP中用於:
- 電位場的空間插值
- 材料參數的溫度/組分插值
- 網格點之間的物理量計算
"""

import numpy as np
from typing import Union, Tuple, Optional, List
from scipy.interpolate import interp1d, RectBivariateSpline
import warnings


class LinearInterpolator:
    """
    一維線性插值器
    
    實現簡單高效的線性插值，特別適用於SEMITIP中的
    電位分佈和材料參數插值。
    """
    
    def __init__(self, x_data: np.ndarray, y_data: np.ndarray, 
                 extrapolate: bool = False):
        """
        初始化線性插值器
        
        參數:
            x_data: x坐標數組（必須單調遞增）
            y_data: y坐標數組
            extrapolate: 是否允許外推
        """
        x_data = np.asarray(x_data)
        y_data = np.asarray(y_data)
        
        if len(x_data) != len(y_data):
            raise ValueError("x_data和y_data長度必須相等")
            
        if len(x_data) < 2:
            raise ValueError("至少需要2個數據點進行插值")
            
        # 檢查x_data是否單調遞增
        if not np.all(np.diff(x_data) > 0):
            # 自動排序
            sort_indices = np.argsort(x_data)
            x_data = x_data[sort_indices]
            y_data = y_data[sort_indices]
            warnings.warn("x_data已自動排序為遞增順序")
            
        self.x_data = x_data
        self.y_data = y_data
        self.extrapolate = extrapolate
        
    def __call__(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        執行線性插值
        
        參數:
            x: 插值點
            
        返回:
            插值結果
        """
        x = np.asarray(x)
        scalar_input = x.ndim == 0
        x = np.atleast_1d(x)
        
        result = np.zeros_like(x, dtype=float)
        
        for i, xi in enumerate(x):
            if xi < self.x_data[0]:
                if self.extrapolate:
                    # 左外推
                    slope = (self.y_data[1] - self.y_data[0]) / (self.x_data[1] - self.x_data[0])
                    result[i] = self.y_data[0] + slope * (xi - self.x_data[0])
                else:
                    result[i] = self.y_data[0]
            elif xi > self.x_data[-1]:
                if self.extrapolate:
                    # 右外推
                    slope = (self.y_data[-1] - self.y_data[-2]) / (self.x_data[-1] - self.x_data[-2])
                    result[i] = self.y_data[-1] + slope * (xi - self.x_data[-1])
                else:
                    result[i] = self.y_data[-1]
            else:
                # 內插
                # 找到插值區間
                idx = np.searchsorted(self.x_data, xi) - 1
                if idx < 0:
                    idx = 0
                elif idx >= len(self.x_data) - 1:
                    idx = len(self.x_data) - 2
                    
                # 線性插值公式
                x0, x1 = self.x_data[idx], self.x_data[idx + 1]
                y0, y1 = self.y_data[idx], self.y_data[idx + 1]
                
                t = (xi - x0) / (x1 - x0)
                result[i] = y0 + t * (y1 - y0)
                
        return result[0] if scalar_input else result


class CubicSplineInterpolator:
    """
    三次樣條插值器
    
    提供平滑的三次多項式插值，適用於需要高階平滑性的場合，
    如SEMITIP中的電位場平滑處理。
    """
    
    def __init__(self, x_data: np.ndarray, y_data: np.ndarray,
                 boundary_condition: str = 'natural',
                 derivative_values: Optional[Tuple[float, float]] = None):
        """
        初始化三次樣條插值器
        
        參數:
            x_data: x坐標數組
            y_data: y坐標數組  
            boundary_condition: 邊界條件 ('natural', 'clamped', 'periodic')
            derivative_values: 當boundary_condition='clamped'時的端點導數值
        """
        x_data = np.asarray(x_data)
        y_data = np.asarray(y_data)
        
        if len(x_data) != len(y_data):
            raise ValueError("x_data和y_data長度必須相等")
            
        if len(x_data) < 3:
            raise ValueError("三次樣條插值至少需要3個數據點")
            
        # 排序數據
        sort_indices = np.argsort(x_data)
        self.x_data = x_data[sort_indices]
        self.y_data = y_data[sort_indices]
        
        self.boundary_condition = boundary_condition
        self.derivative_values = derivative_values
        
        # 計算樣條係數
        self._compute_spline_coefficients()
        
    def _compute_spline_coefficients(self):
        """計算三次樣條的係數"""
        n = len(self.x_data)
        h = np.diff(self.x_data)  # 區間長度
        
        # 建立三對角線系統求解二階導數
        A = np.zeros((n, n))
        b = np.zeros(n)
        
        # 內部節點方程
        for i in range(1, n-1):
            A[i, i-1] = h[i-1]
            A[i, i] = 2 * (h[i-1] + h[i])
            A[i, i+1] = h[i]
            b[i] = 6 * ((self.y_data[i+1] - self.y_data[i]) / h[i] - 
                       (self.y_data[i] - self.y_data[i-1]) / h[i-1])
        
        # 邊界條件
        if self.boundary_condition == 'natural':
            # 自然邊界: S''(x0) = S''(xn) = 0
            A[0, 0] = 1
            A[-1, -1] = 1
            b[0] = b[-1] = 0
        elif self.boundary_condition == 'clamped':
            # 固定邊界: S'(x0) = y0', S'(xn) = yn'
            if self.derivative_values is None:
                raise ValueError("固定邊界條件需要提供導數值")
            dy0, dyn = self.derivative_values
            
            A[0, 0] = 2 * h[0]
            A[0, 1] = h[0]
            b[0] = 6 * ((self.y_data[1] - self.y_data[0]) / h[0] - dy0)
            
            A[-1, -2] = h[-1]
            A[-1, -1] = 2 * h[-1]
            b[-1] = 6 * (dyn - (self.y_data[-1] - self.y_data[-2]) / h[-1])
        
        # 求解二階導數
        self.second_derivatives = np.linalg.solve(A, b)
        
    def __call__(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """執行三次樣條插值"""
        x = np.asarray(x)
        scalar_input = x.ndim == 0
        x = np.atleast_1d(x)
        
        result = np.zeros_like(x, dtype=float)
        
        for i, xi in enumerate(x):
            # 找到所在區間
            if xi <= self.x_data[0]:
                result[i] = self.y_data[0]
            elif xi >= self.x_data[-1]:
                result[i] = self.y_data[-1]
            else:
                # 找到插值區間
                idx = np.searchsorted(self.x_data, xi) - 1
                
                x0, x1 = self.x_data[idx], self.x_data[idx + 1]
                y0, y1 = self.y_data[idx], self.y_data[idx + 1]
                M0, M1 = self.second_derivatives[idx], self.second_derivatives[idx + 1]
                
                h = x1 - x0
                t = xi - x0
                
                # 三次樣條公式
                result[i] = (y0 * (x1 - xi) + y1 * (xi - x0)) / h - \
                           h * ((M0 * (x1 - xi) + M1 * (xi - x0)) / 6) + \
                           (M0 * (x1 - xi)**3 + M1 * (xi - x0)**3) / (6 * h)
                
        return result[0] if scalar_input else result


class BilinearInterpolator:
    """
    雙線性插值器
    
    用於二維矩形網格的插值，在SEMITIP中用於處理
    二維電位場和材料參數的空間分佈。
    """
    
    def __init__(self, x_coords: np.ndarray, y_coords: np.ndarray, 
                 z_values: np.ndarray):
        """
        初始化雙線性插值器
        
        參數:
            x_coords: x坐標數組 (長度為m)
            y_coords: y坐標數組 (長度為n)  
            z_values: z值數組 (形狀為(n, m))
        """
        self.x_coords = np.asarray(x_coords)
        self.y_coords = np.asarray(y_coords)
        self.z_values = np.asarray(z_values)
        
        if self.z_values.shape != (len(y_coords), len(x_coords)):
            raise ValueError("z_values形狀必須為(len(y_coords), len(x_coords))")
            
        # 檢查坐標是否單調
        if not (np.all(np.diff(self.x_coords) > 0) and 
                np.all(np.diff(self.y_coords) > 0)):
            raise ValueError("坐標數組必須嚴格遞增")
            
    def __call__(self, x: Union[float, np.ndarray], 
                 y: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        執行雙線性插值
        
        參數:
            x: x坐標
            y: y坐標
            
        返回:
            插值結果
        """
        x = np.asarray(x)
        y = np.asarray(y)
        
        # 處理標量輸入
        scalar_input = x.ndim == 0 and y.ndim == 0
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        
        if len(x) != len(y):
            raise ValueError("x和y數組長度必須相等")
            
        result = np.zeros_like(x, dtype=float)
        
        for i in range(len(x)):
            xi, yi = x[i], y[i]
            
            # 邊界處理
            if xi <= self.x_coords[0]:
                x_idx = 0
                x_t = 0
            elif xi >= self.x_coords[-1]:
                x_idx = len(self.x_coords) - 2
                x_t = 1
            else:
                x_idx = np.searchsorted(self.x_coords, xi) - 1
                x_t = (xi - self.x_coords[x_idx]) / (self.x_coords[x_idx + 1] - self.x_coords[x_idx])
                
            if yi <= self.y_coords[0]:
                y_idx = 0
                y_t = 0
            elif yi >= self.y_coords[-1]:
                y_idx = len(self.y_coords) - 2
                y_t = 1
            else:
                y_idx = np.searchsorted(self.y_coords, yi) - 1
                y_t = (yi - self.y_coords[y_idx]) / (self.y_coords[y_idx + 1] - self.y_coords[y_idx])
            
            # 雙線性插值
            z00 = self.z_values[y_idx, x_idx]
            z01 = self.z_values[y_idx, x_idx + 1]
            z10 = self.z_values[y_idx + 1, x_idx]
            z11 = self.z_values[y_idx + 1, x_idx + 1]
            
            # 先在x方向插值
            z0 = z00 * (1 - x_t) + z01 * x_t
            z1 = z10 * (1 - x_t) + z11 * x_t
            
            # 再在y方向插值
            result[i] = z0 * (1 - y_t) + z1 * y_t
            
        return result[0] if scalar_input else result


class TrilinearInterpolator:
    """
    三線性插值器
    
    用於三維規則網格的插值，適用於SEMITIP中三維電位場
    和電荷密度分佈的插值計算。
    """
    
    def __init__(self, x_coords: np.ndarray, y_coords: np.ndarray,
                 z_coords: np.ndarray, values: np.ndarray):
        """
        初始化三線性插值器
        
        參數:
            x_coords: x坐標數組 (長度為l)
            y_coords: y坐標數組 (長度為m)
            z_coords: z坐標數組 (長度為n)
            values: 值數組 (形狀為(n, m, l))
        """
        self.x_coords = np.asarray(x_coords)
        self.y_coords = np.asarray(y_coords)
        self.z_coords = np.asarray(z_coords)
        self.values = np.asarray(values)
        
        expected_shape = (len(z_coords), len(y_coords), len(x_coords))
        if self.values.shape != expected_shape:
            raise ValueError(f"values形狀必須為{expected_shape}")
            
    def __call__(self, x: float, y: float, z: float) -> float:
        """執行三線性插值"""
        # 找到插值立方體
        x_idx = np.searchsorted(self.x_coords, x) - 1
        y_idx = np.searchsorted(self.y_coords, y) - 1
        z_idx = np.searchsorted(self.z_coords, z) - 1
        
        # 邊界處理
        x_idx = max(0, min(x_idx, len(self.x_coords) - 2))
        y_idx = max(0, min(y_idx, len(self.y_coords) - 2))
        z_idx = max(0, min(z_idx, len(self.z_coords) - 2))
        
        # 計算插值參數
        x_t = (x - self.x_coords[x_idx]) / (self.x_coords[x_idx + 1] - self.x_coords[x_idx])
        y_t = (y - self.y_coords[y_idx]) / (self.y_coords[y_idx + 1] - self.y_coords[y_idx])
        z_t = (z - self.z_coords[z_idx]) / (self.z_coords[z_idx + 1] - self.z_coords[z_idx])
        
        # 限制在[0,1]範圍內
        x_t = max(0, min(1, x_t))
        y_t = max(0, min(1, y_t))
        z_t = max(0, min(1, z_t))
        
        # 獲取立方體8個頂點的值
        c000 = self.values[z_idx, y_idx, x_idx]
        c001 = self.values[z_idx, y_idx, x_idx + 1]
        c010 = self.values[z_idx, y_idx + 1, x_idx]
        c011 = self.values[z_idx, y_idx + 1, x_idx + 1]
        c100 = self.values[z_idx + 1, y_idx, x_idx]
        c101 = self.values[z_idx + 1, y_idx, x_idx + 1]
        c110 = self.values[z_idx + 1, y_idx + 1, x_idx]
        c111 = self.values[z_idx + 1, y_idx + 1, x_idx + 1]
        
        # 三線性插值
        c00 = c000 * (1 - x_t) + c001 * x_t
        c01 = c010 * (1 - x_t) + c011 * x_t
        c10 = c100 * (1 - x_t) + c101 * x_t
        c11 = c110 * (1 - x_t) + c111 * x_t
        
        c0 = c00 * (1 - y_t) + c01 * y_t
        c1 = c10 * (1 - y_t) + c11 * y_t
        
        return c0 * (1 - z_t) + c1 * z_t


# 便捷函數介面
def linear_interpolation(x_data: np.ndarray, y_data: np.ndarray,
                        x_new: Union[float, np.ndarray],
                        extrapolate: bool = False) -> Union[float, np.ndarray]:
    """
    線性插值便捷函數
    
    參數:
        x_data: 已知x坐標
        y_data: 已知y坐標
        x_new: 待插值x坐標
        extrapolate: 是否允許外推
        
    返回:
        插值結果
    """
    interpolator = LinearInterpolator(x_data, y_data, extrapolate)
    return interpolator(x_new)


def cubic_spline_interpolation(x_data: np.ndarray, y_data: np.ndarray,
                             x_new: Union[float, np.ndarray],
                             boundary_condition: str = 'natural') -> Union[float, np.ndarray]:
    """
    三次樣條插值便捷函數
    
    參數:
        x_data: 已知x坐標
        y_data: 已知y坐標  
        x_new: 待插值x坐標
        boundary_condition: 邊界條件
        
    返回:
        插值結果
    """
    interpolator = CubicSplineInterpolator(x_data, y_data, boundary_condition)
    return interpolator(x_new)


def bilinear_interpolation(x_coords: np.ndarray, y_coords: np.ndarray,
                          z_values: np.ndarray, x_new: Union[float, np.ndarray],
                          y_new: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    雙線性插值便捷函數
    
    參數:
        x_coords: x坐標網格
        y_coords: y坐標網格
        z_values: z值矩陣
        x_new: 待插值x坐標
        y_new: 待插值y坐標
        
    返回:
        插值結果
    """
    interpolator = BilinearInterpolator(x_coords, y_coords, z_values)
    return interpolator(x_new, y_new)


if __name__ == "__main__":
    # 基本測試
    print("Pysemitip插值工具模組測試")
    print("=" * 40)
    
    # 測試線性插值
    print("1. 線性插值測試:")
    x_data = np.array([0, 1, 2, 3, 4])
    y_data = np.array([0, 1, 4, 9, 16])  # y = x²
    x_new = np.array([0.5, 1.5, 2.5])
    y_interp = linear_interpolation(x_data, y_data, x_new)
    print(f"   x_new: {x_new}")
    print(f"   y_interp: {y_interp}")
    print(f"   預期值: [0.5, 2.5, 6.5]")
    
    # 測試雙線性插值
    print("\n2. 雙線性插值測試:")
    x_coords = np.array([0, 1, 2])
    y_coords = np.array([0, 1, 2])
    z_values = np.array([[0, 1, 4],    # y=0: z=x²
                        [1, 2, 5],     # y=1: z=x²+1
                        [4, 5, 8]])    # y=2: z=x²+4
    
    result = bilinear_interpolation(x_coords, y_coords, z_values, 0.5, 0.5)
    print(f"   插值點(0.5, 0.5)的值: {result}")
    print(f"   預期值: 1.25")
    
    print("\n測試完成！")
