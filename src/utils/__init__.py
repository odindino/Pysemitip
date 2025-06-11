"""
Pysemitip基礎工具模組

此模組包含SEMITIP專案的基礎數值計算工具，遵循「物理優先，架構為本」的設計原則。

子模組:
- numerical: 核心數值計算函數（黃金分割優化、費米-狄拉克積分等）
- interpolation: 插值工具函數
"""

__version__ = "1.0.0"
__author__ = "odindino"

# 導入核心數值工具
from .numerical import (
    golden_section_search,
    fermi_dirac_integral,
    fermi_dirac_occupation,
    numerical_derivative,
    adaptive_quadrature,
    trapezoidal_integration
)

from .interpolation import (
    linear_interpolation,
    cubic_spline_interpolation,
    bilinear_interpolation
)

__all__ = [
    'golden_section_search',
    'fermi_dirac_integral',
    'fermi_dirac_occupation', 
    'numerical_derivative',
    'adaptive_quadrature',
    'trapezoidal_integration',
    'linear_interpolation',
    'cubic_spline_interpolation',
    'bilinear_interpolation'
]
