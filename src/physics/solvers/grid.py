import numpy as np
from typing import List, Tuple, Optional, TYPE_CHECKING
from dataclasses import dataclass

if TYPE_CHECKING:
    def set_material_regions(self, regions_config):
        """
        Sets the semiconductor material region configurations for the grid.
        These configurations are used by get_region_id_at_point and to build the region_id_map.

        Args:
            regions_config: A list of semiconductor region configuration objects.
                            Each object is assumed to have 'id', and optionally 
                            'start_depth_nm' and 'end_depth_nm' attributes.
                            For MultInt, if depth attributes are missing, regions are 
                            assigned based on grid geometry.
        """
        self.semiconductor_regions_config = regions_config
        # For MultInt compatibility: only sort if all regions have depth attributes
        if self.semiconductor_regions_config and all(
            hasattr(r, 'start_depth_nm') and hasattr(r, 'end_depth_nm') 
            for r in self.semiconductor_regions_config
        ):
            self.semiconductor_regions_config.sort(key=lambda r: r.start_depth_nm)
        
        self._build_region_id_map() # Build the map after setting new configsonfigSemiconductorRegion will be imported from its actual location
    # For now, using a forward reference or a more generic type hint if not directly available
    from ...core.config_schema import SemiconductorRegion as ConfigSemiconductorRegion

@dataclass
class GridParameters:
    """參數類別用於定義 Grid3D 的配置"""
    nr: int = 32
    nv: int = 32
    ns: int = 32
    nphi: int = 16  # 重命名避免與 numpy 的 np 衝突
    delr: float = 1.0
    delv: float = 1.0
    dels: float = 1.0
    delp: float = np.pi/16
    rmax: float = 31.0
    vmax: float = 31.0
    smax: float = 31.0

class Grid3D:
    """
    3D 網格類別，與測試兼容
    """
    def __init__(self, params: GridParameters):
        """
        初始化 3D 網格
        
        Args:
            params: 網格參數
        """
        self.params = params
        
        # Generate coordinate arrays
        self.r = np.linspace(0, params.rmax, params.nr)
        self.zv = np.linspace(0, params.vmax, params.nv)
        self.zs = np.linspace(0, -params.smax, params.ns)  # Negative for semiconductor region
        self.phi = np.linspace(0, params.delp * (params.nphi - 1), params.nphi)
        
    def __repr__(self):
        return f"Grid3D(nr={self.params.nr}, nv={self.params.nv}, ns={self.params.ns}, nphi={self.params.nphi})"

def create_grid_from_config(*args, **kwargs):
    """創建網格的輔助函數"""
    # Placeholder implementation
    params = GridParameters()
    return Grid3D(params)

def generate_prolate_spheroidal_grid(*args, **kwargs):
    """生成prolate spheroidal網格的輔助函數"""
    # Placeholder implementation  
    params = GridParameters()
    return Grid3D(params)

class HyperbolicGrid:
    """
    生成用於 STM 模擬的雙曲面網格 (prolate spheroidal coordinates)。

    此網格旨在擬合 STM 針尖的雙曲面形狀，從而在針尖-樣品間隙中提供高解析度，
    在遠離該區域時提供較粗的解析度。網格由座標 (eta, nu) 定義。

    Args:
        N_eta (int): 沿 eta 方向 (雙曲面) 的網格點數。
        N_nu (int): 沿 nu 方向 (橢球面) 的網格點數。
        R (float): 針尖的曲率半徑 (nm)。
        Z_TS (float): 針尖-樣品之間的距離 (nm)。
        r_max_factor (float, optional): 用於確定模擬區域徑向範圍的因子。
                                        模擬半徑將是 r_max_factor * R。預設為 5.0。

    Attributes:
        N_eta, N_nu (int): 網格維度。
        R, Z_TS (float): 物理參數。
        f (float): 座標系統的焦點距離。
        eta_tip (float): 描述針尖表面的 eta 座標值。
        eta_grid_max (float): 網格座標 eta_grid 的最大值。
        eta_grid (np.ndarray): 2D 網格座標 eta_grid (從 0 開始，代表針尖)。
        nu (np.ndarray): 2D 網格座標 nu。
        r (np.ndarray): 轉換後的 2D 笛卡爾 r 座標。
        z (np.ndarray): 轉換後的 2D 笛卡爾 z 座標。
        semiconductor_regions_config (List[ConfigSemiconductorRegion], optional): Configuration for semiconductor regions.
    """
    def __init__(self, N_eta, N_nu, R, Z_TS, r_max_factor=5.0):
        self.N_eta = N_eta
        self.N_nu = N_nu
        
        if R <= 0:
            raise ValueError("針尖曲率半徑 R 必須為正。")
        self.R = R
        
        if Z_TS <= 0:
            raise ValueError("針尖-樣品距離 Z_TS 必須為正。")
        self.Z_TS = Z_TS
        
        # 物理限制：此模型要求針尖-樣品距離大於曲率半徑
        # This check was present in the original prompt's context for multint.py,
        # but might be too restrictive for the grid generation itself.
        # For now, I'll keep it as it was implied to be a requirement.
        if Z_TS <= R:
            # Consider if this should be a warning or allow Z_TS = R for specific cases.
            # For hyperbolic model, Z_TS > R is generally assumed for f to be real.
            # However, the formulas for f = sqrt(Z_TS^2 - R^2) implies Z_TS >= R.
            # If Z_TS = R, f = 0, which is a sphere, not a hyperboloid.
            # Let's assume Z_TS > R is a strict requirement for this grid type.
            raise ValueError(f"針尖-樣品距離 Z_TS ({Z_TS} nm) 必須大於針尖半徑 R ({R} nm) "
                             "以形成雙曲面座標。")
            
        self.r_max_factor = r_max_factor

        # 內部參數
        self.f = 0.0
        self.eta_tip = 0.0
        self.eta_grid_max = 0.0
        
        # 網格陣列
        self.eta_grid = None
        self.nu = None
        self.r = None
        self.z = None # Assuming z=0 is sample surface, positive z into sample

        self.semiconductor_regions_config: Optional[List['ConfigSemiconductorRegion']] = None
        self.region_id_map = None # For storing region IDs for each grid point

        self._generate_grid()

    def set_material_regions(self, regions_config: List['ConfigSemiconductorRegion']):
        """
        Sets the semiconductor material region configurations for the grid.
        These configurations are used by get_region_id_at_point and to build the region_id_map.

        Args:
            regions_config: A list of semiconductor region configuration objects.
                            Each object is assumed to have 'id', and optionally 
                            'start_depth_nm' and 'end_depth_nm' attributes.
                            For MultInt, if depth attributes are missing, regions are 
                            assigned based on grid geometry.
        """
        self.semiconductor_regions_config = regions_config
        # For MultInt compatibility: only sort if all regions have depth attributes
        if self.semiconductor_regions_config and all(
            hasattr(r, 'start_depth_nm') and hasattr(r, 'end_depth_nm') 
            for r in self.semiconductor_regions_config
        ):
            self.semiconductor_regions_config.sort(key=lambda r: r.start_depth_nm)
        
        self._build_region_id_map() # Build the map after setting new configs

    def _build_region_id_map(self):
        """
        Pre-computes the region ID for each grid point and stores it in self.region_id_map.
        This map is then used by get_region_id_at_point for quick lookups.
        """
        if self.z is None or self.semiconductor_regions_config is None:
            self.region_id_map = np.full_like(self.z, fill_value=None, dtype=object)
            return

        self.region_id_map = np.full_like(self.z, fill_value=None, dtype=object)
        
        # Check if regions have depth attributes (for MultPlane) or use simple assignment (for MultInt)
        has_depth_attrs = all(
            hasattr(r, 'start_depth_nm') and hasattr(r, 'end_depth_nm') 
            for r in self.semiconductor_regions_config
        )
        
        for eta_idx in range(self.N_eta):
            for nu_idx in range(self.N_nu):
                z_coord = self.z[eta_idx, nu_idx]
                assigned_region_id = None
                
                if has_depth_attrs:
                    # Use depth-based assignment for MultPlane
                    for region_config in self.semiconductor_regions_config:
                        if region_config.start_depth_nm <= z_coord < region_config.end_depth_nm:
                            assigned_region_id = region_config.id
                            break # First match based on sorted order
                else:
                    # For MultInt: use simple assignment based on grid position or region ID
                    # Default behavior: assign regions based on z-coordinate and number of regions
                    if len(self.semiconductor_regions_config) == 1:
                        assigned_region_id = self.semiconductor_regions_config[0].id
                    elif len(self.semiconductor_regions_config) == 2:
                        # Simple 2-region assignment: region 1 for z < 0, region 2 for z >= 0
                        if z_coord < 0:
                            assigned_region_id = self.semiconductor_regions_config[0].id
                        else:
                            assigned_region_id = self.semiconductor_regions_config[1].id
                    else:
                        # For more complex cases, assign first region as default
                        assigned_region_id = self.semiconductor_regions_config[0].id
                        
                self.region_id_map[eta_idx, nu_idx] = assigned_region_id


    def get_region_id_at_point(self, eta_idx: int, nu_idx: int) -> Optional[int]:
        """
        Determines the semiconductor region ID for a given grid point using the pre-computed map.

        Args:
            eta_idx: Index along the eta dimension.
            nu_idx: Index along the nu dimension.

        Returns:
            The ID of the semiconductor region if the point falls within one,
            otherwise None.
        """
        if self.region_id_map is None:
            # This case should ideally not be hit if _build_region_id_map is called appropriately
            # Fallback or raise error, for now, attempt to build it if not present.
            # logger.warning("region_id_map was not built. Attempting to build now.")
            self._build_region_id_map()
            if self.region_id_map is None: # Still None after attempt
                 return None

        # Check bounds to prevent IndexError
        if 0 <= eta_idx < self.N_eta and 0 <= nu_idx < self.N_nu:
            return self.region_id_map[eta_idx, nu_idx]
        else:
            # logger.error(f"Indices ({eta_idx}, {nu_idx}) out of bounds for grid ({self.N_eta}, {self.N_nu}).")
            return None # Indices out of bounds

    def _calculate_physical_parameters(self):
        """
        根據物理輸入 (R, Z_TS) 計算座標系的內部參數 (f, eta_tip)。
        
        推導:
        1. 座標變換:
           r = f * sinh(eta) * sin(nu)
           z = f * cosh(eta) * cos(nu)
        2. 樣品表面在 z=0，對應 nu = pi/2。此變換滿足此條件。
        3. 針尖表面為 eta = eta_tip，其尖端 (nu=0) 在 z = Z_TS。
           => Z_TS = f * cosh(eta_tip)
        4. 針尖在尖端的曲率半徑為 R。
           => R = Z_TS * tanh^2(eta_tip)
        5. 聯立求解 f 和 eta_tip。
        """
        # 求解 eta_tip
        tanh_eta_tip_sq = self.R / self.Z_TS
        self.eta_tip = np.arctanh(np.sqrt(tanh_eta_tip_sq))
        
        # 求解焦點距離 f
        cosh_eta_tip = 1 / np.sqrt(1 - tanh_eta_tip_sq)
        self.f = self.Z_TS / cosh_eta_tip

    def _calculate_eta_grid_max(self):
        """
        計算網格座標 eta_grid 的最大值，以確定網格的外部邊界。
        邊界由最大半徑 r_max = r_max_factor * R 定義。
        在樣品表面 (z=0, nu=pi/2)，r = f * sinh(eta)。
        因此，eta_max = asinh(r_max / f)。
        由於我們的網格從 eta_tip 開始，所以 eta_grid_max = eta_max - eta_tip。
        """
        r_max = self.r_max_factor * self.R
        
        if self.f < 1e-9:
            # 避免除以零的極端情況
            raise ValueError("焦點距離 f 過小，無法建立網格。")

        # 計算 r_max 對應的總 eta 值
        eta_total_max = np.arcsinh(r_max / self.f)
        
        # 從總 eta 中減去針尖的 eta_tip，得到網格的範圍
        self.eta_grid_max = eta_total_max - self.eta_tip

        if self.eta_grid_max <= 0:
            raise ValueError("計算出的網格範圍 (eta_grid_max) 為負或零。請嘗試增加 r_max_factor。")

    def _generate_grid(self):
        """
        生成 eta 和 nu 座標陣列，並將它們轉換為笛卡爾座標 (r, z)。
        """
        # 1. 計算物理參數 f 和 eta_tip
        self._calculate_physical_parameters()
        
        # 2. 計算網格邊界 eta_grid_max
        self._calculate_eta_grid_max()

        # 3. 創建一維線性間隔的網格座標
        # eta_grid 從 0 (針尖表面) 開始
        # nu 從 0 (中心軸) 到 pi/2 (樣品表面)
        eta_1d = np.linspace(0, self.eta_grid_max, self.N_eta)
        nu_1d = np.linspace(0, np.pi / 2, self.N_nu)

        # 4. 創建 2D 網格
        self.eta_grid, self.nu = np.meshgrid(eta_1d, nu_1d, indexing='ij')

        # 5. 執行座標變換
        # 實際的 eta 值需要加上針尖的偏移量 eta_tip
        eta_actual = self.eta_grid + self.eta_tip
        
        # 套用雙曲面座標變換公式
        self.r = self.f * np.sinh(eta_actual) * np.sin(self.nu)
        self.z = self.f * np.cosh(eta_actual) * np.cos(self.nu) - self.Z_TS - self.R # To make z=0 at sample surface on axis

        # 驗證邊界條件 (可選，用於除錯)
        # 針尖尖端 (eta_grid=0, nu=0) 的 z 座標應接近 Z_TS
        z_apex = self.z[0, 0]
        # 樣品表面 (nu = N_nu-1) 的 z 座標應接近 0
        z_sample_max = np.max(np.abs(self.z[:, -1]))
        
        # print(f"驗證: 針尖尖端 z = {z_apex:.4f} (目標: {self.Z_TS})")
        # print(f"驗證: 樣品表面最大 z = {z_sample_max:.4g} (目標: 0)")