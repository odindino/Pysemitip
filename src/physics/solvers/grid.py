import numpy as np

class HyperbolicGrid:
    """
    生成用於 STM 模擬的雙曲面網格 (prolate spheroidal coordinates)。

    此網格旨在擬合 STM 針尖的雙曲面形狀，從而在針尖-樣品間隙中提供高解析度，
    在遠離該區域時提供較粗的解析度。網格由座標 (eta, nu) 定義。

    Args:
        N_eta (int): 沿 eta 方向 (雙曲面) 的網格點數。
        N_nu (int): 沿 nu 方向 (橢球面) 的網格點數。
        R (float): 針尖的曲率半徑 (nm)。
        Z_TS (float): 針尖-樣品間的距離 (nm)。
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
    """
    def __init__(self, N_eta, N_nu, R, Z_TS, r_max_factor=5.0):
        self.N_eta = N_eta
        self.N_nu = N_nu
        
        if R <= 0:
            raise ValueError("針尖曲率半徑 R 必須為正數。")
        self.R = R
        
        if Z_TS <= 0:
            raise ValueError("針尖-樣品距離 Z_TS 必須為正數。")
        self.Z_TS = Z_TS
        
        # 物理限制：此模型要求針尖-樣品距離大於曲率半徑
        if Z_TS <= R:
            raise ValueError(f"此模型要求 Z_TS > R，但 Z_TS={Z_TS}, R={R}。")
            
        self.r_max_factor = r_max_factor

        # 內部參數
        self.f = 0.0
        self.eta_tip = 0.0
        self.eta_grid_max = 0.0
        
        # 網格陣列
        self.eta_grid = None
        self.nu = None
        self.r = None
        self.z = None

        self._generate_grid()

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
        self.z = self.f * np.cosh(eta_actual) * np.cos(self.nu)

        # 驗證邊界條件 (可選，用於除錯)
        # 針尖尖端 (eta_grid=0, nu=0) 的 z 座標應接近 Z_TS
        z_apex = self.z[0, 0]
        # 樣品表面 (nu = N_nu-1) 的 z 座標應接近 0
        z_sample_max = np.max(np.abs(self.z[:, -1]))
        
        # print(f"驗證: 針尖尖端 z = {z_apex:.4f} (目標: {self.Z_TS})")
        # print(f"驗證: 樣品表面最大 z = {z_sample_max:.4g} (目標: 0)")