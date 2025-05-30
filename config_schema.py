"""
SEMITIP 配置資料結構定義

此模組定義了所有 SEMITIP 配置相關的資料類別，包括：
- 探針配置 (TipConfig)
- 半導體區域配置 (SemiconductorRegion)
- 表面態配置 (SurfaceRegion)
- 網格配置 (GridConfig)
- 計算配置 (ComputationConfig)
- 電壓掃描配置 (VoltageScanConfig)
- 模擬類型特有配置 (MultIntConfig, MultPlaneConfig)
- 主配置類別 (SemitipConfig)

作者: Pysemitip 開發團隊
版本: 1.0
日期: 2025-05-30
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class TipConfig:
    """探針配置參數"""
    shank_slope: float = 1.0
    separation: float = 1.0  # nm
    radius: float = 1.0  # nm
    protrusion_radius: float = 0.0  # nm
    contact_potential: float = 0.0  # eV
    x_position: float = 0.0  # nm
    y_position: float = 0.0  # nm
    fermi_energy: float = 8.0  # eV


@dataclass
class EffectiveMass:
    """有效質量參數"""
    conduction_band: float = 0.0635
    heavy_hole: float = 0.643
    light_hole: float = 0.081
    split_off_hole: float = 0.172


@dataclass
class SemiconductorRegion:
    """半導體區域參數"""
    id: int = 1
    donor_concentration: float = 1.0e18  # cm^-3
    acceptor_concentration: float = 0.0  # cm^-3
    band_gap: float = 1.42  # eV
    valence_band_offset: float = 0.0  # eV
    donor_binding_energy: float = 0.006  # eV
    acceptor_binding_energy: float = 0.028  # eV
    effective_mass: EffectiveMass = field(default_factory=EffectiveMass)
    spin_orbit_splitting: float = 0.341  # eV
    degeneracy_indicator: int = 0
    inversion_indicator: int = 0


@dataclass
class SurfaceDistribution:
    """表面態分佈參數"""
    density: float = 0.0  # cm^-2 eV^-1
    neutrality_level: float = 0.0  # eV
    fwhm: float = 0.0  # eV
    centroid_energy: float = 0.0  # eV


@dataclass
class SurfaceRegion:
    """表面區域參數"""
    id: int = 1
    first_distribution: SurfaceDistribution = field(default_factory=SurfaceDistribution)
    second_distribution: SurfaceDistribution = field(default_factory=SurfaceDistribution)


@dataclass
class GridConfig:
    """計算網格配置"""
    mirror_plane: bool = True
    radial_points: int = 32
    vacuum_points: int = 16
    semiconductor_points: int = 32
    angular_points: int = 16
    initial_grid_size: float = 0.5


@dataclass
class ComputationConfig:
    """計算參數配置"""
    scaling_steps: int = 3
    max_iterations: List[int] = field(default_factory=lambda: [80000, 40000, 20000, 10000])
    convergence_parameters: List[float] = field(default_factory=lambda: [1.0e-3, 1.0e-3, 1.0e-4, 1.0e-4])
    charge_density_table_size: int = 20000


@dataclass
class VoltageScanConfig:
    """電壓掃描配置"""
    points: int = 41
    start_voltage: float = -2.0  # V
    end_voltage: float = 2.0  # V
    modulation_voltage: float = 0.050  # V
    negative_ramp: float = 0.0  # nm/V
    positive_ramp: float = 0.0  # nm/V
    custom_voltages: Optional[List[float]] = None


@dataclass
class MultIntConfig:
    """MultInt 特有參數"""
    parallel_wavevectors: int = 20
    energy_points: int = 20
    expansion_factor: int = 20
    semiconductor_depth_fraction: float = 0.75


@dataclass
class MultPlaneConfig:
    """MultPlane 特有參數"""
    vacuum_width: float = 2.0  # nm
    vacuum_spacing: float = 0.05  # nm
    max_energies: Dict[str, float] = field(default_factory=lambda: {
        'light_hole': 6.0,
        'heavy_hole': 1.8,
        'split_off': 4.0,
        'conduction_band': 7.5
    })
    compute_all_bands: bool = False


@dataclass
class SemitipConfig:
    """SEMITIP 主配置類別"""
    version: str = "1.0"
    simulation_type: str = "MultInt"  # "MultInt" or "MultPlane"
    
    # 基本環境
    temperature: float = 300.0  # K
    dielectric_constant: float = 12.9
    
    # 各部分配置
    tip: TipConfig = field(default_factory=TipConfig)
    semiconductor_regions: List[SemiconductorRegion] = field(default_factory=list)
    surface_regions: List[SurfaceRegion] = field(default_factory=list)
    grid: GridConfig = field(default_factory=GridConfig)
    computation: ComputationConfig = field(default_factory=ComputationConfig)
    voltage_scan: VoltageScanConfig = field(default_factory=VoltageScanConfig)
    
    # 表面設定
    electron_affinity: float = 4.07  # eV
    surface_temperature_dependence: bool = False
    
    # 輸出設定
    output_basic: bool = True
    output_contours: bool = False
    output_full_potential: bool = False
    num_contours: int = 6
    contour_spacing: float = 0.0
    contour_angle: float = 0.0
      # 模擬類型特有參數
    multint_config: Optional[MultIntConfig] = None
    multplane_config: Optional[MultPlaneConfig] = None
    
    def __post_init__(self):
        """初始化後的處理"""        # 處理模擬類型特有參數
        if self.simulation_type == "MultInt":
            if not self.multint_config:
                self.multint_config = MultIntConfig()
            elif isinstance(self.multint_config, dict):
                self.multint_config = MultIntConfig(**self.multint_config)
        elif self.simulation_type == "MultPlane":
            if not self.multplane_config:
                self.multplane_config = MultPlaneConfig()
            elif isinstance(self.multplane_config, dict):
                self.multplane_config = MultPlaneConfig(**self.multplane_config)
            
        # 處理字典類型的 tip 參數
        if isinstance(self.tip, dict):
            self.tip = TipConfig(**self.tip)
            
        # 處理字典類型的 grid 參數
        if isinstance(self.grid, dict):
            self.grid = GridConfig(**self.grid)
            
        # 處理字典類型的 computation 參數
        if isinstance(self.computation, dict):
            self.computation = ComputationConfig(**self.computation)
            
        # 處理字典類型的 voltage_scan 參數
        if isinstance(self.voltage_scan, dict):
            self.voltage_scan = VoltageScanConfig(**self.voltage_scan)
            
        # 處理半導體區域
        if self.semiconductor_regions:
            processed_regions = []
            for region in self.semiconductor_regions:
                if isinstance(region, dict):
                    # 處理有效質量
                    if 'effective_mass' in region and isinstance(region['effective_mass'], dict):
                        region['effective_mass'] = EffectiveMass(**region['effective_mass'])
                    processed_regions.append(SemiconductorRegion(**region))
                else:
                    processed_regions.append(region)
            self.semiconductor_regions = processed_regions
            
        # 處理表面區域
        if self.surface_regions:
            processed_surface_regions = []
            for region in self.surface_regions:
                if isinstance(region, dict):
                    # 處理表面分佈
                    if 'first_distribution' in region and isinstance(region['first_distribution'], dict):
                        region['first_distribution'] = SurfaceDistribution(**region['first_distribution'])
                    if 'second_distribution' in region and isinstance(region['second_distribution'], dict):
                        region['second_distribution'] = SurfaceDistribution(**region['second_distribution'])
                    processed_surface_regions.append(SurfaceRegion(**region))
                else:
                    processed_surface_regions.append(region)
            self.surface_regions = processed_surface_regions

    def validate(self):
        """驗證配置參數的有效性"""
        errors = []
        
        # 基本參數驗證
        if self.temperature <= 0:
            errors.append("溫度必須大於 0K")
        
        if self.dielectric_constant <= 0:
            errors.append("介電常數必須大於 0")
        
        # 探針參數驗證
        if self.tip.radius <= 0:
            errors.append("探針半徑必須大於 0")
        
        if self.tip.separation <= 0:
            errors.append("探針分離距離必須大於 0")
        
        # 半導體區域驗證
        if not self.semiconductor_regions:
            errors.append("至少需要一個半導體區域")
        
        for region in self.semiconductor_regions:
            if region.band_gap <= 0:
                errors.append(f"半導體區域 {region.id} 的帶隙必須大於 0")
            
            if region.donor_concentration < 0 or region.acceptor_concentration < 0:
                errors.append(f"半導體區域 {region.id} 的摻雜濃度不能為負值")
        
        # 網格參數驗證
        if self.grid.radial_points <= 0:
            errors.append("徑向網格點數必須大於 0")
        
        if self.grid.vacuum_points <= 0:
            errors.append("真空網格點數必須大於 0")
        
        if self.grid.semiconductor_points <= 0:
            errors.append("半導體網格點數必須大於 0")
        
        if self.grid.initial_grid_size <= 0:
            errors.append("初始網格大小必須大於 0")
        
        # 電壓掃描驗證
        if self.voltage_scan.points <= 0:
            errors.append("電壓點數必須大於 0")
        
        if self.voltage_scan.start_voltage >= self.voltage_scan.end_voltage:
            errors.append("起始電壓必須小於結束電壓")
        
        # 模擬類型特有驗證
        if self.simulation_type == "MultInt":
            if not self.multint_config:
                errors.append("MultInt 模擬類型需要 multint_config")
            else:
                if self.multint_config.parallel_wavevectors <= 0:
                    errors.append("平行波向量數必須大於 0")
                if self.multint_config.energy_points <= 0:
                    errors.append("能量點數必須大於 0")
        
        elif self.simulation_type == "MultPlane":
            if not self.multplane_config:
                errors.append("MultPlane 模擬類型需要 multplane_config")
            else:
                if self.multplane_config.vacuum_width <= 0:
                    errors.append("真空寬度必須大於 0")
                if self.multplane_config.vacuum_spacing <= 0:
                    errors.append("真空間距必須大於 0")
        
        if errors:
            raise ValueError("配置驗證失敗:\n" + "\n".join(f"- {error}" for error in errors))
        
        return True
