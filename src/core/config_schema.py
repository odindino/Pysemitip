"""
SEMITIP 配置資料結構定義

此模組定義了所有 SEMITIP 配置相關的資料類別，包括：
- 環境配置 (EnvironmentConfig)
- 位置配置 (PositionConfig)
- 探針配置 (TipConfig)
- 半導體區域配置 (SemiconductorRegion)
- 表面態配置 (SurfaceRegion)
- 網格配置 (GridConfig)
- 計算配置 (ComputationConfig)
- 電壓掃描配置 (VoltageScanConfig)
- 模擬類型特有配置 (MultIntConfig, MultPlaneConfig)
- 階層配置 (SemiconductorConfig, SurfaceConfig, OutputConfig)
- 主配置類別 (SemitipConfig)

版本 2.0 更新：
- 支援階層化 YAML 結構
- 新增環境、半導體、表面、輸出配置類別
- 更新探針配置使用位置子結構
- 向後相容舊格式

作者: Odindino
版本: 2.0
日期: 2025-05-30
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union


@dataclass
class EnvironmentConfig:
    """環境配置參數"""
    temperature: float = 300.0  # K
    dielectric_constant: float = 12.9


@dataclass
class PositionConfig:
    """位置配置參數"""
    x: float = 0.0  # nm
    y: float = 0.0  # nm


@dataclass
class TipConfig:
    """探針配置參數"""
    shank_slope: float = 1.0
    separation: float = 1.0  # nm
    radius: float = 1.0  # nm
    protrusion_radius: float = 0.0  # nm
    contact_potential: float = 0.0  # eV
    position: PositionConfig = field(default_factory=PositionConfig)
    fermi_energy: float = 8.0  # eV
    work_function: float = 5.3  # eV - tip work function
    
    # 向後相容性屬性
    @property
    def x_position(self) -> float:
        """向後相容：取得 x 位置"""
        return self.position.x
    
    @x_position.setter
    def x_position(self, value: float):
        """向後相容：設定 x 位置"""
        self.position.x = value
    
    @property
    def y_position(self) -> float:
        """向後相容：取得 y 位置"""
        return self.position.y
    
    @y_position.setter
    def y_position(self, value: float):
        """向後相容：設定 y 位置"""
        self.position.y = value
    
    @property
    def slope(self) -> float:
        """向後相容：取得探針斜率"""
        return self.shank_slope
    
    @slope.setter 
    def slope(self, value: float):
        """向後相容：設定探針斜率"""
        self.shank_slope = value


@dataclass
class EffectiveMass:
    """有效質量參數"""
    conduction_band: float = 0.0635
    valence_band_heavy: float = 0.643
    valence_band_light: float = 0.081
    split_off: float = 0.172
    
    # Legacy aliases for backward compatibility
    @property
    def heavy_hole(self) -> float:
        return self.valence_band_heavy
    
    @property
    def light_hole(self) -> float:
        return self.valence_band_light
    
    @property
    def split_off_hole(self) -> float:
        return self.split_off


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
    
    # Additional fields for YAML compatibility
    affinity: float = 4.07  # eV - electron affinity
    permittivity: float = 12.9  # relative permittivity
    allow_degeneracy: bool = True
    allow_inversion: bool = True


@dataclass
class SurfaceDistribution:
    """表面態分佈參數"""
    density: float = 0.0  # cm^-2 eV^-1
    neutrality_level: float = 0.0  # eV
    fwhm: float = 0.0  # eV
    center_energy: float = 0.0  # eV
    
    # Legacy alias
    @property
    def centroid_energy(self) -> float:
        return self.center_energy


@dataclass
class SurfaceRegion:
    """表面區域參數"""
    id: int = 1
    first_distribution: SurfaceDistribution = field(default_factory=SurfaceDistribution)
    second_distribution: Optional[SurfaceDistribution] = None
    position: PositionConfig = field(default_factory=PositionConfig)


@dataclass
class SemiconductorConfig:
    """半導體配置 - 階層結構"""
    regions: List[SemiconductorRegion] = field(default_factory=list)
    electron_affinity: float = 4.07  # eV


@dataclass
class SurfaceConfig:
    """表面配置 - 階層結構"""
    regions: List[SurfaceRegion] = field(default_factory=list)
    temperature_dependence: bool = False


@dataclass
class OutputConfig:
    """輸出配置 - 階層結構"""
    basic_output: bool = True
    equipotential_contours: bool = False
    full_potential: bool = False
    num_contours: int = 6
    contour_spacing: float = 0.0
    contour_angle: float = 0.0


@dataclass
class GridConfig:
    """計算網格配置"""
    mirror_plane: bool = True
    radial_points: int = 32
    vacuum_points: int = 16
    semiconductor_points: int = 32
    angular_points: int = 16
    initial_grid_size: float = 0.5
    
    # Additional extent parameters
    radial_extent: float = 20.0
    vacuum_extent: float = 20.0
    semiconductor_extent: float = 20.0
    energy_points: int = 20


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
    start: float = -2.0  # V (alias for start_voltage)
    end: float = 2.0  # V (alias for end_voltage)
    modulation_voltage: float = 0.050  # V
    negative_ramp: float = 0.0  # nm/V
    positive_ramp: float = 0.0  # nm/V
    custom_voltages: Optional[List[float]] = None
    
    # Legacy aliases
    @property
    def start_voltage(self) -> float:
        return self.start
    
    @property 
    def end_voltage(self) -> float:
        return self.end
    
    @property
    def voltages(self) -> List[float]:
        """生成電壓掃描陣列"""
        if self.custom_voltages:
            return self.custom_voltages
        else:
            # Generate linear space without numpy
            if self.points <= 1:
                return [self.start]
            step = (self.end - self.start) / (self.points - 1)
            return [self.start + i * step for i in range(self.points)]


@dataclass
class MultIntConfig:
    """MultInt 特有參數"""
    parallel_wavevectors: int = 20
    energy_points: int = 20
    expansion_factor: int = 20
    semiconductor_depth_fraction: float = 0.75
    integration_cutoff: float = 0.1


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
    """SEMITIP 主配置類別 - 支援階層化 YAML 結構"""
    version: str = "1.0"
    simulation_type: str = "MultInt"  # "MultInt" or "MultPlane"
    
    # 階層化配置結構 - 對應新 YAML 格式
    environment: EnvironmentConfig = field(default_factory=EnvironmentConfig)
    tip: TipConfig = field(default_factory=TipConfig)
    semiconductor: SemiconductorConfig = field(default_factory=SemiconductorConfig)
    surface: SurfaceConfig = field(default_factory=SurfaceConfig)
    grid: GridConfig = field(default_factory=GridConfig)
    computation: ComputationConfig = field(default_factory=ComputationConfig)
    voltage_scan: VoltageScanConfig = field(default_factory=VoltageScanConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    
    # 模擬類型特有參數 - 對應新 YAML 命名
    multint_specific: Optional[MultIntConfig] = None
    multplane_specific: Optional[MultPlaneConfig] = None
    
    # 向後相容性屬性
    @property
    def temperature(self) -> float:
        """向後相容：取得溫度"""
        return self.environment.temperature
    
    @temperature.setter
    def temperature(self, value: float):
        """向後相容：設定溫度"""
        self.environment.temperature = value
    
    @property
    def dielectric_constant(self) -> float:
        """向後相容：取得介電常數"""
        return self.environment.dielectric_constant
    
    @dielectric_constant.setter
    def dielectric_constant(self, value: float):
        """向後相容：設定介電常數"""
        self.environment.dielectric_constant = value
    
    @property
    def semiconductor_regions(self) -> List[SemiconductorRegion]:
        """向後相容：取得半導體區域"""
        return self.semiconductor.regions
    
    @semiconductor_regions.setter
    def semiconductor_regions(self, value: List[SemiconductorRegion]):
        """向後相容：設定半導體區域"""
        self.semiconductor.regions = value
    
    @property
    def surface_regions(self) -> List[SurfaceRegion]:
        """向後相容：取得表面區域"""
        return self.surface.regions
    
    @surface_regions.setter
    def surface_regions(self, value: List[SurfaceRegion]):
        """向後相容：設定表面區域"""
        self.surface.regions = value
    
    @property
    def contact_potential(self) -> float:
        """向後相容：取得接觸電位"""
        return self.tip.contact_potential
    
    @contact_potential.setter
    def contact_potential(self, value: float):
        """向後相容：設定接觸電位"""
        self.tip.contact_potential = value
    
    @property
    def mirror_symmetry(self) -> bool:
        """向後相容：取得鏡像對稱性"""
        return self.grid.mirror_plane
    
    @mirror_symmetry.setter
    def mirror_symmetry(self, value: bool):
        """向後相容：設定鏡像對稱性"""
        self.grid.mirror_plane = value
    
    @property
    def charge_table_points(self) -> int:
        """向後相容：取得電荷表點數"""
        return self.computation.charge_density_table_size
    
    @charge_table_points.setter
    def charge_table_points(self, value: int):
        """向後相容：設定電荷表點數"""
        self.computation.charge_density_table_size = value
    
    @property
    def max_iterations(self) -> int:
        """向後相容：取得最大迭代次數 (第一個值)"""
        return self.computation.max_iterations[0] if self.computation.max_iterations else 50
    
    @max_iterations.setter
    def max_iterations(self, value: int):
        """向後相容：設定最大迭代次數 (更新第一個值)"""
        if not self.computation.max_iterations:
            self.computation.max_iterations = [value]
        else:
            self.computation.max_iterations[0] = value
    
    @property
    def convergence_tolerance(self) -> float:
        """向後相容：取得收斂容差 (第一個值)"""
        return self.computation.convergence_parameters[0] if self.computation.convergence_parameters else 1e-3
    
    @convergence_tolerance.setter
    def convergence_tolerance(self, value: float):
        """向後相容：設定收斂容差 (更新第一個值)"""
        if not self.computation.convergence_parameters:
            self.computation.convergence_parameters = [value]
        else:
            self.computation.convergence_parameters[0] = value
    
    @property
    def electron_affinity(self) -> float:
        """向後相容：取得電子親和力"""
        return self.semiconductor.electron_affinity
    
    @electron_affinity.setter
    def electron_affinity(self, value: float):
        """向後相容：設定電子親和力"""
        self.semiconductor.electron_affinity = value
    
    @property
    def surface_temperature_dependence(self) -> bool:
        """向後相容：取得表面溫度相依性"""
        return self.surface.temperature_dependence
    
    @surface_temperature_dependence.setter
    def surface_temperature_dependence(self, value: bool):
        """向後相容：設定表面溫度相依性"""
        self.surface.temperature_dependence = value
    
    @property
    def output_basic(self) -> bool:
        """向後相容：取得基本輸出設定"""
        return self.output.basic_output
    
    @output_basic.setter
    def output_basic(self, value: bool):
        """向後相容：設定基本輸出"""
        self.output.basic_output = value
    
    @property
    def output_contours(self) -> bool:
        """向後相容：取得等高線輸出設定"""
        return self.output.equipotential_contours
    
    @output_contours.setter
    def output_contours(self, value: bool):
        """向後相容：設定等高線輸出"""
        self.output.equipotential_contours = value
    
    @property
    def output_full_potential(self) -> bool:
        """向後相容：取得完整電位輸出設定"""
        return self.output.full_potential
    
    @output_full_potential.setter
    def output_full_potential(self, value: bool):
        """向後相容：設定完整電位輸出"""
        self.output.full_potential = value
    
    @property
    def num_contours(self) -> int:
        """向後相容：取得等高線數量"""
        return self.output.num_contours
    
    @num_contours.setter
    def num_contours(self, value: int):
        """向後相容：設定等高線數量"""
        self.output.num_contours = value
    
    @property
    def contour_spacing(self) -> float:
        """向後相容：取得等高線間距"""
        return self.output.contour_spacing
    
    @contour_spacing.setter
    def contour_spacing(self, value: float):
        """向後相容：設定等高線間距"""
        self.output.contour_spacing = value
    
    @property
    def contour_angle(self) -> float:
        """向後相容：取得等高線角度"""
        return self.output.contour_angle
    
    @contour_angle.setter
    def contour_angle(self, value: float):
        """向後相容：設定等高線角度"""
        self.output.contour_angle = value
    
    # 提供舊的參數命名的相容性
    @property
    def multint_config(self) -> Optional[MultIntConfig]:
        """向後相容：取得 MultInt 配置"""
        return self.multint_specific
    
    @multint_config.setter
    def multint_config(self, value: Optional[MultIntConfig]):
        """向後相容：設定 MultInt 配置"""
        self.multint_specific = value
    
    @property
    def multplane_config(self) -> Optional[MultPlaneConfig]:
        """向後相容：取得 MultPlane 配置"""
        return self.multplane_specific
    
    @multplane_config.setter
    def multplane_config(self, value: Optional[MultPlaneConfig]):
        """向後相容：設定 MultPlane 配置"""
        self.multplane_specific = value
    
    def __post_init__(self):
        """初始化後的處理"""
        # 處理模擬類型特有參數
        if self.simulation_type == "MultInt":
            if not self.multint_specific:
                self.multint_specific = MultIntConfig()
        elif self.simulation_type == "MultPlane":
            if not self.multplane_specific:
                self.multplane_specific = MultPlaneConfig()
                  # 處理字典類型的環境參數
        if isinstance(self.environment, dict):
            self.environment = EnvironmentConfig(**self.environment)
            
        # 處理字典類型的 tip 參數
        if isinstance(self.tip, dict):
            # 處理舊格式的 x_position, y_position
            tip_dict = self.tip.copy()
            if 'x_position' in tip_dict or 'y_position' in tip_dict:
                position_dict = {
                    'x': tip_dict.pop('x_position', 0.0),
                    'y': tip_dict.pop('y_position', 0.0)
                }
                tip_dict['position'] = PositionConfig(**position_dict)
            # 處理新格式的巢狀 position 字典
            elif 'position' in tip_dict and isinstance(tip_dict['position'], dict):
                tip_dict['position'] = PositionConfig(**tip_dict['position'])
            self.tip = TipConfig(**tip_dict)
            
        # 處理字典類型的半導體參數
        if isinstance(self.semiconductor, dict):
            self.semiconductor = SemiconductorConfig(**self.semiconductor)
            
        # 處理字典類型的表面參數
        if isinstance(self.surface, dict):
            self.surface = SurfaceConfig(**self.surface)
            
        # 處理字典類型的 grid 參數
        if isinstance(self.grid, dict):
            self.grid = GridConfig(**self.grid)
            
        # 處理字典類型的 computation 參數
        if isinstance(self.computation, dict):
            self.computation = ComputationConfig(**self.computation)
            
        # 處理字典類型的 voltage_scan 參數
        if isinstance(self.voltage_scan, dict):
            self.voltage_scan = VoltageScanConfig(**self.voltage_scan)
            
        # 處理字典類型的輸出參數
        if isinstance(self.output, dict):
            self.output = OutputConfig(**self.output)
            
        # 處理字典類型的模擬特有參數
        if isinstance(self.multint_specific, dict):
            self.multint_specific = MultIntConfig(**self.multint_specific)
            
        if isinstance(self.multplane_specific, dict):
            self.multplane_specific = MultPlaneConfig(**self.multplane_specific)
        
        # 處理半導體區域
        if self.semiconductor.regions:
            for i, region in enumerate(self.semiconductor.regions):
                if isinstance(region, dict):
                    # 處理字典類型的有效質量
                    if 'effective_mass' in region and isinstance(region['effective_mass'], dict):
                        region['effective_mass'] = EffectiveMass(**region['effective_mass'])
                    self.semiconductor.regions[i] = SemiconductorRegion(**region)
        
        # 處理表面區域
        if self.surface.regions:
            for i, region in enumerate(self.surface.regions):
                if isinstance(region, dict):
                    # 處理字典類型的表面分佈
                    if 'first_distribution' in region and isinstance(region['first_distribution'], dict):
                        region['first_distribution'] = SurfaceDistribution(**region['first_distribution'])
                    if 'second_distribution' in region and isinstance(region['second_distribution'], dict):
                        region['second_distribution'] = SurfaceDistribution(**region['second_distribution'])
                    self.surface.regions[i] = SurfaceRegion(**region)

