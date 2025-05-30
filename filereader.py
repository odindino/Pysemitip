"""
SEMITIP YAML 配置檔案讀取器

此模組專門處理 YAML 格式的 SEMITIP 配置檔案，提供：
1. YAML 配置檔案的讀取
2. YAML 配置檔案的寫入
3. 配置驗證和錯誤處理
4. 配置檔案的結構化操作

作者: Odindino
版本: 1.0
日期: 2025-05-30
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional, Union
import logging

from config_schema import SemitipConfig

# 設定日誌
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class YamlConfigReader:
    """YAML 配置檔案讀取器"""
    
    def __init__(self):
        """初始化讀取器"""
        self.config: Optional[SemitipConfig] = None
    
    def load_config(self, file_path: Union[str, Path]) -> SemitipConfig:
        """
        從 YAML 檔案載入配置

        Args:
            file_path: YAML 配置檔案路徑
            
        Returns:
            SemitipConfig: 解析後的配置物件
            
        Raises:
            FileNotFoundError: 檔案不存在
            yaml.YAMLError: YAML 格式錯誤
            ValueError: 配置驗證失敗
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"配置檔案不存在: {file_path}")
        
        logger.info(f"正在載入配置檔案: {file_path}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                yaml_data = yaml.safe_load(f)
              # 將 YAML 資料轉換為配置物件
            config = self._yaml_to_config(yaml_data)
            
            # 執行配置驗證
            config.validate()
            logger.info("配置驗證通過")
            
            self.config = config
            logger.info("配置檔案載入成功")
            return config
            
        except yaml.YAMLError as e:
            logger.error(f"YAML 格式錯誤: {e}")
            raise
        except ValueError as e:
            logger.error(f"配置驗證失敗: {e}")
            raise
        except Exception as e:
            logger.error(f"載入配置檔案時發生錯誤: {e}")
            raise
    def save_config(self, config: SemitipConfig, file_path: Union[str, Path]) -> None:
        """
        將配置儲存為 YAML 檔案
        
        Args:
            config: 要儲存的配置物件
            file_path: 輸出檔案路徑
            
        Raises:
            PermissionError: 沒有寫入權限
            IOError: 檔案寫入錯誤
        """
        file_path = Path(file_path)
        
        # 檢查配置物件是否有效
        if not isinstance(config, SemitipConfig):
            raise ValueError("提供的配置物件無效，必須是 SemitipConfig 類型")
        logger.info(f"正在儲存配置檔案: {file_path}")
        
        
        # 轉換為 YAML 格式
        yaml_data = self._config_to_yaml(config)
        
        try:
            file_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(file_path, 'w', encoding='utf-8') as f:
                yaml.dump(yaml_data, f, default_flow_style=False, 
                         allow_unicode=True, indent=2, sort_keys=False)
            
            logger.info(f"配置檔案已儲存至: {file_path}")
            
        except PermissionError as e:
            logger.error(f"沒有寫入權限: {file_path}")
            raise
        except IOError as e:
            logger.error(f"檔案寫入錯誤: {e}")
            raise
    
    def _yaml_to_config(self, yaml_data: Dict[str, Any]) -> SemitipConfig:
        """
        將 YAML 資料轉換為 SemitipConfig 物件
        
        Args:
            yaml_data: 解析後的 YAML 資料
            
        Returns:
            SemitipConfig: 配置物件
        """
        # 直接使用 dataclass 的 __init__ 方法
        # 這會自動處理巢狀結構
        try:
            return SemitipConfig(**yaml_data)
        except TypeError as e:
            # 如果直接轉換失敗，可能需要更複雜的轉換邏輯
            logger.warning(f"直接轉換失敗，使用手動轉換: {e}")
            return self._manual_yaml_to_config(yaml_data)
    def _manual_yaml_to_config(self, yaml_data: Dict[str, Any]) -> SemitipConfig:
        """
        手動轉換 YAML 資料為配置物件（當自動轉換失敗時使用）
        
        Args:
            yaml_data: 解析後的 YAML 資料
            
        Returns:
            SemitipConfig: 配置物件
        """
        from config_schema import (
            TipConfig, EffectiveMass, SemiconductorRegion, 
            SurfaceDistribution, SurfaceRegion, GridConfig,
            ComputationConfig, VoltageScanConfig, MultIntConfig, MultPlaneConfig
        )
        
        # 處理環境參數
        environment_data = yaml_data.get('environment', {})
        temperature = environment_data.get('temperature', yaml_data.get('temperature', 300.0))
        dielectric_constant = environment_data.get('dielectric_constant', yaml_data.get('dielectric_constant', 12.9))
        
        # 處理 tip 配置 (支援新的 position 結構)
        tip_data = yaml_data.get('tip', {})
        if tip_data:
            tip_copy = tip_data.copy()
            # 處理 position 結構轉換
            if 'position' in tip_copy:
                position = tip_copy.pop('position')
                tip_copy['x_position'] = position.get('x', 0.0)
                tip_copy['y_position'] = position.get('y', 0.0)
            tip_config = TipConfig(**tip_copy)
        else:
            tip_config = TipConfig()        
        # 處理半導體區域 (支援新的階層結構)
        semiconductor_regions = []
        semiconductor_data = yaml_data.get('semiconductor', {})
        regions_data = semiconductor_data.get('regions', yaml_data.get('semiconductor_regions', []))
        
        for region_data in regions_data:
            # 處理區域內的有效質量
            region_copy = region_data.copy()
            if 'effective_mass' in region_copy:
                effective_mass_data = region_copy['effective_mass']
                region_copy['effective_mass'] = EffectiveMass(**effective_mass_data)
            semiconductor_regions.append(SemiconductorRegion(**region_copy))
        
        # 處理電子親和力 (新結構在 semiconductor 中，舊結構在根層級)
        electron_affinity = semiconductor_data.get('electron_affinity', yaml_data.get('electron_affinity', 4.07))
        
        # 處理表面區域 (支援新的階層結構)
        surface_regions = []
        surface_data = yaml_data.get('surface', {})
        surface_regions_data = surface_data.get('regions', yaml_data.get('surface_regions', []))
        
        for surface_region_data in surface_regions_data:
            # 處理 first_distribution 和 second_distribution
            surface_copy = surface_region_data.copy()
            
            if 'first_distribution' in surface_copy:
                first_dist_data = surface_copy['first_distribution']
                surface_copy['first_distribution'] = SurfaceDistribution(**first_dist_data)
            
            if 'second_distribution' in surface_copy:
                second_dist_data = surface_copy['second_distribution']  
                surface_copy['second_distribution'] = SurfaceDistribution(**second_dist_data)
                
            surface_regions.append(SurfaceRegion(**surface_copy))
        
        # 處理表面溫度相依性 (新結構在 surface 中，舊結構在根層級)
        surface_temperature_dependence = surface_data.get('temperature_dependence', 
                                                         yaml_data.get('surface_temperature_dependence', False))
        
        # 處理網格配置
        grid_data = yaml_data.get('grid', {})
        grid_config = GridConfig(**grid_data) if grid_data else GridConfig()
        
        # 處理計算配置
        computation_data = yaml_data.get('computation', {})
        computation_config = ComputationConfig(**computation_data) if computation_data else ComputationConfig()
          # 處理電壓掃描配置
        voltage_scan_data = yaml_data.get('voltage_scan', {})
        voltage_scan_config = VoltageScanConfig(**voltage_scan_data) if voltage_scan_data else VoltageScanConfig()
        
        # 處理特定配置 (支援新的命名)
        multint_data = yaml_data.get('multint_specific', yaml_data.get('multint_config', {}))
        multint_config = MultIntConfig(**multint_data) if multint_data else None
        
        multplane_data = yaml_data.get('multplane_specific', yaml_data.get('multplane_config', {}))
        multplane_config = MultPlaneConfig(**multplane_data) if multplane_data else None
        
        # 處理輸出設定 (支援新的階層結構)
        output_data = yaml_data.get('output', {})
        output_basic = output_data.get('basic_output', yaml_data.get('output_basic', True))
        output_contours = output_data.get('equipotential_contours', yaml_data.get('output_contours', False))
        output_full_potential = output_data.get('full_potential', yaml_data.get('output_full_potential', False))
        num_contours = output_data.get('num_contours', yaml_data.get('num_contours', 6))
        contour_spacing = output_data.get('contour_spacing', yaml_data.get('contour_spacing', 0.0))
        contour_angle = output_data.get('contour_angle', yaml_data.get('contour_angle', 0.0))
        
        # 創建主配置物件
        config_args = {
            'version': yaml_data.get('version', '1.0'),
            'simulation_type': yaml_data.get('simulation_type', 'MultInt'),
            'temperature': temperature,
            'dielectric_constant': dielectric_constant,
            'tip': tip_config,
            'semiconductor_regions': semiconductor_regions,
            'surface_regions': surface_regions,
            'grid': grid_config,
            'computation': computation_config,
            'voltage_scan': voltage_scan_config,
            'electron_affinity': electron_affinity,
            'surface_temperature_dependence': surface_temperature_dependence,
            'output_basic': output_basic,
            'output_contours': output_contours,
            'output_full_potential': output_full_potential,
            'num_contours': num_contours,
            'contour_spacing': contour_spacing,
            'contour_angle': contour_angle,
            'multint_config': multint_config,
            'multplane_config': multplane_config        }
        
        return SemitipConfig(**config_args)
    
    def _config_to_yaml(self, config: SemitipConfig) -> Dict[str, Any]:
        """
        將配置物件轉換為 YAML 格式的字典
        
        Args:
            config: 配置物件
            
        Returns:
            Dict: 可序列化為 YAML 的字典
        """
        result = {}
        
        # 基本資訊
        result['version'] = config.version
        result['simulation_type'] = config.simulation_type
        
        # 階層化結構
        result['environment'] = self._dataclass_to_dict(config.environment)
        result['tip'] = self._dataclass_to_dict(config.tip)
        result['semiconductor'] = self._dataclass_to_dict(config.semiconductor)
        result['surface'] = self._dataclass_to_dict(config.surface)
        result['grid'] = self._dataclass_to_dict(config.grid)
        result['computation'] = self._dataclass_to_dict(config.computation)
        result['voltage_scan'] = self._dataclass_to_dict(config.voltage_scan)
        result['output'] = self._dataclass_to_dict(config.output)
        
        # 模擬類型特有參數
        if config.multint_specific:
            result['multint_specific'] = self._dataclass_to_dict(config.multint_specific)
        
        if config.multplane_specific:
            result['multplane_specific'] = self._dataclass_to_dict(config.multplane_specific)
        
        return result
    
    def _dataclass_to_dict(self, obj) -> Dict[str, Any]:
        """
        將 dataclass 物件轉換為字典
        
        Args:
            obj: dataclass 物件
            
        Returns:
            Dict: 轉換後的字典
        """
        if obj is None:
            return {}
        
        from dataclasses import asdict
        return asdict(obj)
    
    def get_simulation_type(self) -> Optional[str]:
        """
        取得當前配置的模擬類型
        
        Returns:
            str: 模擬類型 ("MultInt" 或 "MultPlane")
        """
        if self.config:
            return self.config.simulation_type
        return None
    
    def validate_config(self) -> bool:
        """
        驗證當前配置的有效性
        
        Returns:
            bool: 驗證是否通過
            
        Raises:
            ValueError: 配置驗證失敗
        """
        if not self.config:
            raise ValueError("尚未載入配置檔案")
        
        return self.config.validate()
        if config.semiconductor_regions:
            semiconductor_regions = []
            for region in config.semiconductor_regions:
                region_dict = self._dataclass_to_dict(region)
                if region.effective_mass:
                    region_dict['effective_mass'] = self._dataclass_to_dict(region.effective_mass)
                semiconductor_regions.append(region_dict)
            
            result['semiconductor'] = {
                'regions': semiconductor_regions,
                'electron_affinity': config.electron_affinity
            }
        
        # 處理表面區域 (新階層結構)
        if config.surface_regions:
            surface_regions = []
            for surface in config.surface_regions:
                surface_dict = self._dataclass_to_dict(surface)
                if surface.first_distribution:
                    surface_dict['first_distribution'] = self._dataclass_to_dict(surface.first_distribution)
                if surface.second_distribution:
                    surface_dict['second_distribution'] = self._dataclass_to_dict(surface.second_distribution)
                surface_regions.append(surface_dict)
            
            result['surface'] = {
                'regions': surface_regions,
                'temperature_dependence': config.surface_temperature_dependence
            }
        
        # 處理網格配置
        if config.grid:
            result['grid'] = self._dataclass_to_dict(config.grid)
        
        # 處理計算配置
        if config.computation:
            result['computation'] = self._dataclass_to_dict(config.computation)
        
        # 處理電壓掃描配置
        if config.voltage_scan:
            result['voltage_scan'] = self._dataclass_to_dict(config.voltage_scan)
        
        # 處理特定配置 (使用新的命名)
        if config.multint_config:
            result['multint_specific'] = self._dataclass_to_dict(config.multint_config)
        
        if config.multplane_config:
            result['multplane_specific'] = self._dataclass_to_dict(config.multplane_config)
        
        # 處理輸出設定 (新階層結構)
        result['output'] = {
            'basic_output': config.output_basic,
            'equipotential_contours': config.output_contours,
            'full_potential': config.output_full_potential,
            'num_contours': config.num_contours,
            'contour_spacing': config.contour_spacing,
            'contour_angle': config.contour_angle
        }
        
        return result
    
    def _dataclass_to_dict(self, obj) -> Dict[str, Any]:
        """
        將 dataclass 物件轉換為字典
        
        Args:
            obj: dataclass 物件
            
        Returns:
            Dict: 字典表示
        """
        import dataclasses
        
        if dataclasses.is_dataclass(obj):
            result = {}
            for field in dataclasses.fields(obj):
                value = getattr(obj, field.name)
                if value is not None:
                    if dataclasses.is_dataclass(value):
                        result[field.name] = self._dataclass_to_dict(value)
                    elif isinstance(value, list):
                        result[field.name] = [
                            self._dataclass_to_dict(item) if dataclasses.is_dataclass(item) else item
                            for item in value
                        ]
                    else:
                        result[field.name] = value
            return result
        else:
            return obj
    
    def get_simulation_type(self) -> Optional[str]:
        """
        取得模擬類型
        
        Returns:
            str: 'multint' 或 'multplane'，如果都沒有則返回 None
        """
        if self.config is None:
            return None
        
        if self.config.multint is not None:
            return 'multint'
        elif self.config.multplane is not None:
            return 'multplane'
        else:
            return None
    
    def validate_config(self) -> bool:
        """
        驗證當前載入的配置
        
        Returns:
            bool: 驗證是否通過
        """
        if self.config is None:
            logger.error("尚未載入任何配置")
            return False
        
        try:
            self.config.validate()
            logger.info("配置驗證通過")
            return True
        except ValueError as e:
            logger.error(f"配置驗證失敗: {e}")
            return False


# 便利函數
def load_yaml_config(file_path: Union[str, Path]) -> SemitipConfig:
    """
    載入 YAML 配置檔案的便利函數
    
    Args:
        file_path: 配置檔案路徑
        
    Returns:
        SemitipConfig: 配置物件
    """
    reader = YamlConfigReader()
    return reader.load_config(file_path)


def save_yaml_config(config: SemitipConfig, file_path: Union[str, Path]) -> None:
    """
    儲存配置為 YAML 檔案的便利函數
    
    Args:
        config: 配置物件
        file_path: 輸出檔案路徑
    """
    reader = YamlConfigReader()
    reader.save_config(config, file_path)


# 範例使用
if __name__ == "__main__":
    # 載入配置檔案範例
    try:
        config = load_yaml_config("config_template.yaml")
        print(f"載入的模擬類型: {config.simulation_type}")
        
        # 儲存配置檔案範例
        save_yaml_config(config, "output_config.yaml")
        print("配置檔案已儲存")
        
    except Exception as e:
        print(f"錯誤: {e}")