"""
使用配置模型讀取 YAML 測試腳本
"""

import yaml
from pathlib import Path
import sys
import logging
import traceback

# 設定日誌
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

# 確保在當前目錄運行
BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

# 導入配置相關類別
try:
    from src.core.config_schema import (
        SemitipConfig, TipConfig, SemiconductorRegion, SurfaceRegion,
        EffectiveMass, SurfaceDistribution, GridConfig, ComputationConfig,
        VoltageScanConfig, MultIntConfig, MultPlaneConfig
    )
    logger.info("成功導入配置模型")
except ImportError as e:
    logger.error(f"導入配置模型失敗: {e}")
    sys.exit(1)

def load_yaml_to_config(file_path):
    """載入 YAML 配置文件並轉換為配置物件"""
    path = Path(file_path)
    logger.info(f"讀取 YAML 配置文件: {path}")
    
    if not path.exists():
        logger.error(f"文件不存在: {path}")
        return None
    
    try:
        # 讀取 YAML 文件
        with open(path, 'r', encoding='utf-8') as f:
            yaml_data = yaml.safe_load(f)
        
        logger.info("YAML 文件讀取成功")
        
        # 手動將 YAML 數據轉換為配置物件
        config = convert_yaml_to_config(yaml_data)
        
        logger.info(f"配置類型: {config.simulation_type}")
        logger.info(f"探針分離: {config.tip.separation} nm")
        logger.info(f"電介質常數: {config.dielectric_constant}")
        logger.info(f"溫度: {config.temperature} K")
        logger.info(f"半導體區域數量: {len(config.semiconductor_regions)}")
        
        # 嘗試驗證配置
        try:
            # 暫時跳過驗證
            # config.validate()
            logger.info("跳過配置驗證")
        except Exception as e:
            logger.error(f"配置驗證失敗: {e}")
            
        return config
        
    except Exception as e:
        logger.error(f"處理 YAML 配置時發生錯誤: {e}")
        logger.error(traceback.format_exc())
        return None

def convert_yaml_to_config(yaml_data):
    """手動轉換 YAML 數據為配置物件"""
    # 處理 tip 配置
    tip_data = yaml_data.get('tip', {})
    tip = TipConfig(**tip_data)
    
    # 處理半導體區域
    semi_regions = []
    for region_data in yaml_data.get('semiconductor_regions', []):
        # 處理有效質量
        mass_data = region_data.get('effective_mass', {})
        effective_mass = EffectiveMass(**mass_data)
        
        # 創建區域，並替換有效質量字典為物件
        region_copy = region_data.copy()
        region_copy['effective_mass'] = effective_mass
        semi_regions.append(SemiconductorRegion(**region_copy))
    
    # 處理表面區域
    surface_regions = []
    for surface_data in yaml_data.get('surface_regions', []):
        # 處理表面分布
        first_dist_data = surface_data.get('first_distribution', {})
        second_dist_data = surface_data.get('second_distribution', {})
        
        first_dist = SurfaceDistribution(**first_dist_data)
        second_dist = SurfaceDistribution(**second_dist_data)
        
        # 創建表面區域
        surface_copy = surface_data.copy()
        surface_copy['first_distribution'] = first_dist
        surface_copy['second_distribution'] = second_dist
        surface_regions.append(SurfaceRegion(**surface_copy))
    
    # 處理網格配置
    grid_data = yaml_data.get('grid', {})
    grid = GridConfig(**grid_data)
    
    # 處理計算配置
    computation_data = yaml_data.get('computation', {})
    computation = ComputationConfig(**computation_data)
    
    # 處理電壓掃描配置
    voltage_scan_data = yaml_data.get('voltage_scan', {})
    voltage_scan = VoltageScanConfig(**voltage_scan_data)
    
    # 處理特有配置
    multint_config = None
    multplane_config = None
    
    if yaml_data.get('simulation_type') == 'MultInt' and 'multint_config' in yaml_data:
        multint_data = yaml_data.get('multint_config', {})
        multint_config = MultIntConfig(**multint_data)
    
    if yaml_data.get('simulation_type') == 'MultPlane' and 'multplane_config' in yaml_data:
        multplane_data = yaml_data.get('multplane_config', {})
        
        # 處理 max_energies 字典
        max_energies = multplane_data.get('max_energies', {})
        multplane_data_copy = multplane_data.copy()
        multplane_data_copy['max_energies'] = max_energies
        
        multplane_config = MultPlaneConfig(**multplane_data_copy)
    
    # 創建主配置物件
    config_args = {
        'version': yaml_data.get('version', '1.0'),
        'simulation_type': yaml_data.get('simulation_type', ''),
        'temperature': yaml_data.get('temperature'),
        'dielectric_constant': yaml_data.get('dielectric_constant'),
        'tip': tip,
        'semiconductor_regions': semi_regions,
        'surface_regions': surface_regions,
        'grid': grid,
        'computation': computation,
        'voltage_scan': voltage_scan,
        'electron_affinity': yaml_data.get('electron_affinity'),
        'surface_temperature_dependence': yaml_data.get('surface_temperature_dependence'),
        'output_basic': yaml_data.get('output_basic', True),
        'output_contours': yaml_data.get('output_contours', False),
        'output_full_potential': yaml_data.get('output_full_potential', False),
        'num_contours': yaml_data.get('num_contours'),
        'contour_spacing': yaml_data.get('contour_spacing'),
        'contour_angle': yaml_data.get('contour_angle'),
        'multint_config': multint_config,
        'multplane_config': multplane_config
    }
    
    return SemitipConfig(**config_args)

def main():
    """主函數"""
    data_dir = BASE_DIR / "data" / "input"
    
    # 測試 MultInt 配置
    multint_path = data_dir / "MultInt_config.yaml"
    multint_config = load_yaml_to_config(multint_path)
    
    print("\n" + "-" * 40 + "\n")
    
    # 測試 MultPlane 配置
    multplane_path = data_dir / "MultPlane_config.yaml"
    multplane_config = load_yaml_to_config(multplane_path)

if __name__ == "__main__":
    main()
