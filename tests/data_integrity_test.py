#!/usr/bin/env python3
"""
數據完整性和邊界情況測試

此測試專門檢查：
1. 數據解析的準確性
2. 邊界值處理
3. 錯誤處理機制
4. 數據類型驗證
5. 缺失欄位處理

作者: Pysemitip 開發團隊
日期: 2024-01-02
"""

import yaml
import sys
import logging
from pathlib import Path
from typing import Dict, Any

# 設定日誌
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 確保能導入專案模組
BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

try:
    from src.core.filereader import YamlConfigReader
    from src.core.config_schema import SemitipConfig
except ImportError as e:
    logger.error(f"模組導入失敗: {e}")
    sys.exit(1)


class DataIntegrityTester:
    """數據完整性測試器"""
    
    def __init__(self):
        self.reader = YamlConfigReader()
    
    def test_data_accuracy(self, file_path: Path):
        """測試數據解析準確性"""
        logger.info(f"\n{'='*60}")
        logger.info(f"測試檔案數據準確性: {file_path.name}")
        logger.info('='*60)
        
        # 讀取原始 YAML
        with open(file_path, 'r', encoding='utf-8') as f:
            original_yaml = yaml.safe_load(f)
        
        # 轉換為配置物件
        config = self.reader.load_config(file_path)
        
        # 詳細比較關鍵數據
        self._compare_environment_data(original_yaml, config)
        self._compare_tip_data(original_yaml, config)
        self._compare_semiconductor_data(original_yaml, config)
        self._compare_surface_data(original_yaml, config)
        self._compare_output_data(original_yaml, config)
        
    def _compare_environment_data(self, yaml_data: Dict, config: SemitipConfig):
        """比較環境數據"""
        logger.info("🔍 檢查環境數據...")
        
        env = yaml_data.get('environment', {})
        yaml_temp = env.get('temperature')
        yaml_dielectric = env.get('dielectric_constant')
        
        logger.info(f"  原始溫度: {yaml_temp} K")
        logger.info(f"  解析溫度: {config.temperature} K")
        logger.info(f"  溫度一致: {'✅' if yaml_temp == config.temperature else '❌'}")
        
        logger.info(f"  原始介電常數: {yaml_dielectric}")
        logger.info(f"  解析介電常數: {config.dielectric_constant}")
        logger.info(f"  介電常數一致: {'✅' if yaml_dielectric == config.dielectric_constant else '❌'}")
    
    def _compare_tip_data(self, yaml_data: Dict, config: SemitipConfig):
        """比較探針數據"""
        logger.info("🔍 檢查探針數據...")
        
        tip = yaml_data.get('tip', {})
        position = tip.get('position', {})
        
        logger.info(f"  原始分離距離: {tip.get('separation')} nm")
        logger.info(f"  解析分離距離: {config.tip.separation} nm")
        logger.info(f"  分離距離一致: {'✅' if tip.get('separation') == config.tip.separation else '❌'}")
        
        logger.info(f"  原始探針半徑: {tip.get('radius')} nm")
        logger.info(f"  解析探針半徑: {config.tip.radius} nm")
        logger.info(f"  探針半徑一致: {'✅' if tip.get('radius') == config.tip.radius else '❌'}")
        
        logger.info(f"  原始 X 位置: {position.get('x', 0.0)} nm")
        logger.info(f"  解析 X 位置: {config.tip.x_position} nm")
        logger.info(f"  X 位置一致: {'✅' if position.get('x', 0.0) == config.tip.x_position else '❌'}")
        
        logger.info(f"  原始 Y 位置: {position.get('y', 0.0)} nm")
        logger.info(f"  解析 Y 位置: {config.tip.y_position} nm")
        logger.info(f"  Y 位置一致: {'✅' if position.get('y', 0.0) == config.tip.y_position else '❌'}")
    
    def _compare_semiconductor_data(self, yaml_data: Dict, config: SemitipConfig):
        """比較半導體數據"""
        logger.info("🔍 檢查半導體數據...")
        
        semiconductor = yaml_data.get('semiconductor', {})
        regions = semiconductor.get('regions', [])
        
        logger.info(f"  原始區域數量: {len(regions)}")
        logger.info(f"  解析區域數量: {len(config.semiconductor_regions)}")
        logger.info(f"  區域數量一致: {'✅' if len(regions) == len(config.semiconductor_regions) else '❌'}")
        
        for i, (yaml_region, config_region) in enumerate(zip(regions, config.semiconductor_regions)):
            logger.info(f"  區域 {i+1}:")
            
            # 檢查摻雜濃度
            yaml_donor = yaml_region.get('donor_concentration')
            config_donor = config_region.donor_concentration
            logger.info(f"    供體濃度: {yaml_donor} -> {config_donor} {'✅' if yaml_donor == config_donor else '❌'}")
            
            # 檢查能隙
            yaml_bandgap = yaml_region.get('band_gap')
            config_bandgap = config_region.band_gap
            logger.info(f"    能隙: {yaml_bandgap} -> {config_bandgap} {'✅' if yaml_bandgap == config_bandgap else '❌'}")
            
            # 檢查有效質量
            if 'effective_mass' in yaml_region and config_region.effective_mass:
                yaml_mass = yaml_region['effective_mass']
                config_mass = config_region.effective_mass
                
                for mass_type in ['conduction_band', 'heavy_hole', 'light_hole', 'split_off_hole']:
                    yaml_val = yaml_mass.get(mass_type)
                    config_val = getattr(config_mass, mass_type)
                    status = '✅' if yaml_val == config_val else '❌'
                    logger.info(f"    {mass_type}: {yaml_val} -> {config_val} {status}")
        
        # 檢查電子親和力
        yaml_affinity = semiconductor.get('electron_affinity')
        config_affinity = config.electron_affinity
        logger.info(f"  電子親和力: {yaml_affinity} -> {config_affinity} {'✅' if yaml_affinity == config_affinity else '❌'}")
    
    def _compare_surface_data(self, yaml_data: Dict, config: SemitipConfig):
        """比較表面數據"""
        logger.info("🔍 檢查表面數據...")
        
        surface = yaml_data.get('surface', {})
        regions = surface.get('regions', [])
        
        logger.info(f"  原始表面區域數量: {len(regions)}")
        logger.info(f"  解析表面區域數量: {len(config.surface_regions)}")
        logger.info(f"  表面區域數量一致: {'✅' if len(regions) == len(config.surface_regions) else '❌'}")
        
        # 檢查溫度相依性
        yaml_temp_dep = surface.get('temperature_dependence')
        config_temp_dep = config.surface_temperature_dependence
        logger.info(f"  溫度相依性: {yaml_temp_dep} -> {config_temp_dep} {'✅' if yaml_temp_dep == config_temp_dep else '❌'}")
    
    def _compare_output_data(self, yaml_data: Dict, config: SemitipConfig):
        """比較輸出數據"""
        logger.info("🔍 檢查輸出數據...")
        
        output = yaml_data.get('output', {})
        
        # 檢查各種輸出設定
        checks = [
            ('basic_output', 'output_basic'),
            ('equipotential_contours', 'output_contours'),
            ('full_potential', 'output_full_potential'),
            ('num_contours', 'num_contours'),
            ('contour_spacing', 'contour_spacing'),
            ('contour_angle', 'contour_angle')
        ]
        
        for yaml_key, config_attr in checks:
            yaml_val = output.get(yaml_key)
            config_val = getattr(config, config_attr)
            status = '✅' if yaml_val == config_val else '❌'
            logger.info(f"  {yaml_key}: {yaml_val} -> {config_val} {status}")
    
    def test_error_handling(self):
        """測試錯誤處理"""
        logger.info(f"\n{'='*60}")
        logger.info("測試錯誤處理機制")
        logger.info('='*60)
        
        # 測試不存在的檔案
        self._test_nonexistent_file()
        
        # 測試格式錯誤的 YAML
        self._test_malformed_yaml()
        
        # 測試缺失必要欄位
        self._test_missing_fields()
        
        # 測試無效數值
        self._test_invalid_values()
    
    def _test_nonexistent_file(self):
        """測試不存在的檔案"""
        logger.info("🔍 測試不存在的檔案...")
        
        try:
            config = self.reader.load_config("nonexistent_file.yaml")
            logger.error("❌ 應該要拋出例外")
        except FileNotFoundError:
            logger.info("✅ 正確處理檔案不存在的情況")
        except Exception as e:
            logger.warning(f"⚠️ 拋出了非預期的例外: {e}")
    
    def _test_malformed_yaml(self):
        """測試格式錯誤的 YAML"""
        logger.info("🔍 測試格式錯誤的 YAML...")
        
        # 創建格式錯誤的 YAML 檔案
        malformed_yaml = """
version: "1.0"
simulation_type: "MultInt"
environment:
  temperature: 300.0
  dielectric_constant: 12.9
tip:
  separation: 1.0
  radius: [1.0  # 缺少右方括號，格式錯誤
"""
        
        temp_file = BASE_DIR / "tests" / "temp_malformed.yaml"
        try:
            with open(temp_file, 'w', encoding='utf-8') as f:
                f.write(malformed_yaml)
            
            config = self.reader.load_config(temp_file)
            logger.error("❌ 應該要拋出 YAML 解析例外")
            
        except yaml.YAMLError:
            logger.info("✅ 正確處理 YAML 格式錯誤")
        except Exception as e:
            logger.warning(f"⚠️ 拋出了非預期的例外: {e}")
        finally:
            if temp_file.exists():
                temp_file.unlink()
    
    def _test_missing_fields(self):
        """測試缺失必要欄位"""
        logger.info("🔍 測試缺失必要欄位...")
        
        # 創建缺少必要欄位的 YAML
        minimal_yaml = """
version: "1.0"
# 缺少 simulation_type
environment:
  temperature: 300.0
  dielectric_constant: 12.9
"""
        
        temp_file = BASE_DIR / "tests" / "temp_minimal.yaml"
        try:
            with open(temp_file, 'w', encoding='utf-8') as f:
                f.write(minimal_yaml)
            
            yaml_data = yaml.safe_load(minimal_yaml)
            config = self.reader._manual_yaml_to_config(yaml_data)
            
            # 檢查是否使用了預設值
            if config.simulation_type == "MultInt":
                logger.info("✅ 正確使用預設模擬類型")
            else:
                logger.warning(f"⚠️ 預設模擬類型不正確: {config.simulation_type}")
            
        except Exception as e:
            logger.warning(f"⚠️ 處理缺失欄位時發生例外: {e}")
        finally:
            if temp_file.exists():
                temp_file.unlink()
    
    def _test_invalid_values(self):
        """測試無效數值"""
        logger.info("🔍 測試無效數值...")
        
        # 創建含有無效數值的 YAML
        invalid_yaml = """
version: "1.0"
simulation_type: "MultInt"
environment:
  temperature: -100.0  # 負溫度，物理上不合理
  dielectric_constant: 12.9
tip:
  separation: 1.0
  radius: -5.0  # 負半徑，物理上不合理
  position:
    x: 0.0
    y: 0.0
"""
        
        temp_file = BASE_DIR / "tests" / "temp_invalid.yaml"
        try:
            with open(temp_file, 'w', encoding='utf-8') as f:
                f.write(invalid_yaml)
            
            yaml_data = yaml.safe_load(invalid_yaml)
            config = self.reader._manual_yaml_to_config(yaml_data)
            
            # 檢查是否成功解析（即使數值不合理）
            logger.info(f"  解析後的溫度: {config.temperature} K")
            logger.info(f"  解析後的探針半徑: {config.tip.radius} nm")
            
            if config.temperature == -100.0:
                logger.warning("⚠️ 系統接受了負溫度值，可能需要添加驗證")
            
            if config.tip.radius == -5.0:
                logger.warning("⚠️ 系統接受了負半徑值，可能需要添加驗證")
            
        except Exception as e:
            logger.info(f"✅ 系統正確拒絕了無效數值: {e}")
        finally:
            if temp_file.exists():
                temp_file.unlink()
    
    def test_performance(self, file_path: Path):
        """測試性能"""
        logger.info(f"\n{'='*60}")
        logger.info(f"測試解析性能: {file_path.name}")
        logger.info('='*60)
        
        import time
        
        # 測試多次載入的時間
        times = []
        num_tests = 10
        
        for i in range(num_tests):
            start_time = time.time()
            config = self.reader.load_config(file_path)
            end_time = time.time()
            times.append(end_time - start_time)
        
        avg_time = sum(times) / len(times)
        min_time = min(times)
        max_time = max(times)
        
        logger.info(f"  平均載入時間: {avg_time:.4f} 秒")
        logger.info(f"  最快載入時間: {min_time:.4f} 秒")
        logger.info(f"  最慢載入時間: {max_time:.4f} 秒")
        
        if avg_time < 0.1:
            logger.info("✅ 載入性能優秀")
        elif avg_time < 0.5:
            logger.info("👍 載入性能良好")
        else:
            logger.warning("⚠️ 載入性能可能需要最佳化")


def main():
    """主測試函數"""
    logger.info("開始執行數據完整性和邊界情況測試")
    
    tester = DataIntegrityTester()
    
    # 測試檔案
    test_files = [
        BASE_DIR / "data" / "input" / "MultInt_config.yaml",
        BASE_DIR / "data" / "input" / "MultPlane_config.yaml"
    ]
    
    for file_path in test_files:
        # 測試數據準確性
        tester.test_data_accuracy(file_path)
        
        # 測試性能
        tester.test_performance(file_path)
    
    # 測試錯誤處理
    tester.test_error_handling()
    
    logger.info(f"\n{'='*60}")
    logger.info("數據完整性測試完成")
    logger.info('='*60)


if __name__ == "__main__":
    main()
