#!/usr/bin/env python3
"""
全面的 YAML 結構測試程序

此測試程序用於驗證 Pysemitip 項目中的新分層 YAML 結構：
1. 測試 YAML 結構讀取準確性
2. 測試 YAML 到配置物件的轉換過程
3. 測試配置物件到 YAML 的回寫功能（往返轉換）
4. 測試新舊格式的向後相容性
5. 測試錯誤處理和邊界情況
6. 對比新舊結構的解析結果

作者: Pysemitip 開發團隊
日期: 2024-01-02
版本: 1.0
"""

import yaml
import sys
import logging
import traceback
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import asdict

# 設定日誌格式
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('test_results.log', encoding='utf-8')
    ]
)
logger = logging.getLogger(__name__)

# 確保能導入專案模組
BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

try:
    from filereader import YamlConfigReader
    from config_schema import SemitipConfig
    logger.info("✅ 成功導入必要的模組")
except ImportError as e:
    logger.error(f"❌ 模組導入失敗: {e}")
    sys.exit(1)


class YamlStructureValidator:
    """YAML 結構驗證器"""
    
    def __init__(self):
        self.reader = YamlConfigReader()
        self.test_results = {
            'total_tests': 0,
            'passed_tests': 0,
            'failed_tests': 0,
            'errors': []
        }
    
    def log_test_result(self, test_name: str, success: bool, error_msg: str = ""):
        """記錄測試結果"""
        self.test_results['total_tests'] += 1
        if success:
            self.test_results['passed_tests'] += 1
            logger.info(f"✅ {test_name}: 通過")
        else:
            self.test_results['failed_tests'] += 1
            self.test_results['errors'].append(f"{test_name}: {error_msg}")
            logger.error(f"❌ {test_name}: 失敗 - {error_msg}")
    
    def test_yaml_file_loading(self, file_path: Path) -> Optional[Dict]:
        """測試 YAML 檔案載入"""
        test_name = f"載入 YAML 檔案: {file_path.name}"
        
        try:
            if not file_path.exists():
                self.log_test_result(test_name, False, f"檔案不存在: {file_path}")
                return None
            
            with open(file_path, 'r', encoding='utf-8') as f:
                yaml_data = yaml.safe_load(f)
            
            if yaml_data is None:
                self.log_test_result(test_name, False, "YAML 檔案為空或格式錯誤")
                return None
            
            # 檢查基本結構
            required_keys = ['version', 'simulation_type']
            missing_keys = [key for key in required_keys if key not in yaml_data]
            
            if missing_keys:
                self.log_test_result(test_name, False, f"缺少必要欄位: {missing_keys}")
                return None
            
            self.log_test_result(test_name, True)
            return yaml_data
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
            return None
    
    def test_hierarchical_structure(self, yaml_data: Dict, file_name: str):
        """測試分層結構的正確性"""
        
        # 測試環境設定塊
        self._test_environment_block(yaml_data, file_name)
        
        # 測試 tip 配置塊
        self._test_tip_block(yaml_data, file_name)
        
        # 測試半導體配置塊
        self._test_semiconductor_block(yaml_data, file_name)
        
        # 測試表面配置塊
        self._test_surface_block(yaml_data, file_name)
        
        # 測試輸出配置塊
        self._test_output_block(yaml_data, file_name)
        
        # 測試模擬特定配置
        self._test_simulation_specific_block(yaml_data, file_name)
    
    def _test_environment_block(self, yaml_data: Dict, file_name: str):
        """測試環境配置塊"""
        test_name = f"{file_name} - 環境配置塊"
        
        try:
            environment = yaml_data.get('environment', {})
            if not environment:
                self.log_test_result(test_name, False, "缺少 environment 塊")
                return
            
            # 檢查必要欄位
            required_fields = ['temperature', 'dielectric_constant']
            missing_fields = [field for field in required_fields if field not in environment]
            
            if missing_fields:
                self.log_test_result(test_name, False, f"environment 塊缺少欄位: {missing_fields}")
                return
            
            # 檢查數值範圍
            temp = environment['temperature']
            dielectric = environment['dielectric_constant']
            
            if not (0 < temp < 1000):
                self.log_test_result(test_name, False, f"溫度值異常: {temp}")
                return
            
            if not (1 < dielectric < 100):
                self.log_test_result(test_name, False, f"介電常數值異常: {dielectric}")
                return
            
            self.log_test_result(test_name, True)
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
    
    def _test_tip_block(self, yaml_data: Dict, file_name: str):
        """測試探針配置塊"""
        test_name = f"{file_name} - Tip 配置塊"
        
        try:
            tip = yaml_data.get('tip', {})
            if not tip:
                self.log_test_result(test_name, False, "缺少 tip 塊")
                return
            
            # 檢查 position 結構
            if 'position' in tip:
                position = tip['position']
                if not isinstance(position, dict):
                    self.log_test_result(test_name, False, "tip.position 應該是字典結構")
                    return
                
                required_pos_fields = ['x', 'y']
                missing_pos_fields = [field for field in required_pos_fields if field not in position]
                
                if missing_pos_fields:
                    self.log_test_result(test_name, False, f"tip.position 缺少欄位: {missing_pos_fields}")
                    return
            
            # 檢查基本 tip 欄位
            required_tip_fields = ['separation', 'radius']
            missing_tip_fields = [field for field in required_tip_fields if field not in tip]
            
            if missing_tip_fields:
                self.log_test_result(test_name, False, f"tip 配置缺少欄位: {missing_tip_fields}")
                return
            
            self.log_test_result(test_name, True)
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
    
    def _test_semiconductor_block(self, yaml_data: Dict, file_name: str):
        """測試半導體配置塊"""
        test_name = f"{file_name} - 半導體配置塊"
        
        try:
            semiconductor = yaml_data.get('semiconductor', {})
            if not semiconductor:
                self.log_test_result(test_name, False, "缺少 semiconductor 塊")
                return
            
            # 檢查 regions
            regions = semiconductor.get('regions', [])
            if not regions:
                self.log_test_result(test_name, False, "semiconductor.regions 為空")
                return
            
            # 檢查每個區域的結構
            for i, region in enumerate(regions):
                if not isinstance(region, dict):
                    self.log_test_result(test_name, False, f"region[{i}] 應該是字典結構")
                    return
                
                # 檢查有效質量結構
                if 'effective_mass' in region:
                    effective_mass = region['effective_mass']
                    if not isinstance(effective_mass, dict):
                        self.log_test_result(test_name, False, f"region[{i}].effective_mass 應該是字典結構")
                        return
                    
                    required_mass_fields = ['conduction_band', 'heavy_hole', 'light_hole', 'split_off_hole']
                    missing_mass_fields = [field for field in required_mass_fields if field not in effective_mass]
                    
                    if missing_mass_fields:
                        self.log_test_result(test_name, False, f"region[{i}].effective_mass 缺少欄位: {missing_mass_fields}")
                        return
            
            # 檢查電子親和力
            if 'electron_affinity' not in semiconductor:
                self.log_test_result(test_name, False, "semiconductor 缺少 electron_affinity")
                return
            
            self.log_test_result(test_name, True)
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
    
    def _test_surface_block(self, yaml_data: Dict, file_name: str):
        """測試表面配置塊"""
        test_name = f"{file_name} - 表面配置塊"
        
        try:
            surface = yaml_data.get('surface', {})
            if not surface:
                # 表面配置可能是可選的
                self.log_test_result(test_name, True, "表面配置為可選項")
                return
            
            # 檢查 regions
            regions = surface.get('regions', [])
            
            # 檢查每個表面區域的結構
            for i, region in enumerate(regions):
                if not isinstance(region, dict):
                    self.log_test_result(test_name, False, f"surface region[{i}] 應該是字典結構")
                    return
            
            # 檢查溫度相依性
            if 'temperature_dependence' not in surface:
                self.log_test_result(test_name, False, "surface 缺少 temperature_dependence")
                return
            
            self.log_test_result(test_name, True)
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
    
    def _test_output_block(self, yaml_data: Dict, file_name: str):
        """測試輸出配置塊"""
        test_name = f"{file_name} - 輸出配置塊"
        
        try:
            output = yaml_data.get('output', {})
            if not output:
                self.log_test_result(test_name, False, "缺少 output 塊")
                return
            
            # 檢查基本輸出設定
            required_output_fields = ['basic_output', 'equipotential_contours', 'full_potential']
            missing_output_fields = [field for field in required_output_fields if field not in output]
            
            if missing_output_fields:
                self.log_test_result(test_name, False, f"output 配置缺少欄位: {missing_output_fields}")
                return
            
            self.log_test_result(test_name, True)
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
    
    def _test_simulation_specific_block(self, yaml_data: Dict, file_name: str):
        """測試模擬特定配置塊"""
        test_name = f"{file_name} - 模擬特定配置塊"
        
        try:
            simulation_type = yaml_data.get('simulation_type', '')
            
            if simulation_type == 'MultInt':
                multint_specific = yaml_data.get('multint_specific', {})
                if not multint_specific:
                    self.log_test_result(test_name, False, "MultInt 模擬缺少 multint_specific 塊")
                    return
                    
            elif simulation_type == 'MultPlane':
                multplane_specific = yaml_data.get('multplane_specific', {})
                if not multplane_specific:
                    self.log_test_result(test_name, False, "MultPlane 模擬缺少 multplane_specific 塊")
                    return
            
            self.log_test_result(test_name, True)
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
    
    def test_yaml_to_config_conversion(self, file_path: Path):
        """測試 YAML 到配置物件的轉換"""
        test_name = f"YAML→Config 轉換: {file_path.name}"
        
        try:
            config = self.reader.load_config(file_path)
            
            if config is None:
                self.log_test_result(test_name, False, "轉換結果為 None")
                return None
            
            # 驗證基本屬性
            if not hasattr(config, 'simulation_type'):
                self.log_test_result(test_name, False, "配置物件缺少 simulation_type")
                return None
            
            if not hasattr(config, 'tip'):
                self.log_test_result(test_name, False, "配置物件缺少 tip")
                return None
            
            # 檢查分層結構是否正確轉換
            if not config.semiconductor_regions:
                self.log_test_result(test_name, False, "半導體區域為空")
                return None
            
            self.log_test_result(test_name, True)
            return config
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
            return None
    
    def test_config_to_yaml_conversion(self, config: SemitipConfig, original_file: str):
        """測試配置物件到 YAML 的轉換"""
        test_name = f"Config→YAML 轉換: {original_file}"
        
        try:
            # 使用內部方法進行轉換
            yaml_dict = self.reader._config_to_yaml(config)
            
            if not yaml_dict:
                self.log_test_result(test_name, False, "轉換結果為空")
                return None
            
            # 檢查必要的分層結構
            required_blocks = ['environment', 'tip', 'semiconductor', 'output']
            missing_blocks = [block for block in required_blocks if block not in yaml_dict]
            
            if missing_blocks:
                self.log_test_result(test_name, False, f"轉換結果缺少塊: {missing_blocks}")
                return None
            
            # 檢查 tip 位置結構
            if 'position' not in yaml_dict['tip']:
                self.log_test_result(test_name, False, "tip 配置缺少 position 結構")
                return None
            
            self.log_test_result(test_name, True)
            return yaml_dict
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
            return None
    
    def test_round_trip_conversion(self, file_path: Path):
        """測試往返轉換 (YAML → Config → YAML)"""
        test_name = f"往返轉換測試: {file_path.name}"
        
        try:
            # 第一次轉換: YAML → Config
            original_config = self.reader.load_config(file_path)
            if original_config is None:
                self.log_test_result(test_name, False, "第一次轉換失敗")
                return
            
            # 第二次轉換: Config → YAML
            yaml_dict = self.reader._config_to_yaml(original_config)
            if yaml_dict is None:
                self.log_test_result(test_name, False, "第二次轉換失敗")
                return
            
            # 第三次轉換: YAML → Config (再次)
            yaml_str = yaml.dump(yaml_dict, default_flow_style=False, allow_unicode=True)
            yaml_data_roundtrip = yaml.safe_load(yaml_str)
            final_config = self.reader._manual_yaml_to_config(yaml_data_roundtrip)
            
            # 比較關鍵屬性
            if original_config.simulation_type != final_config.simulation_type:
                self.log_test_result(test_name, False, "模擬類型不一致")
                return
            
            if original_config.temperature != final_config.temperature:
                self.log_test_result(test_name, False, "溫度值不一致")
                return
            
            if original_config.tip.radius != final_config.tip.radius:
                self.log_test_result(test_name, False, "探針半徑不一致")
                return
            
            self.log_test_result(test_name, True)
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
    
    def test_backward_compatibility(self, file_path: Path):
        """測試向後相容性"""
        test_name = f"向後相容性測試: {file_path.name}"
        
        try:
            # 讀取現有檔案
            with open(file_path, 'r', encoding='utf-8') as f:
                original_yaml = yaml.safe_load(f)
            
            # 創建舊格式的測試資料
            old_format_yaml = self._convert_to_old_format(original_yaml)
            
            # 測試舊格式是否能正確解析
            old_config = self.reader._manual_yaml_to_config(old_format_yaml)
            new_config = self.reader._manual_yaml_to_config(original_yaml)
            
            # 比較關鍵值是否一致
            if old_config.simulation_type != new_config.simulation_type:
                self.log_test_result(test_name, False, "模擬類型解析不一致")
                return
            
            if old_config.temperature != new_config.temperature:
                self.log_test_result(test_name, False, "溫度解析不一致")
                return
            
            self.log_test_result(test_name, True)
            
        except Exception as e:
            self.log_test_result(test_name, False, str(e))
    
    def _convert_to_old_format(self, new_format_yaml: Dict) -> Dict:
        """將新格式轉換為舊格式用於測試"""
        old_format = {}
        
        # 基本資訊
        old_format['version'] = new_format_yaml.get('version', '1.0')
        old_format['simulation_type'] = new_format_yaml.get('simulation_type', 'MultInt')
        
        # 環境設定 (展平)
        environment = new_format_yaml.get('environment', {})
        old_format['temperature'] = environment.get('temperature', 300.0)
        old_format['dielectric_constant'] = environment.get('dielectric_constant', 12.9)
        
        # Tip 配置 (展平 position)
        tip = new_format_yaml.get('tip', {})
        old_format['tip'] = tip.copy()
        if 'position' in tip:
            position = tip['position']
            old_format['tip']['x_position'] = position.get('x', 0.0)
            old_format['tip']['y_position'] = position.get('y', 0.0)
            del old_format['tip']['position']
        
        # 半導體區域 (展平)
        semiconductor = new_format_yaml.get('semiconductor', {})
        old_format['semiconductor_regions'] = semiconductor.get('regions', [])
        old_format['electron_affinity'] = semiconductor.get('electron_affinity', 4.07)
        
        # 表面區域 (展平)
        surface = new_format_yaml.get('surface', {})
        old_format['surface_regions'] = surface.get('regions', [])
        old_format['surface_temperature_dependence'] = surface.get('temperature_dependence', False)
        
        # 輸出設定 (展平)
        output = new_format_yaml.get('output', {})
        old_format['output_basic'] = output.get('basic_output', True)
        old_format['output_contours'] = output.get('equipotential_contours', False)
        old_format['output_full_potential'] = output.get('full_potential', False)
        old_format['num_contours'] = output.get('num_contours', 6)
        old_format['contour_spacing'] = output.get('contour_spacing', 0.0)
        old_format['contour_angle'] = output.get('contour_angle', 0.0)
        
        # 特定配置
        if 'multint_specific' in new_format_yaml:
            old_format['multint_config'] = new_format_yaml['multint_specific']
        if 'multplane_specific' in new_format_yaml:
            old_format['multplane_config'] = new_format_yaml['multplane_specific']
        
        # 其他配置
        for key in ['grid', 'computation', 'voltage_scan']:
            if key in new_format_yaml:
                old_format[key] = new_format_yaml[key]
        
        return old_format
    
    def generate_test_report(self):
        """生成測試報告"""
        logger.info("\n" + "="*80)
        logger.info("YAML 結構測試報告")
        logger.info("="*80)
        logger.info(f"總測試數: {self.test_results['total_tests']}")
        logger.info(f"通過測試: {self.test_results['passed_tests']}")
        logger.info(f"失敗測試: {self.test_results['failed_tests']}")
        
        if self.test_results['failed_tests'] > 0:
            logger.info("\n失敗的測試:")
            for error in self.test_results['errors']:
                logger.error(f"  - {error}")
        
        success_rate = (self.test_results['passed_tests'] / self.test_results['total_tests']) * 100
        logger.info(f"\n成功率: {success_rate:.1f}%")
        
        if success_rate >= 90:
            logger.info("🎉 測試結果優秀！")
        elif success_rate >= 70:
            logger.info("👍 測試結果良好，但還有改進空間")
        else:
            logger.warning("⚠️  測試結果需要改進")
        
        logger.info("="*80)


def main():
    """主測試函數"""
    logger.info("開始執行 YAML 結構全面測試")
    
    validator = YamlStructureValidator()
    
    # 定義測試檔案
    test_files = [
        BASE_DIR / "Import_Files" / "MultInt_config.yaml",
        BASE_DIR / "Import_Files" / "MultPlane_config.yaml"
    ]
    
    for file_path in test_files:
        logger.info(f"\n{'='*60}")
        logger.info(f"測試檔案: {file_path.name}")
        logger.info('='*60)
        
        # 1. 測試 YAML 檔案載入
        yaml_data = validator.test_yaml_file_loading(file_path)
        if yaml_data is None:
            continue
        
        # 2. 測試分層結構
        validator.test_hierarchical_structure(yaml_data, file_path.name)
        
        # 3. 測試 YAML 到配置轉換
        config = validator.test_yaml_to_config_conversion(file_path)
        if config is None:
            continue
        
        # 4. 測試配置到 YAML 轉換
        yaml_dict = validator.test_config_to_yaml_conversion(config, file_path.name)
        
        # 5. 測試往返轉換
        validator.test_round_trip_conversion(file_path)
        
        # 6. 測試向後相容性
        validator.test_backward_compatibility(file_path)
    
    # 生成測試報告
    validator.generate_test_report()
    
    # 返回成功碼
    return 0 if validator.test_results['failed_tests'] == 0 else 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
