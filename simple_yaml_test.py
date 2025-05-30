"""
簡化版 YAML 測試腳本
"""

import yaml
from pathlib import Path
import logging
import sys

# 設定日誌
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def test_yaml_file(file_path):
    """測試 YAML 文件是否可讀取"""
    path = Path(file_path)
    logger.info(f"測試 YAML 文件: {path}")
    
    if not path.exists():
        logger.error(f"文件不存在: {path}")
        return False
        
    try:
        with open(path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)
            
        logger.info("YAML 文件讀取成功")
        logger.info(f"文件內容概要:")
        logger.info(f"- 模擬類型: {data.get('simulation_type')}")
        logger.info(f"- 溫度: {data.get('temperature')} K")
        logger.info(f"- 半導體區域數量: {len(data.get('semiconductor_regions', []))}")
        logger.info(f"- 表面區域數量: {len(data.get('surface_regions', []))}")
        
        return True
    except Exception as e:
        logger.error(f"讀取 YAML 文件時發生錯誤: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def main():
    """主函數"""
    base_dir = Path(__file__).resolve().parent
    conversion_dir = base_dir / "converted_configs"
    
    # 測試 MultInt 配置文件
    multint_path = conversion_dir / "MultInt_config.yaml"
    test_yaml_file(multint_path)
    
    print("\n" + "-" * 40 + "\n")
    
    # 測試 MultPlane 配置文件
    multplane_path = conversion_dir / "MultPlane_config.yaml"
    test_yaml_file(multplane_path)

if __name__ == "__main__":
    main()
