"""
測試 YAML 配置文件讀取腳本

此腳本用於測試 filereader.py 對 YAML 配置文件的讀取功能。
"""

import sys
import os
import traceback
from pathlib import Path
import logging

# 設定日誌
logging.basicConfig(level=logging.DEBUG, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 確保在腳本所在目錄中運行
BASE_DIR = Path(__file__).resolve().parent
os.chdir(BASE_DIR)
logger.info(f"工作目錄設置為: {BASE_DIR}")

# 確保導入路徑正確
sys.path.insert(0, str(BASE_DIR))
logger.info(f"導入路徑: {sys.path}")

# 檢查文件是否存在
filereader_path = BASE_DIR / "filereader.py"
config_schema_path = BASE_DIR / "config_schema.py"

if not filereader_path.exists():
    logger.error(f"找不到 filereader.py: {filereader_path}")
    sys.exit(1)
else:
    logger.info(f"找到 filereader.py: {filereader_path}")

if not config_schema_path.exists():
    logger.error(f"找不到 config_schema.py: {config_schema_path}")
    sys.exit(1)
else:
    logger.info(f"找到 config_schema.py: {config_schema_path}")

# 導入 filereader 模組
try:
    from filereader import YamlConfigReader
    logger.info("成功導入 YamlConfigReader")
except ImportError as e:
    logger.error(f"無法導入 YamlConfigReader: {e}")
    logger.error(f"錯誤追蹤: {traceback.format_exc()}")
    sys.exit(1)
except Exception as e:
    logger.error(f"導入 YamlConfigReader 時發生未預期的錯誤: {e}")
    logger.error(f"錯誤追蹤: {traceback.format_exc()}")
    sys.exit(1)

def test_read_yaml(yaml_path):
    """測試讀取 YAML 文件"""
    reader = YamlConfigReader()
    try:
        full_path = Path(yaml_path).resolve()
        logger.info(f"正在讀取 YAML 文件: {full_path}")
        
        # 檢查文件是否存在
        if not full_path.exists():
            logger.error(f"找不到 YAML 文件: {full_path}")
            return None
            
        config = reader.load_config(full_path)
        logger.info(f"成功讀取 YAML 文件: {full_path}")
        
        # 檢查配置對象是否正確創建
        if config is None:
            logger.error("讀取配置時返回 None")
            return None
            
        logger.info(f"模擬類型: {config.simulation_type}")
        logger.info(f"溫度: {config.temperature} K")
        logger.info(f"半導體區域數量: {len(config.semiconductor_regions)}")
        logger.info(f"表面區域數量: {len(config.surface_regions)}")
        
        # 檢查特有配置
        if config.simulation_type == "MultInt":
            if config.multint_config:
                logger.info(f"MultInt 配置: {config.multint_config}")
            else:
                logger.warning("缺少 MultInt 特有配置")
        elif config.simulation_type == "MultPlane":
            if config.multplane_config:
                logger.info(f"MultPlane 配置: {config.multplane_config}")
            else:
                logger.warning("缺少 MultPlane 特有配置")
        
        # 驗證配置 (暫時跳過以方便測試)
        try:
            # 暫時注釋掉驗證，如果需要，可以取消注釋
            config.validate()
        except ValueError as e:
            logger.error(f"配置驗證失敗: {e}")
        except Exception as e:
            logger.error(f"驗證過程中發生未預期錯誤: {e}")
            logger.error(f"錯誤追蹤: {traceback.format_exc()}")
        
        return config
    except Exception as e:
        logger.error(f"讀取 YAML 文件時發生錯誤: {e}")
        logger.error(f"錯誤追蹤: {traceback.format_exc()}")
        return None

def main():
    """主函數"""
    try:
        # 測試 MultInt YAML 文件
        multint_path = BASE_DIR / "converted_configs" / "MultInt_config.yaml"
        logger.info(f"測試 MultInt 配置文件: {multint_path}")
        multint_config = test_read_yaml(multint_path)
        
        # 測試 MultPlane YAML 文件
        multplane_path = BASE_DIR / "converted_configs" / "MultPlane_config.yaml"
        logger.info(f"測試 MultPlane 配置文件: {multplane_path}")
        multplane_config = test_read_yaml(multplane_path)
        
        logger.info("測試完成")
    except Exception as e:
        logger.error(f"測試過程中發生未預期錯誤: {e}")
        logger.error(f"錯誤追蹤: {traceback.format_exc()}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.critical(f"程序執行失敗: {e}")
        logger.critical(f"錯誤追蹤: {traceback.format_exc()}")
        sys.exit(1)
