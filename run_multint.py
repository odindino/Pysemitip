import argparse
import logging
import sys
import yaml

from src.simulation.multint import MultInt
# 【註】在您的版本中，您將 validate_config 改為了 SemitipConfig.validate
# 我將採用您更新後的版本。
from src.core.config_schema import SemitipConfig 

def setup_logging():
    """設定日誌記錄，將資訊同時輸出到控制台和檔案。"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler("logs/simulation.log", mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )

def main():
    """
    主執行函式：解析參數、讀取設定、執行模擬。
    """
    setup_logging()
    
    # --- 【修改】更新參數解析邏輯 ---
    parser = argparse.ArgumentParser(description="執行 Pysemitip 自洽模擬。")
    # 1. 將設定檔路徑改為 "位置參數"
    parser.add_argument(
        'config_path',  # 參數名稱，不再有 '--'
        type=str,
        help='設定檔的路徑 (YAML 格式)。'
    )
    # 2. 新增對 "--plot" 旗標的定義
    parser.add_argument(
        '--plot',
        action='store_true',  # 當出現此旗標時，其值為 True
        help='執行完畢後繪製結果圖。'
    )
    args = parser.parse_args()
    
    # --- 讀取與驗證設定檔 ---
    logger = logging.getLogger(__name__)
    # 【修改】使用新的參數名稱 args.config_path
    logger.info(f"正在從 {args.config_path} 讀取設定檔...")
    try:
        with open(args.config_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        
        # 驗證設定檔結構
        validated_config = SemitipConfig(**config_dict)
        logger.info("設定檔讀取並驗證成功。")
        
    except FileNotFoundError:
        logger.error(f"錯誤：找不到設定檔 {args.config_path}")
        return
    except Exception as e:
        logger.error(f"讀取或驗證設定檔時發生錯誤: {e}")
        return

    # --- 執行模擬 ---
    try:
        simulation = MultInt(validated_config)
        simulation.run_self_consistent_loop()
            
        logger.info("模擬執行完畢。")
        
        # 我們可以在這裡使用 args.plot 的值
        if args.plot:
            logger.info("需要繪圖... (繪圖功能待實現)")
            # 在此處呼叫繪圖函式
            # plot_results(simulation.results)

    except Exception as e:
        logger.error("模擬過程中發生未預期的錯誤。", exc_info=True)


if __name__ == '__main__':
    main()