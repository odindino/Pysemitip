import argparse
import logging
import sys
import yaml
import os # 新增導入
from datetime import datetime # 新增導入

from src.simulation.multint import MultInt
# 【註】在您的版本中，您將 validate_config 改為了 SemitipConfig.validate
# 我將採用您更新後的版本。
from src.core.config_schema import SemitipConfig 

def setup_logging(config_path: str): # 修改：增加 config_path 參數
    """設定日誌記錄，將資訊同時輸出到控制台和檔案。

    日誌將儲存到 data/output/results/<config_name>_<timestamp>/logs/simulation.log
    """
    # 從設定檔路徑中提取不含副檔名的檔案名稱
    config_filename = os.path.basename(config_path)
    config_name_without_ext = os.path.splitext(config_filename)[0]
    
    # 產生基於時間戳和設定檔名的唯一執行目錄
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_name = f"{config_name_without_ext}_{timestamp}"
    
    base_output_dir = "data/output/results"
    run_output_dir = os.path.join(base_output_dir, run_name)
    log_dir = os.path.join(run_output_dir, "logs") # 日誌儲存在執行專屬目錄下的 logs 子目錄

    if not os.path.exists(log_dir):
        os.makedirs(log_dir) # 創建 logs 目錄以及其父目錄 run_output_dir
        
    log_file_path = os.path.join(log_dir, "simulation.log")

    # 移除現有的 handlers，以避免重複記錄 (如果此函數可能被多次呼叫或在 Jupyter 環境中)
    # 這對於確保每次執行 setup_logging 時都能正確設定 handlers 很重要
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        handler.close() # 關閉 handler 以釋放檔案資源

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return run_output_dir # 返回此次執行的輸出目錄路徑

def main():
    """
    主執行函式：解析參數、讀取設定、執行模擬。
    """
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
    
    # --- 設定日誌 ---
    # 修改：將設定檔路徑傳遞給 setup_logging，並獲取執行輸出目錄
    run_output_dir = setup_logging(args.config_path)
    
    logger = logging.getLogger(__name__) # 在 setup_logging 之後獲取 logger
    logger.info(f"模擬執行專屬輸出目錄: {run_output_dir}")
    
    # --- 讀取與驗證設定檔 ---
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
        # 修改：將 run_output_dir 傳遞給 MultInt 的建構函數
        simulation = MultInt(validated_config, run_output_dir)
        simulation.run_self_consistent_loop()
            
        logger.info("模擬執行完畢。")
        
        # 我們可以在這裡使用 args.plot 的值
        if args.plot:
            logger.info("需要繪圖... (繪圖功能待實現)")
            # 在此處呼叫繪圖函式
            # plot_results(simulation.results, run_output_dir) # 繪圖結果也應儲存到 run_output_dir

    except Exception as e:
        logger.error("模擬過程中發生未預期的錯誤。", exc_info=True)


if __name__ == '__main__':
    main()