#!/usr/bin/env python3
"""
測試新的 config_schema.py 對階層化 YAML 結構的支援

此腳本測試：
1. 載入 MultInt 和 MultPlane 配置檔案
2. 檢查是否能夠直接轉換（不需要手動轉換）
3. 驗證配置物件的屬性是否正確
"""

import sys
from pathlib import Path

# 添加項目路徑
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from filereader import YamlConfigReader
import logging

# 設定日誌
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_config_loading():
    """測試配置檔案載入"""
    reader = YamlConfigReader()
    
    # 測試檔案路徑
    multint_file = project_root / "Import_Files" / "MultInt_config.yaml"
    multplane_file = project_root / "Import_Files" / "MultPlane_config.yaml"
    
    print("=" * 60)
    print("🧪 測試新版 config_schema.py 對階層化 YAML 的支援")
    print("=" * 60)
    
    # 測試 MultInt 配置
    print("\n📂 測試 MultInt 配置檔案...")
    try:
        multint_config = reader.load_config(multint_file)
        print("✅ MultInt 配置載入成功！")
        
        # 檢查階層結構
        print(f"   ├─ 模擬類型: {multint_config.simulation_type}")
        print(f"   ├─ 溫度: {multint_config.environment.temperature} K")
        print(f"   ├─ 介電常數: {multint_config.environment.dielectric_constant}")
        print(f"   ├─ 探針位置: ({multint_config.tip.position.x}, {multint_config.tip.position.y}) nm")
        print(f"   ├─ 半導體區域數: {len(multint_config.semiconductor.regions)}")
        print(f"   ├─ 表面區域數: {len(multint_config.surface.regions)}")
        print(f"   └─ MultInt 特有參數: {multint_config.multint_specific is not None}")
        
        # 測試向後相容性
        print(f"\n🔄 向後相容性測試:")
        print(f"   ├─ config.temperature: {multint_config.temperature}")
        print(f"   ├─ config.tip.x_position: {multint_config.tip.x_position}")
        print(f"   └─ config.semiconductor_regions: {len(multint_config.semiconductor_regions)}")
        
    except Exception as e:
        print(f"❌ MultInt 配置載入失敗: {e}")
        return False
    
    # 測試 MultPlane 配置
    print("\n📂 測試 MultPlane 配置檔案...")
    try:
        multplane_config = reader.load_config(multplane_file)
        print("✅ MultPlane 配置載入成功！")
        
        # 檢查階層結構
        print(f"   ├─ 模擬類型: {multplane_config.simulation_type}")
        print(f"   ├─ 溫度: {multplane_config.environment.temperature} K")
        print(f"   ├─ 介電常數: {multplane_config.environment.dielectric_constant}")
        print(f"   ├─ 探針位置: ({multplane_config.tip.position.x}, {multplane_config.tip.position.y}) nm")
        print(f"   ├─ 半導體區域數: {len(multplane_config.semiconductor.regions)}")
        print(f"   ├─ 表面區域數: {len(multplane_config.surface.regions)}")
        print(f"   └─ MultPlane 特有參數: {multplane_config.multplane_specific is not None}")
        
        # 測試向後相容性
        print(f"\n🔄 向後相容性測試:")
        print(f"   ├─ config.temperature: {multplane_config.temperature}")
        print(f"   ├─ config.tip.x_position: {multplane_config.tip.x_position}")
        print(f"   └─ config.semiconductor_regions: {len(multplane_config.semiconductor_regions)}")
        
    except Exception as e:
        print(f"❌ MultPlane 配置載入失敗: {e}")
        return False
    
    print("\n✨ 所有測試通過！新版 config_schema.py 成功支援階層化 YAML 結構")
    return True

def test_config_saving():
    """測試配置檔案儲存"""
    print("\n💾 測試配置檔案儲存功能...")
    reader = YamlConfigReader()
    
    try:
        # 載入配置
        multint_file = project_root / "Import_Files" / "MultInt_config.yaml"
        config = reader.load_config(multint_file)
        
        # 儲存到測試檔案
        output_file = project_root / "test_output_config.yaml"
        reader.save_config(config, output_file)
        
        print(f"✅ 配置檔案已儲存至: {output_file}")
        
        # 重新載入驗證
        reloaded_config = reader.load_config(output_file)
        print("✅ 重新載入測試通過")
        
        # 清理測試檔案
        output_file.unlink()
        print("🧹 測試檔案已清理")
        
        return True
        
    except Exception as e:
        print(f"❌ 配置儲存測試失敗: {e}")
        return False

if __name__ == "__main__":
    print("🚀 開始測試 Pysemitip 配置系統")
    
    success = True
    success &= test_config_loading()
    success &= test_config_saving()
    
    if success:
        print("\n🎉 所有測試通過！配置系統運作正常")
        sys.exit(0)
    else:
        print("\n💥 部分測試失敗，請檢查錯誤訊息")
        sys.exit(1)
