#!/usr/bin/env python3
"""
System readiness test - comprehensive check before full simulation.
"""

import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def test_system_readiness():
    """Comprehensive test to verify system is ready for simulation."""
    print("=== Pysemitip 系統準備狀況檢查 ===\n")
    
    all_passed = True
    
    # 1. 基本匯入測試
    try:
        from src.core.filereader import YamlConfigReader, load_yaml_config
        from src.core.config_schema import SemitipConfig
        print("✓ 核心模組匯入成功")
    except Exception as e:
        print(f"✗ 核心模組匯入失敗: {e}")
        all_passed = False
    
    # 2. 配置檔案測試
    try:
        config_path = project_root / "data/input/examples/quick_test.yaml"
        config = load_yaml_config(str(config_path))
        print("✓ YAML 配置載入成功")
    except Exception as e:
        print(f"✗ YAML 配置載入失敗: {e}")
        all_passed = False
        return False
    
    # 3. 關鍵屬性存取測試
    critical_properties = [
        ('tip.radius', lambda: config.tip.radius),
        ('tip.slope', lambda: config.tip.slope),
        ('contact_potential', lambda: config.contact_potential),
        ('semiconductor_regions', lambda: len(config.semiconductor_regions)),
        ('surface_regions', lambda: len(config.surface_regions)),
        ('voltage_scan.voltages', lambda: len(config.voltage_scan.voltages)),
        ('mirror_symmetry', lambda: config.mirror_symmetry),
        ('temperature', lambda: config.temperature),
        ('charge_table_points', lambda: config.charge_table_points),
    ]
    
    for prop_name, prop_getter in critical_properties:
        try:
            value = prop_getter()
            print(f"✓ {prop_name}: {value}")
        except Exception as e:
            print(f"✗ {prop_name}: 錯誤 - {e}")
            all_passed = False
    
    # 4. 數值計算相關模組測試（如果有 numpy）
    try:
        import numpy as np
        print("✓ NumPy 可用")
        
        try:
            import scipy
            print("✓ SciPy 可用")
        except ImportError:
            print("! SciPy 未安裝 - 某些功能可能受限")
            
        try:
            import matplotlib
            print("✓ Matplotlib 可用")
        except ImportError:
            print("! Matplotlib 未安裝 - 無法生成圖表")
            
    except ImportError:
        print("! NumPy 未安裝 - 數值計算功能不可用")
        print("  請執行: pip install numpy scipy matplotlib pyyaml")
    
    # 5. 資料夾結構檢查
    required_dirs = [
        'data/input/examples',
        'data/output/results',
        'src/core',
        'src/physics',
        'src/simulation',
        'tests',
        'docs'
    ]
    
    for dir_path in required_dirs:
        full_path = project_root / dir_path
        if full_path.exists():
            print(f"✓ 目錄存在: {dir_path}")
        else:
            print(f"✗ 目錄缺失: {dir_path}")
            all_passed = False
    
    # 6. 主執行檔案檢查
    main_files = ['run_multint.py']
    for file_name in main_files:
        file_path = project_root / file_name
        if file_path.exists():
            print(f"✓ 主執行檔案存在: {file_name}")
        else:
            print(f"✗ 主執行檔案缺失: {file_name}")
            all_passed = False
    
    # 7. 總結
    print(f"\n=== 檢查結果 ===")
    if all_passed:
        print("🎉 系統已準備就緒！")
        print("\n下一步:")
        print("1. 安裝科學計算套件: pip install numpy scipy matplotlib pyyaml")
        print("2. 執行模擬: python run_multint.py data/input/examples/quick_test.yaml --plot")
        return True
    else:
        print("❌ 系統尚未完全準備好")
        print("請解決上述問題後再次執行")
        return False

if __name__ == "__main__":
    success = test_system_readiness()
    sys.exit(0 if success else 1)