#!/usr/bin/env python3
"""
測試配置驗證功能
"""

import os
import sys
from pathlib import Path

# 確保能夠導入專案模組
current_dir = Path(os.getcwd())
project_root = current_dir.parent if 'tests' in current_dir.name else current_dir
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.core.filereader import YamlConfigReader
from src.core.config_schema import SemitipConfig

# 創建一個正常配置
reader = YamlConfigReader()
valid_config = reader._yaml_to_config({
    'simulation_type': 'MultInt', 
    'environment': {'temperature': 300}
})

# 創建一個錯誤的配置 (溫度為負)
invalid_config = reader._yaml_to_config({
    'simulation_type': 'MultInt', 
    'environment': {'temperature': -100}
})

# 測試正常配置
print("\n=== 測試正常配置 ===")
try:
    valid_config.validate()
    print("✓ 正常配置驗證通過")
except ValueError as e:
    print(f"✗ 正常配置驗證失敗 (不應發生): {e}")

# 測試錯誤配置
print("\n=== 測試錯誤配置 ===")
try:
    invalid_config.validate()
    print("✗ 錯誤配置驗證通過 (不應發生)")
except ValueError as e:
    print(f"✓ 錯誤配置驗證失敗 (預期結果): {e}")
