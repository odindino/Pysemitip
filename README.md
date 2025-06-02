# Pysemitip - Modern Python SEMITIP Implementation

🎯 **狀態**: 系統已完全準備就緒，配置載入和所有屬性存取正常運作

A modern Python implementation of the SEMITIP Fortran package for STM (Scanning Tunneling Microscopy) simulations on semiconductor surfaces.

## ✅ 目前狀態 (2025-June-01)

### 已完成功能
- ✅ **配置系統**: YAML 載入、驗證、型別轉換完全正常
- ✅ **向後相容性**: 支援平面和階層 YAML 結構
- ✅ **模組匯入**: 所有相對匯入路徑已修正
- ✅ **屬性存取**: MultInt 模擬所需的所有屬性可用
- ✅ **專案結構**: 檔案整理完成，結構清晰
- ✅ **測試覆蓋**: 全面的測試確保穩定性

### 測試狀態
```bash
# 所有測試通過
python tests/test_system_readiness.py  # 系統準備狀況檢查
python tests/test_simulation_properties.py  # 模擬屬性測試
python tests/test_voltage_scan.py  # 電壓掃描測試
```

## 🚀 快速開始

### 1. 安裝依賴
```bash
pip install numpy scipy matplotlib pyyaml
```

### 2. 執行模擬
```bash
python run_multint.py data/input/examples/quick_test.yaml --plot
```

### 3. 檢查結果
結果將保存在 `data/output/results/` 目錄中，包含：
- `multint_results.pkl` - 主要結果檔案
- `output_MultInt.log` - 詳細計算日誌
- 圖表檔案 (如果使用 `--plot` 選項)

## 📁 專案結構

```
Pysemitip/
├── run_multint.py          # 主執行檔案
├── src/                    # 原始碼
│   ├── core/              # 配置和檔案處理
│   ├── physics/           # 物理模擬模組
│   ├── simulation/        # 主要模擬程式
│   └── utils/             # 工具函數
├── data/                  # 資料檔案
│   ├── input/             # 輸入配置
│   └── output/            # 輸出結果
├── tests/                 # 測試檔案
├── docs/                  # 文件檔案
├── requirements.txt       # Python 依賴
└── environment.yml        # Conda 環境檔案
```

## ⚙️ 配置檔案

使用 YAML 格式的配置檔案：

```yaml
version: "1.0"
simulation_type: "MultInt"

# 環境參數
environment:
  temperature: 300.0
  dielectric_constant: 12.9

# STM 探針參數
tip:
  radius: 1.0
  separation: 1.0
  work_function: 5.3

# 電壓掃描
voltage_scan:
  start: -2.0
  end: 2.0
  points: 41

# 半導體區域
semiconductor_regions:
  - id: 1
    donor_concentration: 1.0e18
    band_gap: 1.42
    affinity: 4.07
```

完整範例請參見 `data/input/examples/` 目錄。

## 🧪 測試

### 系統檢查
```bash
python tests/test_system_readiness.py
```

### 特定功能測試
```bash
python tests/test_basic.py              # 基本匯入
python tests/test_simulation_properties.py  # 模擬屬性
python tests/test_voltage_scan.py       # 電壓掃描
```

## 📊 與 Fortran 版本比較

Python 實作的結果應該與原始 Fortran 版本的輸出 (`fort_MultInt.16`) 一致。比較方法：

1. 執行 Python 模擬
2. 檢查 `output_MultInt.log` 檔案
3. 將關鍵數值與 Fortran 輸出比較

## 📚 文件

- **專案狀態**: `docs/PROJECT_STATUS.md`
- **修正紀錄**: `docs/FIXES_APPLIED.md`
- **快速指南**: `docs/QUICK_START.md`
- **環境設定**: `docs/ENVIRONMENT_SETUP.md`

## 🛠️ 開發

### 代碼風格
- 遵循 PEP 8 標準
- 使用 Black 格式化工具
- 類型提示 (Type Hints)

### 測試覆蓋率
- 單元測試覆蓋所有核心功能
- 整合測試確保端到端功能
- 配置測試驗證 YAML 處理

## 🎯 核心特色

1. **現代化架構**: 模組化設計，清晰的關注點分離
2. **靈活配置**: 人類可讀的 YAML 配置檔案
3. **強大的驗證**: 自動型別轉換和配置驗證
4. **向後相容**: 支援現有的 YAML 檔案結構
5. **全面測試**: 高覆蓋率的測試確保穩定性

## 📈 技術亮點

- **無 NumPy 依賴的配置系統**: 配置載入不需要科學計算套件
- **智能型別轉換**: 自動處理 YAML 中的科學記號字串
- **階層到平面映射**: 通過屬性提供向後相容的存取方式
- **模組化物理引擎**: 可擴展的物理模擬框架

---

**準備好開始模擬了！** 🚀