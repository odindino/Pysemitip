# Pysemitip

## 專案概述

Pysemitip 是一個基於 Python 的掃描隧道顯微鏡 (STM) 模擬工具，源自卡內基美隆大學 (CMU) 開發的 SEMITIP 程式。本專案旨在將原始的 Fortran 程式重新實作為 Python 版本，提供更好的可維護性、擴展性和易用性。

### 主要目標

- **現代化重寫**：將 CMU 的 Fortran SEMITIP 程式轉換為 Python 實作
- **提高可維護性**：使用現代程式設計實踐和清晰的程式架構
- **增強易用性**：提供友善的配置格式和更直觀的使用介面
- **促進協作**：降低學習門檻，讓更多研究人員能參與開發和使用
- **模型改進**：為未來的實驗數據適配和模型優化提供基礎

## 專案結構

```
Pysemitip/
├── README.md                    # 專案說明文件
├── environment.yml              # Conda 環境配置
├── requirements.txt             # Python 依賴套件
├── config_schema.py             # 配置資料結構定義
├── filereader.py               # YAML 配置檔讀取器
├── fileconverter.py            # Fort.9 格式轉換器 (開發中)
├── test/                       # 測試腳本
│   ├── yaml_config_test.py     # YAML 配置測試
│   └── simple_yaml_test.py     # 簡單 YAML 讀取測試
├── converted_configs/          # 轉換後的配置檔案
│   ├── MultInt_config.yaml     # MultInt 模擬配置
│   └── MultPlane_config.yaml   # MultPlane 模擬配置
├── examples/                   # 範例配置檔案
│   ├── example_multint_fixed.yaml
│   └── example_multplane_fixed.yaml
└── Fortran/                    # 原始 Fortran 參考檔案
    ├── MultInt/
    │   └── fort_MultInt.9
    └── MultPlane/
        └── fort_MultPlane.9
```

## 核心功能

### 1. 配置管理系統

- **統一的資料結構**：使用 Python dataclass 定義所有配置參數
- **YAML 格式配置**：人類可讀的配置檔案格式，取代難以維護的 fort.9 格式
- **參數驗證**：自動驗證配置參數的合理性和一致性
- **巢狀結構支援**：支援複雜的配置結構，如多個半導體區域和表面態

### 2. 支援的模擬類型

#### MultInt 模擬
- 適用於多重積分方法的 STM 模擬
- 支援平行波向量和能量點設定
- 可配置的半導體深度分數

#### MultPlane 模擬
- 適用於多平面方法的 STM 模擬
- 支援真空寬度和間距設定
- 可配置的最大能量設定
- 支援所有能帶計算選項

### 3. 配置參數類別

#### 探針配置 (TipConfig)
- 探針幾何參數：錐度、半徑、突出半徑
- 位置參數：x, y 座標
- 電性參數：接觸電位、費米能級

#### 半導體配置 (SemiconductorRegion)
- 載子濃度：施子/受子濃度
- 能帶參數：能隙、價帶偏移
- 有效質量：導帶、重電洞、輕電洞、分裂軌道電洞
- 結合能：施子/受子結合能

#### 表面態配置 (SurfaceRegion)
- 雙分佈模型：支援兩個表面態分佈
- 每個分佈包含：密度、中性能級、半高寬、重心能量

#### 計算網格配置 (GridConfig)
- 網格點數：徑向、真空、半導體、角度方向
- 鏡像平面設定
- 初始網格大小

#### 計算參數配置 (ComputationConfig)
- 縮放步驟數
- 最大迭代次數陣列
- 收斂參數陣列
- 電荷密度表格大小

#### 電壓掃描配置 (VoltageScanConfig)
- 掃描範圍：起始/結束電壓
- 掃描點數
- 調變電壓
- 斜坡參數：正/負斜坡

## 安裝與環境設置

### 1. 建立 Conda 環境

```bash
# 使用 environment.yml 建立環境
conda env create -f environment.yml

# 啟動環境
conda activate pysemitip
```

### 2. 手動安裝 (可選)

```bash
# 建立新環境
conda create -n pysemitip python=3.12

# 啟動環境
conda activate pysemitip

# 安裝依賴套件
conda install pyyaml numpy scipy matplotlib pandas pytest
```

## 使用方法

### 1. 讀取 YAML 配置檔

```python
from filereader import YamlConfigReader

# 建立讀取器
reader = YamlConfigReader()

# 讀取配置檔案
config = reader.load_config("MultInt_config.yaml")

# 檢視配置資訊
print(f"模擬類型: {config.simulation_type}")
print(f"探針半徑: {config.tip.radius} nm")
print(f"半導體區域數: {len(config.semiconductor_regions)}")
```

### 2. 建立新的配置

```python
from config_schema import SemitipConfig, TipConfig, SemiconductorRegion

# 建立基本配置
config = SemitipConfig(
    simulation_type="MultInt",
    temperature=300.0,
    dielectric_constant=12.9
)

# 設定探針參數
config.tip = TipConfig(
    radius=1.0,
    separation=1.0,
    contact_potential=0.0
)

# 新增半導體區域
region = SemiconductorRegion(
    donor_concentration=1e18,
    band_gap=1.42
)
config.semiconductor_regions.append(region)

# 儲存配置
reader = YamlConfigReader()
reader.save_config(config, "my_config.yaml")
```

### 3. 驗證配置

```python
# 配置會自動驗證
try:
    config.validate()
    print("配置驗證通過")
except ValueError as e:
    print(f"配置錯誤: {e}")
```

## 配置檔案範例

### MultInt 配置範例

```yaml
version: "1.0"
simulation_type: "MultInt"
temperature: 300.0
dielectric_constant: 12.9

tip:
  shank_slope: 1.0
  separation: 1.0
  radius: 1.0
  protrusion_radius: 0.0
  contact_potential: 0.0
  x_position: 0.0
  y_position: 0.0
  fermi_energy: 8.0

semiconductor_regions:
  - id: 1
    donor_concentration: 1.0e18
    acceptor_concentration: 0.0
    band_gap: 1.42
    effective_mass:
      conduction_band: 0.0635
      heavy_hole: 0.643
      light_hole: 0.081
      split_off_hole: 0.172
```

## 開發狀態

### ✅ 已完成功能
- 配置資料結構定義 (`config_schema.py`)
- YAML 配置檔讀取器 (`filereader.py`)
- 基本的配置驗證
- 測試腳本和範例配置檔

### 🚧 開發中功能
- Fort.9 格式轉換器 (`fileconverter.py`)
- STM 模擬核心引擎
- 結果分析和視覺化工具

### 📋 計劃功能
- 圖形化使用者介面 (GUI)
- 參數掃描和最佳化
- 實驗數據比較工具
- 進階模型和物理效應

## 測試

```bash
# 執行 YAML 配置測試
python test/yaml_config_test.py

# 執行簡單 YAML 讀取測試
python test/simple_yaml_test.py
```

## 從 SEMITIP Fortran 移植

本專案保持與原始 SEMITIP 程式的物理模型一致性，同時提供以下改進：

1. **配置格式**：從難以維護的 fort.9 格式改為 YAML 格式
2. **程式結構**：模組化設計，易於擴展和修改
3. **錯誤處理**：完整的錯誤檢查和使用者友善的錯誤訊息
4. **文件化**：完整的程式碼註釋和使用說明

## 貢獻指南

歡迎對本專案的貢獻！請遵循以下指南：

1. **程式碼風格**：使用 Black 進行程式碼格式化
2. **型別註解**：使用 mypy 進行型別檢查
3. **測試**：為新功能添加相應的測試
4. **文件**：更新相關的文件和註釋

## 依賴套件

### 核心依賴
- Python >= 3.9
- PyYAML >= 6.0
- NumPy >= 1.21.0
- SciPy >= 1.7.0

### 開發依賴
- pytest >= 7.0.0
- black >= 22.0.0
- mypy >= 0.991
- flake8 >= 5.0.0

### 可選依賴
- matplotlib >= 3.5.0 (繪圖)
- pandas >= 1.3.0 (數據分析)
- jupyter >= 1.0.0 (互動式分析)

## 授權

本專案採用 MIT 授權條款。詳見 LICENSE 檔案。

## 致謝

- 感謝卡內基美隆大學開發原始的 SEMITIP 程式
- 感謝所有對本專案做出貢獻的研究人員和開發者
- This program is the python version of Semitip originally developed by Prof. Feenstra group in CMU(https://www.andrew.cmu.edu/user/feenstra/semitip_v6/).

## 聯絡資訊

如有問題或建議，請透過 GitHub Issues 聯絡我們。

---

**注意**：本專案目前處於早期開發階段，API 可能會有變更。建議在正式使用前查看最新的文件和範例。

---

**Joke:**
有一天，小美走在路上踢到一個小東西，蹲下來看發現是一根針，她就向周圍問「這是誰的針呀？」  
小明聽到就及著急地跑過來說「拍謝 My Tip」「謝謝你，我找了好久，終於出現了」  
沒錯，Pysemitip問世了，幫助你量測STM後，分析數據和以理論模擬輔助實驗解釋的好幫手！  
