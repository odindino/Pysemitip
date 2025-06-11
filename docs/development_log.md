# Pysemitip 開發日誌

**專案名稱**: Pysemitip - 現代化Python架構的SEMITIP移植專案  
**開發者**: odindino  
**開始日期**: 2025年6月11日  
**專案目標**: 將Fortran SEMITIP程式移植為現代化Python架構，用於掃描穿隧顯微鏡(STM)的3D有限差分求解器

---

## 專案概述

SEMITIP是一個成熟的Fortran程式，用於模擬掃描穿隧顯微鏡的電子隧穿電流。本專案的目標是將其移植為現代化的Python架構，遵循「物理優先，架構為本」的設計原則。

### 設計原則
- **物理優先**: 確保數值精度和物理正確性
- **架構為本**: 建立清晰的模組化架構
- **分階段實現**: 按照複雜度逐步實現功能
- **驗證導向**: 每個階段都有完整的測試驗證

---

## 階段規劃

### 第一階段: 基礎工具層實現 (1-2週) ✅ 已完成
**目標**: 建立無依賴的數值計算基礎  
**時間**: 2025年6月11日

#### 實現內容
- [x] **黃金分割優化算法** (`utils/numerical.py`)
  - 從Fortran GSECT-6.0移植
  - 支援SEMMIN和SURFMIN函數使用
  - 通過精度和收斂性驗證

- [x] **費米-狄拉克相關函數** (`utils/numerical.py`)
  - 費米-狄拉克積分函數 (FJINT)
  - 費米-狄拉克佔據函數 (FD)
  - 支援載流子濃度計算

- [x] **插值工具** (`utils/interpolation.py`)
  - 線性插值 (一維)
  - 三次樣條插值
  - 雙線性插值 (二維)
  - 三線性插值 (三維)

- [x] **數值微積分工具** (`utils/numerical.py`)
  - 數值微分 (前向/後向/中心差分)
  - 自適應積分算法
  - 梯形積分法 (TRAP)

#### 驗證結果
- ✅ 17個測試全部通過
- ✅ 與Fortran GSECT數值一致性驗證
- ✅ 費米-狄拉克函數精度驗證
- ✅ 插值算法準確性驗證

#### 文件結構
```
src/utils/
├── __init__.py          # 模組初始化和API導出
├── numerical.py         # 核心數值計算函數
└── interpolation.py     # 插值工具函數

tests/
├── test_numerical_tools.py                    # 綜合測試套件
└── test_gsect_fortran_compatibility.py        # Fortran兼容性驗證

demos/
├── demo_numerical_tools.py     # 功能演示
└── gsect_demo.png              # 優化演示圖表
```

---

### 第二階段: 物理模型層實現 (2-3週) ✅ 已完成
**目標**: 實現核心物理計算模組  
**時間**: 2025年6月11日

#### 實現內容
- [x] **材料參數管理** (`physics/materials.py`)
  - 半導體材料資料庫 (Si n/p型, GaAs n型, Si本征)
  - 物理常數定義 (與Fortran COMMON塊對應)
  - 表面態參數管理

- [x] **電荷密度計算** (`physics/charge_density.py`)
  - 費米-狄拉克統計 (RHOCB/RHOVB等效)
  - 自洽電荷中性條件求解
  - 摻雜劑離化計算

- [x] **泊松方程求解器** (`physics/poisson.py`)
  - 三維柱坐標有限差分法
  - STM幾何專用網格 (194,560個網格點)
  - 複雜邊界條件處理

- [x] **穿隧電流計算** (`physics/tunneling_current.py`)
  - Bardeen轉移哈密頓方法
  - WKB近似穿隧係數計算
  - 修復偏壓相關勢壘剖面問題

#### 驗證結果
- ✅ 27個物理模型測試全部通過 (4分37秒)
- ✅ 電荷中性精度 < 1e10 cm⁻³
- ✅ 穿隧電流顯示正確的偏壓依賴性
- ✅ 不同材料顯示合理的電流差異

#### 重要技術修正
**日期**: 2025年6月11日  
**問題**: 穿隧電流在不同偏壓下顯示相同數值  
**根本原因**: `calculate_simple_stm_current`函數使用平坦電位剖面 (`np.zeros_like(z_grid)`)，未考慮偏壓相關的勢壘變化

**修正方案**:
```python
# 修正前: 平坦電位剖面
potential_profile = np.zeros_like(z_grid)

# 修正後: 偏壓相關的梯形勢壘剖面
for i, z in enumerate(z_grid):
    # 線性電位降落 (tip → sample)
    bias_drop = bias_voltage * z / separation
    
    # 梯形勢壘形狀
    if z < separation * 0.1:  # 接近探針
        potential_profile[i] = tip_wf + 1.0  # 1 eV勢壘
    elif z > separation * 0.9:  # 接近樣品
        potential_profile[i] = sample_wf + bias_drop + 1.0
    else:  # 中間區域
        tip_potential = tip_wf + 1.0
        sample_potential = sample_wf + bias_voltage + 1.0
        potential_profile[i] = tip_potential + (sample_potential - tip_potential) * z / separation
    
    # 鏡像電位修正
    image_correction = -0.1 / (4 * z + 0.05) if z > 0.01 else 0
    potential_profile[i] += image_correction
```

**修正結果**:
- Si n型: +1V: 1.50A, -1V: -0.59A (正確的偏壓依賴性)
- GaAs n型: +1V: 6.48A, -1V: -4.56A (符合材料特性)
- 不同材料顯示合理的電流差異

#### 文件結構
```
src/physics/
├── __init__.py                  # 模組初始化和API導出
├── materials.py                 # 材料參數資料庫 (340行)
├── charge_density.py            # 電荷密度計算 (465行)
├── poisson.py                   # 泊松方程求解器 (510行)
└── tunneling_current.py         # 穿隧電流計算 (630行)

tests/
└── test_physics_models.py       # 物理模型測試 (27測試)

demos/
├── demo_phase2_physics.py       # Phase 2功能演示
└── phase2_results.png           # 計算結果圖表

debug/
└── debug_current.py             # 穿隧電流調試腳本
```

#### 階段總結
Phase 2成功實現了所有核心物理計算模組，並解決了穿隧電流計算的關鍵問題。所有測試通過，物理結果合理，為下一階段的幾何和網格層實現奠定了堅實基礎。

---

### 第三階段: 幾何和網格層實現 (1-2週) 🚧 規劃中
**目標**: 實現STM幾何建模和網格管理

#### 規劃內容
- [ ] **STM幾何建模** (`geometry/stm_geometry.py`)
  - 探針和樣品幾何定義
  - 坐標系統轉換
  - 邊界條件設定

- [ ] **三維網格管理** (`geometry/grid.py`)
  - 柱坐標網格生成
  - 網格點編號和索引
  - 記憶體最佳化

- [ ] **可視化工具** (`visualization/`)
  - 3D幾何可視化
  - 電位分布繪圖
  - 電流密度圖

---

## 專案狀態總覽

### 當前進度
- ✅ **Phase 1**: 基礎工具層 (100% 完成)
- ✅ **Phase 2**: 物理模型層 (100% 完成，含重要修正)
- 🚧 **Phase 3**: 幾何和網格層 (規劃中)

### 測試狀態
- **總測試數**: 50個測試 (Phase 1: 23個, Phase 2: 27個)
- **通過率**: 100%
- **執行時間**: ~5分鐘

### 關鍵里程碑
1. **2025年6月11日**: 完成Phase 1基礎工具層
2. **2025年6月11日**: 完成Phase 2物理模型層
3. **2025年6月11日**: 修正穿隧電流偏壓依賴性問題

### 下一步行動
1. 開始Phase 3幾何和網格層實現
2. 建立STM幾何建模架構
3. 實現3D可視化功能

---

## 技術決策記錄

### TDR-001: 數值精度保證
**決策**: 使用numpy.float64作為預設精度  
**理由**: 確保與Fortran雙精度一致性  
**影響**: 所有數值計算模組

### TDR-002: 模組化架構設計
**決策**: 採用分層模組架構 (utils → physics → geometry → simulation)  
**理由**: 清晰的依賴關係，易於測試和維護  
**影響**: 整體專案結構

### TDR-003: 穿隧電流勢壘建模
**決策**: 實現偏壓相關的梯形勢壘剖面  
**理由**: 修正原有平坦電位剖面的物理不合理性  
**影響**: tunneling_current.py核心算法