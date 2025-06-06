# Python 程式碼結構分析

## 📁 Python 模組總覽

### 主要模組結構
```
src/
├── __init__.py
├── simulation/
│   ├── __init__.py
│   └── multint.py              # 主要模擬控制器
├── physics/
│   ├── __init__.py
│   ├── core/
│   │   ├── __init__.py
│   │   ├── poisson.py          # Poisson 方程求解器
│   │   ├── charge_density.py   # 電荷密度計算
│   │   ├── potential.py        # 電位相關計算
│   │   └── schrodinger.py      # Schrödinger 方程求解
│   ├── materials/
│   │   ├── __init__.py
│   │   ├── semiconductor.py    # 半導體材料
│   │   ├── tip.py              # STM 探針
│   │   └── surface_states.py   # 表面態
│   ├── solvers/
│   │   ├── __init__.py
│   │   └── grid.py             # 網格求解器
│   └── visualization/
│       ├── __init__.py
│       ├── plotter.py          # 結果繪圖
│       └── contour.py          # 等高線圖
├── core/
│   ├── __init__.py
│   ├── config.py               # 配置管理
│   └── file_handler.py         # 檔案處理
└── utils/
    ├── __init__.py
    └── constants.py            # 物理常數
```

## 🔍 各模組詳細分析

### 1. simulation/multint.py (主控制器)

#### 基本資訊
- **對應 Fortran**: MultInt3-6.4.f
- **行數**: ~500 行
- **設計模式**: 物件導向，單一職責

#### 類別結構
```python
class MultIntSimulation:
    ├── __init__()              # 初始化
    ├── run()                   # 主執行函數
    ├── _setup_simulation()     # 模擬設定
    ├── _voltage_scan_loop()    # 電壓掃描迴圈
    ├── _compute_charge_tables() # 電荷密度表格
    ├── _solve_poisson()        # Poisson 求解
    ├── _calculate_current()    # 電流計算
    └── _save_results()         # 結果儲存
```

#### 關鍵方法對應
| Python 方法 | Fortran 對應 | 狀態 |
|-------------|-------------|------|
| `_setup_simulation()` | 初始化部分 | ✅ |
| `_voltage_scan_loop()` | 主迴圈 | ✅ |
| `_solve_poisson()` | CALL SEMITIP3 | ⚠️ |
| `_calculate_current()` | CALL INTCURR | ❌ |

### 2. physics/core/poisson.py (Poisson 求解器)

#### 基本資訊
- **對應 Fortran**: semitip3-6.1.f
- **行數**: ~400 行
- **數值方法**: SOR 迭代法

#### 類別結構
```python
class PoissonSolver:
    ├── __init__()                    # 初始化網格
    ├── solve()                       # 主求解函數
    ├── _setup_grids()               # 網格設定
    ├── _set_boundary_conditions()   # 邊界條件
    ├── _sor_iteration()             # SOR 迭代
    ├── _check_convergence()         # 收斂檢查
    ├── _multigrid_solve()           # 多網格求解
    └── _update_potential()          # 電位更新
```

#### 目前問題
- **收斂過快**: 迭代次數固定為 200，應該是動態收斂
- **非線性求解**: 缺少類似 GSECT 的非線性處理
- **邊界條件**: 需要更精確的邊界處理

### 3. physics/core/charge_density.py (電荷密度)

#### 基本資訊
- **對應 Fortran**: semirhomult-6.0.f, surfrhomult-6.2.f
- **行數**: ~300 行
- **功能**: 電荷密度插值表建立

#### 類別結構
```python
class ChargeDensityCalculator:
    ├── __init__()                      # 初始化
    ├── create_bulk_charge_table()      # 體電荷密度表
    ├── create_surface_charge_table()   # 表面電荷密度表
    ├── get_bulk_charge_density()       # 體電荷密度插值
    ├── get_surface_charge_density()    # 表面電荷密度插值
    └── _setup_energy_range()           # 能量範圍設定
```

#### 實現狀態
- ✅ **插值表建立**: 正確實現
- ✅ **能量範圍**: 與 Fortran 一致
- ✅ **插值函數**: 工作正常

### 4. physics/materials/semiconductor.py (半導體)

#### 基本資訊
- **對應 Fortran**: COMMON/SEMI/ 區塊
- **行數**: ~200 行
- **功能**: 半導體材料屬性

#### 類別結構
```python
class SemiconductorRegion:
    ├── __init__()                  # 初始化參數
    ├── calculate_fermi_level()     # 費米能階計算
    ├── get_carrier_density()       # 載流子濃度
    ├── get_charge_density()        # 電荷密度
    └── get_band_properties()       # 能帶屬性
```

### 5. physics/materials/tip.py (STM 探針)

#### 基本資訊
- **對應 Fortran**: 探針相關變數
- **行數**: ~100 行
- **關鍵修復**: Tip potential 動態更新

#### 類別結構
```python
class Tip:
    ├── __init__()                  # 初始化
    ├── tip_potential (property)    # 探針電位 (動態)
    ├── update_bias()              # 更新偏壓
    └── get_geometry()             # 幾何參數
```

#### 重要修復
```python
@property
def tip_potential(self):
    """Dynamic tip potential calculation"""
    return self.bias_voltage + self.contact_potential
```

### 6. physics/materials/surface_states.py (表面態)

#### 基本資訊
- **對應 Fortran**: surfrhomult-6.2.f 部分
- **行數**: ~150 行
- **功能**: 表面態分佈計算

### 7. physics/core/schrodinger.py (Schrödinger 求解)

#### 基本資訊
- **對應 Fortran**: intcurr-6.2.f
- **行數**: ~200 行
- **狀態**: ⚠️ 部分實現

#### 目前問題
- **電流計算**: 出現 NaN 值
- **局域化態**: 搜尋不完整
- **積分方法**: 需要完善

## 📊 架構比較

### Fortran vs Python 設計模式

| 特性 | Fortran | Python |
|------|---------|--------|
| **結構** | 程序式 | 物件導向 |
| **資料共享** | COMMON 區塊 | 類別屬性 |
| **模組化** | 子程序 | 類別/方法 |
| **陣列索引** | 1-based | 0-based |
| **記憶體管理** | 靜態陣列 | 動態陣列 |

### 模組對應關係

```
Fortran 檔案              Python 模組
├── MultInt3-6.4.f    →   simulation/multint.py
├── semitip3-6.1.f    →   physics/core/poisson.py
├── semirhomult-6.0.f →   physics/core/charge_density.py
├── surfrhomult-6.2.f →   physics/materials/surface_states.py
├── intcurr-6.2.f     →   physics/core/schrodinger.py
├── potcut3-6.0.f     →   physics/core/potential.py
├── potexpand-6.1.f   →   ❌ 未實現
├── gsect-6.0.f       →   ⚠️ 分散在多個模組
└── contr3-6.0.f      →   ❌ 未實現
```

## 🎯 實現狀態分析

### ✅ 已完成模組
1. **simulation/multint.py**: 主控邏輯 (80%)
2. **physics/core/charge_density.py**: 電荷密度 (95%)
3. **physics/materials/semiconductor.py**: 半導體 (90%)
4. **physics/materials/tip.py**: 探針 (85%)

### ⚠️ 部分完成模組
1. **physics/core/poisson.py**: Poisson 求解 (60%)
   - 缺少非線性求解
   - 收斂邏輯需改進
2. **physics/core/schrodinger.py**: 電流計算 (30%)
   - 出現 NaN 問題
   - 積分方法不完整

### ❌ 未實現功能
1. **potexpand-6.1.f**: 電位多極展開
2. **contr3-6.0.f**: 輔助控制函數
3. **gsect-6.0.f**: 非線性求根方法

## 🔧 設計模式分析

### 優點
1. **模組化設計**: 清晰的關注點分離
2. **物件導向**: 容易維護和擴展
3. **型別提示**: 增強代碼可讀性
4. **配置管理**: 靈活的 YAML 配置系統

### 需要改進
1. **數值精度**: 確保與 Fortran 一致
2. **效能最佳化**: 關鍵迴圈的最佳化
3. **錯誤處理**: 更完善的異常處理
4. **測試覆蓋**: 增加單元測試

## 📝 關鍵發現

### 成功翻譯的特點
1. **配置系統**: YAML vs Fortran 輸入檔案
2. **資料結構**: 物件屬性 vs COMMON 區塊
3. **模組邊界**: 清晰的介面定義

### 翻譯挑戰
1. **數值方法**: 複雜算法的精確移植
2. **記憶體佈局**: 多維陣列的索引轉換
3. **控制流程**: GOTO 語句的結構化重寫

---

**分析完成日期**: 2025-06-06  
**下一步**: 建立詳細的函數對應關係
