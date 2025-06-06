# Fortran-Python 程式碼對照表

## 概述
本文檔詳細比較 Fortran SEMITIP 與 Python 實現，確保所有功能都有正確對應。

## 主要檔案對照

### 1. 主程序結構

| Fortran 檔案 | Python 檔案 | 狀態 | 備註 |
|-------------|-------------|------|------|
| MultInt3-6.4.f | src/simulation/multint.py | ✅ 完成 | 主要模擬流程 |
| semitip3-6.1.f | src/physics/core/poisson.py | ⚠️ 需檢查 | 核心 Poisson 求解器 |

### 2. 核心子程序對照

#### A. 主要模擬流程 (MultInt3-6.4.f)

| Fortran 子程序/區塊 | 對應 Python 函數 | 檔案位置 | 狀態 | 關鍵差異 |
|-------------------|-----------------|----------|------|----------|
| **主程式初始化 (lines 1-100)** | `MultIntSimulation.__init__()` | multint.py:50-85 | ✅ | 無 |
| **讀取參數 (lines 52-65)** | `YamlConfigReader.load_config()` | filereader.py | ✅ | YAML vs Fortran格式 |
| **網格初始化 (lines 90-120)** | `Grid3D.__init__()` | grid.py:80-150 | ⚠️ | 需檢查NV參數 |
| **能量範圍計算 (lines 342-351)** | `_calculate_energy_range()` | multint.py:463-530 | ✅ | 邏輯匹配 |
| **電荷密度表格 (lines 355-380)** | `ChargeDensityTables.create_charge_tables()` | charge_density.py | ⚠️ | 需詳細檢查 |
| **多網格循環 (lines 400-450)** | `_solve_single_bias()` | multint.py:242-380 | ✅ | 邏輯匹配 |
| **電流計算 (lines 500-600)** | 未完整實現 | - | ❌ 缺失 | 需要實現 |

#### B. Poisson 求解器 (semitip3-6.1.f)

| Fortran 子程序 | 對應 Python 函數 | 檔案位置 | 狀態 | 關鍵參數檢查 |
|---------------|-----------------|----------|------|-------------|
| **SEMITIP3 主程序** | `PoissonSolver.solve()` | poisson.py:176-287 | ✅ | - |
| **真空電位更新 (lines 446-567)** | `_update_vacuum_fortran_style()` | poisson.py:289-350 | ⚠️ | 需檢查邊界條件 |
| **介面電位更新 (lines 578-622)** | `_update_interface_fortran_style()` | poisson.py:352-409 | ⚠️ | EEP常數需驗證 |
| **半導體電位更新 (lines 633-731)** | `_update_semiconductor_fortran_style()` | poisson.py:411-480 | ⚠️ | 單位轉換需檢查 |
| **初始電位設定 (lines 170-206)** | `_set_initial_potential()` | poisson.py:482-523 | ✅ | IINIT=1邏輯正確 |
| **Band bending計算 (PCENT)** | `get_band_bending()` | poisson.py:749-771 | ⚠️ | 權重計算需驗證 |

#### C. 金色分割搜索 (GSECT)

| Fortran | Python | 檔案位置 | 狀態 | 備註 |
|---------|--------|----------|------|------|
| **GSECT 子程序** | `golden_section_search()` | poisson.py:20-79 | ✅ | 算法匹配 |
| **GS常數 0.3819660** | `gs = 0.3819660` | poisson.py:45 | ✅ | 數值匹配 |

### 3. 關鍵常數對照

| 項目 | Fortran 值 | Python 值 | 檔案位置 | 狀態 |
|------|-----------|-----------|----------|------|
| **EEP常數** | `1.80943E-20` | `1.80943e-20` | poisson.py:225 | ✅ |
| **黃金分割常數** | `0.3819660` | `0.3819660` | poisson.py:45 | ✅ |
| **GaAs介電常數** | `12.9` | `12.9` | 多處 | ✅ |
| **溫度** | `300.0` | `300.0` | config | ✅ |

### 4. 數值方法對照

#### A. 有限差分格式

| 區域 | Fortran 實現 | Python 實現 | 狀態 | 備註 |
|------|-------------|-------------|------|------|
| **真空區拉普拉斯算子** | semitip3 lines 500-520 | `_update_vacuum_fortran_style()` | ⚠️ | 需檢查係數 |
| **半導體區泊松方程** | semitip3 lines 650-700 | `_update_semiconductor_fortran_style()` | ⚠️ | 需檢查係數 |
| **介面邊界條件** | semitip3 lines 585-590 | `_update_interface_fortran_style()` | ⚠️ | 3階差分需驗證 |

#### B. 非線性求解

| 項目 | Fortran | Python | 狀態 | 需檢查 |
|------|---------|--------|------|--------|
| **表面電位搜索範圍** | `DELSURF = BIAS/1E6` | `delta_surf = bias/1e3` | ❌ | 範圍不同 |
| **體電位搜索範圍** | `DELSEM = BIAS/1E6` | `delta_sem = bias/1e3` | ❌ | 範圍不同 |
| **收斂容忍度** | 需從代碼確認 | `1e-4` | ⚠️ | 需確認 |

## 發現的問題

### ✅ 已修正問題

1. **搜索範圍**: ✅ 已恢復到 Fortran 原始值 `BIAS/1E6`
2. **EEP 單位因子**: ✅ 已確認正確 (`eep * 1e7` for surface, `eep / eps` for bulk)
3. **收斂條件**: ✅ 已匹配 Fortran 邏輯 (`change1 < EP`, `change2 < 2*EP`)
4. **Fermi-Dirac 積分**: ✅ 已實現精確數值積分匹配 Fortran FJINT

### 🚨 根本問題 - Band Bending 量級差異

**現狀**:
- Fortran: ~0.07V (70mV)  
- Python: ~0.00003V (0.03mV)
- **差異: 約2300倍**

**已排除的原因**:
- ❌ 搜索範圍 (已修正)
- ❌ 物理常數 (6.815E21 正確)  
- ❌ Fermi積分精度 (已實現精確積分)
- ❌ 收斂邏輯 (已匹配Fortran)
- ❌ 網格參數 (完全匹配)

### 🔍 待調查的可能原因

1. **電荷密度表格插值方法**:
   - Fortran 可能使用不同的插值方式
   - 需要檢查 table lookup 實現

2. **初始條件差異**:
   - Fortran 初始 band bending 設定
   - 可能的前置計算步驟

3. **隱藏的單位轉換**:
   - 電位單位 (V vs eV)
   - 長度單位 (nm vs cm)
   - 密度單位轉換因子

### ⚠️ 中優先級問題

1. **Band bending 計算權重**:
   - Fortran: `(9.*VSINT(1,I,K)-VSINT(1,I+1,K))/8.`
   - Python: 實現相同但需要驗證索引

2. **收斂標準**:
   - Fortran 具體標準需要從代碼中提取
   - Python 可能過於嚴格

## 下一步行動計劃

### 立即修正 (今天)
1. ✅ 恢復搜索範圍到 Fortran 原始值
2. ✅ 檢查並修正 EEP 單位因子
3. ✅ 驗證介面電位計算

### 短期改進 (本週)
1. 詳細分析 Fortran 收斂標準
2. 實現電流計算部分
3. 驗證所有數值常數

### 驗證方法
1. 逐步對比 Fortran 輸出的每100次迭代
2. 檢查中間變量值
3. 確認物理單位一致性

---

**更新日期**: 2025-06-06
**狀態**: 初始分析完成，發現關鍵問題待修正