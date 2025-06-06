# Python-Fortran 數值差異除錯報告

## 📊 問題識別

**測試日期**: 2025-01-27  
**測試目的**: 分析 Python 版本與 Fortran 版本的關鍵數值差異  
**測試階段**: Phase 1 - 初始問題診斷

## 🔍 主要差異分析

### 1. Poisson 求解器收斂行為差異

#### 📈 Band Bending 數值比較
| 項目 | Fortran 結果 | Python 結果 | 差異倍數 |
|------|-------------|-------------|----------|
| Solution #1 | ~0.07 V | 0.00001319 V | **5,000倍** |
| Solution #2 | ~0.07 V | 0.00003296 V | **2,100倍** |
| Solution #3 | ~0.07 V | 0.00005273 V | **1,300倍** |

**❌ 關鍵問題**: Python 版本的 band bending 值比 Fortran 版本小了 3-4 個數量級

#### 🔄 迭代行為差異
| 項目 | Fortran | Python |
|------|---------|--------|
| 迭代次數 | 動態 (100-3500) | 固定 200 次 |
| 收斂檢查 | ✅ 每 100 次輸出 | ❌ 未實現 |
| 中途輸出 | 詳細 ITER,Pot0 | 無 |
| 最終收斂 | 達到穩定值 | 固定次數停止 |

### 2. 電流計算差異

#### 📊 電流結果比較
| 項目 | Fortran | Python | 狀態 |
|------|---------|--------|------|
| VB 電流 | -1.81937781E-14 A | **NaN** | ❌ 失敗 |
| CB 電流 | 0.0 A | 0.0 A | ✅ 一致 |
| 局域化態數量 | 0 (各能帶) | 198 (VB) | ❌ 差異巨大 |

**❌ 關鍵問題**: Python 版本電流計算返回 NaN，局域化態數量計算錯誤

### 3. 載流子密度驗證

#### ✅ 成功項目
| 項目 | Fortran | Python | 相對誤差 |
|------|---------|--------|----------|
| CB 載流子密度 | 2.94679424E+17 | 2.950922e+17 | 0.14% |
| VB 載流子密度 | 57.446033 | 57.44 | 0.01% |
| Fermi 能階 | 1.4186435 | 1.418687 | 0.003% |

**✅ 結論**: 載流子密度計算已正確實現

## 🎯 根本原因分析

### 問題 1: Poisson 求解器非線性收斂失敗

#### Fortran 正確行為
```fortran
! 詳細迭代過程，動態收斂
ITER,Pot0 =         100 -8.27837288E-02
ITER,Pot0 =         200 -8.84749368E-02
...逐步收斂到穩定值...
ITER,Pot0 =        3500  7.00571761E-02
BAND BENDING AT MIDPOINT = 6.98396191E-02
```

#### Python 問題行為
```python
# 固定 200 次迭代，無動態收斂
for iteration in range(200):
    # SOR 迭代但未整合 Golden Section Search
    # 無中途收斂檢查和輸出
BAND BENDING AT MIDPOINT = 0.00005273  # 錯誤的小值
```

**根本原因**: 
1. 未實現 Golden Section Search 與 SOR 的整合
2. 固定迭代次數而非動態收斂
3. 缺少非線性求解的正確邏輯

### 問題 2: 電流計算積分算法錯誤

#### Fortran 正確行為
```fortran
number of VB light-hole localized states = 0
number of VB heavy-hole localized states = 0
number of CB localized states = 0
valence band current ext,loc = -1.81937781E-14 0.0000000
```

#### Python 問題行為
```python
number of VB light-hole localized states = 198  # 錯誤：應該是 0
valence band current ext,loc = nan nan          # 積分失敗
```

**根本原因**:
1. 局域化態搜尋算法錯誤
2. 波函數積分邊界條件處理錯誤
3. 電流密度積分範圍設定問題

## 🔧 解決方案規劃

### Phase 1: Poisson 求解器修復 (優先級: 🔥 Critical)

#### 1.1 實現動態收斂準則
```python
# 需要實現的邏輯
max_iterations = 5000
convergence_tolerance = 1e-6
for iteration in range(max_iterations):
    # SOR 迭代
    residual = calculate_residual()
    if iteration % 100 == 0:
        print(f"ITER,Pot0 = {iteration:8d} {potential_at_center:.8E}")
    if residual < convergence_tolerance:
        break
```

#### 1.2 整合 Golden Section Search
```python
# 需要整合 GSECT 到主迭代循環
def poisson_solve_with_gsect():
    for iteration in range(max_iterations):
        # SOR 更新
        # 每隔 N 次迭代調用 Golden Section
        if iteration % golden_section_interval == 0:
            optimal_potential = golden_section_search()
            update_boundary_potential(optimal_potential)
```

### Phase 2: 電流計算修復 (優先級: 🔥 Critical)

#### 2.1 修復局域化態搜尋
```python
# 需要修復的邏輯
def find_localized_states():
    # 正確的邊界條件檢查
    # 能量範圍和積分限制
    # 波函數歸一化
    pass
```

#### 2.2 修復電流積分
```python
# 需要實現正確的電流計算
def calculate_current():
    # 修復積分範圍
    # 正確的透射係數計算
    # 避免 NaN 值的數值穩定性改進
    pass
```

## 📋 測試計畫

### 測試階段 1: Poisson 求解器驗證
- [ ] 實現動態收斂邏輯
- [ ] 整合 Golden Section Search
- [ ] 驗證 band bending 值達到 ~0.07 V
- [ ] 檢查迭代行為匹配 Fortran

### 測試階段 2: 電流計算驗證
- [ ] 修復局域化態計算 (目標: 0 個局域化態)
- [ ] 修復電流積分 (目標: ~-1.8E-14 A)
- [ ] 消除 NaN 值問題
- [ ] 驗證完整電壓掃描結果

### 測試階段 3: 端到端驗證
- [ ] 完整 multint 運行
- [ ] 多個偏壓點驗證
- [ ] 性能和穩定性測試
- [ ] 文件更新和維護指南

## 📊 測試基準和成功標準

### 數值精度目標
| 項目 | 目標精度 | 當前狀態 |
|------|----------|----------|
| Band bending | < 5% 相對誤差 | ❌ 99.9% 偏差 |
| 電流值 | < 10% 相對誤差 | ❌ NaN 錯誤 |
| 載流子密度 | < 1% 相對誤差 | ✅ 0.14% |

### 行為一致性目標
- [ ] 迭代收斂行為匹配
- [ ] 中途輸出格式一致
- [ ] 電流計算無 NaN 值
- [ ] 完整電壓掃描成功

## 🚀 下一步行動計畫

### 立即行動 (今日)
1. 檢查和修復 Poisson 求解器的動態收斂邏輯
2. 實現正確的 Golden Section Search 整合
3. 建立詳細的除錯輸出機制

### 短期目標 (3-5 天)
1. 修復電流計算的積分算法
2. 解決局域化態搜尋問題
3. 完成端到端數值驗證

### 中期目標 (1-2 週)
1. 性能優化和穩定性改進
2. 完整的回歸測試套件
3. 更新所有映射文件

---

**測試負責人**: Claude AI Assistant  
**下次更新**: 修復第一個關鍵問題後  
**緊急聯絡**: 如 band bending 值未達到 0.01 V 以上需立即檢查收斂邏輯
