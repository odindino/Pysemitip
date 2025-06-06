# Python-Fortran 比較與修復文檔

## 問題分析階段

### 關鍵差異發現

#### 1. 迭代收斂行為差異
**Fortran 版本**:
- 顯示詳細的迭代過程 (`ITER,Pot0 = 100 -8.27837288E-02`)
- 經過3500次迭代達到收斂
- 最終 Band Bending = 6.98396191E-02

**Python 版本**:
- 只顯示最終迭代次數 (200次)
- 最終 Band Bending = 0.00005040
- 差異巨大：0.0698 vs 0.00005

#### 實際調試發現的問題

**問題1: 半導體區域電位變化極小**
```
Semiconductor potential max change: 7.64e-07  (極小變化)
Interface potential max change: 7.64e-07      (極小變化)
```
對比 Fortran 的大幅度電位變化，Python版本幾乎沒有更新。

**問題2: Golden Section Search 未生效**
界面和半導體電位更新使用的非線性求解器沒有產生有效的電位變化。

**問題3: 電荷密度耦合問題**
測試表明 dummy 電荷函數能正常回傳值，但電位更新邏輯有問題。

#### 2. 電流計算結果差異
**Fortran 版本**:
```
valence band current ext,loc = -1.81937781E-14   0.0000000
```

**Python 版本**:
```
valence band current ext,loc = nan nan
```

#### 3. 局域態計算差異
**Fortran 版本**:
```
number of VB light-hole localized states = 0
number of VB heavy-hole localized states = 0
number of VB split-off localized states = 0
```

**Python 版本**:
```
number of VB light-hole localized states = 198
```

## 根本原因分析

### 問題1: Poisson 方程求解器收斂問題
Python版本的Poisson求解器未能正確收斂，導致band bending計算錯誤。

### 問題2: 電流計算中的NaN錯誤
這通常是由於數值積分中的除零或無效數學運算導致的。

### 問題3: 局域態計算邏輯錯誤
Python版本錯誤地計算了大量的局域態，而Fortran版本正確地得到0個。

## 修復計劃

### 階段1: Poisson求解器修復
- [ ] 檢查relaxation parameter設定
- [ ] 驗證邊界條件實現
- [ ] 修復迭代收斂判定條件

### 階段2: 電流計算修復
- [ ] 檢查數值積分邊界
- [ ] 修復除零保護
- [ ] 驗證物理常數

### 階段3: 局域態計算修復
- [ ] 檢查局域態判定條件
- [ ] 驗證能帶計算
- [ ] 修復能級比較邏輯

---

**建立時間**: 2025-06-06  
**維護者**: GitHub Copilot  
**下次更新**: 修復完成後
