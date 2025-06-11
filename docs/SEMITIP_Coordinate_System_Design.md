# SEMITIP 座標系統設計規範

**文件版本**: 1.0  
**創建日期**: 2025年6月11日  
**作者**: odindino  
**目的**: 為 Pysemitip Phase 3 幾何層實現提供完整的座標系統設計規範

---

## 1. 概述

SEMITIP (Scanning tunneling microscopy Electronic and Magnetic Interaction with Tips) 使用一個精心設計的三維座標系統來模擬 STM 幾何配置。本文件詳細描述了該座標系統的數學定義、物理意義、網格結構以及在 Python 實現中的具體要求。

### 1.1 設計目標

- **物理直觀性**: 座標系統直接對應 STM 實驗幾何
- **數值效率**: 支援自適應網格細化和對稱性優化
- **計算精度**: 在關鍵區域（探針尖端、表面附近）提供高解析度
- **Fortran 相容性**: 與原始 SEMITIP Fortran 程式碼數值一致

---

## 2. 座標系統定義

### 2.1 主座標系統：柱座標系 (r, φ, z)

**基本定義**:
- **r**: 徑向距離 [nm]，範圍 [0, R_max]
- **φ**: 方位角 [弧度]，範圍 [0, 2π] 或 [0, π]（鏡像對稱時）
- **z**: 軸向位置 [nm]，垂直於半導體表面

**座標原點**:
- **z = 0**: 半導體表面（參考面）
- **r = 0**: 柱座標軸，通常對應探針中心軸
- **φ = 0**: 參考方位，可選擇任意水平方向

**z 軸方向定義**（基於 Fortran 程式碼分析）:
```
z < 0: 真空區域（探針側）
z = 0: 半導體表面  
z > 0: 半導體內部
```

**依據**: `semitip3-6.1.f` 第 736-741 行輸出格式
```fortran
DO 487 J=NV,1,-1
  WRITE(11,*) -J*DELV(1),VAC(1,1,J,1),VAC(1,NR,J,1)  ! 負 z (真空)
487 CONTINUE
WRITE(11,*) 0.,VSINT(1,1,1),VSINT(1,NR,1)           ! z=0 (表面)
DO 490 J=1,NS
   WRITE(11,*) S(J),SEM(1,1,J,1),SEM(1,NR,J,1)       ! 正 z (半導體)
490 CONTINUE
```

### 2.2 探針幾何：雙曲線座標系 (ξ, η, φ)

SEMITIP 使用雙曲線座標系精確描述探針幾何，這對於準確計算探針表面的電位和電場至關重要。

**核心參數** (`semitip3-6.1.f` 第 97-101 行):
```fortran
ETAT=1./SQRT(1.+1./SLOPE**2)    ! 探針錐角參數
A=RAD*SLOPE**2/ETAT             ! 雙曲線焦點位置 [nm]
SPRIME=A*ETAT                   ! 探針特徵長度 [nm]
Z0=SEP-SPRIME                   ! 雙曲座標原點位置 [nm]  
C=Z0/SPRIME                     ! 幾何參數（無量綱）
```

**物理意義**:
- **SLOPE**: 探針錐角的正切值
- **RAD**: 探針曲率半徑 [nm]
- **SEP**: 探針末端到半導體表面的分離距離 [nm]
- **Z0**: 雙曲座標系原點在 z 軸上的位置

**探針形狀函數** (`semitip3-6.1.f` 第 600-605 行):
```fortran
real function p(r)
p=0.
if (r.lt.rad2) p=sqrt(rad2**2-r**2)  ! 球形突起
return
end
```

---

## 3. 網格結構設計

### 3.1 網格維度參數

**Fortran 參數定義** (`MultInt3-6.4.f` 第 27-28 行):
```fortran
PARAMETER(NRDIM=512,NVDIM=64,NSDIM=512,NPDIM=64)
```

**Python 實現對應**:
```python
class GridDimensions:
    MAX_R_POINTS = 512    # NRDIM: 徑向網格點數上限
    MAX_V_POINTS = 64     # NVDIM: 真空區軸向網格點數上限  
    MAX_S_POINTS = 512    # NSDIM: 半導體區軸向網格點數上限
    MAX_P_POINTS = 64     # NPDIM: 角向網格點數上限
```

### 3.2 網格分佈函數

**徑向網格（r 方向）**:

**Fortran 實現** (`semitip3-6.1.f` 第 116 行):
```fortran
R(I)=(2*NR*DELR0/PI)*TAN(PI*(I-0.5)/(2.*NR))
```

**數學表達式**:
```
r(i) = (2 × NR × Δr₀ / π) × tan(π × (i - 0.5) / (2 × NR))
```

**特點**:
- 小 r 處密集分佈（探針軸附近）
- 大 r 處稀疏分佈（遠場區域）
- 保證計算關鍵區域的高精度

**角向網格（φ 方向）**:

**Fortran 實現** (`semitip3-6.1.f` 第 397-400 行):
```fortran
IF (MIRROR.EQ.1) THEN
   DELP=PI/FLOAT(NP)        ! 鏡像對稱：0 到 π
ELSE
   DELP=2.*PI/FLOAT(NP)     ! 完整圓：0 到 2π
END IF
```

**軸向網格（z 方向）**:

**半導體區域** (`semitip3-6.1.f` 第 164 行):
```fortran
S(J)=(2*NS*DELS0/PI)*TAN(PI*(J-0.5)/(2.*NS))
```

**真空區域**: 
- 依賴於徑向位置的變網格間距
- 考慮探針幾何的影響
- 在表面附近密集分佈

### 3.3 三層自適應細化策略

**細化過程** (`semitip3-6.1.f` 第 227-240 行):

```fortran
NR=NR*2     ! 徑向網格點數翻倍
NS=NS*2     ! 半導體軸向網格點數翻倍  
NV=NV*2     ! 真空軸向網格點數翻倍
NP=NP*2     ! 角向網格點數翻倍
DELR0=DELR0/2.   ! 基礎網格間距減半
DELS0=DELS0/2.   ! 基礎網格間距減半
```

**實現策略**:

1. **第一層（粗網格）**:
   - 快速收斂到近似解
   - 提供全局電位分佈
   
2. **第二層（中等網格）**:
   - 使用第一層結果作為初值
   - 細化關鍵區域的計算
   
3. **第三層（細網格）**:
   - 達到最終精度要求
   - 確保數值收斂性

**收斂控制**:
```fortran
EP(IP)      ! 各細化層級的收斂容差
ITMAX(IP)   ! 各層級的最大迭代次數
```

---

## 4. 邊界條件與幾何配置

### 4.1 探針邊界條件

**探針電位設定**:
```
V_tip = BIAS + CPot
```
其中：
- **BIAS**: 施加的偏壓 [V]
- **CPot**: 接觸電位差 [V]

**探針識別**: 
使用邏輯陣列 `TIP(I,J,K)` 標識探針內部的網格點

### 4.2 半導體表面邊界

**界面連續性條件**:
- 電位連續: `V_vacuum(surface) = V_semiconductor(surface)`
- 電場法向分量連續: `ε_vac × E_n,vac = ε_sem × E_n,sem`

**表面電荷密度**:
由 `RHOSURF` 函數計算，包含表面態貢獻

### 4.3 外邊界條件

**徑向邊界**:
- `IBC=0`: Dirichlet 邊界條件（固定電位）
- `IBC=1`: Neumann 邊界條件（固定電場）

**軸向邊界**:
- 深入半導體：接地條件或體電位
- 遠離探針：零場條件

---

## 5. 對稱性利用

### 5.1 鏡像對稱（MIRROR=1）

**適用條件**:
- 探針位於對稱軸上（Y0=0）
- 幾何配置關於某個平面對稱

**實現方式**:
- 只計算 φ ∈ [0, π] 範圍
- 計算量減少約 50%
- 結果通過鏡像映射獲得完整解

**邊界處理**:
```fortran
IF (MIRROR.EQ.1) THEN
   ! 在 φ=0 和 φ=π 處實施鏡像邊界條件
   ! 確保解的對稱性
END IF
```

### 5.2 其他對稱性

**柱對稱（默認）**:
- 電位不依賴於 φ 角度
- 進一步簡化計算

**平移對稱**:
- 在某些配置下可利用 z 方向的平移對稱性

---

## 6. 實現要點

### 6.1 數值精度要求

**浮點精度**:
- 使用雙精度浮點數（`numpy.float64`）
- 確保與 Fortran 計算的一致性

**收斂判據**:
```python
def convergence_check(pot_old, pot_new, tolerance):
    """檢查電位收斂性"""
    relative_change = np.max(np.abs(pot_new - pot_old) / 
                           (np.abs(pot_new) + 1e-12))
    return relative_change < tolerance
```

### 6.2 記憶體管理

**陣列分配策略**:
```python
# 真空區電位陣列
VAC = np.zeros((2, NRDIM, NVDIM, NPDIM), dtype=np.float64)

# 半導體區電位陣列  
SEM = np.zeros((2, NRDIM, NSDIM, NPDIM), dtype=np.float64)

# 表面電位陣列
VSINT = np.zeros((2, NRDIM, NPDIM), dtype=np.float64)
```

**記憶體優化**:
- 利用對稱性減少陣列大小
- 實施稀疏存儲（如適用）
- 及時釋放中間計算陣列

### 6.3 與現有模組的介面

**與 physics 模組整合**:
```python
from physics import PoissonSolver, ChargeDensityCalculator
from geometry import STMGeometry, Grid3D

# 座標系統應與現有的物理模型無縫整合
```

**向後相容性**:
- 保持與現有 API 的相容性
- 提供座標轉換工具
- 支援舊格式的輸入/輸出

---

## 7. 驗證與測試策略

### 7.1 與 Fortran 的數值比較

**測試案例**:
1. **基礎網格生成**: 比較網格點位置的數值精度
2. **探針幾何**: 驗證雙曲座標轉換的正確性
3. **邊界條件**: 確保邊界處理的一致性
4. **自適應細化**: 比較各層級網格的收斂性

**精度要求**:
- 網格點位置: 相對誤差 < 1e-12
- 幾何參數: 相對誤差 < 1e-10
- 邊界條件: 完全一致

### 7.2 物理合理性檢查

**幾何合理性**:
- 探針形狀的連續性
- 網格間距的單調性
- 邊界條件的物理意義

**數值穩定性**:
- 細化過程的單調收斂
- 對稱性的保持
- 邊界條件的穩定性

---

## 8. 實現順序

### 8.1 第一階段：基礎幾何類

```python
# 優先實現
class STMGeometry:
    """STM 基礎幾何配置"""
    
class Grid3D:  
    """三維柱座標網格管理"""
    
class TipGeometry:
    """探針雙曲座標幾何"""
```

### 8.2 第二階段：網格細化

```python
class AdaptiveGridRefinement:
    """三層自適應網格細化"""
    
class BoundaryConditions:
    """邊界條件處理"""
```

### 8.3 第三階段：對稱性與優化

```python
class SymmetryHandler:
    """鏡像對稱和其他對稱性處理"""
    
class GridOptimization:
    """網格優化和記憶體管理"""
```

---

## 9. 參考資料

### 9.1 Fortran 程式碼參考

- `src/fortran/MultInt/semitip3-6.1.f`: 主要座標系統實現
- `src/fortran/MultInt/MultInt3-6.4.f`: 主程式和參數定義
- `src/fortran/MultInt/potcut3-6.0.f`: 電位輸出格式

### 9.2 技術文件參考

- `docs/Fortran-semitip/NewCoords_diagram.pdf`: 座標系統圖表
- `docs/Fortran-semitip/NewCoords_doc.pdf`: 座標系統文檔
- `docs/Fortran-semitip/SEMITIP V6, Technical Manual.pdf`: 技術手冊

### 9.3 物理模型參考

- 現有的 `src/physics/poisson.py`: 泊松方程求解器
- 現有的 `src/physics/charge_density.py`: 電荷密度計算
- 現有的 `src/utils/numerical.py`: 數值計算工具

---

## 10. 更新日誌

**版本 1.0** (2025年6月11日):
- 初始版本
- 完整的座標系統定義
- 基於 Fortran 程式碼分析的具體依據
- 詳細的實現規範和測試策略

---

*本文件將隨著 Phase 3 實現的進展持續更新，確保座標系統設計的準確性和完整性。*