# 詳細映射：intcurr-6.2.f ↔ physics/core/schrodinger.py

## 📁 檔案資訊

**Fortran 原始檔**: `src/fortran/MultInt/intcurr-6.2.f`  
**Python 對應檔**: `src/physics/core/schrodinger.py`  
**映射完成度**: 35% ⚠️  
**優先級**: **HIGH** (電流計算是模擬的最終目標)

## 📝 檔案描述

### Fortran 檔案功能
INTCURR 模組負責計算穿隧電流，包含：
- 使用 Bardeen 公式和 T&H 近似
- 求解 1D 薛丁格方程進行數值積分
- 計算延伸態 (extended states) 和局域態 (localized states) 的電流貢獻
- 自洽電荷密度計算 (CDESEM, CDESURF, CDLSEM, CDLSURF, CDEVAC, CDLVAC)

### Python 檔案功能
SchrodingerSolver 類別實現：
- 基於 transfer matrix 方法的量子力學穿隧電流計算
- WKB 近似求解
- 物件導向設計的能帶結構處理

## 🔄 主要函數對應關係

### 1. 主控函數
| Fortran | Python | 狀態 |
|---------|--------|------|
| `SUBROUTINE INTCURR(...)` | `SchrodingerSolver.solve_tunneling_current()` | ⚠️ 部分完成 |

**關鍵對應**:
```fortran
! Fortran 主控邏輯 (lines 37-178)
SUBROUTINE INTCURR(IMPOT,BARR,PROF,NBARR1,NV,NS,NSP,...)
```
```python
# Python 對應方法 (lines 87-139)
def solve_tunneling_current(self, potential_profile: PotentialProfile,
                          band_params: Dict, bias_voltage: float,
                          temperature: float = 300.0) -> TunnelCurrent:
```

### 2. 價帶電流計算
| Fortran | Python | 狀態 |
|---------|--------|------|
| `SUBROUTINE VBCURR1(...)` | `_calculate_extended_current()` | ❌ 不完整 |
| `CALL VBwf(...)` | `_solve_1d_schrodinger()` | ⚠️ 基礎實現 |
| `CALL VBloc(...)` | `_find_localized_states()` | ❌ 未正確實現 |

### 3. 導帶電流計算
| Fortran | Python | 狀態 |
|---------|--------|------|
| `SUBROUTINE CBCURR1(...)` | `_calculate_extended_current()` | ❌ 不完整 |
| `CALL CBwf(...)` | `_solve_1d_schrodinger()` | ⚠️ 基礎實現 |
| `CALL CBloc(...)` | `_find_localized_states()` | ❌ 未正確實現 |

## 📊 詳細行對行映射

### 主要結構對應

#### A. 初始化和參數設定
```fortran
! Fortran (lines 37-75)
SUBROUTINE INTCURR(IMPOT,BARR,PROF,NBARR1,NV,NS,NSP,...)
DIMENSION S(NSDIM),BARR(NVDIM1),PROF(NSDIM),NLOC(4),...
DATA RQUANT/12900./
PI=4.*ATAN(1.)
tk1=tk
tk2=tk
```

↔

```python
# Python (lines 60-86)
class SchrodingerSolver:
    def __init__(self, mass_electron: float = PC.M0):
        self.m0 = mass_electron
        self.hbar = PC.HBAR
        self.energy_tolerance = 1e-6  # eV
        self.max_iterations = 100
```

#### B. 能帶剖面創建
```fortran
! Fortran VB profile (lines 78-104)
DO 100 J=1,NSP
   SZ=S(J)
   VBPROF(J)=PROF(J)+VBEDGE(SZ)
   IF (J.EQ.1) THEN
      PMAX=VBPROF(J)
   ELSE
      PMAX=AMAX1(PMAX,VBPROF(J))
   END IF
100 CONTINUE
EV=PROF(NS)+VBEDGE(S(NS))
```

↔

```python
# Python (lines 141-166)
def _create_band_profiles(self, potential_profile: PotentialProfile,
                        band_gap: float) -> BandProfile:
    z_combined, pot_combined = potential_profile.get_combined_profile()
    vb_profile = pot_combined - band_gap
    cb_profile = pot_combined
    vb_max = np.max(vb_profile)
    cb_min = np.min(cb_profile)
    vb_bulk = vb_profile[-1]
    cb_bulk = cb_profile[-1]
```

#### C. 價帶電流計算主迴圈
```fortran
! Fortran VBCURR1 延伸態 (lines 209-340)
emax=EV
emin=amin1(ef-10.*tk1,ef+bias-10.*tk2)
dele=(emax-emin)/ne
do 120 iwky=0,nwk-1
   wky=iwky*delwk
   do 115 iwkx=0,nwk-1
      wkx=iwkx*delwk
      wkparr=sqrt(wkx**2+wky**2)
      do 110 ie=1,ne
         ener=emax-(ie-0.5)*dele
         call VBwf(...)
         trans=2.*nwkdeg*(2.*wf)**2*WKFTIP/(WKSEM/EFFM)
         sum=sum+trans*occdiff
```

↔

```python
# Python (lines 301-350) - 需要實現
def _calculate_extended_current(self, band_profile: BandProfile,
                              band_type: str, effective_masses: List[float],
                              fermi_level: float, tip_fermi: float,
                              temperature: float) -> float:
    # 目前實現不完整 - 返回 NaN
    return float('nan')
```

## ⚠️ 關鍵問題識別

### 1. 電流計算返回 NaN
**Fortran 行為**: 
- 精確數值積分計算穿隧概率
- 費米分佈權重正確計算
- 波矢量空間積分 (k-space integration)

**Python 問題**: 
- `_calculate_extended_current()` 返回 `NaN`
- 缺少正確的 k-space 積分
- 費米函數差值 (`occdiff`) 未正確計算

### 2. 局域態搜尋不完整
**Fortran 實現**:
```fortran
! VB 局域態搜尋 (lines 341-420)
call VBloc(IMPOT,n,wf,wfderiv,ener,wkparr,sep,bias,...)
if (n.eq.nsav) go to 310
IF (PSISEM(1).NE.0) THEN
   NLOC=NLOC+1
   write(6,*) 'VB localized state at energy ',ener
```

**Python 問題**:
- `_find_bound_states_1d()` 使用簡化的 shooting method
- 未考慮平行動量 (k_parallel) 的影響
- 邊界條件檢查過於簡化

### 3. 波函數求解方法差異
**Fortran**: 調用專門的 `VBwf`, `CBwf`, `VBloc`, `CBloc` 函數  
**Python**: 使用通用的 `_solve_1d_schrodinger()` 有限差分法

## 🔧 修復計劃

### Phase 1: 修復電流計算 (HIGH)
1. **實現正確的 k-space 積分**
   ```python
   def _calculate_extended_current(self, ...):
       current = 0.0
       for iwky in range(nwk):
           for iwkx in range(nwk):
               k_parallel = np.sqrt(wkx**2 + wky**2)
               for energy in energy_grid:
                   # 費米分佈差值
                   occ_diff = self._fermi_difference(energy, bias, ...)
                   # 穿隧概率
                   transmission = self._calculate_transmission(...)
                   current += transmission * occ_diff * degeneracy
       return current
   ```

2. **添加正確的費米函數計算**
   ```python
   def _fermi_difference(self, energy, bias, fermi_level, temperature):
       occ_tip = 1.0 / (1.0 + np.exp((energy - bias - fermi_level) / (PC.KB * temperature)))
       occ_sem = 1.0 / (1.0 + np.exp((energy - fermi_level) / (PC.KB * temperature)))
       return occ_tip - occ_sem
   ```

### Phase 2: 修復局域態搜尋 (MEDIUM)
1. **實現 Fortran 風格的能量掃描**
   ```python
   def _find_localized_states_fortran_style(self, ...):
       localized_count = 0
       for energy in energy_range:
           n, wavefunction = self._solve_bound_state(energy, k_parallel)
           if n != n_previous and wavefunction[0] != 0:
               localized_count += 1
               # 記錄局域態
   ```

### Phase 3: 整合電荷密度計算 (MEDIUM)
1. **添加自洽電荷密度陣列**
   ```python
   class ChargeAccumulator:
       def __init__(self):
           self.cde_sem = np.zeros(...)  # CDESEM
           self.cde_surf = 0.0           # CDESURF
           self.cdl_sem = np.zeros(...)  # CDLSEM
           # ... 其他陣列
   ```

## 📈 當前完成狀態

### ✅ 已實現功能
- 基本類別結構 (`SchrodingerSolver`)
- 能帶剖面創建 (`_create_band_profiles`)
- 1D 薛丁格方程求解框架 (`_solve_1d_schrodinger`)
- 波函數正規化 (`_normalize_wavefunction`)

### ⚠️ 部分完成功能
- 局域態搜尋 (過於簡化)
- 有限差分法實現 (缺少邊界條件)
- 數據結構定義 (缺少電荷密度)

### ❌ 未實現功能
- **延伸態電流計算** (關鍵缺失)
- **局域態電流計算** (關鍵缺失)
- **k-space 積分** (關鍵缺失)
- **費米分佈權重** (關鍵缺失)
- **穿隧概率計算** (關鍵缺失)
- **自洽電荷密度累積** (關鍵缺失)

## 🎯 驗證標準

修復完成後，應該能夠：
1. **電流數值匹配**: Python 計算的電流值與 Fortran 結果在 5% 內
2. **局域態數量匹配**: 找到的局域態數量與 Fortran 一致
3. **電荷密度分佈匹配**: 自洽電荷密度陣列與 Fortran 結果一致
4. **無 NaN 結果**: 所有電流計算返回有效數值

## 💡 實施建議

1. **優先修復延伸態電流計算** - 這是導致當前 NaN 結果的主因
2. **參考 POTEXPAND 實現** - 理解 Fortran 的網格擴展邏輯
3. **逐步對比數值結果** - 在每個修復階段驗證中間結果
4. **保持數值精度** - 使用 double precision 避免累積誤差

**修復這個模組是實現完整 SEMITIP 功能的關鍵一步！** 🚀
