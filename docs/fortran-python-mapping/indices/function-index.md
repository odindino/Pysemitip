# 函數級別交叉引用索引

## 📚 概述

本文件提供 Fortran 和 Python 實現之間的詳細函數對應關係，便於程式維護和調試時快速定位對應代碼。

## 🔍 主要子程序對應表

### A. 主控制流程

| Fortran 子程序 | 檔案位置 | Python 對應函數 | 檔案位置 | 映射狀態 |
|---------------|----------|----------------|----------|----------|
| `PROGRAM MultInt3` | MultInt3-6.4.f:1 | `MultIntSimulation.run()` | simulation/multint.py:45 | ✅ 完成 |
| 參數讀取區塊 | MultInt3-6.4.f:52-65 | `YamlConfigReader.load_config()` | core/filereader.py:30 | ✅ 完成 |
| 電壓掃描迴圈 | MultInt3-6.4.f:200-600 | `MultIntSimulation._voltage_sweep()` | simulation/multint.py:150 | ✅ 完成 |

### B. Poisson 求解器

| Fortran 子程序 | 檔案位置 | Python 對應函數 | 檔案位置 | 映射狀態 |
|---------------|----------|----------------|----------|----------|
| `SUBROUTINE SEMITIP3` | semitip3-6.1.f:80 | `PoissonSolver.solve()` | physics/core/poisson.py:85 | ⚠️ 部分完成 |
| `FUNCTION SEMMIN` | semitip3-6.1.f:45 | `PoissonSolver._semiconductor_residual()` | physics/core/poisson.py:352 | ✅ 完成 |
| `FUNCTION SURFMIN` | semitip3-6.1.f:62 | `PoissonSolver._surface_residual()` | physics/core/poisson.py:380 | ✅ 完成 |
| SOR 迭代核心 | semitip3-6.1.f:300-500 | `PoissonSolver._sor_iteration()` | physics/core/poisson.py:180 | ⚠️ 需檢查 |
| Band bending 計算 | semitip3-6.1.f:750-800 | `PoissonSolver.get_band_bending()` | physics/core/poisson.py:749 | ⚠️ 需驗證 |

### C. 數值方法

| Fortran 子程序 | 檔案位置 | Python 對應函數 | 檔案位置 | 映射狀態 |
|---------------|----------|----------------|----------|----------|
| `SUBROUTINE GSECT` | gsect-6.0.f:10 | `golden_section_search()` | physics/core/poisson.py:520 | ✅ 完成 |
| 黃金分割初始化 | gsect-6.0.f:15-25 | `golden_section_search()` 初始化 | physics/core/poisson.py:525 | ✅ 完成 |
| 收斂檢查 | gsect-6.0.f:45-60 | `golden_section_search()` 收斂邏輯 | physics/core/poisson.py:550 | ✅ 完成 |

### D. 電荷密度計算

| Fortran 子程序 | 檔案位置 | Python 對應函數 | 檔案位置 | 映射狀態 |
|---------------|----------|----------------|----------|----------|
| `SUBROUTINE SEMIRHOMULT` | semirhomult-6.0.f:10 | `ChargeDensityCalculator.calculate_bulk_density()` | physics/core/charge_density.py:85 | ✅ 完成 |
| `FUNCTION RHOBULK` | MultInt3-6.4.f:650 | `ChargeDensityCalculator.get_bulk_density()` | physics/core/charge_density.py:190 | ✅ 完成 |
| `FUNCTION RHOCB` | semirhomult-6.0.f:350 | `ChargeDensityCalculator._electron_density()` | physics/core/charge_density.py:220 | ✅ 完成 |
| `FUNCTION RHOVB` | semirhomult-6.0.f:370 | `ChargeDensityCalculator._hole_density()` | physics/core/charge_density.py:250 | ✅ 完成 |
| `FUNCTION RHOA` | semirhomult-6.0.f:390 | `ChargeDensityCalculator._acceptor_density()` | physics/core/charge_density.py:280 | ✅ 完成 |
| `FUNCTION RHOD` | semirhomult-6.0.f:410 | `ChargeDensityCalculator._donor_density()` | physics/core/charge_density.py:310 | ✅ 完成 |

### E. 表面態計算

| Fortran 子程序 | 檔案位置 | Python 對應函數 | 檔案位置 | 映射狀態 |
|---------------|----------|----------------|----------|----------|
| `SUBROUTINE SURFRHOMULT` | surfrhomult-6.2.f:10 | `SurfaceStatesCalculator.calculate_surface_density()` | physics/materials/surface_states.py:75 | ✅ 完成 |
| `FUNCTION RHOSURF` | MultInt3-6.4.f:720 | `SurfaceStatesCalculator.get_surface_density()` | physics/materials/surface_states.py:200 | ✅ 完成 |
| `FUNCTION RHOS1` | surfrhomult-6.2.f:200 | `SurfaceStatesCalculator._distribution_1()` | physics/materials/surface_states.py:230 | ✅ 完成 |
| `FUNCTION RHOS2` | surfrhomult-6.2.f:250 | `SurfaceStatesCalculator._distribution_2()` | physics/materials/surface_states.py:260 | ✅ 完成 |
| `FUNCTION RHOS` | surfrhomult-6.2.f:300 | `SurfaceStatesCalculator._mixed_distribution()` | physics/materials/surface_states.py:290 | ✅ 完成 |
| `FUNCTION SIGSUM` | surfrhomult-6.2.f:350 | `SurfaceStatesCalculator._sigma_sum()` | physics/materials/surface_states.py:320 | ✅ 完成 |

### F. 電位處理

| Fortran 子程序 | 檔案位置 | Python 對應函數 | 檔案位置 | 映射狀態 |
|---------------|----------|----------------|----------|----------|
| `SUBROUTINE POTCUT3` | potcut3-6.0.f:20 | `PotentialProcessor.extract_potential_profile()` | physics/core/potential.py:250 | ✅ 完成 |
| 真空區域插值 | potcut3-6.0.f:80-120 | `PotentialProcessor._extract_vacuum_potential()` | physics/core/potential.py:350 | ✅ 完成 |
| 半導體區域插值 | potcut3-6.0.f:130-170 | `PotentialProcessor._extract_semiconductor_potential()` | physics/core/potential.py:400 | ✅ 完成 |
| `SUBROUTINE POTEXPAND` | potexpand-6.1.f:15 | `PotentialProcessor.expand_potential()` | physics/core/potential.py:500 | ⚠️ 部分實現 |

### G. 電流計算

| Fortran 子程序 | 檔案位置 | Python 對應函數 | 檔案位置 | 映射狀態 |
|---------------|----------|----------------|----------|----------|
| `SUBROUTINE INTCURR` | intcurr-6.2.f:20 | `SchrodingerSolver.calculate_current()` | physics/core/schrodinger.py:50 | ❌ 部分實現 |
| Schrödinger 積分 | intcurr-6.2.f:100-200 | `SchrodingerSolver._integrate_schrodinger()` | physics/core/schrodinger.py:120 | ❌ 問題存在 |
| 傳輸係數計算 | intcurr-6.2.f:250-300 | `SchrodingerSolver._transmission_coefficient()` | physics/core/schrodinger.py:200 | ❌ NaN 問題 |

### H. 輔助功能

| Fortran 子程序 | 檔案位置 | Python 對應函數 | 檔案位置 | 映射狀態 |
|---------------|----------|----------------|----------|----------|
| `SUBROUTINE CONTR3` | contr3-6.0.f:10 | `ContourPlotter.plot_potential_contours()` | visualization/contour_plots.py:30 | ❌ 未實現 |
| `INTEGER FUNCTION IGETAR` | MultInt3-6.4.f:630 | `GridGeometry._get_surface_area()` | physics/core/grid.py:450 | ✅ 完成 |
| `INTEGER FUNCTION IGETREG` | MultInt3-6.4.f:640 | `GridGeometry._get_material_region()` | physics/core/grid.py:470 | ✅ 完成 |

## 🔧 核心算法對應

### 1. SOR 迭代 (Successive Over-Relaxation)

**Fortran 核心邏輯**:
```fortran
! semitip3-6.1.f 第300-350行
DO 40 I=1,NR
   DO 40 J=1,NS
      CALL GSECT(SEMMIN,Pot1,Pot2,DELPOT1,ITER1,DELSEMP)
      SEM(1,I,J,K)=(1.-OMEGA)*SEM(1,I,J,K)+OMEGA*Pot1
40 CONTINUE
```

**Python 對應**:
```python
# physics/core/poisson.py 第180-220行
def _sor_iteration(self):
    for i in range(self.grid.nr):
        for j in range(self.grid.ns):
            optimal_potential = self.golden_section_search(
                self._semiconductor_residual, pot_min, pot_max, tolerance)
            self.semiconductor_grid[i, j, k] = ((1 - self.omega) * 
                self.semiconductor_grid[i, j, k] + self.omega * optimal_potential)
```

### 2. 電荷密度插值

**Fortran 核心邏輯**:
```fortran
! MultInt3-6.4.f 第650-680行
FUNCTION RHOBULK(IREG,POTEN,X,Y,S,I,J,K,NR,NS,NP)
J=NINT((EF+POTEN-ESTART)/DELE)+1
IF (J.LT.1) J=1
IF (J.GT.NE) J=NE
RHOBULK=RHOBTAB(IREG,J)
```

**Python 對應**:
```python
# physics/core/charge_density.py 第190-210行
def get_bulk_density(self, region_idx, potential, x, y, z, i, j, k):
    energy = self.fermi_energy + potential
    energy_index = int(round((energy - self.energy_start) / self.energy_step))
    energy_index = max(0, min(energy_index, len(self.energy_grid) - 1))
    return self.charge_density_tables[region_idx][energy_index]
```

### 3. 黃金分割搜索

**Fortran 核心邏輯**:
```fortran
! gsect-6.0.f 第25-45行
GS = 0.3819660
X1 = XMIN + GS * (XMAX - XMIN)
X2 = XMAX - GS * (XMAX - XMIN)
F1 = F(X1)
F2 = F(X2)
```

**Python 對應**:
```python
# physics/core/poisson.py 第520-540行
def golden_section_search(self, func, xmin, xmax, tolerance, max_iter=100):
    golden_ratio = 0.3819660
    x1 = xmin + golden_ratio * (xmax - xmin)
    x2 = xmax - golden_ratio * (xmax - xmin)
    f1 = func(x1)
    f2 = func(x2)
```

## 📊 映射統計

### 完成度統計
- **完全實現**: 15 個函數 (62.5%)
- **部分實現**: 6 個函數 (25.0%)
- **未實現**: 3 個函數 (12.5%)

### 關鍵問題函數
1. **intcurr-6.2.f 相關**: Schrödinger 求解器有 NaN 問題
2. **potexpand-6.1.f**: 電位展開算法不完整
3. **contr3-6.0.f**: 等高線繪圖未實現

## 🔍 快速查找索引

### 按 Python 檔案查找
- **simulation/multint.py**: MultInt3-6.4.f 主程式
- **physics/core/poisson.py**: semitip3-6.1.f + gsect-6.0.f
- **physics/core/charge_density.py**: semirhomult-6.0.f
- **physics/materials/surface_states.py**: surfrhomult-6.2.f
- **physics/core/potential.py**: potcut3-6.0.f + potexpand-6.1.f
- **physics/core/schrodinger.py**: intcurr-6.2.f
- **visualization/contour_plots.py**: contr3-6.0.f

### 按功能領域查找
- **主控制流程**: MultInt3-6.4.f ↔ simulation/multint.py
- **數值求解**: semitip3-6.1.f + gsect-6.0.f ↔ physics/core/poisson.py
- **物理模型**: semirhomult-6.0.f + surfrhomult-6.2.f ↔ physics/core/charge_density.py + physics/materials/surface_states.py
- **資料處理**: potcut3-6.0.f + potexpand-6.1.f ↔ physics/core/potential.py
- **量子計算**: intcurr-6.2.f ↔ physics/core/schrodinger.py
- **可視化**: contr3-6.0.f ↔ visualization/contour_plots.py

---

**索引建立**: 2025-06-06  
**涵蓋函數**: 24 個主要子程序  
**最後更新**: 完成核心映射後  
**維護頻率**: 每次重大修改後更新
