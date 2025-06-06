# 變數名稱翻譯對照表

## 📚 概述

本文件提供 Fortran 和 Python 實現之間的變數名稱對應關係，包括 COMMON 區塊變數、陣列名稱、物理常數等的完整對照。

## 🏗️ COMMON 區塊對應

### A. /SEMI/ 區塊 - 半導體參數

| Fortran 變數 | 型態 | Python 對應 | 檔案位置 | 說明 |
|-------------|------|-------------|----------|------|
| `TK` | REAL | `temperature` | physics/materials/semiconductor.py:15 | 溫度 (K) |
| `EGAP(NREGDIM)` | REAL | `band_gap` | physics/materials/semiconductor.py:20 | 能隙 (eV) |
| `ED(NREGDIM)` | REAL | `donor_level` | physics/materials/semiconductor.py:25 | 施主能階 (eV) |
| `EA(NREGDIM)` | REAL | `acceptor_level` | physics/materials/semiconductor.py:30 | 受體能階 (eV) |
| `ACB(NREGDIM)` | REAL | `conduction_band_offset` | physics/materials/semiconductor.py:35 | 導帶偏移 (eV) |
| `AVB(NREGDIM)` | REAL | `valence_band_offset` | physics/materials/semiconductor.py:40 | 價帶偏移 (eV) |
| `CD(NREGDIM)` | REAL | `electron_concentration` | physics/materials/semiconductor.py:45 | 電子有效密度 (cm⁻³) |
| `CA(NREGDIM)` | REAL | `hole_concentration` | physics/materials/semiconductor.py:50 | 電洞有效密度 (cm⁻³) |
| `IDEG(NREGDIM)` | INTEGER | `degeneracy_flag` | physics/materials/semiconductor.py:55 | 簡併標記 |
| `IINV(NREGDIM)` | INTEGER | `inversion_flag` | physics/materials/semiconductor.py:60 | 反轉標記 |
| `DELVB(NREGDIM)` | REAL | `valence_band_shift` | physics/materials/semiconductor.py:65 | 價帶位移 (eV) |

### B. /SURF/ 區塊 - 表面態參數

| Fortran 變數 | 型態 | Python 對應 | 檔案位置 | 說明 |
|-------------|------|-------------|----------|------|
| `ISTK` | INTEGER | `temperature_dependent` | physics/materials/surface_states.py:15 | 溫度依賴標記 |
| `TK1` | REAL | `surface_temperature` | physics/materials/surface_states.py:20 | 表面溫度 (K) |
| `EN0(NARDIM)` | REAL | `charge_neutrality_level` | physics/materials/surface_states.py:25 | 電荷中性點 (eV) |
| `EN(NARDIM,2)` | REAL | `energy_levels` | physics/materials/surface_states.py:30 | 表面態能階 (eV) |
| `DENS(NARDIM,2)` | REAL | `densities` | physics/materials/surface_states.py:35 | 表面態密度 (cm⁻²) |
| `FWHM(NARDIM,2)` | REAL | `fwhm` | physics/materials/surface_states.py:40 | 高斯分佈半寬 (eV) |
| `ECENT(NARDIM,2)` | REAL | `centroids` | physics/materials/surface_states.py:45 | 分佈中心能量 (eV) |

### C. /CD/ 區塊 - 電荷密度表格

| Fortran 變數 | 型態 | Python 對應 | 檔案位置 | 說明 |
|-------------|------|-------------|----------|------|
| `EF` | REAL | `fermi_energy` | physics/core/charge_density.py:15 | 費米能階 (eV) |
| `ESTART` | REAL | `energy_start` | physics/core/charge_density.py:20 | 能量範圍起點 (eV) |
| `DELE` | REAL | `energy_step` | physics/core/charge_density.py:25 | 能量步長 (eV) |
| `NE` | INTEGER | `num_energies` | physics/core/charge_density.py:30 | 能量點數 |
| `RHOBTAB(NREGDIM,NEDIM)` | REAL | `charge_density_tables` | physics/core/charge_density.py:35 | 體電荷密度表格 |
| `RHOSTAB(NARDIM,NEDIM)` | REAL | `surface_density_tables` | physics/materials/surface_states.py:50 | 表面電荷密度表格 |

## 🔢 網格和陣列對應

### A. 主要網格陣列

| Fortran 陣列 | 維度 | Python 對應 | 檔案位置 | 說明 |
|-------------|------|-------------|----------|------|
| `VAC(NRDIM,NVDIM,NPDIM)` | 3D | `vacuum_grid` | physics/core/grid.py:100 | 真空區域電位 |
| `TIP(NRDIM,NVDIM,NPDIM)` | 3D | `tip_grid` | physics/core/grid.py:105 | 探針電位 |
| `SEM(NRDIM,NSDIM,NPDIM)` | 3D | `semiconductor_grid` | physics/core/grid.py:110 | 半導體電位 |
| `VSINT(NRDIM,NPDIM)` | 2D | `interface_grid` | physics/core/grid.py:115 | 界面電位 |
| `R(NRDIM)` | 1D | `r_coordinates` | physics/core/grid.py:120 | 徑向坐標 (nm) |
| `S(NSDIM)` | 1D | `z_coordinates_sem` | physics/core/grid.py:125 | 半導體 z 坐標 (nm) |
| `DELV(NRDIM)` | 1D | `z_spacing_vac` | physics/core/grid.py:130 | 真空 z 間距 (nm) |

### B. 電位剖面陣列

| Fortran 陣列 | 維度 | Python 對應 | 檔案位置 | 說明 |
|-------------|------|-------------|----------|------|
| `BARR(NVDIM1)` | 1D | `barrier_potential` | physics/core/potential.py:200 | 真空勢壘電位 |
| `PROF(NSDIM2)` | 1D | `profile_potential` | physics/core/potential.py:205 | 半導體電位剖面 |
| `BARR2(NVDIM2)` | 1D | `expanded_barrier` | physics/core/potential.py:210 | 展開真空勢壘 |
| `PROF2(NSDIM2)` | 1D | `expanded_profile` | physics/core/potential.py:215 | 展開半導體剖面 |

## 🔧 物理參數對應

### A. 探針參數

| Fortran 變數 | Python 對應 | 檔案位置 | 說明 |
|-------------|-------------|----------|------|
| `SEP` | `tip_sample_separation` | physics/materials/tip.py:15 | 探針-樣品距離 (nm) |
| `RAD` | `tip_radius` | physics/materials/tip.py:20 | 探針半徑 (nm) |
| `SLOPE` | `tip_cone_angle` | physics/materials/tip.py:25 | 探針錐角 |
| `BIAS` | `bias_voltage` | physics/materials/tip.py:30 | 偏壓 (V) |
| `CPot` | `contact_potential` | physics/materials/tip.py:35 | 接觸電位差 (V) |
| `PotTIP` | `tip_potential` | physics/materials/tip.py:40 | 探針電位 (V) |

### B. 計算控制參數

| Fortran 變數 | Python 對應 | 檔案位置 | 說明 |
|-------------|-------------|----------|------|
| `ITMAX` | `max_iterations` | physics/core/poisson.py:70 | 最大迭代次數 |
| `EP` | `convergence_tolerance` | physics/core/poisson.py:75 | 收斂容差 |
| `OMEGA` | `relaxation_parameter` | physics/core/poisson.py:80 | SOR 鬆弛參數 |
| `CONV` | `convergence_reached` | physics/core/poisson.py:85 | 收斂標記 |
| `IERR` | `error_code` | physics/core/poisson.py:90 | 錯誤代碼 |

### C. 維度參數

| Fortran 常數 | 值 | Python 對應 | 檔案位置 | 說明 |
|-------------|---|-------------|----------|------|
| `NRDIM` | 512 | `MAX_RADIAL_POINTS` | physics/core/grid.py:10 | 最大徑向點數 |
| `NVDIM` | 64 | `MAX_VACUUM_POINTS` | physics/core/grid.py:15 | 最大真空點數 |
| `NSDIM` | 512 | `MAX_SEMICONDUCTOR_POINTS` | physics/core/grid.py:20 | 最大半導體點數 |
| `NPDIM` | 64 | `MAX_AZIMUTHAL_POINTS` | physics/core/grid.py:25 | 最大方位角點數 |
| `NEDIM` | 50000 | `MAX_ENERGY_POINTS` | physics/core/charge_density.py:10 | 最大能量點數 |
| `NREGDIM` | 2 | `MAX_REGIONS` | physics/materials/semiconductor.py:10 | 最大區域數 |
| `NARDIM` | 2 | `MAX_SURFACE_AREAS` | physics/materials/surface_states.py:10 | 最大表面區域數 |

## 📐 數值常數對應

### A. 物理常數

| Fortran 常數 | 值 | Python 對應 | 檔案位置 | 說明 |
|-------------|---|-------------|----------|------|
| `8.617E-5` | 8.617×10⁻⁵ | `BOLTZMANN_EV` | physics/constants.py:10 | 波茲曼常數 (eV/K) |
| `1.80943E-20` | 1.80943×10⁻²⁰ | `EEP_CONSTANT` | physics/constants.py:15 | 電荷密度常數 |
| `0.02585` | 0.02585 | `THERMAL_ENERGY_300K` | physics/constants.py:20 | 300K 熱能 (eV) |
| `2.355` | 2.355 | `FWHM_TO_SIGMA` | physics/constants.py:25 | FWHM 轉標準差係數 |
| `2.507` | 2.507 | `GAUSSIAN_NORM` | physics/constants.py:30 | 高斯分佈歸一化常數 |

### B. 數值方法常數

| Fortran 常數 | 值 | Python 對應 | 檔案位置 | 說明 |
|-------------|---|-------------|----------|------|
| `0.3819660` | 0.3819660 | `GOLDEN_RATIO` | physics/core/poisson.py:15 | 黃金分割比例 |
| `1.E-6` | 1×10⁻⁶ | `DEFAULT_TOLERANCE` | physics/core/poisson.py:20 | 預設收斂容差 |
| `100` | 100 | `MAX_SOR_ITERATIONS` | physics/core/poisson.py:25 | 最大 SOR 迭代 |

## 🔄 命名規則轉換

### A. Fortran → Python 轉換規則

1. **大寫 → 小寫**: `TEMP` → `temp`
2. **底線分隔**: `PotTIP` → `tip_potential`
3. **描述性命名**: `BARR` → `barrier_potential`
4. **陣列 → 複數**: `VAC` → `vacuum_grid`
5. **去除維度後綴**: `RHOSTAB` → `surface_density_tables`

### B. 特殊轉換案例

| Fortran 原名 | 轉換邏輯 | Python 新名 | 原因 |
|-------------|----------|-------------|------|
| `BBIAS` | B + BIAS | `bias_voltages` | 去除匈牙利命名法 |
| `ICUT` | I + CUT | `cut_position` | 更描述性 |
| `DELR` | DEL + R | `radial_spacing` | 明確物理意義 |
| `DELS` | DEL + S | `semiconductor_spacing` | 明確空間維度 |
| `PCENT` | P + CENT | `band_bending_center` | 完整描述功能 |

## 🔍 快速查找表

### 按字母順序 (Fortran → Python)

| A-E | F-O | P-Z |
|-----|-----|-----|
| `ACB` → `conduction_band_offset` | `FWHM` → `fwhm` | `PCENT` → `band_bending_center` |
| `AVB` → `valence_band_offset` | `IERR` → `error_code` | `PotTIP` → `tip_potential` |
| `BARR` → `barrier_potential` | `ITMAX` → `max_iterations` | `PROF` → `profile_potential` |
| `BBIAS` → `bias_voltages` | `NE` → `num_energies` | `R` → `r_coordinates` |
| `CD` → `electron_concentration` | `OMEGA` → `relaxation_parameter` | `RHOBTAB` → `charge_density_tables` |
| `DELE` → `energy_step` | | `S` → `z_coordinates_sem` |
| `EF` → `fermi_energy` | | `SEP` → `tip_sample_separation` |
| `EGAP` → `band_gap` | | `TK` → `temperature` |

### 按功能分類查找

- **網格相關**: VAC→vacuum_grid, SEM→semiconductor_grid, R→r_coordinates
- **電荷密度**: RHOBTAB→charge_density_tables, RHOSTAB→surface_density_tables
- **電位相關**: BARR→barrier_potential, PROF→profile_potential, PotTIP→tip_potential
- **材料參數**: EGAP→band_gap, ACB→conduction_band_offset, AVB→valence_band_offset
- **數值方法**: OMEGA→relaxation_parameter, EP→convergence_tolerance, ITMAX→max_iterations

---

**變數總數**: 80+ 個主要變數  
**涵蓋範圍**: 所有 COMMON 區塊和主要陣列  
**最後更新**: 2025-06-06  
**維護週期**: 每次添加新變數後更新
