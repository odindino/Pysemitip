# semitip3-6.1.f ↔ physics/core/poisson.py 詳細對應

## 📋 檔案基本資訊

| 項目 | Fortran | Python |
|------|---------|--------|
| **檔案名** | semitip3-6.1.f | physics/core/poisson.py |
| **行數** | 758 行 | 810 行 |
| **版本** | 6.1 (2011年2月) | 對應實現 |
| **主要功能** | 3D Poisson 方程求解器 | 同左 |

## 🔗 整體結構對應

### Fortran 程式結構
```fortran
SUBROUTINE SEMITIP3
├── 網格初始化 (構建 TIP 和 VACUUM 網格)
├── 坐標系設定 (雙曲坐標系)
├── 邊界條件設定
├── 多網格求解迴圈 (IP=1,IPMAX)
│   ├── SOR 迭代主迴圈
│   │   ├── 真空區域更新
│   │   ├── 半導體區域更新 (呼叫 GSECT)
│   │   ├── 表面區域更新 (呼叫 GSECT)
│   │   └── 收斂檢查
│   ├── 網格細化 (加倍解析度)
│   └── 電位插值到新網格
└── 返回結果
```

### Python 對應結構
```python
class PoissonSolver:
├── __init__()                       # 對應網格初始化
├── solve()                          # 對應 SUBROUTINE SEMITIP3
├── _setup_coordinate_system()       # 對應坐標系設定
├── _set_boundary_conditions()       # 對應邊界條件
├── _multigrid_solve()              # 對應多網格迴圈
├── _sor_iteration()                # 對應 SOR 迭代
├── _update_vacuum_region()         # 對應真空區域更新
├── _update_semiconductor_region()  # 對應半導體區域更新
├── _update_surface_region()        # 對應表面區域更新
├── _check_convergence()            # 對應收斂檢查
└── _refine_grid()                  # 對應網格細化
```

## 📊 詳細程式碼對應

### 1. 主函數簽名和初始化

#### Fortran: semitip3-6.1.f 第84-95行
```fortran
SUBROUTINE SEMITIP3(SEP,RAD,SLOPE,ETAT,A,Z0,C,VAC,TIP,SEM,
     &VSINT,R,S,DELV,DELR0,DELS0,DELP,NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,
     &NS,NP,BIAS,IWRIT,ITMAX,EP,IPMAX,Pot0,IERR,IINIT,MIRROR,EPSIL)

DIMENSION VAC(2,NRDIM,NVDIM,NPDIM),SEM(2,NRDIM,NSDIM,NPDIM),
     &VSINT(2,NRDIM,NPDIM),R(NRDIM),DELR(NRDIM),DELV(NRDIM),S(NSDIM),
     &DELS(NSDIM),ITMAX(10),EP(10),DELXSI(NRDIM)
LOGICAL TIP(NRDIM,NVDIM,NPDIM)
```

#### Python: physics/core/poisson.py 第120-140行
```python
class PoissonSolver:
    def __init__(self, grid: Grid3D, tip: TipModel, 
                 dielectric_constant: float = 12.9):
        """
        Initialize Poisson solver.
        
        Args:
            grid: 3D computational grid (corresponds to VAC/SEM arrays)
            tip: Tip model with geometry (corresponds to TIP array)
            dielectric_constant: Semiconductor permittivity (EPSIL)
        """
        self.grid = grid                    # VAC, SEM, VSINT arrays
        self.tip = tip                      # TIP geometry
        self.dielectric_constant = dielectric_constant  # EPSIL
        
        # Grid parameters (corresponds to NR, NV, NS, NP)
        self.nr = grid.nr
        self.nv = grid.nv  
        self.ns = grid.ns
        self.np = grid.np
        
        # Coordinate arrays (corresponds to R, S, DELV arrays)
        self.r_array = grid.r_coords
        self.s_array = grid.s_coords
        self.delv_array = grid.delv
```

**狀態**: ✅ **完成** - 物件導向設計完全對應函數參數

### 2. 坐標系設定

#### Fortran: semitip3-6.1.f 第105-120行
```fortran
ETAT=1./SQRT(1.+1./SLOPE**2)
A=RAD*SLOPE**2/ETAT
SPRIME=A*ETAT
Z0=SEP-SPRIME
C=Z0/SPRIME
PI=4.*ATAN(1.)

WRITE(6,*) 'ETAT, A, Z0, C =',ETAT,A,Z0,C
WRITE(16,*) 'ETAT, A, Z0, C =',ETAT,A,Z0,C
```

#### Python: physics/core/poisson.py 第180-200行
```python
def _setup_coordinate_system(self, separation: float, 
                           tip_radius: float, slope: float) -> dict:
    """Set up hyperbolic coordinate system parameters."""
    
    # Calculate coordinate system parameters (exact Fortran match)
    etat = 1.0 / np.sqrt(1.0 + 1.0 / slope**2)
    a = tip_radius * slope**2 / etat
    sprime = a * etat
    z0 = separation - sprime
    c = z0 / sprime
    
    # Log parameters matching Fortran output
    coord_params = {
        'etat': etat, 'a': a, 'z0': z0, 'c': c
    }
    
    # This corresponds to Fortran WRITE statements
    print(f"ETAT, A, Z0, C = {etat:.8f} {a:.8f} {z0:.8e} {c:.8e}")
    
    return coord_params
```

**狀態**: ✅ **完成** - 數學公式完全一致

### 3. 多網格求解主迴圈

#### Fortran: semitip3-6.1.f 第200-250行
```fortran
C   START LOOP ON SUCCESSIVELY DOUBLING GRID DENSITY
      DO 20 IP=1,IPMAX
         WRITE(6,*) ' '
         WRITE(16,*) ' '
         WRITE(6,*) 'NR,NS,NV,NP =',NR,NS,NV,NP
         WRITE(16,*) 'NR,NS,NV,NP =',NR,NS,NV,NP
         WRITE(6,*) 'DELR,DELS,DELV,DELP =',DELR,DELS,DELV(2),DELP
         WRITE(16,*) 'DELR,DELS,DELV,DELP =',DELR,DELS,DELV(2),DELP
         
C   SOR ITERATION LOOP
         DO 60 ITER=1,ITMAX(IP)
            ...SOR iterations...
         60 CONTINUE
         
         WRITE(6,*) 'SOLUTION #',IP
         WRITE(16,*) 'SOLUTION #',IP
         WRITE(6,*) 'NUMBER OF ITERATIONS =',ITER-1
         WRITE(16,*) 'NUMBER OF ITERATIONS =',ITER-1
```

#### Python: physics/core/poisson.py 第250-300行
```python
def _multigrid_solve(self, parameters: PoissonSolverParameters,
                    charge_tables: ChargeDensityTables) -> Tuple[np.ndarray, dict]:
    """Multi-grid solving loop corresponding to DO 20 IP=1,IPMAX."""
    
    convergence_info = {'solutions': []}
    
    # Multi-grid refinement loop (corresponds to IP=1,IPMAX)
    for grid_level in range(self.max_grid_levels):
        
        # Log grid parameters (matches Fortran WRITE statements)
        print(f"NR,NS,NV,NP = {self.nr} {self.ns} {self.nv} {self.np}")
        print(f"DELR,DELS,DELV,DELP = {self.delr:.5f} {self.dels:.5f} "
              f"{self.delv:.5f} {self.delp:.5f}")
        
        # SOR iteration (corresponds to DO 60 ITER=1,ITMAX(IP))
        iteration_info = self._sor_iteration(parameters, charge_tables)
        
        # Log solution info (matches Fortran)  
        print(f"SOLUTION # {grid_level + 1}")
        print(f"NUMBER OF ITERATIONS = {iteration_info['iterations']}")
        
        # Calculate band bending at midpoint
        band_bending = self._calculate_band_bending()
        print(f"BAND BENDING AT MIDPOINT = {band_bending:.8f}")
        
        convergence_info['solutions'].append({
            'grid_level': grid_level + 1,
            'iterations': iteration_info['iterations'],
            'band_bending': band_bending
        })
        
        # Refine grid for next level (if not final)
        if grid_level < self.max_grid_levels - 1:
            self._refine_grid()
    
    return self.potential_3d, convergence_info
```

**狀態**: ✅ **完成** - 多網格邏輯完全對應

### 4. SOR 迭代核心

#### Fortran: semitip3-6.1.f 第300-400行
```fortran
C   SOR ITERATION LOOP
DO 60 ITER=1,ITMAX(IP)
   C   SOLVE IN VACUUM REGION
   DO 30 K=1,NP
   DO 30 I=1,NR
   DO 30 J=1,NV
      IF(.NOT.TIP(I,J,K)) THEN
         ...vacuum update equations...
      END IF
   30 CONTINUE
   
   C   SOLVE IN SEMICONDUCTOR REGION  
   DO 40 K=1,NP
   DO 40 I=1,NR
   DO 40 J=1,NS
      ...semiconductor update with GSECT...
      CALL GSECT(SEMMIN,Pot1,Pot2,DELPOT1,ITER1,DELSEMP)
   40 CONTINUE
```

#### Python: physics/core/poisson.py 第350-450行
```python
def _sor_iteration(self, parameters: PoissonSolverParameters,
                  charge_tables: ChargeDensityTables) -> dict:
    """SOR iteration corresponding to DO 60 ITER=1,ITMAX."""
    
    iteration_count = 0
    converged = False
    
    # Main SOR loop (corresponds to DO 60 ITER=1,ITMAX(IP))
    for iteration in range(parameters.max_iterations):
        old_potential = self.potential_3d.copy()
        
        # Update vacuum region (corresponds to vacuum DO loops)
        self._update_vacuum_region(parameters)
        
        # Update semiconductor region (corresponds to semiconductor DO loops)  
        self._update_semiconductor_region(parameters, charge_tables)
        
        # Update surface region (corresponds to surface updates)
        self._update_surface_region(parameters, charge_tables)
        
        iteration_count += 1
        
        # Check convergence every N iterations (Fortran style)
        if iteration % parameters.convergence_check_interval == 0:
            if self._check_convergence(old_potential, parameters):
                converged = True
                break
    
    return {
        'iterations': iteration_count,
        'converged': converged,
        'final_residual': self._calculate_residual()
    }
```

**狀態**: ⚠️ **部分完成** - 結構正確，但固定迭代數問題

### 5. 非線性求解 (GSECT 呼叫)

#### Fortran: semitip3-6.1.f 第420-450行
```fortran
C   SEMICONDUCTOR REGION WITH NONLINEAR CHARGE
TEMP1=...
TEMP2=...
DENOM=...
Pot1=TEMP1/DENOM
Pot2=TEMP2/DENOM

CALL GSECT(SEMMIN,Pot1,Pot2,DELPOT1,ITER1,DELSEMP)
SEM(1,I,J,K)=(1.-OMEGA)*SEM(1,I,J,K)+OMEGA*Pot1
```

#### Python: physics/core/poisson.py 第500-550行
```python
def _update_semiconductor_region(self, parameters: PoissonSolverParameters,
                               charge_tables: ChargeDensityTables):
    """Update semiconductor region with nonlinear charge, using GSECT."""
    
    for k in range(self.np):
        for i in range(self.nr):
            for j in range(self.ns):
                
                # Calculate linear terms (corresponds to TEMP1/TEMP2/DENOM)
                temp1, temp2, denom = self._calculate_linear_terms(i, j, k)
                pot1 = temp1 / denom
                pot2 = temp2 / denom
                
                # Define function for golden section search (SEMMIN equivalent)
                def semiconductor_objective(potential):
                    # Get charge density at this potential
                    x, y, s = self._get_coordinates(i, j, k)
                    rho = charge_tables.get_bulk_charge_density(
                        potential, x, y, s, i, j, k
                    )
                    
                    # Calculate residual (corresponds to SEMMIN function)
                    temp = temp1 - rho * PC.elementary_charge / self.dielectric_constant
                    return abs(potential - temp / denom)
                
                # Golden section search (corresponds to CALL GSECT)
                optimal_potential = golden_section_search(
                    semiconductor_objective, 
                    pot1, pot2,
                    tolerance=parameters.golden_section_tolerance
                )
                
                # SOR update (corresponds to OMEGA update)
                omega = parameters.omega
                self.potential_3d[i, j, k] = (
                    (1.0 - omega) * self.potential_3d[i, j, k] + 
                    omega * optimal_potential
                )
```

**狀態**: ❌ **未完成** - GSECT 對應已實現但未整合到主迴圈

### 6. 收斂檢查

#### Fortran: semitip3-6.1.f 第600-650行
```fortran
C   CHECK FOR CONVERGENCE
IF(MOD(ITER,100).EQ.0.OR.ITER.EQ.ITMAX(IP)) THEN
   DIFF=0.
   DO 50 K=1,NP
   DO 50 I=1,NR
   DO 50 J=1,NV+NS
      DIFF=AMAX1(DIFF,ABS(VAC(2,I,J,K)-VAC(1,I,J,K)))
   50 CONTINUE
   
   IF(DIFF.LT.EP(IP)) GOTO 70
END IF
```

#### Python: physics/core/poisson.py 第600-650行
```python
def _check_convergence(self, old_potential: np.ndarray, 
                      parameters: PoissonSolverParameters) -> bool:
    """Check convergence corresponding to Fortran convergence check."""
    
    # Calculate maximum difference (corresponds to Fortran DIFF calculation)
    diff = np.max(np.abs(self.potential_3d - old_potential))
    
    # Log convergence info (matching Fortran style)
    if hasattr(self, 'iteration_count'):
        if (self.iteration_count % 100 == 0 or 
            self.iteration_count >= parameters.max_iterations):
            print(f"Iteration {self.iteration_count}: Max diff = {diff:.6e}")
    
    # Check convergence (corresponds to IF(DIFF.LT.EP(IP)))
    if diff < parameters.tolerance:
        print(f"Converged at iteration {self.iteration_count}")
        return True
    
    return False
```

**狀態**: ⚠️ **部分完成** - 邏輯正確但需要整合

## 📈 關鍵變數對應表

| Fortran 變數 | Python 對應 | 狀態 | 說明 |
|-------------|-------------|------|------|
| `VAC(2,I,J,K)` | `self.potential_3d[i,j,k]` | ✅ | 真空區電位 |
| `SEM(2,I,J,K)` | `self.potential_3d[i,j,k]` | ✅ | 半導體區電位 |
| `VSINT(2,I,K)` | `self.surface_potential[i,k]` | ✅ | 表面電位 |
| `ETAT,A,Z0,C` | `coord_params['etat']` 等 | ✅ | 坐標系參數 |
| `OMEGA` | `parameters.omega` | ✅ | SOR 鬆弛參數 |
| `ITMAX(IP)` | `parameters.max_iterations` | ⚠️ | 最大迭代數 |
| `EP(IP)` | `parameters.tolerance` | ✅ | 收斂容許誤差 |
| `DIFF` | `np.max(np.abs(diff))` | ✅ | 收斂判據 |

## 📊 關鍵函數對應

| Fortran 函數 | Python 對應 | 狀態 | 說明 |
|-------------|-------------|------|------|
| `SEMMIN` | `semiconductor_objective` | ✅ | 半導體最佳化目標函數 |
| `SURFMIN` | `surface_objective` | ✅ | 表面最佳化目標函數 |
| `GSECT` | `golden_section_search` | ✅ | 黃金分割搜尋 |
| `RHOBULK` | `charge_tables.get_bulk_charge_density` | ✅ | 體電荷密度 |
| `RHOSURF` | `charge_tables.get_surface_charge_density` | ✅ | 表面電荷密度 |

## 🔧 已修復的關鍵問題

### 1. 坐標系參數計算 ✅
精確對應 Fortran 的雙曲坐標系設定公式。

### 2. 多網格結構 ✅  
正確實現了從粗網格到細網格的求解策略。

### 3. 邊界條件設定 ✅
邊界條件邏輯與 Fortran 一致。

## ⚠️ 當前問題分析

### 1. 迭代次數固定問題 (關鍵問題)
**現象**: 每個網格層級都固定執行 200 次迭代
**原因**: 收斂檢查沒有正確整合到主迴圈
**位置**: `_sor_iteration` 方法
**Fortran 對應**: 
```fortran
DO 60 ITER=1,ITMAX(IP)
   ...iterations...
   IF(DIFF.LT.EP(IP)) GOTO 70  ! Early exit on convergence
60 CONTINUE
70 CONTINUE
```

**需要修復**: 在每次迭代中檢查收斂，而不是固定迭代數

### 2. 非線性求解未整合 (關鍵問題)
**現象**: 半導體區域更新沒有使用 GSECT 非線性求解
**原因**: `_update_semiconductor_region` 方法未被正確呼叫
**需要**: 將 GSECT 邏輯完全整合到 SOR 迭代中

### 3. Band Bending 計算錯誤
**現象**: Band bending 值過小 (~1e-5)
**應該**: 有明顯的 band bending (~0.1V 數量級)
**原因**: Poisson 求解未收斂到正確解

## 📋 修復優先順序

### 高優先級 (立即修復)
1. **整合非線性 GSECT 求解** - 對應 Fortran 第420-450行
2. **修復動態收斂檢查** - 對應 Fortran 第600-650行  
3. **調整 SOR 參數** - omega 值和收斂條件

### 中優先級
1. **優化網格細化策略**
2. **改善邊界條件處理**
3. **增加詳細的除錯輸出**

## 📊 完成度統計

| 功能模組 | 完成度 | 狀態 |
|---------|--------|------|
| 座標系設定 | 100% | ✅ |
| 網格初始化 | 95% | ✅ |
| 多網格結構 | 85% | ✅ |
| SOR 迭代框架 | 70% | ⚠️ |
| 非線性求解 | 40% | ❌ |
| 收斂檢查 | 60% | ⚠️ |
| 邊界條件 | 80% | ✅ |

**總體完成度**: **70%**

---

**更新日期**: 2025-06-06  
**下一步重點**: 修復 SOR 迭代中的非線性求解和動態收斂檢查
