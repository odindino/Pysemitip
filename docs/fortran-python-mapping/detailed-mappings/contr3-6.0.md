# 詳細映射：contr3-6.0.f ↔ visualization/contour_plots.py

## 📁 檔案資訊

**Fortran 原始檔**: `src/fortran/MultInt/contr3-6.0.f`  
**Python 對應模組**: `src/visualization/contour_plots.py::ContourPlotter`  
**映射完成度**: 15% ❌  
**優先級**: **LOW** (可視化輔助功能)

## 📝 檔案描述

### Fortran 檔案功能
CONTR3 負責生成電位等高線圖：
- 計算指定電位值的等電位線
- 在 2D 平面上繪製電位輪廓
- 輸出適合繪圖軟體的資料格式
- 支援鏡像對稱處理

### Python 檔案功能
`ContourPlotter` 類別（**部分實現**）：
- 使用 matplotlib 生成現代化等高線圖
- 支援互動式可視化
- 3D 等電位面視覺化
- 可匯出多種圖形格式

## 🔄 函數對應關係

### 主要函數映射
| Fortran | Python | 狀態 |
|---------|--------|------|
| `SUBROUTINE CONTR3(ETA1,VAC,TIP,SEM,VSINT,...)` | `ContourPlotter.plot_potential_contours()` | ⚠️ 部分實現 |
| 等電位線追蹤算法 | `ContourPlotter._trace_contour_lines()` | ❌ 未實現 |
| 鏡像對稱處理 | `ContourPlotter._apply_mirror_symmetry()` | ❌ 未實現 |
| 輸出格式化 | `ContourPlotter._format_output()` | ⚠️ 現代化實現 |

## 📊 詳細行對行映射

### A. 主子程序簽名

#### Fortran: contr3-6.0.f 第10-25行
```fortran
SUBROUTINE CONTR3(ETA1,VAC,TIP,SEM,VSINT,R,S,DELV,NRDIM,NVDIM,NSDIM,
 &NPDIM,NR,NV,NS,NP,NUMC,DELPOT,MIRROR,KPLOT1,KPLOT2)

DIMENSION VAC(NRDIM,NVDIM,NPDIM),TIP(NRDIM,NVDIM,NPDIM),
 &SEM(NRDIM,NSDIM,NPDIM),VSINT(NRDIM,NPDIM),R(NRDIM),S(NSDIM),
 &DELV(NRDIM)

LOGICAL MIRROR
```

↔

#### Python: visualization/contour_plots.py 第30-50行
```python
class ContourPlotter:
    def plot_potential_contours(self, eta1, vacuum_grid, tip_grid, 
                               semiconductor_grid, interface_grid,
                               r_coordinates, s_coordinates, z_spacing,
                               num_contours, contour_spacing, 
                               mirror_symmetry, plot_planes):
        """
        對應 Fortran CONTR3 主子程序
        
        Parameters:
        -----------
        eta1 : float
            基準電位值
        vacuum_grid, tip_grid, semiconductor_grid : ndarray
            電位網格陣列
        num_contours : int
            等高線數量 (對應 NUMC)
        contour_spacing : float
            等高線間距 (對應 DELPOT)
        mirror_symmetry : bool
            鏡像對稱標記 (對應 MIRROR)
        plot_planes : tuple
            繪圖平面索引 (對應 KPLOT1, KPLOT2)
        """
```

**狀態**: ⚠️ **部分實現** - 參數對應，但實現不完整

### B. 等電位線計算（**未完整實現**）

#### Fortran: contr3-6.0.f 第40-80行
```fortran
C   CALCULATE CONTOUR LEVELS
DO 100 IC=1,NUMC
   CONLEV=ETA1+(IC-1)*DELPOT
   
   C   TRACE CONTOUR LINES IN VACUUM REGION
   DO 200 I=1,NR-1
      DO 200 J=1,NV-1
         CALL FINDCONT(VAC,TIP,I,J,CONLEV,R,DELV,...)
200   CONTINUE

   C   TRACE CONTOUR LINES IN SEMICONDUCTOR REGION  
   DO 300 I=1,NR-1
      DO 300 J=1,NS-1
         CALL FINDCONT(SEM,NONE,I,J,CONLEV,R,S,...)
300   CONTINUE
100 CONTINUE
```

↔

#### Python: visualization/contour_plots.py 第70-110行
```python
def _calculate_contour_levels(self, eta1, num_contours, contour_spacing):
    """計算等高線電位值"""
    contour_levels = []
    for ic in range(num_contours):
        level = eta1 + ic * contour_spacing
        contour_levels.append(level)
    return np.array(contour_levels)

def _trace_contour_lines(self, grid, contour_level, coordinates):
    """追蹤等電位線 - 需要實現"""
    # TODO: 實現等電位線追蹤算法
    # 對應 Fortran 的 FINDCONT 呼叫
    pass
```

**狀態**: ❌ **未實現** - 核心等電位線追蹤算法缺失

### C. 現有的 matplotlib 實現

#### Python: 當前實現（不對應 Fortran）
```python
def plot_2d_contours(self, potential_data, x_coords, y_coords, 
                    num_levels=20, colormap='viridis'):
    """使用 matplotlib 的現代化等高線圖"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 使用 matplotlib 的內建等高線算法
    contour = ax.contour(x_coords, y_coords, potential_data, 
                        levels=num_levels, cmap=colormap)
    ax.clabel(contour, inline=True, fontsize=10)
    
    ax.set_xlabel('Position (nm)')
    ax.set_ylabel('Position (nm)')
    ax.set_title('Potential Contours')
    
    return fig, ax
```

**狀態**: ✅ **現代化實現** - 功能等效但技術不同

### D. 鏡像對稱處理（**未實現**）

#### Fortran: contr3-6.0.f 第120-150行
```fortran
C   APPLY MIRROR SYMMETRY IF REQUESTED
IF (MIRROR) THEN
   DO 400 IC=1,NUMC
      C   DUPLICATE CONTOUR LINES FOR NEGATIVE R VALUES
      DO 500 IPOINT=1,NPOINTS(IC)
         RMIRR=-RPLOT(IC,IPOINT)
         ZMIRR=ZPLOT(IC,IPOINT)
         C   ADD MIRRORED POINT TO PLOT
500   CONTINUE
400 CONTINUE
END IF
```

↔

#### Python: visualization/contour_plots.py （**需要實現**）
```python
def _apply_mirror_symmetry(self, contour_data):
    """應用鏡像對稱 - 需要實現"""
    # TODO: 實現鏡像對稱邏輯
    # 將正 r 值的等高線複製到負 r 值
    pass
```

**狀態**: ❌ **未實現** - 需要完整實現

## 🔧 實現差異分析

### 1. 算法差異
**Fortran**: 自定義等電位線追蹤算法
**Python 現狀**: 使用 matplotlib 內建 contour 函數
**影響**: 功能等效，但細節控制不同

### 2. 輸出格式
**Fortran**: 文字檔案，適合傳統繪圖軟體
**Python**: 現代化圖形，支援互動和多格式匯出

### 3. 性能考量
**Fortran**: 針對大型網格優化
**Python**: matplotlib 對中小型網格效率高

## ❌ 映射缺失項目

### 關鍵缺失功能
1. **精確等電位線追蹤**: 需要實現與 Fortran 相同的追蹤算法
2. **鏡像對稱處理**: 軸對稱系統的對稱繪圖
3. **資料格式相容**: 與 Fortran 輸出格式的相容性
4. **大網格處理**: 針對大型 3D 網格的優化

### 需要創建的函數
```python
# 需要實現的核心函數
def _find_contour_intersections(self, grid, level, i, j):
    """對應 Fortran FINDCONT 子程序"""
    pass

def _interpolate_contour_point(self, v1, v2, level, r1, r2):
    """等電位線與網格邊界的交點插值"""
    pass

def _export_fortran_format(self, contour_data, filename):
    """匯出與 Fortran 相容的資料格式"""
    pass
```

## 📋 實現計劃

### 短期目標（優先級低）
1. **基本等電位線追蹤**: 實現 2D 等電位線追蹤算法
2. **鏡像對稱**: 添加軸對稱處理功能
3. **格式相容**: 支援 Fortran 輸出格式

### 長期目標
1. **3D 等電位面**: 完整的 3D 可視化功能
2. **互動式探索**: 即時調整等高線參數
3. **高性能處理**: GPU 加速大型網格處理

## ⚠️ 注意事項

### 實現建議
- **優先級**: 由於是輔助可視化功能，建議最後實現
- **技術選擇**: 可考慮使用現代 Python 可視化庫（plotly, mayavi）
- **相容性**: 保留與 Fortran 輸出的格式相容性

### 替代方案
- 使用現有的 matplotlib contour 功能進行基本可視化
- 針對特殊需求再實現自定義追蹤算法
- 考慮使用專業科學可視化軟體（ParaView, VisIt）

---

**映射完成度**: 15% ❌  
**實現狀態**: 現代化替代方案存在，精確對應需要額外開發  
**建議優先級**: 低（輔助功能）  
**最後更新**: 2025-06-06  
**下次檢查**: 核心功能完成後考慮實現
