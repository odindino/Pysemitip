# Pysemitip 專案 - 第一階段：MultInt 架構分析

## 文檔說明

本文檔為 Pysemitip 專案第一階段架構分析的核心產出，旨在為 MultInt 主程式的 Python 移植工作提供完整的架構指導。此分析遵循「物理優先，架構為本」的核心原則，拋棄逐行翻譯的思維，專注於理解物理模型與程式結構的深層邏輯。

**分析日期**: 2025年6月11日  
**分析範圍**: MultInt3-6.4.f 及其相關子程式  
**目標**: 建立現代化 Python 架構的理論基礎

---

## 1. MultInt 程式高級執行流程

MultInt 是一個專門用於計算掃描隧道顯微鏡 (STM) 穿隧電流的 3D 有限差分求解器。其執行流程遵循「先解泊松方程求靜電勢，再解薛丁格方程算穿隧電流」的物理邏輯：

### 1.1 初始化階段
1. **讀取輸入參數** (`fort.9`) - 從輸入檔案讀取所有物理與數值參數
2. **設定半導體材料參數** - 包括能隙、摻雜濃度、有效質量等物理常數
3. **配置表面態分布** - 設定表面態密度、能量分布等界面物理參數
4. **計算費米能階** (`EFFIND`) - 確定半導體的平衡費米能階位置
5. **建構網格系統** - 設定真空、半導體的 3D 空間網格

### 1.2 主要計算迴圈（偏壓掃描）
6. **偏壓迴圈開始** - 對每個設定的偏壓值進行以下計算
7. **估算空乏層寬度** - 使用 1D 近似計算空乏層特徵尺度
8. **建構電荷密度查找表** - 預計算體電荷 (`SEMIRHO`) 和表面電荷密度 (`SURFRHO`)
9. **求解泊松方程** (`SEMITIP3`) - 疊代求解 3D 靜電勢分布直到收斂
10. **提取勢能剖面** (`POTCUT3`) - 沿針尖軸向提取一維勢能剖面
11. **計算穿隧電流** (`INTCURR`) - 求解薛丁格方程並積分得到穿隧電流

### 1.3 後處理與輸出
12. **電流分量分析** - 分離價帶、導帶電流貢獻
13. **反轉層修正** - 根據物理約束調整反轉層電流
14. **電導計算** - 計算微分電導
15. **結果輸出** - 輸出電流-電壓特性與相關物理量

---

## 2. 全域資料結構分析 (COMMON Blocks)

### 2.1 /SEMI/ - 半導體物理參數區塊
**物理意義**: 儲存半導體的本質物理特性，是材料模型的核心

**Fortran 變數定義**:
```fortran
COMMON/SEMI/TK,EGAP(NREGDIM),ED(NREGDIM),EA(NREGDIM),ACB(NREGDIM),
&AVB(NREGDIM),CD(NREGDIM),CA(NREGDIM),IDEG(NREGDIM),IINV(NREGDIM),
&DELVB(NREGDIM)
```

**Python Dataclass 建議**:
```python
@dataclass
class SemiconductorPhysics:
    temperature_kT: float  # TK - 溫度 (eV)
    band_gaps: List[float]  # EGAP - 能隙 (eV)
    donor_binding: List[float]  # ED - 施體結合能 (eV)
    acceptor_binding: List[float]  # EA - 受體結合能 (eV)
    cb_effective_mass: List[float]  # ACB - 導帶有效質量
    vb_effective_mass: List[float]  # AVB - 價帶有效質量
    donor_concentration: List[float]  # CD - 施體濃度 (cm⁻³)
    acceptor_concentration: List[float]  # CA - 受體濃度 (cm⁻³)
    degeneracy_flags: List[int]  # IDEG - 簡併性指標
    inversion_flags: List[int]  # IINV - 反轉層指標
    vb_offsets: List[float]  # DELVB - 價帶偏移 (eV)
```

### 2.2 /SURF/ - 表面態物理參數區塊
**物理意義**: 描述半導體表面的量子態分布，決定表面電荷行為

**Fortran 變數定義**:
```fortran
COMMON/SURF/ISTK,TK1,EN0(NARDIM),EN(NARDIM,2),DENS(NARDIM,2),
&FWHM(NARDIM,2),ECENT(NARDIM,2)
```

**Python Dataclass 建議**:
```python
@dataclass
class SurfaceStates:
    temp_dependent: int  # ISTK - 溫度相依性
    temperature: float  # TK1 - 溫度
    neutrality_levels: List[float]  # EN0 - 電荷中性能階
    distribution_peaks: List[List[float]]  # EN - 分布峰值能量
    state_densities: List[List[float]]  # DENS - 態密度
    gaussian_widths: List[List[float]]  # FWHM - 高斯分布寬度
    energy_centroids: List[List[float]]  # ECENT - 能量重心
```

### 2.3 /CD/ - 電荷密度查找表區塊
**物理意義**: 預計算的電荷密度數據，用於快速查表而非即時計算

**Fortran 變數定義**:
```fortran
COMMON/CD/EF,ESTART,DELE,NE,RHOBTAB(NREGDIM,NEDIM),
&RHOSTAB(NARDIM,NEDIM)
```

**Python Dataclass 建議**:
```python
@dataclass
class ChargeDensityTables:
    fermi_level: float  # EF - 費米能階 (eV)
    energy_start: float  # ESTART - 能量表起始值
    energy_step: float  # DELE - 能量步長
    num_points: int  # NE - 表格點數
    bulk_charge_table: np.ndarray  # RHOBTAB - 體電荷密度表
    surface_charge_table: np.ndarray  # RHOSTAB - 表面電荷密度表
```

### 2.4 /PROTRU/ - 針尖幾何參數區塊
**物理意義**: 描述 STM 針尖的幾何形狀，影響電場分布

**Python Dataclass 建議**:
```python
@dataclass
class TipGeometry:
    protrusion_radius: float  # RAD2 - 半球突起半徑 (nm)
```

### 2.5 /TIPPOS/ - 針尖位置參數區塊
**物理意義**: 定義針尖在空間中的位置，用於區域判斷

**Python Dataclass 建議**:
```python
@dataclass
class TipPosition:
    x_position: float  # X0 - X 方向位置 (nm)
    y_position: float  # Y0 - Y 方向位置 (nm)
```

---

## 3. Fortran 子程式功能分群

### A. 數學與數值工具 (Math & Numerical Utilities)
**rseispackALL.f** - EISPACK 線性代數程式庫，矩陣特徵值計算
**gsect-6.0.f** - 黃金分割法數值優化，用於勢能最佳化
**dgsect-6.0.f** - 黃金分割法的雙精度版本

### B. 物理參數與模型設定 (Physics & Model Setup)
**semirhomult-6.0.f** - 多區域半導體體電荷密度計算，實現費米-狄拉克統計
**surfrhomult-6.2.f** - 多區域表面電荷密度計算，處理表面態分布
**MultInt3-6.4.f** (EFFIND, ENFIND 子程式) - 費米能階與電荷中性能階搜尋

### C. 靜電勢計算 (Electrostatic Potential Calculation)
**semitip3-6.1.f** - 3D 有限差分泊松求解器，程式的心臟
**MultInt3-6.4.f** (SEMMIN, SURFMIN, RHOSURF, RHOBULK) - 電荷密度回呼函式
**MultInt3-6.4.f** (IGETREG, IGETAR) - 空間區域判斷函式

### D. 勢能處理 (Potential Processing)
**potcut3-6.0.f** - 從 3D 勢能場提取 1D 剖面，連接靜電勢與量子計算
**potexpand-6.1.f** - 勢能曲線插值與擴展，提高數值精度

### E. 穿隧電流積分 (Tunneling Current Integration)
**intcurr-6.2.f** - 薛丁格方程求解與電流積分，實現 Bardeen 理論
**MultInt3-6.4.f** (主迴圈中的電流計算部分) - 電流分量分析與修正

### F. 輔助與視覺化 (Auxiliary & Visualization)
**contr3-6.0.f** - 等勢線繪製，用於結果視覺化
**MultInt3-6.4.f** (P 函式) - 針尖幾何定義函式
**MultInt3-6.4.f** (SIG 函式) - 表面態頻譜函式
**MultInt3-6.4.f** (VBEDGE, CBEDGE) - 能帶邊界定義函式

---

## 4. 現代化 Python 模組結構提案

基於上述分析，提出以下 Python 架構，嚴格遵循關注點分離原則：

```
pysemitip/
├── core/                     # 現有的檔案處理核心
│   ├── filereader.py        # fort.9 讀取器
│   ├── fileconverter.py     # 格式轉換器
│   └── config_schema.py     # 配置模式定義
├── physics/                  # 物理模型與參數管理
│   ├── materials.py         # 半導體材料參數 (來自 B 組)
│   ├── surface_states.py    # 表面態模型 (來自 B 組)
│   └── band_structure.py    # 能帶結構計算 (來自 B 組)
├── geometry/                 # 幾何與空間管理
│   ├── grid.py              # 3D 網格系統
│   ├── tip_geometry.py      # 針尖幾何定義
│   └── regions.py           # 空間區域判斷
├── solvers/                  # 核心求解器
│   ├── poisson.py           # 泊松方程求解器 (來自 C 組)
│   ├── schrodinger.py       # 薛丁格方程求解器 (來自 E 組)
│   └── charge_density.py    # 電荷密度計算器 (來自 B 組)
├── simulation/               # 模擬狀態與流程控制
│   ├── state.py             # 全域狀態管理 (COMMON 區塊)
│   ├── workflow.py          # 模擬工作流程
│   └── convergence.py       # 收斂性判斷
├── processing/               # 後處理與分析
│   ├── potential.py         # 勢能處理 (來自 D 組)
│   ├── current.py           # 電流分析 (來自 E 組)
│   └── visualization.py     # 視覺化 (來自 F 組)
├── utils/                    # 數值與數學工具
│   ├── numerical.py         # 數值方法 (來自 A 組)
│   ├── interpolation.py     # 插值算法
│   └── optimization.py      # 優化算法
└── cli/                      # 命令列界面
    ├── multint_cli.py       # MultInt 主程式入口
    └── config_generator.py  # 配置檔案生成器
```

### 4.1 關鍵設計原則
1. **物理邏輯分離**: `physics/` 專注純物理模型，與數值實現解耦
2. **求解器獨立**: `solvers/` 中的各求解器可獨立測試與驗證
3. **狀態集中管理**: `simulation/state.py` 統一管理所有 COMMON 區塊
4. **工具函式重用**: `utils/` 提供可重用的數值算法
5. **漸進式遷移**: 保持與現有 `core/` 模組的兼容性

---

## 5. 優先級翻譯路線圖 (Translation Roadmap)

基於依賴關係分析，制定自下而上的翻譯序列：

### 第一優先級：基礎工具層 (1-2 週)
**目標**: 建立無依賴的數值積木
- **A 組**: `utils/numerical.py` - 數學函式、插值、優化算法
- **理由**: 這些是純數學工具，無物理依賴，最容易驗證正確性
- **驗證**: 與 Fortran 原始函式進行數值比較測試

### 第二優先級：物理模型層 (2-3 週)
**目標**: 建立物理參數管理與基本計算
- **B 組核心**: `physics/materials.py`, `simulation/state.py`
- **包含**: 費米-狄拉克統計、電荷密度計算、材料參數管理
- **理由**: 物理模型是後續計算的基礎，需要先建立正確的物理框架
- **驗證**: 對比 SEMIRHO 和 SURFRHO 的輸出結果

### 第三優先級：幾何與網格層 (1-2 週)
**目標**: 建立空間離散化框架
- **幾何模組**: `geometry/grid.py`, `geometry/tip_geometry.py`
- **包含**: 3D 網格生成、針尖幾何、區域判斷
- **理由**: 為求解器提供空間基礎，相對獨立且可測試
- **驗證**: 確保網格與 Fortran 版本一致

### 第四優先級：泊松求解器 (3-4 週)
**目標**: 實現靜電勢計算核心
- **C 組**: `solvers/poisson.py`
- **包含**: 3D 有限差分、疊代收斂、邊界條件
- **理由**: 這是程式的計算核心，複雜度最高
- **驗證**: 對比完整的勢能分布結果

### 第五優先級：勢能處理層 (1-2 週)
**目標**: 連接靜電勢與量子計算
- **D 組**: `processing/potential.py`
- **包含**: 勢能剖面提取、插值擴展
- **理由**: 依賴泊松求解器結果，為電流計算做準備
- **驗證**: 對比 1D 勢能曲線

### 第六優先級：薛丁格求解器 (4-5 週)
**目標**: 實現量子穿隧電流計算
- **E 組**: `solvers/schrodinger.py`, `processing/current.py`
- **包含**: 薛丁格方程積分、電流分量計算
- **理由**: 依賴所有前述模組，是最終的物理目標
- **驗證**: 對比完整的 I-V 特性曲線

### 第七優先級：整合與優化 (2-3 週)
**目標**: 系統整合與性能優化
- **整合**: `simulation/workflow.py`, `cli/multint_cli.py`
- **優化**: 性能調優、記憶體管理、並行化
- **理由**: 在功能完整後進行整體優化
- **驗證**: 端到端系統測試

### 第八優先級：現代化與API設計 (3-4 週)
**目標**: 建立現代化用戶界面
- **高級API**: 簡化的物理建模界面
- **視覺化**: 現代化結果展示
- **文檔與範例**: 完整的使用說明
- **理由**: 在核心功能穩定後，提升易用性

### 總預估時間: 17-25 週 (約 4-6 個月)

---

## 結論

本架構分析為 Pysemitip 專案提供了清晰的技術路線圖。通過深入理解 MultInt 的物理邏輯與程式結構，我們建立了一個現代化、模組化的 Python 架構。這個架構不僅保持了原始程式的物理正確性，更為後續的維護、擴展和現代化提供了堅實的基礎。

遵循此架構進行開發，將確保：
1. **物理正確性**: 每個模組都有明確的物理意義
2. **可測試性**: 自下而上的依賴結構便於單元測試
3. **可維護性**: 清晰的關注點分離降低複雜度
4. **可擴展性**: 模組化設計支持功能擴展
5. **現代化**: 充分利用 Python 生態系統優勢

此文檔將作為後續所有開發工作的指導原則，確保專案朝向既定目標穩步推進。
