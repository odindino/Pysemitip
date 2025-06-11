# Pysemitip 專案 - 第 1.5 階段：物理與計算模型深度解析## 文檔說明本文檔為 Pysemitip 專案第 1.5 階段的核心產出，旨在深入解析 MultInt 程式背後的物理模型、核心方程式，以及數值計算方法。此分析基於 `docs/analysis/Phase1_MultInt_Architecture.md` 的架構框架，通過交叉引用 Fortran 原始碼和學術文獻，確保 Python 版本在物理上與 Fortran 版本完全等效。**分析日期**: 2025年6月11日  **基礎框架**: Phase1_MultInt_Architecture.md  **分析目標**: 建立權威的物理計算模型參考文檔---## 1. 核心自洽計算迴圈 (The Self-Consistent Loop)### 1.1 物理邏輯Semitip 的核心計算基於一個根本的物理原理：**半導體中的靜電勢和電荷分佈是相互耦合的**。這種耦合表現為：1. **靜電勢決定電荷分佈**: 靜電勢 `V(r)` 影響能帶彎曲，進而決定載子濃度 `n(r)` 和 `p(r)`2. **電荷分佈決定靜電勢**: 電荷密度 `ρ(r)` 通過泊松方程決定靜電勢分佈3. **表面態電荷耦合**: 表面態的佔據依賴於表面勢，而表面電荷又影響整體勢分佈### 1.2 自洽迭代的必要性由於這種強耦合，不能單獨求解任一物理量。MultInt 採用自洽迭代策略：```fortran! 在 MultInt3-6.4.f 中的主迭代迴圈DO 590 IP=1,IPMAX   ! 1. 使用當前勢分佈計算電荷密度   ! 2. 用新的電荷密度求解泊松方程得到新勢   ! 3. 檢查收斂性   IF (convergence_met) EXIT590 CONTINUE```### 1.3 收斂條件與檢查程式在 `semitip3-6.1.f` 中實現收斂性檢查：**收斂判據**:- **勢分佈變化**: `|V_new - V_old| < tolerance`- **電荷守恆**: 總電荷平衡檢查- **最大迭代次數**: `IPMAX` 限制防止無限迭代**具體實現位置**:```fortran! 在 semitip3-6.1.f 第 590 行附近DO 590 IP=1,IPMAX   IF (IWRIT.NE.0) WRITE(6,*) 'SOLUTION #',IP   ! ... 迭代計算過程 ...   ! 收斂性檢查在網格加密和勢分佈更新後進行590 CONTINUE```**網格自適應策略**:程式採用網格加密 (grid doubling) 來提高精度：```fortran! 檢查是否需要網格加密IF (NR*2.GT.NRDIM.OR.NV*2.GT.NVDIM.OR.NS*2.GT.NSDIM.OR.     NP*2.GT.NPDIM) GO TO 500```---## 2. 靜電勢模型 (Electrostatic Potential Model)### 2.1 控制方程 (Governing Equation)MultInt 求解的核心偏微分方程是 **泊松方程**：```∇·[ε(r)∇V(r)] = -ρ(r)```其中：- `V(r)`: 靜電勢 (V)- `ε(r)`: 介電常數分佈- `ρ(r)`: 電荷密度 (cm⁻³)**在不同區域的形式**:- **真空**: `ε = ε₀`, `ρ = 0`- **半導體**: `ε = ε₀εᵣ`, `ρ = ρ_bulk + ρ_surface`- **表面**: 包含表面電荷密度的邊界條件### 2.2 邊界條件 (Boundary Conditions)程式在模擬區域邊界設定以下條件：**1. 針尖表面 (Dirichlet 邊界)**:```fortran! 針尖偏壓條件BARR(J+1) = chi + egap + (BIAS + CPot)```**2. 半導體表面 (Mixed 邊界)**:在 `semitip3-6.1.f` 中實現複雜的表面邊界條件：```fortran! 表面處的邊界條件結合真空和半導體的電場連續性STEMP = (3.*VAC(1,I,1,K)-(9./6.)*VAC(1,I,2,K)+         (1./3.)*VAC(1,I,3,K))/DELV(I)+        EPSIL*(3.75*SEM(1,I,1,K)-(5./6.)*SEM(1,I,2,K)+               0.15*SEM(1,I,3,K))/DELS0```**3. 模擬區域邊界**:- **遠場邊界**: 使用週期性或固定勢邊界條件- **對稱邊界**: 利用系統對稱性減少計算量### 2.3 數值方法 (Numerical Method)**有限差分法 (Finite-Difference Method)**:MultInt 使用三維有限差分法，在柱坐標系統中：```fortran! 在 semitip3-6.1.f 中的有限差分方程! 徑向方向二階差分! 軸向方向二階差分  ! 角向方向二階差分```**網格系統**:- **針尖幾何**: 使用雙曲坐標系統適配針尖形狀- **自適應網格**: 動態網格加密提高精度- **不等間距**: 在界面附近加密網格**數值穩定性**:使用 `GSECT` 優化方法確保迭代穩定：```fortranCALL GSECT(SEMMIN,SEMOLD,SEMNEW,DELSEM)```### 2.4 輸入與輸出**輸入**:- `ρ(r)`: 來自 `RHOBULK` 和 `RHOSURF` 的電荷密度- 幾何參數: 針尖半徑、分離距離、半導體厚度- 材料參數: 介電常數、摻雜濃度**輸出**:- `V(r)`: 三維靜電勢分佈 (`VAC`, `SEM`, `VSINT` 陣列)- `Pot0`: 表面最大能帶彎曲---## 3. 電荷密度模型 (Charge Density Model)### 3.1 總電荷密度組成總電荷密度 `ρ(r)` 由以下部分組成：```ρ(r) = ρ_bulk(r) + ρ_surface(r)```其中 **體電荷密度** `ρ_bulk` 包括：```ρ_bulk = -n(r) + p(r) + N_D⁺(r) - N_A⁻(r)```- `n(r)`: 電子濃度- `p(r)`: 電洞濃度  - `N_D⁺(r)`: 離化施體濃度- `N_A⁻(r)`: 離化受體濃度### 3.2 自由載子濃度計算**費米-狄拉克統計 (Fermi-Dirac Statistics)**:在 `semirhomult-6.0.f` 中實現：**導帶電子濃度**:```fortranFUNCTION RHOCB(IREG,EF,Pot)! 常數 C = (2/√π)*2*(m/(2πℏ²))^1.5 in eV^-1.5 cm^-3DATA C/6.815E21/! 零溫情況IF (TK.NE.0.) GO TO 200IF ((EF-EGAP(IREG)-DELVB(IREG)-Pot).LE.0.) GO TO 150RHOCB=(2.*C/3.)*SQRT((ACB(IREG)*(EF-EGAP(IREG)-DELVB(IREG)-Pot))**3)RETURN! 有限溫度情況 - 費米-狄拉克積分200 RHOCB=C*SQRT((ACB(IREG)*TK)**3)*FJINT(1,(EF-EGAP(IREG)-DELVB(IREG)-Pot)/TK)```**價帶電洞濃度**:```fortranFUNCTION RHOVB(IREG,EF,Pot)! 類似導帶，但計算電洞濃度```**費米-狄拉克積分**:```fortranREAL FUNCTION FJINT(J,ETA)! J=1 對應 1/2 積分，J=3 對應 3/2 積分! 使用梯形積分實現數值積分```### 3.3 表面態電荷計算**表面電荷模型** (基於 `surfrhomult-6.2.f`):表面態電荷通過積分表面態密度獲得：```fortranFUNCTION RHOS(IAR,EF1,DELE)RHOS=RHOS1(IAR,EF1,DELE)+RHOS2(IAR,EF1,DELE)```**積分實現**:分為兩種情況：1. **溫度依賴** (`ISTK=1`):```fortran! 積分費米函數乘以表面態密度IF (EF1.LT.EN(IAR,1)) GO TO 200100 E=E+DELE    IF (E.GT.(EF1+10.*TK)) GO TO 900    SUM=SUM+SIG(IAR,1,E)*fd(e,ef1,tk)*DELE    GO TO 100```2. **零溫近似** (`ISTK=0`):```fortran! 簡單的階梯函數積分300 IF (EF1.LT.EN(IAR,1)) GO TO 500400 E=E+DELE    IF (E.GT.EF1) GO TO 900    SUM=SUM+SIG(IAR,1,E)*DELE    GO TO 400```**表面態分佈函數**:```fortranREAL*8 FUNCTION SIG(IAR,ID,ENER)! 高斯分佈或均勻分佈表面態sig=exp(-(ENER-ECENT(IAR,ID))**2/(2.*width**2))sig=sig*dens(IAR,ID)/(SQRT(2.*pi)*width)```---## 4. 態密度模型 (Density of States Model)
### 4.1 物理模型基礎

MultInt 採用 **有效質量近似** 描述半導體能帶結構，基於以下假設：

- 拋物線能帶近似
- 各向異性有效質量
- 能谷簡併度考慮

**基本態密度公式**:

```
g(E) = (1/2π²) * (2m*/ℏ²)^(3/2) * √(E-E_edge)
```

### 4.2 能帶結構處理

**價帶分裂**:

在 `intcurr-6.2.f` 中分別處理三個價帶：

```fortran
! 輕電洞帶
CALL VBCURR1(...,AVBL,...)
! 重電洞帶  
CALL VBCURR1(...,AVBH,...)
! 分裂帶
CALL VBCURR1(...,AVBSO,...)
```

**能帶邊界定義**:

```fortran
FUNCTION VBEDGE(SZ)
! 價帶邊界隨位置變化
FUNCTION CBEDGE(SZ)  
! 導帶邊界隨位置變化
```

### 4.3 在穿隧電流計算中的作用

態密度直接影響穿隧電流積分：

```fortran
! 在 intcurr-6.2.f 中的電流積分
trans=2.*nwkdeg*(2.*wf)**2*WKFTIP/(WKSEM/EFFM)
sum=sum+trans*occdiff
```

其中 `WKSEM/EFFM` 項包含半導體態密度信息。

---

## 5. 穿隧電流計算 (Tunneling Current Calculation)

### 5.1 理論基礎

**Bardeen 轉移哈密頓方法 (Bardeen's Transfer Hamiltonian Approach)**:

穿隧電流基於 Bardeen 理論，表達式為：

```
I = (e/ℏ) ∫ |M|² ρ_tip(E) ρ_sample(E) [f_tip(E) - f_sample(E)] dE
```

其中：

- `|M|²`: 穿隧矩陣元素
- `ρ_tip(E)`, `ρ_sample(E)`: 針尖和樣品態密度
- `f_tip(E)`, `f_sample(E)`: 費米分佈函數

### 5.2 數值實現

**在 `intcurr-6.2.f` 中的核心積分迴圈**:

```fortran
! 能量積分
do 110 ie=1,ne
   ener=emax-(ie-0.5)*dele
   ! 費米函數差
   occtip=fd(ener-bias,ef,tk2)
   occsem=fd(ener,ef,tk1)  
   occdiff=occtip-occsem
   
   ! 計算穿隧波函數
   call VBwf(IMPOT,wf,wfderiv,WKSEM,ener,wkparr,sep,bias,...)
   
   ! 穿隧機率
   EPERP=ener-WKPARR**2/C
   KAPPA=SQRT(C*(BARR2(NBARR2)-EPERP))
   trans=2.*nwkdeg*(2.*wf)**2*WKFTIP/(WKSEM/EFFM)
   
   ! 累積電流
   sum=sum+trans*occdiff
110 continue
```

**平行動量積分**:

```fortran
! k 空間積分
do 120 iwky=0,nwk-1
   wky=iwky*delwk
   do 115 iwkx=0,nwk-1
      wkx=iwkx*delwk
      wkparr=sqrt(wkx**2+wky**2)
      ! 簡併度因子
      nwkdeg=8
      if (iwkx.eq.0) nwkdeg=nwkdeg/2
      if (iwky.eq.0) nwkdeg=nwkdeg/2
      if (iwkx.eq.iwky) nwkdeg=nwkdeg/2
```

### 5.3 局域態和擴展態

**擴展態電流** (Extended States):

直接通過上述積分計算，對應 bulk 能帶狀態的電流貢獻。

**局域態電流** (Localized States):

```fortran
! 尋找局域態
call VBloc(IMPOT,n,wf,wfderiv,ener,wkparr,sep,bias,...)
if (n.eq.nsav) go to 310  ! 沒有新的局域態
! 找到局域態，計算其對電流的貢獻
```

**電流分量整合**:

```fortran
! 總價帶電流
currv=currvL+currvH+currvSO
! 總電流
curr=currv+currc  ! 價帶 + 導帶
```

### 5.4 像勢效應

程式包含像勢 (Image Potential) 修正：

```fortran
! 在 potexpand-6.1.f 中
IF (IMPOT.EQ.1) THEN
   do 200 j=2,NBARR2-1
      barr2(j)=barr2(j)-1.15*lambda*(NBARR2-1.)**2/
               ((j-1.)*(float(NBARR2)-j))
200 continue
END IF
```

---

## 6. 關鍵物理參數與其來源

### 6.1 核心物理參數對照表

| 物理意義 | Fortran 變數名 | 定義位置 | 單位 | 備註 |
|---------|--------------|---------|------|------|
| 半導體介電常數 | `EPSIL` | `/SEMI/` COMMON 塊 | 相對介電常數 | 影響泊松方程求解 |
| 能隙 | `EGAP(NREGDIM)` | `/SEMI/` COMMON 塊 | eV | 多區域能隙陣列 |
| 導帶有效質量 | `ACB(NREGDIM)` | `/SEMI/` COMMON 塊 | m₀ 單位 | 決定載子濃度 |
| 價帶有效質量 | `AVB(NREGDIM)` | `/SEMI/` COMMON 塊 | m₀ 單位 | 分輕重電洞質量 |
| 重電洞有效質量 | `AVBH(NREGDIM)` | MultInt3-6.4.f | m₀ 單位 | 穿隧電流計算專用 |
| 輕電洞有效質量 | `AVBL(NREGDIM)` | MultInt3-6.4.f | m₀ 單位 | 穿隧電流計算專用 |
| 分裂帶有效質量 | `AVBSO(NREGDIM)` | MultInt3-6.4.f | m₀ 單位 | 自旋軌道分裂帶 |
| 自旋軌道分裂能 | `ESO(NREGDIM)` | MultInt3-6.4.f | eV | 價帶分裂參數 |
| 施體濃度 | `CD(NREGDIM)` | `/SEMI/` COMMON 塊 | cm⁻³ | 摻雜濃度 |
| 受體濃度 | `CA(NREGDIM)` | `/SEMI/` COMMON 塊 | cm⁻³ | 摻雜濃度 |
| 施體束縛能 | `ED(NREGDIM)` | `/SEMI/` COMMON 塊 | eV | 離化能 |
| 受體束縛能 | `EA(NREGDIM)` | `/SEMI/` COMMON 塊 | eV | 離化能 |
| 價帶偏移 | `DELVB(NREGDIM)` | `/SEMI/` COMMON 塊 | eV | 異質結構用 |
| 溫度 | `TK` | `/SEMI/` COMMON 塊 | eV | kT 以電子伏特為單位 |
| 費米能階 | `EF` | `/CD/` COMMON 塊 | eV | 相對於價帶頂 |
| 針尖費米能階 | `EFTIP` | MultInt3-6.4.f | eV | 針尖材料費米能階 |
| 電子親和力 | `CHI` | MultInt3-6.4.f | eV | 半導體電子親和力 |
| 接觸電勢 | `CPot` | MultInt3-6.4.f | eV | 針尖-樣品功函數差 |

### 6.2 表面態參數表

| 物理意義 | Fortran 變數名 | 定義位置 | 單位 | 備註 |
|---------|--------------|---------|------|------|
| 表面態密度 | `DENS(NARDIM,2)` | `/SURF/` COMMON 塊 | cm⁻²eV⁻¹ | 每個表面的態密度 |
| 表面態能階 | `EN(NARDIM,2)` | `/SURF/` COMMON 塊 | eV | 表面態特徵能量 |
| 電荷中性能級 | `EN0(NARDIM)` | `/SURF/` COMMON 塊 | eV | 綜合中性能級 |
| 高斯分佈中心 | `ECENT(NARDIM,2)` | `/SURF/` COMMON 塊 | eV | 高斯分佈表面態中心 |
| 高斯分佈寬度 | `FWHM(NARDIM,2)` | `/SURF/` COMMON 塊 | eV | 半高全寬 |
| 溫度依賴標誌 | `ISTK` | `/SURF/` COMMON 塊 | 無量綱 | 0=不計溫度,1=計溫度 |

### 6.3 幾何結構參數表

| 物理意義 | Fortran 變數名 | 定義位置 | 單位 | 備註 |
|---------|--------------|---------|------|------|
| 針尖半徑 | `RAD` | MultInt3-6.4.f | nm | 主針尖半徑 |
| 半球突起半徑 | `RAD2` | `/PROTRU/` COMMON 塊 | nm | 針尖端部半球半徑 |
| 針尖-樣品距離 | `SEP` | MultInt3-6.4.f | nm | 最近距離 |
| 針尖位置 | `X0, Y0` | `/TIPPOS/` COMMON 塊 | nm | 針尖橫向位置 |
| 針尖斜率 | `SLOPE` | MultInt3-6.4.f | 無量綱 | tan(90°-θ/2) |
| 針尖角度 | `THETA` | MultInt3-6.4.f | 度 | 針尖開口角 |
| 樣品厚度 | `W` | MultInt3-6.4.f | nm | 半導體樣品厚度 |

### 6.4 網格維度參數表

| 參數名稱 | 預設值 | 物理意義 | 程式中位置 |
|---------|--------|---------|-----------|
| `NRDIM` | 512 | 徑向網格最大維度 | PARAMETER 語句 |
| `NVDIM` | 64 | 真空軸向網格最大維度 | PARAMETER 語句 |
| `NSDIM` | 512 | 半導體軸向網格最大維度 | PARAMETER 語句 |
| `NPDIM` | 64 | 角向網格最大維度 | PARAMETER 語句 |
| `NEDIM` | 50000 | 電荷密度查找表大小 | PARAMETER 語句 |
| `NREGDIM` | 2 | 最大半導體區域數 | PARAMETER 語句 |
| `NARDIM` | 2 | 最大表面區域數 | PARAMETER 語句 |
| `NVDIM1` | NVDIM+1 | 真空勢壘陣列維度 | PARAMETER 語句 |
| `NVDIM2` | 2048 | 擴展真空陣列維度 | PARAMETER 語句 |
| `NSDIM2` | 20000 | 擴展半導體陣列維度 | PARAMETER 語句 |

### 6.5 數值計算控制參數表

| 參數名稱 | 典型值 | 物理意義 | 程式中位置 |
|---------|--------|---------|-----------|
| `IPMAX` | 3 | 網格加密最大步數 | 輸入檔案 |
| `ITMAX(IP)` | 20000,10000,5000 | 各步最大迭代次數 | 輸入檔案 |
| `EP(IP)` | 1e-3,1e-3,1e-4 | 各步收斂精度 | 輸入檔案 |
| `NWK` | 20 | 平行動量積分點數 | 輸入檔案 |
| `NEE` | 20 | 能量積分點數 | 輸入檔案 |
| `EXPANI` | 20 | 薛丁格方程積分擴展因子 | 輸入檔案 |
| `FRACZ` | 0.75 | 半導體深度積分比例 | 輸入檔案 |
| `BMOD` | 0.050 | 調制電壓 | 輸入檔案 |

### 6.6 重要物理常數表

| 常數名稱 | 數值 | 單位 | 位置 | 物理意義 |
|---------|------|------|------|---------|
| `EEP` | 1.80943E-20 | V·cm | semitip3-6.1.f | e/ε₀ × 10⁻¹⁴ cm²/nm² |
| `C` | 26.254 | nm⁻²·eV⁻¹ | intcurr-6.2.f | 2m/ℏ² 動能常數 |
| `RQUANT` | 12900 | Ω | intcurr-6.2.f | 電阻量子 h/e² |
| `EPSIL0` | 8.854185E-12 | F/m | MultInt3-6.4.f | 真空介電常數 |
| `E` | 1.60210E-19 | C | MultInt3-6.4.f | 元電荷 |
| `PI` | 4.*ATAN(1.) | 無量綱 | 各程式 | 圓周率 |

### 6.7 電荷密度查找表參數

| 參數名稱 | 物理意義 | 典型範圍 |
|---------|---------|---------|
| `ESTART` | 查找表起始能量 | 相對於費米能階 |
| `DELE` | 能量步長 | 通常為 kT/10 量級 |
| `NE` | 能量點數 | 通常為 NEDIM |
| `RHOBTAB(IREG,IE)` | 體電荷密度表 | 按區域和能量索引 |
| `RHOSTAB(IAR,IE)` | 表面電荷密度表 | 按表面區域和能量索引 |

---

## 結論

本深度解析文檔建立了 MultInt 程式物理模型與數值實現的完整對映關係。通過詳細分析自洽計算迴圈、泊松方程求解、電荷密度計算、態密度模型和穿隧電流積分，我們建立了以下關鍵理解：

### 關鍵物理洞察

1. **多尺度耦合**: 程式巧妙處理了從原子尺度（表面態）到微米尺度（半導體層）的多尺度物理
2. **自洽循環**: 靜電勢-電荷耦合的自洽處理是程式的核心，確保了物理自洽性
3. **量子輸運**: Bardeen 理論的數值實現準確捕捉了穿隧電流的量子本質
4. **數值穩定性**: 網格自適應和優化算法確保了計算的穩定性和精度

### Python 移植指導原則

1. **物理正確性優先**: 每個模組必須保持與 Fortran 版本完全相同的物理內容
2. **數值等效性**: 關鍵數值常數和積分方法必須精確複製
3. **模組化設計**: 利用 Python 的優勢實現更好的代碼組織
4. **測試驅動**: 每個物理模組都應有對應的數值驗證測試

### 參數移植要點

1. **參數層次結構**: 將 COMMON 塊轉換為 Python 類或配置對象
2. **維度參數**: 轉換為動態陣列或配置參數，避免硬編碼限制
3. **物理常數**: 建立統一的常數模組，確保數值精度
4. **輸入檔案格式**: 設計現代化的配置檔案格式（如 YAML 或 JSON）

此文檔將作為 Python 版本開發的第一性原理參考，確保移植過程中的物理準確性和數值等效性。每個表格中的參數都代表了 Semitip 程式中的關鍵物理量，必須在 Python 版本中得到完整和精確的實現。
