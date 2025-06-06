# Fortran-Python 完整對應關係文件

## 📁 文件組織結構

```
fortran-python-mapping/
├── README.md                     # 本文件 - 總覽和導航
├── analysis/                     # 程式碼分析文件
│   ├── fortran-structure.md      # Fortran 程式碼結構分析
│   ├── python-structure.md       # Python 程式碼結構分析
│   └── comparison-summary.md     # 比較摘要
├── detailed-mappings/            # 詳細對應關係
│   ├── MultInt3-6.4.md          # 主程式對應
│   ├── semitip3-6.1.md          # Poisson 求解器對應
│   ├── intcurr-6.2.md           # 電流計算對應
│   ├── semirhomult-6.0.md       # 體電荷密度對應
│   ├── surfrhomult-6.2.md       # 表面電荷密度對應
│   ├── potcut3-6.0.md           # 電位剖面對應
│   ├── potexpand-6.1.md         # 電位展開對應
│   ├── gsect-6.0.md             # 數值方法對應
│   └── contr3-6.0.md            # 輔助函數對應
├── reports/                      # 分析報告
│   ├── translation-status.md    # 翻譯狀態報告
│   ├── missing-features.md      # 缺失功能清單
│   ├── differences.md           # 已知差異報告
│   └── validation-results.md    # 驗證結果
└── index/                       # 索引和查找表
    ├── function-index.md        # 函數對應索引
    ├── variable-index.md        # 變數對應索引
    └── cross-reference.md       # 交叉引用表
```

## 🎯 使用指南

### 快速查找
- **函數對應**: 查看 `index/function-index.md`
- **變數對應**: 查看 `index/variable-index.md`
- **特定檔案**: 查看 `detailed-mappings/` 下對應的 `.md` 檔案

### 狀態檢查
- **翻譯完成度**: 查看 `reports/translation-status.md`
- **已知問題**: 查看 `reports/differences.md`
- **缺失功能**: 查看 `reports/missing-features.md`

## 📊 總覽統計

| Fortran 檔案 | Python 對應 | 狀態 | 完成度 |
|-------------|-------------|------|--------|
| MultInt3-6.4.f | simulation/multint.py | 🔄 進行中 | 80% |
| semitip3-6.1.f | physics/core/poisson.py | ⚠️ 部分完成 | 60% |
| intcurr-6.2.f | physics/core/schrodinger.py | ❌ 未完成 | 30% |
| semirhomult-6.0.f | physics/core/charge_density.py | ✅ 完成 | 95% |
| surfrhomult-6.2.f | physics/materials/surface_states.py | ✅ 完成 | 90% |
| potcut3-6.0.f | physics/core/potential.py | ✅ 完成 | 85% |
| potexpand-6.1.f | 待確認 | ❌ 未找到對應 | 0% |
| gsect-6.0.f | 分散在多個模組 | ⚠️ 部分完成 | 40% |
| contr3-6.0.f | 待確認 | ❌ 未找到對應 | 0% |

**圖例**: ✅ 完成 | ⚠️ 部分完成 | 🔄 進行中 | ❌ 未完成

## 🔍 關鍵發現

### 已完成的核心功能
1. **Tip Potential 更新機制**: 已正確實現
2. **載流子密度計算**: 與 Fortran 結果一致
3. **電荷密度函數**: 正常工作
4. **能量範圍計算**: 邏輯正確

### 需要關注的問題
1. **Poisson 求解器收斂**: 需要實現非線性求解
2. **電流計算**: 出現 NaN 值問題
3. **表面態積分**: 需要完善

### 未翻譯的功能
1. **potexpand-6.1.f**: 電位展開相關函數
2. **contr3-6.0.f**: 輔助計算函數
3. **gsect-6.0.f**: 部分數值方法

## 📝 更新記錄

- **2025-01-06**: 創建文件結構，開始 Fortran 程式碼分析
- **待更新**: 完成各檔案詳細對應關係
- **待更新**: 完成驗證和測試報告

---

**維護說明**: 每次修改 Python 程式碼時，請同步更新對應的文件以保持一致性。
