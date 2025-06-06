# Fortran-Python 映射文件維護指南

## 📋 維護概覽

**目的**: 確保 Fortran-Python 對應關係文件與代碼實現保持同步  
**責任**: 所有參與 Pysemitip 開發的團隊成員  
**頻率**: 每次代碼變更時同步更新  

## 🔄 維護工作流程

### 1. 代碼變更觸發文件更新

#### 當修改 Python 代碼時
```bash
# 步驟 1: 識別影響範圍
git diff --name-only HEAD~1 HEAD | grep -E '\.(py)$'

# 步驟 2: 確定對應的映射文件
# 參考 docs/fortran-python-mapping/README.md 中的對應表

# 步驟 3: 更新相關文件
# 見下文具體更新指南
```

#### 當發現新的 Fortran-Python 差異時
1. 更新對應的詳細映射文件
2. 在 `reports/differences.md` 中記錄
3. 更新 `reports/translation-status.md` 中的完成度
4. 如需要，更新索引文件

### 2. 文件更新檢查清單

#### ✅ 必須更新的情況
- [ ] 新增 Python 函數或類
- [ ] 修改核心算法實現
- [ ] 發現數值精度問題
- [ ] 修復已知的實現差異
- [ ] 更改函數簽名或接口

#### ⚠️ 建議更新的情況
- [ ] 性能優化改進
- [ ] 代碼重構（不影響功能）
- [ ] 添加註釋或文檔字符串
- [ ] 修復小的數值舍入差異

#### ❌ 無需更新的情況
- 純粹的代碼格式化
- 變數重命名（不影響邏輯）
- 日誌輸出格式調整
- 測試代碼修改

## 📝 具體更新指南

### 1. 詳細映射文件更新

#### 文件位置
```
docs/fortran-python-mapping/detailed-mappings/
├── MultInt3-6.4.md
├── semitip3-6.1.md
├── intcurr-6.2.md
├── semirhomult-6.0.md
├── surfrhomult-6.2.md
├── potcut3-6.0.md
├── potexpand-6.1.md
├── gsect-6.0.md
└── contr3-6.0.md
```

#### 更新模板
當修改某個 Python 函數時：

```markdown
## 函數對應: [function_name]

### Fortran 原始碼
```fortran
[更新的 Fortran 代碼片段]
```

### Python 實現
```python
[更新的 Python 代碼片段]
```

### 🔄 變更記錄
**日期**: YYYY-MM-DD  
**變更類型**: [修復/新增/優化]  
**變更描述**: [簡短描述]  
**數值影響**: [是否影響計算結果]  
**測試狀態**: [是否已驗證]
```

### 2. 索引文件更新

#### 函數索引更新 (`indices/function-index.md`)
```markdown
| Fortran 函數 | Python 對應 | 狀態 | 最後更新 |
|-------------|-------------|------|----------|
| FUNCTION_NAME | module.function_name | ✅/⚠️/❌ | YYYY-MM-DD |
```

#### 變數索引更新 (`indices/variable-index.md`)
```markdown
| Fortran 變數 | COMMON 區塊 | Python 對應 | 類型 | 最後驗證 |
|-------------|-------------|-------------|------|----------|
| VAR_NAME | /BLOCK/ | class.attr_name | type | YYYY-MM-DD |
```

### 3. 狀態報告更新

#### 完成度更新 (`reports/translation-status.md`)
每次修復問題或新增功能時：

```markdown
#### [模組名稱]
- **完成度**: [新的百分比]% (原: [舊百分比]%)
- **變更**: [YYYY-MM-DD] [簡短描述]
- **驗證狀態**: [是否重新驗證]
```

#### 驗證結果更新 (`reports/validation-results.md`)
當數值結果發生變化時：

```markdown
| 物理量 | Fortran 值 | Python 值 | 相對誤差 | 狀態 | 更新日期 |
|--------|------------|-----------|----------|------|----------|
| [量名] | [值] | [新值] | [新誤差] | ✅/❌ | YYYY-MM-DD |
```

## 🔍 質量保證檢查

### 1. 文件一致性檢查

#### 自動檢查腳本 (建議實現)
```bash
#!/bin/bash
# check_doc_consistency.sh

echo "檢查文件一致性..."

# 檢查所有映射文件是否存在
for f in MultInt3-6.4 semitip3-6.1 intcurr-6.2 semirhomult-6.0 \
         surfrhomult-6.2 potcut3-6.0 potexpand-6.1 gsect-6.0 contr3-6.0; do
    if [[ ! -f "docs/fortran-python-mapping/detailed-mappings/${f}.md" ]]; then
        echo "❌ 缺少映射文件: ${f}.md"
    else
        echo "✅ 映射文件存在: ${f}.md"
    fi
done

# 檢查索引文件
for idx in function-index variable-index; do
    if [[ ! -f "docs/fortran-python-mapping/indices/${idx}.md" ]]; then
        echo "❌ 缺少索引文件: ${idx}.md"
    else
        echo "✅ 索引文件存在: ${idx}.md"
    fi
done
```

#### 手動檢查項目
- [ ] 所有詳細映射文件更新日期是否與代碼變更匹配
- [ ] 函數和變數索引是否包含新增項目
- [ ] 完成度統計是否反映實際情況
- [ ] 已知問題列表是否及時更新

### 2. 數值驗證檢查

#### 當核心算法變更時
```python
# validation_check.py
def validate_core_modules():
    """驗證核心模組的數值一致性"""
    modules_to_check = [
        'physics.core.charge_density',
        'physics.materials.surface_states', 
        'physics.core.potential',
        'physics.core.poisson'
    ]
    
    for module in modules_to_check:
        result = run_numerical_validation(module)
        update_validation_report(module, result)
```

## 📊 維護統計和追蹤

### 1. 維護指標

#### 每月維護報告
```markdown
# 維護統計報告 - YYYY年MM月

## 文件更新統計
- 詳細映射文件更新: X 次
- 索引文件更新: X 次
- 驗證報告更新: X 次
- 新增問題記錄: X 個
- 解決問題記錄: X 個

## 質量指標
- 文件一致性: XX%
- 數值驗證通過率: XX%
- 過期文件數量: X 個

## 下月計畫
- [ ] 待更新項目 1
- [ ] 待更新項目 2
```

### 2. Git 提交規範

#### 提交訊息格式
```
type(scope): description

[body]

Docs-Updated: 
- detailed-mappings/[filename].md
- indices/[index-type]-index.md
- reports/[report-name].md
```

#### 範例
```
fix(poisson): implement dynamic convergence criteria

Fixed the Poisson solver to use dynamic iteration stopping
instead of fixed 200 iterations, matching Fortran behavior.

Docs-Updated:
- detailed-mappings/semitip3-6.1.md
- reports/translation-status.md
- reports/validation-results.md
```

## 🚨 緊急維護程序

### 當發現重大數值差異時
1. **立即行動** (1小時內)
   - 在 `reports/differences.md` 中記錄問題
   - 標記相關模組狀態為 ⚠️ 或 ❌
   - 通知團隊成員

2. **詳細分析** (24小時內)
   - 更新詳細映射文件，說明差異原因
   - 在驗證報告中記錄具體數值
   - 提供修復時間估計

3. **問題解決** (根據嚴重程度)
   - Critical: 立即修復
   - High: 1週內修復
   - Medium: 2週內修復
   - Low: 下次維護週期修復

### 當代碼結構大幅變更時
1. **規劃階段**
   - 評估對映射文件的影響範圍
   - 制定文件更新計畫
   - 預留足夠的文件更新時間

2. **實施階段**
   - 優先更新影響的詳細映射文件
   - 同步更新索引系統
   - 重新運行數值驗證

3. **驗證階段**
   - 檢查所有相關文件一致性
   - 確認數值結果無回歸
   - 更新完成度統計

## 📚 培訓和知識傳承

### 新團隊成員培訓
1. **文件系統介紹** (30分鐘)
   - 文件組織結構
   - 查找和導航方法
   - 更新工作流程

2. **實際操作練習** (1小時)
   - 修改一個簡單函數
   - 更新對應的映射文件
   - 運行驗證檢查

3. **質量標準培訓** (30分鐘)
   - 數值精度要求
   - 文件更新標準
   - 常見錯誤避免

### 知識文件化
- 維護過程中的常見問題和解決方案
- 數值差異的典型原因和修復方法
- 文件更新的最佳實踐案例

## 🔮 未來改進方向

### 自動化改進
1. **自動文件一致性檢查**: CI/CD 整合
2. **數值回歸測試**: 自動運行和報告
3. **文件更新提醒**: 基於代碼變更的自動通知

### 工具改進
1. **文件生成工具**: 從代碼自動生成基本映射框架
2. **差異檢測工具**: 自動識別實現差異
3. **統計分析工具**: 維護指標的自動化統計

### 流程優化
1. **簡化更新流程**: 減少手動工作量
2. **提高檢查效率**: 更快的質量保證流程
3. **改進協作機制**: 團隊間的更好協調

---

**維護原則**: 「代碼即文件，文件即代碼」- 保持兩者的完全同步是項目成功的關鍵。

**責任**: 每位開發者都是文件維護者，優質的文件是團隊的共同財產。
