# Pysemitip Environment Setup Guide

## 使用 Conda 環境 (推薦)

### 1. 創建並啟動環境
```bash
# 從 environment.yml 創建環境
conda env create -f environment.yml

# 啟動環境
conda activate pysemitip
```

### 2. 驗證安裝
```bash
# 檢查 Python 版本
python --version

# 檢查主要套件
python -c "import yaml; print('PyYAML:', yaml.__version__)"
python -c "import numpy; print('NumPy:', numpy.__version__)"
```

## 使用 pip (替代方案)

### 1. 創建虛擬環境
```bash
# 創建虛擬環境
python -m venv pysemitip_env

# 啟動環境 (Windows)
pysemitip_env\Scripts\activate

# 啟動環境 (macOS/Linux)
source pysemitip_env/bin/activate
```

### 2. 安裝依賴
```bash
pip install -r requirements.txt
```

## 環境管理

### 更新環境
```bash
# 更新所有套件
conda env update -f environment.yml

# 或使用 pip
pip install -r requirements.txt --upgrade
```

### 匯出環境
```bash
# 匯出當前環境 (包含確切版本)
conda env export > environment_lock.yml
pip freeze > requirements_lock.txt
```

### 移除環境
```bash
conda env remove -n pysemitip
```

## 開發流程

1. 啟動環境: `conda activate pysemitip`
2. 開發程式碼
3. 執行測試: `pytest`
4. 程式碼格式化: `black .`
5. 程式碼檢查: `flake8`
6. 型別檢查: `mypy`

## 套件說明

### 核心套件
- **PyYAML**: YAML 檔案處理
- **numpy**: 數值計算
- **scipy**: 科學計算

### 開發工具
- **pytest**: 單元測試
- **black**: 程式碼格式化
- **flake8**: 程式碼風格檢查
- **mypy**: 靜態型別檢查

### 分析工具
- **matplotlib**: 資料視覺化
- **pandas**: 資料分析
- **jupyter**: 互動式開發

### 文件工具
- **sphinx**: 自動化文件生成
