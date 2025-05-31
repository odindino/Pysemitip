#!/usr/bin/env python3
"""
測試 filereader.py 讀取由 fileconverter.py 產生的 YAML 檔案

這個腳本會：
1. 讀取轉換後的 YAML 檔案
2. 顯示讀取結果
3. 驗證配置是否有效
4. 顯示關鍵參數以確認轉換正確性
"""

import sys
import logging
from pathlib import Path
from filereader import YamlConfigReader, load_yaml_config

# 設定日誌格式
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_yaml_file(file_path: str):
    """測試讀取單個 YAML 檔案"""
    print(f"\n{'='*60}")
    print(f"測試檔案: {file_path}")
    print(f"{'='*60}")
    
    try:
        # 使用 filereader 載入配置
        reader = YamlConfigReader()
        config = reader.load_config(file_path)
        
        print(f"✓ 成功載入配置檔案")
        print(f"\n基本資訊:")
        print(f"  - 版本: {config.version}")
        print(f"  - 模擬類型: {config.simulation_type}")
        
        # 顯示環境參數
        if hasattr(config, 'environment') and config.environment:
            print(f"\n環境參數:")
            print(f"  - 溫度: {config.environment.temperature} K")
            print(f"  - 介電常數: {config.environment.dielectric_constant}")
        else:
            # 舊格式兼容
            print(f"\n環境參數:")
            print(f"  - 溫度: {config.temperature} K")
            print(f"  - 介電常數: {config.dielectric_constant}")
        
        # 顯示探針參數
        print(f"\n探針參數:")
        print(f"  - 錐角參數: {config.tip.shank_slope}")
        print(f"  - 探針-樣品距離: {config.tip.separation} nm")
        print(f"  - 探針半徑: {config.tip.radius} nm")
        if hasattr(config.tip, 'position'):
            print(f"  - 位置: ({config.tip.position.x}, {config.tip.position.y}) nm")
        else:
            print(f"  - 位置: ({config.tip.x_position}, {config.tip.y_position}) nm")
        print(f"  - 費米能: {config.tip.fermi_energy} eV")
        
        # 顯示半導體區域資訊
        if hasattr(config, 'semiconductor') and config.semiconductor:
            regions = config.semiconductor.regions
            print(f"\n半導體區域: {len(regions)} 個")
        else:
            regions = config.semiconductor_regions
            print(f"\n半導體區域: {len(regions)} 個")
            
        for i, region in enumerate(regions):
            print(f"\n  區域 {region.id}:")
            print(f"    - 施體濃度: {region.donor_concentration:.2e} cm^-3")
            print(f"    - 受體濃度: {region.acceptor_concentration:.2e} cm^-3")
            print(f"    - 能隙: {region.band_gap} eV")
            if hasattr(region, 'effective_mass'):
                print(f"    - 導帶有效質量: {region.effective_mass.conduction_band}")
        
        # 顯示表面區域資訊
        if hasattr(config, 'surface') and config.surface:
            surface_regions = config.surface.regions
            print(f"\n表面區域: {len(surface_regions)} 個")
        else:
            surface_regions = config.surface_regions
            print(f"\n表面區域: {len(surface_regions)} 個")
            
        for i, surface in enumerate(surface_regions):
            print(f"\n  表面區域 {surface.id}:")
            if surface.first_distribution:
                print(f"    第一分佈:")
                print(f"      - 密度: {surface.first_distribution.density:.2e} cm^-2 eV^-1")
                print(f"      - 中性能級: {surface.first_distribution.neutrality_level} eV")
        
        # 顯示網格參數
        print(f"\n網格參數:")
        print(f"  - 鏡面對稱: {config.grid.mirror_plane}")
        print(f"  - 徑向點數: {config.grid.radial_points}")
        print(f"  - 真空點數: {config.grid.vacuum_points}")
        print(f"  - 半導體點數: {config.grid.semiconductor_points}")
        
        # 顯示計算參數
        print(f"\n計算參數:")
        print(f"  - 縮放步驟數: {config.computation.scaling_steps}")
        print(f"  - 最大迭代次數: {config.computation.max_iterations}")
        print(f"  - 收斂參數: {config.computation.convergence_parameters}")
        
        # 顯示電壓掃描參數
        print(f"\n電壓掃描參數:")
        print(f"  - 點數: {config.voltage_scan.points}")
        print(f"  - 起始電壓: {config.voltage_scan.start_voltage} V")
        print(f"  - 結束電壓: {config.voltage_scan.end_voltage} V")
        
        # 顯示模擬特有參數
        if config.simulation_type == "MultInt":
            if hasattr(config, 'multint_specific') and config.multint_specific:
                multint = config.multint_specific
            elif hasattr(config, 'multint_config') and config.multint_config:
                multint = config.multint_config
            else:
                multint = None
                
            if multint:
                print(f"\nMultInt 特有參數:")
                print(f"  - 平行波向量數: {multint.parallel_wavevectors}")
                print(f"  - 能量點數: {multint.energy_points}")
                print(f"  - 擴展因子: {multint.expansion_factor}")
                
        elif config.simulation_type == "MultPlane":
            if hasattr(config, 'multplane_specific') and config.multplane_specific:
                multplane = config.multplane_specific
            elif hasattr(config, 'multplane_config') and config.multplane_config:
                multplane = config.multplane_config
            else:
                multplane = None
                
            if multplane:
                print(f"\nMultPlane 特有參數:")
                print(f"  - 真空寬度: {multplane.vacuum_width} nm")
                print(f"  - 真空間距: {multplane.vacuum_spacing} nm")
                if hasattr(multplane, 'max_energies'):
                    print(f"  - 最大能量展開:")
                    print(f"    - 輕電洞: {multplane.max_energies.get('light_hole', 'N/A')} eV")
                    print(f"    - 重電洞: {multplane.max_energies.get('heavy_hole', 'N/A')} eV")
        
        # 驗證配置
        print(f"\n驗證配置...")
        if reader.validate_config():
            print(f"✓ 配置驗證通過")
        else:
            print(f"✗ 配置驗證失敗")
            
        return True
        
    except Exception as e:
        print(f"\n✗ 讀取失敗: {str(e)}")
        logger.error(f"詳細錯誤資訊:", exc_info=True)
        return False


def main():
    """主函數"""
    print("開始測試 filereader.py 讀取轉換後的 YAML 檔案")
    
    # 測試檔案列表
    test_files = [
        "test_multint.yaml",
        "test_multplane.yaml",
        "converted_MultInt_config.yaml", 
        "converted_MultPlane_config.yaml"
    ]
    
    # 檢查並測試存在的檔案
    success_count = 0
    tested_count = 0
    
    for file_name in test_files:
        file_path = Path(file_name)
        if file_path.exists():
            tested_count += 1
            if test_yaml_file(file_name):
                success_count += 1
        else:
            print(f"\n跳過不存在的檔案: {file_name}")
    
    # 總結
    print(f"\n{'='*60}")
    print(f"測試總結:")
    print(f"  - 測試檔案數: {tested_count}")
    print(f"  - 成功數: {success_count}")
    print(f"  - 失敗數: {tested_count - success_count}")
    print(f"{'='*60}")
    
    # 返回退出碼
    return 0 if success_count == tested_count else 1


if __name__ == "__main__":
    sys.exit(main())