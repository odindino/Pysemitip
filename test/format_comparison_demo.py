#!/usr/bin/env python3
"""
新舊 YAML 格式對比展示程序

此程序展示：
1. 新舊 YAML 格式的對比
2. 分層結構的優勢
3. 往返轉換示例
4. 可讀性改進

作者: Pysemitip 開發團隊
日期: 2024-01-02
"""

import yaml
import sys
import json
from pathlib import Path
from typing import Dict, Any

# 設定路徑
BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

try:
    from filereader import YamlConfigReader
    from config_schema import SemitipConfig
except ImportError as e:
    print(f"模組導入失敗: {e}")
    sys.exit(1)


class FormatComparisonDemo:
    """格式對比展示器"""
    
    def __init__(self):
        self.reader = YamlConfigReader()
    
    def demonstrate_structure_comparison(self, file_path: Path):
        """展示結構對比"""
        print(f"\n{'='*80}")
        print(f"📋 YAML 結構對比展示: {file_path.name}")
        print('='*80)
        
        # 讀取新格式
        with open(file_path, 'r', encoding='utf-8') as f:
            new_format = yaml.safe_load(f)
        
        # 生成對應的舊格式
        old_format = self._convert_to_old_format(new_format)
        
        # 展示關鍵結構差異
        self._show_environment_comparison(new_format, old_format)
        self._show_tip_comparison(new_format, old_format)
        self._show_semiconductor_comparison(new_format, old_format)
        self._show_output_comparison(new_format, old_format)
    
    def _show_environment_comparison(self, new_format: Dict, old_format: Dict):
        """展示環境設定對比"""
        print("\n🌡️  環境設定對比:")
        print("-" * 50)
        
        print("新格式 (分層結構):")
        env_yaml = yaml.dump({'environment': new_format['environment']}, 
                           default_flow_style=False, allow_unicode=True)
        print(f"  {env_yaml.replace(chr(10), chr(10) + '  ')}")
        
        print("舊格式 (平面結構):")
        old_env = {
            'temperature': old_format['temperature'],
            'dielectric_constant': old_format['dielectric_constant']
        }
        old_env_yaml = yaml.dump(old_env, default_flow_style=False, allow_unicode=True)
        print(f"  {old_env_yaml.replace(chr(10), chr(10) + '  ')}")
        
        print("✨ 新格式優勢: 環境參數邏輯分組，更易理解和維護")
    
    def _show_tip_comparison(self, new_format: Dict, old_format: Dict):
        """展示探針設定對比"""
        print("\n🔍 探針設定對比:")
        print("-" * 50)
        
        print("新格式 (分層 position 結構):")
        tip_excerpt = {
            'tip': {
                'separation': new_format['tip']['separation'],
                'radius': new_format['tip']['radius'],
                'position': new_format['tip']['position']
            }
        }
        tip_yaml = yaml.dump(tip_excerpt, default_flow_style=False, allow_unicode=True)
        print(f"  {tip_yaml.replace(chr(10), chr(10) + '  ')}")
        
        print("舊格式 (平面 x_position/y_position):")
        old_tip = {
            'tip': {
                'separation': old_format['tip']['separation'],
                'radius': old_format['tip']['radius'],
                'x_position': old_format['tip']['x_position'],
                'y_position': old_format['tip']['y_position']
            }
        }
        old_tip_yaml = yaml.dump(old_tip, default_flow_style=False, allow_unicode=True)
        print(f"  {old_tip_yaml.replace(chr(10), chr(10) + '  ')}")
        
        print("✨ 新格式優勢: position 子結構更直觀，支援未來擴展 (如 z 座標)")
    
    def _show_semiconductor_comparison(self, new_format: Dict, old_format: Dict):
        """展示半導體設定對比"""
        print("\n🔬 半導體設定對比:")
        print("-" * 50)
        
        print("新格式 (邏輯分組):")
        semi_excerpt = {
            'semiconductor': {
                'regions': new_format['semiconductor']['regions'][:1],  # 只顯示第一個區域
                'electron_affinity': new_format['semiconductor']['electron_affinity']
            }
        }
        semi_yaml = yaml.dump(semi_excerpt, default_flow_style=False, allow_unicode=True)
        print(f"  {semi_yaml.replace(chr(10), chr(10) + '  ')}")
        
        print("舊格式 (平面結構):")
        old_semi = {
            'semiconductor_regions': old_format['semiconductor_regions'][:1],
            'electron_affinity': old_format['electron_affinity']
        }
        old_semi_yaml = yaml.dump(old_semi, default_flow_style=False, allow_unicode=True)
        print(f"  {old_semi_yaml.replace(chr(10), chr(10) + '  ')}")
        
        print("✨ 新格式優勢: 相關參數分組，electron_affinity 與 regions 邏輯關聯")
    
    def _show_output_comparison(self, new_format: Dict, old_format: Dict):
        """展示輸出設定對比"""
        print("\n📊 輸出設定對比:")
        print("-" * 50)
        
        print("新格式 (統一 output 塊):")
        output_yaml = yaml.dump({'output': new_format['output']}, 
                              default_flow_style=False, allow_unicode=True)
        print(f"  {output_yaml.replace(chr(10), chr(10) + '  ')}")
        
        print("舊格式 (分散的 output_* 參數):")
        old_output = {k: v for k, v in old_format.items() if k.startswith('output_') or k.startswith('num_contours') or k.startswith('contour_')}
        old_output_yaml = yaml.dump(old_output, default_flow_style=False, allow_unicode=True)
        print(f"  {old_output_yaml.replace(chr(10), chr(10) + '  ')}")
        
        print("✨ 新格式優勢: 所有輸出相關參數統一管理，避免參數散落")
    
    def demonstrate_round_trip_conversion(self, file_path: Path):
        """展示往返轉換"""
        print(f"\n{'='*80}")
        print(f"🔄 往返轉換展示: {file_path.name}")
        print('='*80)
        
        print("步驟 1: 讀取原始 YAML 檔案")
        with open(file_path, 'r', encoding='utf-8') as f:
            original_yaml = yaml.safe_load(f)
        print(f"  ✅ 原始檔案包含 {len(original_yaml)} 個頂層參數")
        
        print("\n步驟 2: 轉換為配置物件")
        config = self.reader.load_config(file_path)
        print(f"  ✅ 成功轉換為 {config.simulation_type} 配置物件")
        print(f"  📊 包含 {len(config.semiconductor_regions)} 個半導體區域")
        print(f"  📊 包含 {len(config.surface_regions)} 個表面區域")
        
        print("\n步驟 3: 轉換回 YAML 格式")
        converted_yaml = self.reader._config_to_yaml(config)
        print(f"  ✅ 成功轉換回 YAML，包含 {len(converted_yaml)} 個頂層塊")
        
        print("\n步驟 4: 關鍵數值對比")
        comparisons = [
            ("模擬類型", original_yaml['simulation_type'], config.simulation_type),
            ("溫度", original_yaml['environment']['temperature'], config.temperature),
            ("探針半徑", original_yaml['tip']['radius'], config.tip.radius),
            ("探針分離", original_yaml['tip']['separation'], config.tip.separation),
        ]
        
        for name, original, converted in comparisons:
            status = "✅" if original == converted else "❌"
            print(f"  {status} {name}: {original} → {converted}")
    
    def demonstrate_readability_improvements(self, file_path: Path):
        """展示可讀性改進"""
        print(f"\n{'='*80}")
        print(f"📖 可讀性改進展示: {file_path.name}")
        print('='*80)
        
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        print("🔍 新格式的可讀性特點:")
        print("\n1. 清晰的邏輯區塊:")
        
        # 找出主要區塊
        blocks = []
        current_block = None
        for i, line in enumerate(lines, 1):
            stripped = line.strip()
            if stripped and not stripped.startswith('#') and not stripped.startswith(' ') and ':' in stripped:
                if current_block:
                    blocks.append(current_block)
                current_block = {'name': stripped.split(':')[0], 'line': i}
        if current_block:
            blocks.append(current_block)
        
        for block in blocks:
            print(f"   📁 {block['name']} (第 {block['line']} 行)")
        
        print("\n2. 分層結構範例:")
        print("   環境設定:")
        print("   ├── 溫度")
        print("   └── 介電常數")
        print("   ")
        print("   探針配置:")
        print("   ├── 基本參數 (半徑、分離等)")
        print("   └── 位置")
        print("       ├── x")
        print("       └── y")
        print("   ")
        print("   半導體:")
        print("   ├── 區域列表")
        print("   │   ├── 摻雜濃度")
        print("   │   ├── 能隙")
        print("   │   └── 有效質量")
        print("   └── 電子親和力")
        
        print("\n3. 中文註解說明:")
        comment_lines = [line.strip() for line in lines if line.strip().startswith('#')]
        print(f"   📝 包含 {len(comment_lines)} 行中文註解")
        print("   📝 每個參數都有詳細的物理意義說明")
        print("   📝 單位標註清楚 (nm, eV, K, cm⁻³)")
        
        print("\n✨ 整體改進:")
        print("   • 從平面結構升級為階層結構")
        print("   • 邏輯相關的參數組織在一起")
        print("   • 參數命名更加語義化")
        print("   • 支援未來功能擴展")
        print("   • 維護和修改更加容易")
    
    def _convert_to_old_format(self, new_format: Dict) -> Dict:
        """將新格式轉換為舊格式用於對比"""
        old_format = {}
        
        # 基本資訊
        old_format['version'] = new_format.get('version', '1.0')
        old_format['simulation_type'] = new_format.get('simulation_type', 'MultInt')
        
        # 環境設定 (展平)
        environment = new_format.get('environment', {})
        old_format['temperature'] = environment.get('temperature', 300.0)
        old_format['dielectric_constant'] = environment.get('dielectric_constant', 12.9)
        
        # Tip 配置 (展平 position)
        tip = new_format.get('tip', {})
        old_format['tip'] = tip.copy()
        if 'position' in tip:
            position = tip['position']
            old_format['tip']['x_position'] = position.get('x', 0.0)
            old_format['tip']['y_position'] = position.get('y', 0.0)
            del old_format['tip']['position']
        
        # 半導體區域 (展平)
        semiconductor = new_format.get('semiconductor', {})
        old_format['semiconductor_regions'] = semiconductor.get('regions', [])
        old_format['electron_affinity'] = semiconductor.get('electron_affinity', 4.07)
        
        # 表面區域 (展平)
        surface = new_format.get('surface', {})
        old_format['surface_regions'] = surface.get('regions', [])
        old_format['surface_temperature_dependence'] = surface.get('temperature_dependence', False)
        
        # 輸出設定 (展平)
        output = new_format.get('output', {})
        old_format['output_basic'] = output.get('basic_output', True)
        old_format['output_contours'] = output.get('equipotential_contours', False)
        old_format['output_full_potential'] = output.get('full_potential', False)
        old_format['num_contours'] = output.get('num_contours', 6)
        old_format['contour_spacing'] = output.get('contour_spacing', 0.0)
        old_format['contour_angle'] = output.get('contour_angle', 0.0)
        
        # 特定配置
        if 'multint_specific' in new_format:
            old_format['multint_config'] = new_format['multint_specific']
        if 'multplane_specific' in new_format:
            old_format['multplane_config'] = new_format['multplane_specific']
        
        # 其他配置
        for key in ['grid', 'computation', 'voltage_scan']:
            if key in new_format:
                old_format[key] = new_format[key]
        
        return old_format


def main():
    """主展示程序"""
    print("🎯 Pysemitip YAML 結構現代化展示")
    print("=" * 80)
    print("本程序展示從平面 YAML 結構升級到階層結構的改進")
    
    demo = FormatComparisonDemo()
    
    # 測試檔案
    test_files = [
        BASE_DIR / "Import_Files" / "MultInt_config.yaml",
        BASE_DIR / "Import_Files" / "MultPlane_config.yaml"
    ]
    
    for file_path in test_files:
        if not file_path.exists():
            print(f"⚠️ 檔案不存在: {file_path}")
            continue
        
        # 展示結構對比
        demo.demonstrate_structure_comparison(file_path)
        
        # 展示往返轉換
        demo.demonstrate_round_trip_conversion(file_path)
        
        # 展示可讀性改進
        demo.demonstrate_readability_improvements(file_path)
    
    print(f"\n{'='*80}")
    print("🎉 展示完成！")
    print("新的階層 YAML 結構提供了更好的：")
    print("  • 📖 可讀性 - 邏輯分組和清晰結構")
    print("  • 🔧 可維護性 - 相關參數組織在一起")
    print("  • 🚀 擴展性 - 易於添加新功能")
    print("  • 🛡️ 相容性 - 向後相容舊格式")
    print("=" * 80)


if __name__ == "__main__":
    main()
