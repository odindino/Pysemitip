#!/usr/bin/env python3
"""
快速測試完整的 MultInt 模擬，驗證 Pot0 修復效果
"""
import yaml
import sys
import os
import logging
from datetime import datetime

from src.simulation.multint import MultInt
from src.core.config_schema import SemitipConfig

# 設置簡單日誌
logging.basicConfig(level=logging.INFO, format='%(name)s - %(message)s')
logger = logging.getLogger(__name__)

def test_multint_with_fixes():
    """測試修復後的 MultInt 模擬"""
    print("Testing MultInt simulation with nonlinear Poisson fixes...")
    
    try:
        # 載入配置
        config_path = "data/input/examples/test/quick_test.yaml"
        print(f"Loading config from {config_path}")
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = yaml.safe_load(f)
        
        # 修改配置以加快測試
        if 'computation' not in config_data:
            config_data['computation'] = {}
        
        # 大幅減少迭代次數
        config_data['computation']['max_iterations'] = [2, 2, 2, 2]  # 只做2次 SCF
        config_data['computation']['convergence_parameters'] = [1e-2, 1e-2, 1e-2, 1e-2]  # 放寬容差
        
        # 減少網格大小
        if 'grid' not in config_data:
            config_data['grid'] = {}
        config_data['grid']['radial_points'] = 8   # 減少到 8 個徑向點
        config_data['grid']['angular_points'] = 6  # 減少到 6 個角度點
        
        config = SemitipConfig(**config_data)
        
        # 創建輸出目錄
        output_dir = f"data/output/results/quick_test_fix_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        os.makedirs(output_dir, exist_ok=True)
        
        # 初始化 MultInt
        multint = MultInt(config, output_dir)
        
        print("Running quick SCF loop test...")
        print(f"Grid size: {multint.grid.N_eta} x {multint.grid.N_nu}")
        print(f"Max SCF iterations: 2")
        
        # 運行自洽迴圈
        multint.run_self_consistent_loop()
        
        # 檢查結果
        if multint.results:
            for bias_V, result in multint.results.items():
                pot0_equivalent = result.get('adjusted_bias_V', 'N/A')
                scf_iters = result.get('scf_iterations', 0)
                tip_potential = result.get('tip_potential_V', 'N/A')
                
                print(f"\n✅ Results for bias {bias_V}V:")
                print(f"  Adjusted bias: {pot0_equivalent}V")
                print(f"  Tip potential: {tip_potential}V") 
                print(f"  SCF iterations: {scf_iters}")
                
                # 檢查電位合理性
                potential = result.get('potential_Volts')
                if potential is not None:
                    pot_min, pot_max = float(potential.min()), float(potential.max())
                    pot_range = pot_max - pot_min
                    print(f"  Potential range: [{pot_min:.3f}, {pot_max:.3f}]V (span: {pot_range:.3f}V)")
                    
                    if pot_range > 0.1 and pot_range < 20.0:
                        print(f"  ✅ Reasonable potential variation")
                    else:
                        print(f"  ⚠️  Unusual potential variation")
                
        print(f"\n✅ MultInt test completed successfully!")
        print(f"Output saved to: {output_dir}")
        return True
        
    except Exception as e:
        print(f"❌ MultInt test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_multint_with_fixes()
    print(f"\n{'='*60}")
    if success:
        print("✅ MultInt nonlinear Poisson test PASSED")
    else:
        print("❌ MultInt nonlinear Poisson test FAILED")
    print(f"{'='*60}")
    sys.exit(0 if success else 1)