# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
識別並定義對接的結合位點
"""

import json
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser
import logging

logger = logging.getLogger(__name__)

class BindingSiteIdentifier:
    def __init__(self, receptor_file="structures/processed/6OAS_VP2_chain_A.pdb"):
        self.receptor_file = Path(receptor_file)
        self.analysis_file = Path("analysis/binding_site_analysis.json")
        
    def load_interface_analysis(self):
        """載入之前的界面分析結果"""
        if self.analysis_file.exists():
            with open(self.analysis_file, 'r') as f:
                return json.load(f)
        else:
            logger.warning("⚠️ 未找到界面分析結果，將使用預設位點")
            return None
    
    def calculate_pocket_center(self, key_residues):
        """基於關鍵殘基計算口袋中心"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", self.receptor_file)
        
        coordinates = []
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    # 檢查是否為關鍵殘基
                    res_id = f"{chain.id}:{residue.resname}{residue.id[1]}"
                    
                    if any(res_id in key_res for key_res in key_residues.keys()):
                        # 使用 CA 原子坐標
                        if 'CA' in residue:
                            coordinates.append(residue['CA'].coord)
        
        if coordinates:
            center = np.mean(coordinates, axis=0)
            return center
        else:
            logger.warning("⚠️ 未找到關鍵殘基，使用幾何中心")
            return self._calculate_geometric_center()
    
    def _calculate_geometric_center(self):
        """計算蛋白質的幾何中心"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", self.receptor_file)
        
        coordinates = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if 'CA' in residue:
                        coordinates.append(residue['CA'].coord)
        
        return np.mean(coordinates, axis=0)
    
    def define_grid_box(self, center=None, size=20.0):
        """定義對接格點盒子"""
        if center is None:
            # 嘗試從界面分析獲取中心
            analysis = self.load_interface_analysis()
            if analysis and "key_residues" in analysis:
                center = self.calculate_pocket_center(analysis["key_residues"])
            else:
                center = self._calculate_geometric_center()
        
        # 定義盒子參數
        box_config = {
            "center_x": float(center[0]),
            "center_y": float(center[1]), 
            "center_z": float(center[2]),
            "size_x": size,
            "size_y": size,
            "size_z": size
        }
        
        logger.info(f"📦 對接盒子定義:")
        logger.info(f"   中心: ({box_config['center_x']:.2f}, "
                   f"{box_config['center_y']:.2f}, "
                   f"{box_config['center_z']:.2f})")
        logger.info(f"   尺寸: {size} x {size} x {size} Å")
        
        return box_config
    
    def generate_vina_config(self, box_config, output_file="docking_results/vina_config.txt"):
        """生成 AutoDock Vina 配置文件"""
        config_content = f"""# AutoDock Vina configuration file
# Generated for CPV VP2 protein docking

# Receptor
receptor = structures/receptors/receptor.pdbqt

# Grid box definition
center_x = {box_config['center_x']:.3f}
center_y = {box_config['center_y']:.3f}
center_z = {box_config['center_z']:.3f}

size_x = {box_config['size_x']}
size_y = {box_config['size_y']}
size_z = {box_config['size_z']}

# Search parameters
exhaustiveness = 8
num_modes = 9
energy_range = 3

# Output options
out = docking_results/result.pdbqt
log = docking_results/result.log
"""
        
        # 確保輸出目錄存在
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            f.write(config_content)
        
        logger.info(f"📝 Vina 配置文件已創建: {output_path}")
        return output_path
    
    def visualize_binding_site(self, box_config):
        """生成 PyMOL 腳本用於視覺化結合位點"""
        pymol_script = f"""# PyMOL script to visualize binding site
# Load receptor
load {self.receptor_file}

# Show cartoon representation
show cartoon
color lightblue

# Create binding site selection (within 8 Å of box center)
select binding_site, br. all within 8 of ({box_config['center_x']}, {box_config['center_y']}, {box_config['center_z']})

# Highlight binding site
show sticks, binding_site
color orange, binding_site

# Create pseudoatom at box center
pseudoatom box_center, pos=[{box_config['center_x']}, {box_config['center_y']}, {box_config['center_z']}]
show spheres, box_center
color red, box_center

# Create box outline (approximate)
draw_box_outline({box_config['center_x']}, {box_config['center_y']}, {box_config['center_z']}, {box_config['size_x']})

print("Binding site visualization complete!")
print("Key residues are shown in orange")
print("Red sphere shows box center")
"""
        
        script_file = Path("analysis/visualizations/binding_site.pml")
        script_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(script_file, 'w') as f:
            f.write(pymol_script)
        
        logger.info(f"🎨 PyMOL 視覺化腳本已創建: {script_file}")
        logger.info("   執行: pymol analysis/visualizations/binding_site.pml")
        
        return script_file
    
    def run_analysis(self):
        """執行完整的結合位點分析"""
        logger.info("🎯 開始結合位點識別")
        
        try:
            # 定義格點盒子
            box_config = self.define_grid_box()
            
            # 生成 Vina 配置
            config_file = self.generate_vina_config(box_config)
            
            # 生成視覺化腳本
            vis_script = self.visualize_binding_site(box_config)
            
            # 保存結果
            results = {
                "binding_site_analysis": {
                    "box_center": [box_config['center_x'], box_config['center_y'], box_config['center_z']],
                    "box_size": [box_config['size_x'], box_config['size_y'], box_config['size_z']],
                    "config_file": str(config_file),
                    "visualization_script": str(vis_script)
                }
            }
            
            results_file = Path("analysis/binding_site_definition.json")
            with open(results_file, 'w') as f:
                json.dump(results, f, indent=2)
            
            logger.info("✅ 結合位點分析完成")
            logger.info(f"📁 結果保存至: {results_file}")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ 結合位點分析失敗: {e}")
            return False

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    identifier = BindingSiteIdentifier()
    identifier.run_analysis()