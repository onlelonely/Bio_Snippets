# ---------------------------------------------
# Title: 操作型指南_CPV
# Description: From: Source/3. Efforts/Res_Drug_&_Transcriptome/待整合/操作型指南_CPV.md
# ---------------------------------------------

#!/usr/bin/env python3
"""
分析 6OAS 結構中的 TfR-VP2 結合界面
"""

import requests
import json
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from pathlib import Path

class BindingSiteAnalyzer:
    def __init__(self, pdb_id="6OAS"):
        self.pdb_id = pdb_id
        self.pdb_file = f"structures/original/{pdb_id}.pdb"
        
    def download_structure(self):
        """從 RCSB PDB 下載結構"""
        url = f"https://files.rcsb.org/download/{self.pdb_id}.pdb"
        
        print(f"📥 下載 {self.pdb_id} 結構...")
        response = requests.get(url)
        
        if response.status_code == 200:
            Path("structures/original").mkdir(parents=True, exist_ok=True)
            with open(self.pdb_file, 'w') as f:
                f.write(response.text)
            print(f"✅ 結構已保存至 {self.pdb_file}")
        else:
            raise Exception(f"❌ 下載失敗: {response.status_code}")
    
    def analyze_interface(self):
        """分析 TfR-VP2 界面殘基"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(self.pdb_id, self.pdb_file)
        
        # 定義鏈 (需要根據實際 PDB 結構調整)
        vp2_chains = ['A', 'B', 'C']  # VP2 capsid protein chains
        tfr_chains = ['D', 'E']       # Transferrin receptor chains
        
        interface_residues = []
        cutoff_distance = 5.0  # Å
        
        print(f"🔍 分析結合界面 (距離閾值: {cutoff_distance} Å)")
        
        for model in structure:
            vp2_atoms = []
            tfr_atoms = []
            
            # 收集原子座標
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if chain.id in vp2_chains:
                            vp2_atoms.append((atom, chain.id, residue))
                        elif chain.id in tfr_chains:
                            tfr_atoms.append((atom, chain.id, residue))
            
            # 計算界面殘基
            for vp2_atom, vp2_chain, vp2_res in vp2_atoms:
                for tfr_atom, tfr_chain, tfr_res in tfr_atoms:
                    distance = np.linalg.norm(
                        vp2_atom.coord - tfr_atom.coord
                    )
                    
                    if distance <= cutoff_distance:
                        interface_residues.append({
                            'vp2_chain': vp2_chain,
                            'vp2_residue': f"{vp2_res.resname}{vp2_res.id[1]}",
                            'tfr_chain': tfr_chain,
                            'tfr_residue': f"{tfr_res.resname}{tfr_res.id[1]}",
                            'distance': distance
                        })
        
        return interface_residues
    
    def identify_binding_pocket(self):
        """識別藥物結合口袋"""
        interface_data = self.analyze_interface()
        
        # 統計 VP2 界面殘基
        vp2_residues = {}
        for contact in interface_data:
            res_key = f"{contact['vp2_chain']}:{contact['vp2_residue']}"
            if res_key not in vp2_residues:
                vp2_residues[res_key] = 0
            vp2_residues[res_key] += 1
        
        # 排序並顯示關鍵殘基
        sorted_residues = sorted(vp2_residues.items(), 
                               key=lambda x: x[1], reverse=True)
        
        print("\n🎯 VP2 關鍵界面殘基 (按接觸頻率排序):")
        for residue, count in sorted_residues[:20]:
            print(f"   {residue:<15} {count:>3} 個接觸")
        
        return sorted_residues
    
    def generate_report(self):
        """生成分析報告"""
        print(f"\n📊 {self.pdb_id} 結構分析報告")
        print("=" * 50)
        
        try:
            # 下載結構 (如果不存在)
            if not Path(self.pdb_file).exists():
                self.download_structure()
            
            # 分析結合界面
            key_residues = self.identify_binding_pocket()
            
            # 保存結果
            report_file = "analysis/binding_site_analysis.json"
            Path("analysis").mkdir(exist_ok=True)
            
            report_data = {
                "pdb_id": self.pdb_id,
                "analysis_date": str(Path().cwd()),
                "key_residues": dict(key_residues[:10]),
                "total_interface_residues": len(key_residues)
            }
            
            with open(report_file, 'w') as f:
                json.dump(report_data, f, indent=2)
            
            print(f"\n💾 分析結果已保存至 {report_file}")
            
        except Exception as e:
            print(f"❌ 分析過程發生錯誤: {e}")

if __name__ == "__main__":
    analyzer = BindingSiteAnalyzer()
    analyzer.generate_report()