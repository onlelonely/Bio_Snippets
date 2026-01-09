# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
配體品質控制和過濾
"""

import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import logging
import json

logger = logging.getLogger(__name__)

class LigandValidator:
    def __init__(self, ligand_dir="ligands/processed"):
        self.ligand_dir = Path(ligand_dir)
        self.results = []
        
    def calculate_drug_properties(self, mol):
        """計算藥物化學性質"""
        try:
            properties = {
                "molecular_weight": Descriptors.MolWt(mol),
                "logp": Crippen.MolLogP(mol),
                "hbd": rdMolDescriptors.CalcNumHBD(mol),
                "hba": rdMolDescriptors.CalcNumHBA(mol),
                "tpsa": rdMolDescriptors.CalcTPSA(mol),
                "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
                "heavy_atoms": mol.GetNumHeavyAtoms()
            }
            
            # Lipinski's Rule of Five 違反次數
            properties["lipinski_violations"] = sum([
                properties["molecular_weight"] > 500,
                properties["logp"] > 5,
                properties["hbd"] > 5,
                properties["hba"] > 10
            ])
            
            # QED (Drug-likeness)
            from rdkit.Chem import QED
            properties["qed"] = QED.qed(mol)
            
            return properties
            
        except Exception as e:
            logger.warning(f"⚠️ 計算性質失敗: {e}")
            return None
    
    def check_pains(self, mol):
        """檢查 PAINS (Pan-Assay Interference Compounds)"""
        try:
            # 創建 PAINS 過濾器
            params = FilterCatalogParams()
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
            catalog = FilterCatalog(params)
            
            # 檢查分子
            matches = catalog.GetMatches(mol)
            
            return {
                "is_pains": len(matches) > 0,
                "pains_alerts": len(matches),
                "pains_descriptions": [match.GetDescription() for match in matches]
            }
            
        except Exception as e:
            logger.warning(f"⚠️ PAINS 檢查失敗: {e}")
            return {"is_pains": False, "pains_alerts": 0, "pains_descriptions": []}
    
    def validate_structure(self, mol):
        """驗證分子結構"""
        issues = []
        
        try:
            # 檢查分子是否有效
            if mol is None:
                issues.append("Invalid molecule")
                return issues
            
            # 檢查原子數
            if mol.GetNumAtoms() == 0:
                issues.append("No atoms")
            elif mol.GetNumAtoms() > 100:
                issues.append("Too many atoms (>100)")
            
            # 檢查連通性
            from rdkit.Chem import rdmolops
            if rdmolops.GetMolFrags(mol).__len__() > 1:
                issues.append("Multiple fragments")
            
            # 檢查價態
            try:
                Chem.SanitizeMol(mol)
            except:
                issues.append("Sanitization failed")
            
            # 檢查立體化學
            if mol.GetNumConformers() == 0:
                issues.append("No 3D coordinates")
            
        except Exception as e:
            issues.append(f"Structure validation error: {e}")
        
        return issues
    
    def process_pdbqt_file(self, pdbqt_file):
        """處理單個 PDBQT 文件"""
        try:
            # 嘗試從 PDBQT 文件讀取分子
            # 注意：RDKit 不直接支援 PDBQT，需要轉換
            
            # 首先轉換為 SDF 格式
            temp_sdf = pdbqt_file.with_suffix('.temp.sdf')
            
            import subprocess
            cmd = ["obabel", str(pdbqt_file), "-O", str(temp_sdf)]
            result = subprocess.run(cmd, capture_output=True)
            
            if result.returncode != 0 or not temp_sdf.exists():
                return None
            
            # 從 SDF 讀取分子
            mol = Chem.MolFromMolFile(str(temp_sdf), removeHs=False)
            
            # 清理臨時文件
            temp_sdf.unlink()
            
            if mol is None:
                return None
            
            # 分析分子
            analysis = {
                "file_name": pdbqt_file.name,
                "file_size": pdbqt_file.stat().st_size,
            }
            
            # 計算性質
            properties = self.calculate_drug_properties(mol)
            if properties:
                analysis.update(properties)
            
            # PAINS 檢查
            pains_result = self.check_pains(mol)
            analysis.update(pains_result)
            
            # 結構驗證
            structure_issues = self.validate_structure(mol)
            analysis["structure_issues"] = structure_issues
            analysis["has_structure_issues"] = len(structure_issues) > 0
            
            # 整體評分
            analysis["overall_score"] = self.calculate_overall_score(analysis)
            
            return analysis
            
        except Exception as e:
            logger.warning(f"⚠️ 處理文件失敗 {pdbqt_file}: {e}")
            return None
    
    def calculate_overall_score(self, analysis):
        """計算整體品質評分 (0-1)"""
        score = 1.0
        
        # Lipinski violations penalty
        if "lipinski_violations" in analysis:
            score -= analysis["lipinski_violations"] * 0.2
        
        # PAINS penalty
        if analysis.get("is_pains", False):
            score -= 0.5
        
        # Structure issues penalty
        if analysis.get("has_structure_issues", False):
            score -= 0.3
        
        # QED bonus
        if "qed" in analysis:
            score += analysis["qed"] * 0.2
        
        return max(0.0, min(1.0, score))
    
    def run_validation(self):
        """執行完整的驗證流程"""
        logger.info("🔍 開始配體品質控制")
        
        pdbqt_files = list(self.ligand_dir.glob("*.pdbqt"))
        
        if not pdbqt_files:
            logger.error("❌ 未找到 PDBQT 文件")
            return False
        
        logger.info(f"📁 找到 {len(pdbqt_files)} 個 PDBQT 文件")
        
        # 處理所有文件
        for i, pdbqt_file in enumerate(pdbqt_files):
            if i % 1000 == 0:
                logger.info(f"   處理進度: {i}/{len(pdbqt_files)}")
            
            analysis = self.process_pdbqt_file(pdbqt_file)
            if analysis:
                self.results.append(analysis)
        
        # 創建 DataFrame
        df = pd.DataFrame(self.results)
        
        if df.empty:
            logger.error("❌ 沒有有效的分析結果")
            return False
        
        # 生成統計報告
        self.generate_statistics_report(df)
        
        # 保存結果
        self.save_results(df)
        
        # 過濾高品質分子
        self.filter_high_quality_ligands(df)
        
        return True
    
    def generate_statistics_report(self, df):
        """生成統計報告"""
        logger.info("📊 生成統計報告")
        
        stats = {
            "total_ligands": len(df),
            "lipinski_compliant": len(df[df["lipinski_violations"] == 0]),
            "pains_free": len(df[~df["is_pains"]]),
            "structure_valid": len(df[~df["has_structure_issues"]]),
            "high_quality": len(df[df["overall_score"] >= 0.7]),
            "medium_quality": len(df[(df["overall_score"] >= 0.5) & (df["overall_score"] < 0.7)]),
            "low_quality": len(df[df["overall_score"] < 0.5])
        }
        
        # 性質統計
        if "molecular_weight" in df.columns:
            stats["mw_stats"] = {
                "mean": df["molecular_weight"].mean(),
                "std": df["molecular_weight"].std(),
                "min": df["molecular_weight"].min(),
                "max": df["molecular_weight"].max()
            }
        
        if "qed" in df.columns:
            stats["qed_stats"] = {
                "mean": df["qed"].mean(),
                "std": df["qed"].std(),
                "min": df["qed"].min(),
                "max": df["qed"].max()
            }
        
        # 打印統計信息
        logger.info(f"   總配體數: {stats['total_ligands']}")
        logger.info(f"   Lipinski 合規: {stats['lipinski_compliant']} ({stats['lipinski_compliant']/stats['total_ligands']:.1%})")
        logger.info(f"   無 PAINS: {stats['pains_free']} ({stats['pains_free']/stats['total_ligands']:.1%})")
        logger.info(f"   結構有效: {stats['structure_valid']} ({stats['structure_valid']/stats['total_ligands']:.1%})")
        logger.info(f"   高品質: {stats['high_quality']} ({stats['high_quality']/stats['total_ligands']:.1%})")
        
        return stats
    
    def save_results(self, df):
        """保存分析結果"""
        # 保存完整結果
        output_file = self.ligand_dir / "validation_results.csv"
        df.to_csv(output_file, index=False)
        
        # 保存統計摘要
        stats_file = self.ligand_dir / "validation_statistics.json"
        stats = self.generate_statistics_report(df)
        
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2, default=str)
        
        logger.info(f"💾 結果已保存:")
        logger.info(f"   📊 {output_file}")
        logger.info(f"   📈 {stats_file}")
    
    def filter_high_quality_ligands(self, df):
        """過濾並保存高品質配體"""
        # 定義過濾標準
        high_quality = df[
            (df["overall_score"] >= 0.7) &
            (df["lipinski_violations"] <= 1) &
            (~df["is_pains"]) &
            (~df["has_structure_issues"])
        ]
        
        # 創建高品質配體目錄
        hq_dir = self.ligand_dir / "high_quality"
        hq_dir.mkdir(exist_ok=True)
        
        # 複製高品質配體文件
        for _, row in high_quality.iterrows():
            source_file = self.ligand_dir / row["file_name"]
            target_file = hq_dir / row["file_name"]
            
            if source_file.exists():
                import shutil
                shutil.copy2(source_file, target_file)
        
        # 保存高品質配體列表
        hq_list_file = hq_dir / "high_quality_ligands.csv"
        high_quality.to_csv(hq_list_file, index=False)
        
        logger.info(f"✨ 高品質配體篩選完成:")
        logger.info(f"   篩選出 {len(high_quality)} 個高品質配體")
        logger.info(f"   保存至 {hq_dir}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    validator = LigandValidator()
    validator.run_validation()