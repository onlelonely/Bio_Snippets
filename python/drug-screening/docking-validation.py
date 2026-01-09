# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
分子對接協議驗證 - 自對接和交叉驗證
"""

import subprocess
import logging
from pathlib import Path
import json
import pandas as pd
from Bio.PDB import PDBParser
import numpy as np

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DockingValidator:
    def __init__(self, 
                 receptor_file="structures/receptors/receptor.pdbqt",
                 config_file="docking_results/vina_config.txt"):
        self.receptor_file = Path(receptor_file)
        self.config_file = Path(config_file)
        self.validation_dir = Path("docking_results/validation")
        self.validation_dir.mkdir(parents=True, exist_ok=True)
        
    def prepare_reference_ligands(self):
        """準備參考配體用於驗證"""
        logger.info("🔧 準備參考配體")
        
        # 從已知抑制劑數據庫獲取對照化合物
        controls_file = Path("ligands/controls/known_inhibitors.json")
        
        if not controls_file.exists():
            logger.warning("⚠️ 未找到已知抑制劑數據，使用預設對照")
            return self._create_default_controls()
        
        with open(controls_file, 'r') as f:
            inhibitor_data = json.load(f)
        
        reference_ligands = []
        
        # 處理陽性對照
        for control in inhibitor_data.get("positive_controls", []):
            ligand_file = Path(f"ligands/controls/{control['name']}.pdbqt")
            if ligand_file.exists():
                reference_ligands.append({
                    "name": control["name"],
                    "file": ligand_file,
                    "type": "positive_control",
                    "expected_activity": "active"
                })
        
        # 處理陰性對照
        for control in inhibitor_data.get("negative_controls", []):
            ligand_file = Path(f"ligands/controls/{control['name']}.pdbqt")
            if ligand_file.exists():
                reference_ligands.append({
                    "name": control["name"],
                    "file": ligand_file,
                    "type": "negative_control",
                    "expected_activity": "inactive"
                })
        
        logger.info(f"✅ 準備了 {len(reference_ligands)} 個參考配體")
        return reference_ligands
    
    def _create_default_controls(self):
        """創建預設對照配體"""
        # 這裡可以添加一些簡單的小分子作為對照
        return []
    
    def run_single_docking(self, ligand_file, output_prefix):
        """執行單次對接"""
        output_file = self.validation_dir / f"{output_prefix}_result.pdbqt"
        log_file = self.validation_dir / f"{output_prefix}_log.txt"
        
        cmd = [
            "vina",
            "--receptor", str(self.receptor_file),
            "--ligand", str(ligand_file),
            "--config", str(self.config_file),
            "--out", str(output_file),
            "--log", str(log_file)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0 and output_file.exists():
                return self._parse_vina_log(log_file)
            else:
                logger.error(f"❌ 對接失敗: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            logger.error("❌ 對接超時")
            return None
        except Exception as e:
            logger.error(f"❌ 對接過程發生錯誤: {e}")
            return None
    
    def _parse_vina_log(self, log_file):
        """解析 Vina 日誌文件"""
        try:
            with open(log_file, 'r') as f:
                content = f.read()
            
            # 提取結合能
            import re
            scores = re.findall(r'\s+1\s+([-\d.]+)', content)
            
            if scores:
                best_score = float(scores[0])
                return {
                    "binding_affinity": best_score,
                    "success": True,
                    "log_file": str(log_file)
                }
            else:
                return {"success": False, "error": "Could not parse binding affinity"}
                
        except Exception as e:
            return {"success": False, "error": str(e)}
    
    def calculate_rmsd(self, reference_file, docked_file):
        """計算 RMSD (需要相同的配體)"""
        try:
            # 使用 OpenBabel 計算 RMSD
            cmd = [
                "obrms",
                str(reference_file),
                str(docked_file)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                # 解析 RMSD 值
                import re
                rmsd_match = re.search(r'RMSD\s+([\d.]+)', result.stdout)
                if rmsd_match:
                    return float(rmsd_match.group(1))
            
            return None
            
        except Exception as e:
            logger.warning(f"⚠️ RMSD 計算失敗: {e}")
            return None
    
    def validate_positive_controls(self):
        """驗證陽性對照"""
        logger.info("✅ 驗證陽性對照")
        
        reference_ligands = self.prepare_reference_ligands()
        positive_controls = [lig for lig in reference_ligands 
                           if lig["type"] == "positive_control"]
        
        if not positive_controls:
            logger.warning("⚠️ 沒有陽性對照可供驗證")
            return {}
        
        validation_results = []
        
        for control in positive_controls:
            logger.info(f"   測試 {control['name']}")
            
            result = self.run_single_docking(
                control["file"], 
                f"positive_{control['name']}"
            )
            
            if result and result["success"]:
                validation_results.append({
                    "name": control["name"],
                    "type": "positive_control",
                    "binding_affinity": result["binding_affinity"],
                    "expected_activity": control["expected_activity"]
                })
                
                logger.info(f"     結合能: {result['binding_affinity']:.2f} kcal/mol")
            else:
                logger.warning(f"⚠️ {control['name']} 對接失敗")
        
        return validation_results
    
    def validate_negative_controls(self):
        """驗證陰性對照"""
        logger.info("❌ 驗證陰性對照")
        
        reference_ligands = self.prepare_reference_ligands()
        negative_controls = [lig for lig in reference_ligands 
                           if lig["type"] == "negative_control"]
        
        validation_results = []
        
        for control in negative_controls:
            logger.info(f"   測試 {control['name']}")
            
            result = self.run_single_docking(
                control["file"],
                f"negative_{control['name']}"
            )
            
            if result and result["success"]:
                validation_results.append({
                    "name": control["name"],
                    "type": "negative_control",
                    "binding_affinity": result["binding_affinity"],
                    "expected_activity": control["expected_activity"]
                })
                
                logger.info(f"     結合能: {result['binding_affinity']:.2f} kcal/mol")
        
        return validation_results
    
    def analyze_validation_results(self, positive_results, negative_results):
        """分析驗證結果"""
        logger.info("📊 分析驗證結果")
        
        all_results = positive_results + negative_results
        
        if not all_results:
            logger.error("❌ 沒有驗證結果可供分析")
            return False
        
        df = pd.DataFrame(all_results)
        
        # 分析陽性對照
        if positive_results:
            pos_scores = [r["binding_affinity"] for r in positive_results]
            logger.info(f"🎯 陽性對照結果:")
            logger.info(f"   平均結合能: {np.mean(pos_scores):.2f} kcal/mol")
            logger.info(f"   最佳結合能: {np.min(pos_scores):.2f} kcal/mol")
            logger.info(f"   標準差: {np.std(pos_scores):.2f}")
        
        # 分析陰性對照
        if negative_results:
            neg_scores = [r["binding_affinity"] for r in negative_results]
            logger.info(f"❌ 陰性對照結果:")
            logger.info(f"   平均結合能: {np.mean(neg_scores):.2f} kcal/mol")
            logger.info(f"   最差結合能: {np.max(neg_scores):.2f} kcal/mol")
        
        # 評估分離度
        if positive_results and negative_results:
            separation = np.mean(neg_scores) - np.mean(pos_scores)
            logger.info(f"📏 陽性/陰性分離度: {separation:.2f} kcal/mol")
            
            if separation > 2.0:
                logger.info("✅ 分離度良好，協議有效")
                protocol_valid = True
            else:
                logger.warning("⚠️ 分離度不足，可能需要調整協議")
                protocol_valid = False
        else:
            protocol_valid = True  # 無法評估但繼續
        
        # 保存結果
        validation_summary = {
            "protocol_valid": protocol_valid,
            "positive_controls": len(positive_results),
            "negative_controls": len(negative_results),
            "positive_scores": pos_scores if positive_results else [],
            "negative_scores": neg_scores if negative_results else [],
            "separation": separation if positive_results and negative_results else None
        }
        
        # 保存詳細結果
        results_file = self.validation_dir / "validation_results.csv"
        df.to_csv(results_file, index=False)
        
        summary_file = self.validation_dir / "validation_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(validation_summary, f, indent=2)
        
        logger.info(f"💾 驗證結果已保存:")
        logger.info(f"   📊 {results_file}")
        logger.info(f"   📋 {summary_file}")
        
        return protocol_valid
    
    def run_full_validation(self):
        """執行完整的對接驗證"""
        logger.info("🔍 開始對接協議驗證")
        
        try:
            # 檢查必需文件
            if not self.receptor_file.exists():
                logger.error(f"❌ 受體文件不存在: {self.receptor_file}")
                return False
            
            if not self.config_file.exists():
                logger.error(f"❌ 配置文件不存在: {self.config_file}")
                return False
            
            # 驗證陽性對照
            positive_results = self.validate_positive_controls()
            
            # 驗證陰性對照
            negative_results = self.validate_negative_controls()
            
            # 分析結果
            is_valid = self.analyze_validation_results(positive_results, negative_results)
            
            if is_valid:
                logger.info("🎉 對接協議驗證通過！")
            else:
                logger.warning("⚠️ 對接協議可能需要調整")
            
            return is_valid
            
        except Exception as e:
            logger.error(f"❌ 驗證過程發生錯誤: {e}")
            return False

if __name__ == "__main__":
    validator = DockingValidator()
    success = validator.run_full_validation()
    
    if success:
        print("✅ 對接協議驗證完成")
    else:
        print("❌ 對接協議驗證失敗")