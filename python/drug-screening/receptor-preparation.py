# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
自動化蛋白質結構預處理管道
"""

import os
import subprocess
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
import logging

# 設置日誌
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class VP2Selector(Select):
    """選擇器類 - 只保留 VP2 蛋白鏈"""
    
    def __init__(self, chain_ids=['A']):
        self.chain_ids = chain_ids
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids
    
    def accept_residue(self, residue):
        # 只接受標準氨基酸殘基
        return is_aa(residue, standard=True)
    
    def accept_atom(self, atom):
        # 排除氫原子 (稍後會重新添加)
        return atom.element != 'H'

class ReceptorPreparator:
    def __init__(self, input_pdb="structures/original/6OAS.pdb"):
        self.input_pdb = Path(input_pdb)
        self.output_dir = Path("structures/processed")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def extract_vp2_chain(self, chain_id='A'):
        """提取 VP2 蛋白鏈"""
        logger.info(f"🧬 提取 VP2 鏈 {chain_id}")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("6OAS", self.input_pdb)
        
        # 使用選擇器提取特定鏈
        io = PDBIO()
        io.set_structure(structure)
        
        output_file = self.output_dir / f"6OAS_VP2_chain_{chain_id}.pdb"
        io.save(str(output_file), VP2Selector([chain_id]))
        
        logger.info(f"✅ VP2 鏈已保存至 {output_file}")
        return output_file
    
    def add_hydrogens_chimera(self, input_file):
        """使用 Chimera 添加氫原子"""
        logger.info("⚛️ 使用 Chimera 添加氫原子")
        
        output_file = input_file.with_suffix('.H.pdb')
        
        # Chimera 腳本
        chimera_script = f"""
import chimera
from chimera import runCommand

# 打開結構
runCommand("open {input_file}")

# 添加氫原子
runCommand("addh")

# 優化氫原子位置
runCommand("minimize spec #0 nsteps 100")

# 保存結果
runCommand("write format pdb #0 {output_file}")

# 退出
runCommand("stop")
"""
        
        script_file = self.output_dir / "add_hydrogens.py"
        with open(script_file, 'w') as f:
            f.write(chimera_script)
        
        try:
            # 執行 Chimera (需要安裝)
            subprocess.run([
                "chimera", "--nogui", "--script", str(script_file)
            ], check=True, capture_output=True)
            
            logger.info(f"✅ 氫原子已添加: {output_file}")
            return output_file
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("⚠️ Chimera 不可用，跳過氫原子添加")
            return input_file
    
    def convert_to_pdbqt(self, input_file):
        """轉換為 AutoDock PDBQT 格式"""
        logger.info("🔄 轉換為 PDBQT 格式")
        
        output_file = input_file.with_suffix('.pdbqt')
        
        try:
            # 使用 AutoDock Tools 的 prepare_receptor4.py
            cmd = [
                "python", "-m", "AutoDockTools.Utilities24.prepare_receptor4",
                "-R", str(input_file),
                "-o", str(output_file),
                "-A", "hydrogens",  # 保留極性氫
                "-U", "nphs"        # 移除非極性氫
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info(f"✅ PDBQT 文件已創建: {output_file}")
                return output_file
            else:
                logger.error(f"❌ PDBQT 轉換失敗: {result.stderr}")
                return None
                
        except FileNotFoundError:
            logger.warning("⚠️ AutoDock Tools 不可用，嘗試 OpenBabel")
            return self._convert_with_babel(input_file)
    
    def _convert_with_babel(self, input_file):
        """使用 OpenBabel 作為備用轉換方法"""
        output_file = input_file.with_suffix('.pdbqt')
        
        try:
            cmd = [
                "obabel", 
                str(input_file),
                "-O", str(output_file),
                "-p", "7.4"  # pH 7.4
            ]
            
            subprocess.run(cmd, check=True)
            logger.info(f"✅ OpenBabel 轉換成功: {output_file}")
            return output_file
            
        except subprocess.CalledProcessError:
            logger.error("❌ OpenBabel 轉換也失敗")
            return None
    
    def validate_structure(self, pdb_file):
        """驗證結構完整性"""
        logger.info("🔍 驗證結構完整性")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", pdb_file)
        
        stats = {
            "total_residues": 0,
            "missing_atoms": 0,
            "non_standard_residues": 0
        }
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    stats["total_residues"] += 1
                    
                    # 檢查標準氨基酸
                    if not is_aa(residue, standard=True):
                        stats["non_standard_residues"] += 1
                    
                    # 檢查缺失的主鏈原子
                    required_atoms = ['N', 'CA', 'C', 'O']
                    present_atoms = [atom.name for atom in residue]
                    
                    for atom_name in required_atoms:
                        if atom_name not in present_atoms:
                            stats["missing_atoms"] += 1
        
        # 打印統計信息
        logger.info(f"📊 結構統計:")
        logger.info(f"   總殘基數: {stats['total_residues']}")
        logger.info(f"   缺失原子: {stats['missing_atoms']}")
        logger.info(f"   非標準殘基: {stats['non_standard_residues']}")
        
        # 判斷結構是否可用
        is_valid = (stats["missing_atoms"] < stats["total_residues"] * 0.1 and
                   stats["non_standard_residues"] == 0)
        
        if is_valid:
            logger.info("✅ 結構驗證通過")
        else:
            logger.warning("⚠️ 結構可能有問題，請檢查")
        
        return is_valid, stats
    
    def process_receptor(self):
        """完整的受體處理流程"""
        logger.info("🚀 開始受體處理流程")
        
        try:
            # 步驟 1: 提取 VP2 鏈
            vp2_file = self.extract_vp2_chain()
            
            # 步驟 2: 驗證結構
            is_valid, stats = self.validate_structure(vp2_file)
            
            if not is_valid:
                logger.warning("⚠️ 結構驗證未通過，但繼續處理")
            
            # 步驟 3: 添加氫原子 (可選)
            vp2_h_file = self.add_hydrogens_chimera(vp2_file)
            
            # 步驟 4: 轉換為 PDBQT
            receptor_pdbqt = self.convert_to_pdbqt(vp2_h_file)
            
            if receptor_pdbqt:
                logger.info("🎉 受體處理完成!")
                logger.info(f"📁 最終文件: {receptor_pdbqt}")
                
                # 創建符號鏈接便於使用
                receptor_link = Path("structures/receptors/receptor.pdbqt")
                receptor_link.parent.mkdir(exist_ok=True)
                
                if receptor_link.exists():
                    receptor_link.unlink()
                receptor_link.symlink_to(receptor_pdbqt.resolve())
                
                logger.info(f"🔗 創建符號鏈接: {receptor_link}")
                return receptor_pdbqt
            else:
                logger.error("❌ 受體處理失敗")
                return None
                
        except Exception as e:
            logger.error(f"❌ 處理過程發生錯誤: {e}")
            return None

if __name__ == "__main__":
    preparator = ReceptorPreparator()
    result = preparator.process_receptor()
    
    if result:
        print("✅ 受體準備完成，可以進行下一步")
    else:
        print("❌ 受體準備失敗，請檢查錯誤信息")