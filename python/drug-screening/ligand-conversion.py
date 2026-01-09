# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
批次轉換配體庫為 AutoDock Vina 格式
"""

import os
import subprocess
import multiprocessing
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import time
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class LigandConverter:
    def __init__(self, input_dir="ligands/raw", output_dir="ligands/processed"):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.stats = {
            "total_input": 0,
            "successful_conversions": 0,
            "failed_conversions": 0,
            "start_time": None,
            "end_time": None
        }
    
    def split_sdf_library(self, sdf_file, chunk_size=1000):
        """將大型 SDF 文件分割為小塊便於並行處理"""
        logger.info(f"📦 分割 SDF 文件: {sdf_file}")
        
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
        chunk_files = []
        current_chunk = []
        chunk_num = 0
        
        for i, mol in enumerate(supplier):
            if mol is None:
                continue
                
            current_chunk.append(mol)
            
            if len(current_chunk) >= chunk_size:
                # 保存當前塊
                chunk_file = self.output_dir / f"chunk_{chunk_num:04d}.sdf"
                writer = Chem.SDWriter(str(chunk_file))
                
                for chunk_mol in current_chunk:
                    writer.write(chunk_mol)
                writer.close()
                
                chunk_files.append(chunk_file)
                chunk_num += 1
                current_chunk = []
                
                logger.info(f"   創建塊 {chunk_num}: {chunk_file}")
        
        # 處理最後一塊
        if current_chunk:
            chunk_file = self.output_dir / f"chunk_{chunk_num:04d}.sdf"
            writer = Chem.SDWriter(str(chunk_file))
            
            for chunk_mol in current_chunk:
                writer.write(chunk_mol)
            writer.close()
            
            chunk_files.append(chunk_file)
        
        logger.info(f"✅ 分割完成，共創建 {len(chunk_files)} 個塊")
        return chunk_files
    
    def convert_single_molecule(self, mol_data):
        """轉換單個分子"""
        mol, mol_id, output_path = mol_data
        
        try:
            # 基本分子驗證
            if mol is None:
                return False, f"Invalid molecule {mol_id}"
            
            # 檢查分子是否符合基本標準
            mw = rdMolDescriptors.CalcExactMolWt(mol)
            if mw < 100 or mw > 600:
                return False, f"Molecular weight out of range: {mw}"
            
            # 使用 RDKit 生成 PDB 格式
            mol_with_h = Chem.AddHs(mol)
            
            # 嘗試生成 3D 坐標
            from rdkit.Chem import AllChem
            if AllChem.EmbedMolecule(mol_with_h) == 0:
                AllChem.MMFFOptimizeMolecule(mol_with_h)
            
            # 保存為臨時 PDB 文件
            temp_pdb = output_path.with_suffix('.pdb')
            Chem.MolToPDBFile(mol_with_h, str(temp_pdb))
            
            # 使用 OpenBabel 轉換為 PDBQT
            cmd = [
                "obabel",
                str(temp_pdb),
                "-O", str(output_path),
                "-p", "7.4",  # pH 7.4
                "--gen3d"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            # 清理臨時文件
            if temp_pdb.exists():
                temp_pdb.unlink()
            
            if result.returncode == 0 and output_path.exists():
                return True, f"Successfully converted {mol_id}"
            else:
                return False, f"OpenBabel failed for {mol_id}: {result.stderr}"
                
        except Exception as e:
            return False, f"Exception converting {mol_id}: {str(e)}"
    
    def process_sdf_chunk(self, chunk_file):
        """處理單個 SDF 塊"""
        logger.info(f"🔄 處理塊: {chunk_file}")
        
        supplier = Chem.SDMolSupplier(str(chunk_file), removeHs=False)
        results = []
        
        for i, mol in enumerate(supplier):
            if mol is None:
                continue
            
            # 生成輸出文件名
            mol_id = f"{chunk_file.stem}_{i:04d}"
            output_file = self.output_dir / f"{mol_id}.pdbqt"
            
            # 轉換分子
            success, message = self.convert_single_molecule((mol, mol_id, output_file))
            results.append((success, message))
            
            if success:
                self.stats["successful_conversions"] += 1
            else:
                self.stats["failed_conversions"] += 1
        
        # 清理塊文件
        chunk_file.unlink()
        
        return results
    
    def parallel_conversion(self, sdf_files, max_workers=None):
        """並行轉換配體文件"""
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), 8)
        
        logger.info(f"🚀 開始並行轉換，使用 {max_workers} 個進程")
        
        # 首先分割所有 SDF 文件
        all_chunks = []
        for sdf_file in sdf_files:
            chunks = self.split_sdf_library(sdf_file)
            all_chunks.extend(chunks)
        
        self.stats["total_input"] = len(all_chunks) * 1000  # 估算
        
        # 並行處理所有塊
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(self.process_sdf_chunk, chunk): chunk 
                      for chunk in all_chunks}
            
            for i, future in enumerate(as_completed(futures)):
                chunk = futures[future]
                try:
                    results = future.result()
                    logger.info(f"✅ 完成塊 {i+1}/{len(all_chunks)}: {chunk}")
                except Exception as e:
                    logger.error(f"❌ 處理塊失敗 {chunk}: {e}")
    
    def validate_conversions(self):
        """驗證轉換結果"""
        logger.info("🔍 驗證轉換結果...")
        
        pdbqt_files = list(self.output_dir.glob("*.pdbqt"))
        valid_files = []
        
        for pdbqt_file in pdbqt_files:
            try:
                # 檢查文件大小
                if pdbqt_file.stat().st_size < 100:
                    logger.warning(f"⚠️ 文件過小: {pdbqt_file}")
                    continue
                
                # 檢查 PDBQT 格式
                with open(pdbqt_file, 'r') as f:
                    content = f.read()
                    
                # 基本格式檢查
                if "ATOM" in content and "TORSDOF" in content:
                    valid_files.append(pdbqt_file)
                else:
                    logger.warning(f"⚠️ 格式無效: {pdbqt_file}")
                    
            except Exception as e:
                logger.warning(f"⚠️ 檢查文件失敗 {pdbqt_file}: {e}")
        
        validation_stats = {
            "total_pdbqt_files": len(pdbqt_files),
            "valid_files": len(valid_files),
            "invalid_files": len(pdbqt_files) - len(valid_files),
            "validation_rate": len(valid_files) / len(pdbqt_files) if pdbqt_files else 0
        }
        
        logger.info(f"📊 驗證統計:")
        logger.info(f"   總文件數: {validation_stats['total_pdbqt_files']}")
        logger.info(f"   有效文件: {validation_stats['valid_files']}")
        logger.info(f"   無效文件: {validation_stats['invalid_files']}")
        logger.info(f"   有效率: {validation_stats['validation_rate']:.2%}")
        
        return validation_stats
    
    def run_conversion(self):
        """執行完整的轉換流程"""
        logger.info("🚀 開始配體庫轉換")
        self.stats["start_time"] = time.time()
        
        try:
            # 查找輸入 SDF 文件
            sdf_files = list(self.input_dir.glob("*.sdf"))
            
            if not sdf_files:
                logger.error("❌ 未找到 SDF 文件")
                return False
            
            logger.info(f"📁 找到 {len(sdf_files)} 個 SDF 文件")
            
            # 執行並行轉換
            self.parallel_conversion(sdf_files)
            
            # 驗證結果
            validation_stats = self.validate_conversions()
            
            self.stats["end_time"] = time.time()
            self.stats["duration"] = self.stats["end_time"] - self.stats["start_time"]
            
            # 保存統計信息
            final_stats = {**self.stats, **validation_stats}
            
            stats_file = self.output_dir / "conversion_stats.json"
            with open(stats_file, 'w') as f:
                json.dump(final_stats, f, indent=2)
            
            logger.info("🎉 配體轉換完成!")
            logger.info(f"📊 最終統計: {final_stats}")
            logger.info(f"💾 統計信息保存至: {stats_file}")
            
            return validation_stats["validation_rate"] > 0.8
            
        except Exception as e:
            logger.error(f"❌ 轉換過程發生錯誤: {e}")
            return False

if __name__ == "__main__":
    converter = LigandConverter()
    success = converter.run_conversion()
    
    if success:
        print("✅ 配體庫轉換成功")
    else:
        print("❌ 配體庫轉換失敗，請檢查日誌")