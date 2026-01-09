# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
大規模平行分子對接執行
"""

import os
import subprocess
import multiprocessing
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import time
import json
import shutil

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ParallelDocker:
    def __init__(self, 
                 ligand_dir="ligands/processed/high_quality",
                 receptor_file="structures/receptors/receptor.pdbqt",
                 config_file="docking_results/vina_config.txt",
                 output_dir="docking_results/individual"):
        
        self.ligand_dir = Path(ligand_dir)
        self.receptor_file = Path(receptor_file)
        self.config_file = Path(config_file)
        self.output_dir = Path(output_dir)
        self.failed_dir = Path("docking_results/failed")
        
        # 創建輸出目錄
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.failed_dir.mkdir(parents=True, exist_ok=True)
        
        # 統計信息
        self.stats = {
            "total_ligands": 0,
            "successful_dockings": 0,
            "failed_dockings": 0,
            "start_time": None,
            "end_time": None
        }
    
    def find_ligand_files(self):
        """查找所有配體文件"""
        ligand_files = list(self.ligand_dir.glob("*.pdbqt"))
        
        if not ligand_files:
            logger.error(f"❌ 在 {self.ligand_dir} 中未找到 PDBQT 文件")
            return []
        
        logger.info(f"📁 找到 {len(ligand_files)} 個配體文件")
        self.stats["total_ligands"] = len(ligand_files)
        
        return ligand_files
    
    def dock_single_ligand(self, ligand_file):
        """對接單個配體"""
        ligand_name = ligand_file.stem
        output_file = self.output_dir / f"{ligand_name}_result.pdbqt"
        log_file = self.output_dir / f"{ligand_name}_log.txt"
        
        # 如果結果已存在，跳過
        if output_file.exists() and log_file.exists():
            return self._parse_existing_result(log_file, ligand_name)
        
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
                return self._parse_docking_result(log_file, ligand_name, "success")
            else:
                # 移動失敗的文件
                self._handle_failed_docking(ligand_file, ligand_name, result.stderr)
                return self._parse_docking_result(None, ligand_name, "failed", result.stderr)
                
        except subprocess.TimeoutExpired:
            error_msg = "Docking timeout"
            self._handle_failed_docking(ligand_file, ligand_name, error_msg)
            return self._parse_docking_result(None, ligand_name, "failed", error_msg)
        except Exception as e:
            error_msg = str(e)
            self._handle_failed_docking(ligand_file, ligand_name, error_msg)
            return self._parse_docking_result(None, ligand_name, "failed", error_msg)
    
    def _parse_existing_result(self, log_file, ligand_name):
        """解析已存在的結果"""
        result = self._parse_docking_result(log_file, ligand_name, "success")
        if result["binding_affinity"] is not None:
            return result
        else:
            # 如果無法解析，視為失敗
            return self._parse_docking_result(None, ligand_name, "failed", "Could not parse existing result")
    
    def _parse_docking_result(self, log_file, ligand_name, status, error_msg=None):
        """解析對接結果"""
        result = {
            "ligand_name": ligand_name,
            "status": status,
            "binding_affinity": None,
            "error_message": error_msg
        }
        
        if status == "success" and log_file and Path(log_file).exists():
            try:
                with open(log_file, 'r') as f:
                    content = f.read()
                
                # 提取最佳結合能
                import re
                scores = re.findall(r'\s+1\s+([-\d.]+)', content)
                
                if scores:
                    result["binding_affinity"] = float(scores[0])
                else:
                    result["status"] = "failed"
                    result["error_message"] = "Could not parse binding affinity"
                    
            except Exception as e:
                result["status"] = "failed"
                result["error_message"] = f"Log parsing error: {e}"
        
        # 更新統計
        if result["status"] == "success":
            self.stats["successful_dockings"] += 1
        else:
            self.stats["failed_dockings"] += 1
        
        return result
    
    def _handle_failed_docking(self, ligand_file, ligand_name, error_msg):
        """處理失敗的對接"""
        # 複製失敗的配體文件到失敗目錄
        failed_ligand = self.failed_dir / ligand_file.name
        try:
            shutil.copy2(ligand_file, failed_ligand)
        except Exception:
            pass
        
        # 記錄錯誤
        error_log = self.failed_dir / f"{ligand_name}_error.txt"
        with open(error_log, 'w') as f:
            f.write(f"Ligand: {ligand_name}\n")
            f.write(f"Error: {error_msg}\n")
            f.write(f"Time: {time.ctime()}\n")
    
    def create_batches(self, ligand_files, batch_size=1000):
        """將配體文件分批處理"""
        batches = []
        for i in range(0, len(ligand_files), batch_size):
            batch = ligand_files[i:i+batch_size]
            batches.append(batch)
        
        logger.info(f"📦 將 {len(ligand_files)} 個配體分為 {len(batches)} 批")
        return batches
    
    def run_parallel_docking(self, max_workers=None, batch_size=1000):
        """執行平行對接"""
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), 8)
        
        logger.info(f"🚀 開始平行對接，使用 {max_workers} 個進程")
        
        # 查找配體文件
        ligand_files = self.find_ligand_files()
        if not ligand_files:
            return False
        
        # 分批處理
        batches = self.create_batches(ligand_files, batch_size)
        
        all_results = []
        self.stats["start_time"] = time.time()
        
        try:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # 提交所有任務
                future_to_ligand = {
                    executor.submit(self.dock_single_ligand, ligand_file): ligand_file 
                    for ligand_file in ligand_files
                }
                
                # 收集結果
                completed = 0
                for future in as_completed(future_to_ligand):
                    ligand_file = future_to_ligand[future]
                    
                    try:
                        result = future.result()
                        all_results.append(result)
                        completed += 1
                        
                        # 每 100 個完成時報告進度
                        if completed % 100 == 0:
                            logger.info(f"📈 進度: {completed}/{len(ligand_files)} "
                                      f"({completed/len(ligand_files)*100:.1f}%)")
                            
                    except Exception as e:
                        logger.error(f"❌ 處理 {ligand_file} 時發生錯誤: {e}")
        
        except KeyboardInterrupt:
            logger.warning("⚠️ 收到中斷信號，正在保存當前結果...")
        
        self.stats["end_time"] = time.time()
        self.stats["duration"] = self.stats["end_time"] - self.stats["start_time"]
        
        # 保存結果
        self.save_results(all_results)
        
        return len(all_results) > 0
    
    def save_results(self, results):
        """保存對接結果"""
        logger.info("💾 保存對接結果")
        
        # 保存詳細結果
        import pandas as pd
        df = pd.DataFrame(results)
        
        results_file = Path("docking_results/summaries/all_docking_results.csv")
        results_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(results_file, index=False)
        
        # 保存成功的結果（排序）
        successful_results = df[df["status"] == "success"].copy()
        if not successful_results.empty:
            successful_results = successful_results.sort_values("binding_affinity")
            
            success_file = Path("docking_results/summaries/successful_results.csv")
            successful_results.to_csv(success_file, index=False)
            
            # 保存前 1000 名
            top_results = successful_results.head(1000)
            top_file = Path("docking_results/summaries/top_1000_results.csv")
            top_results.to_csv(top_file, index=False)
        
        # 保存統計信息
        final_stats = {
            **self.stats,
            "success_rate": self.stats["successful_dockings"] / self.stats["total_ligands"] if self.stats["total_ligands"] > 0 else 0,
            "dockings_per_second": self.stats["successful_dockings"] / self.stats["duration"] if self.stats["duration"] > 0 else 0
        }
        
        stats_file = Path("docking_results/summaries/docking_statistics.json")
        with open(stats_file, 'w') as f:
            json.dump(final_stats, f, indent=2)
        
        # 打印統計信息
        logger.info("📊 對接統計:")
        logger.info(f"   總配體數: {final_stats['total_ligands']}")
        logger.info(f"   成功對接: {final_stats['successful_dockings']}")
        logger.info(f"   失敗對接: {final_stats['failed_dockings']}")
        logger.info(f"   成功率: {final_stats['success_rate']:.1%}")
        logger.info(f"   總耗時: {final_stats['duration']/3600:.1f} 小時")
        logger.info(f"   速度: {final_stats['dockings_per_second']:.1f} 對接/秒")
        
        logger.info(f"📁 結果文件:")
        logger.info(f"   📊 {results_file}")
        if not successful_results.empty:
            logger.info(f"   ✅ {success_file}")
            logger.info(f"   🏆 {top_file}")
        logger.info(f"   📈 {stats_file}")
    
    def resume_docking(self):
        """恢復中斷的對接"""
        logger.info("🔄 恢復對接任務")
        
        # 檢查已完成的對接
        completed_files = set()
        if self.output_dir.exists():
            for log_file in self.output_dir.glob("*_log.txt"):
                ligand_name = log_file.stem.replace("_log", "")
                completed_files.add(f"{ligand_name}.pdbqt")
        
        logger.info(f"📋 已完成 {len(completed_files)} 個對接")
        
        # 找到剩餘的配體
        all_ligands = self.find_ligand_files()
        remaining_ligands = [lig for lig in all_ligands 
                           if lig.name not in completed_files]
        
        if not remaining_ligands:
            logger.info("✅ 所有對接已完成")
            return True
        
        logger.info(f"🔄 還需對接 {len(remaining_ligands)} 個配體")
        
        # 繼續對接剩餘的配體
        return self.run_parallel_docking()

if __name__ == "__main__":
    Docker = ParallelDocker()
    
    # 檢查是否為恢復模式
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--resume":
        success = docker.resume_docking()
    else:
        success = docker.run_parallel_docking()
    
    if success:
        print("✅ 對接完成")
    else:
        print("❌ 對接失敗")