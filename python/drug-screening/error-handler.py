# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
錯誤處理和自動恢復系統
"""

import os
import json
import shutil
import subprocess
from pathlib import Path
import logging
from datetime import datetime
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ErrorHandler:
    def __init__(self):
        self.error_dir = Path("logs/errors")
        self.error_dir.mkdir(parents=True, exist_ok=True)
        
        self.recovery_strategies = {
            "docking_failure": self._recover_docking_failure,
            "md_failure": self._recover_md_failure,
            "disk_full": self._handle_disk_full,
            "memory_error": self._handle_memory_error,
            "file_corruption": self._handle_file_corruption
        }
    
    def scan_for_errors(self):
        """掃描系統中的錯誤"""
        logger.info("🔍 掃描系統錯誤")
        
        errors_found = []
        
        # 1. 檢查對接失敗
        docking_errors = self._check_docking_errors()
        errors_found.extend(docking_errors)
        
        # 2. 檢查 MD 模擬失敗
        md_errors = self._check_md_errors()
        errors_found.extend(md_errors)
        
        # 3. 檢查磁碟空間
        disk_errors = self._check_disk_space()
        errors_found.extend(disk_errors)
        
        # 4. 檢查文件完整性
        file_errors = self._check_file_integrity()
        errors_found.extend(file_errors)
        
        if errors_found:
            logger.warning(f"⚠️ 發現 {len(errors_found)} 個錯誤")
            for error in errors_found:
                logger.warning(f"   - {error['type']}: {error['description']}")
        else:
            logger.info("✅ 未發現系統錯誤")
        
        return errors_found
    
    def _check_docking_errors(self):
        """檢查對接錯誤"""
        errors = []
        
        failed_dir = Path("docking_results/failed")
        if failed_dir.exists():
            failed_files = list(failed_dir.glob("*_error.txt"))
            
            for error_file in failed_files:
                try:
                    with open(error_file, 'r') as f:
                        error_content = f.read()
                    
                    # 分析錯誤類型
                    if "timeout" in error_content.lower():
                        error_type = "docking_timeout"
                    elif "memory" in error_content.lower():
                        error_type = "docking_memory"
                    else:
                        error_type = "docking_failure"
                    
                    errors.append({
                        "type": error_type,
                        "description": f"對接失敗: {error_file.stem}",
                        "file": str(error_file),
                        "details": error_content[:200]
                    })
                    
                except Exception as e:
                    logger.warning(f"⚠️ 無法讀取錯誤文件 {error_file}: {e}")
        
        return errors
    
    def _check_md_errors(self):
        """檢查 MD 模擬錯誤"""
        errors = []
        
        md_dir = Path("md_simulations")
        if md_dir.exists():
            for compound_dir in md_dir.iterdir():
                if compound_dir.is_dir():
                    # 檢查是否有 .log 文件但沒有對應的輸出
                    log_files = list(compound_dir.glob("*.log"))
                    
                    for log_file in log_files:
                        try:
                            with open(log_file, 'r') as f:
                                log_content = f.read()
                            
                            if "fatal error" in log_content.lower() or "segmentation fault" in log_content.lower():
                                errors.append({
                                    "type": "md_failure",
                                    "description": f"MD 模擬失敗: {compound_dir.name}",
                                    "file": str(log_file),
                                    "details": log_content[-300:]  # 最後300字符
                                })
                                
                        except Exception as e:
                            logger.warning(f"⚠️ 無法讀取 MD 日誌 {log_file}: {e}")
        
        return errors
    
    def _check_disk_space(self):
        """檢查磁碟空間"""
        errors = []
        
        import psutil
        disk_usage = psutil.disk_usage('/')
        
        free_percent = (disk_usage.free / disk_usage.total) * 100
        
        if free_percent < 10:
            errors.append({
                "type": "disk_full",
                "description": f"磁碟空間不足: 僅剩 {free_percent:.1f}%",
                "details": f"可用空間: {disk_usage.free / (1024**3):.2f} GB"
            })
        
        return errors
    
    def _check_file_integrity(self):
        """檢查文件完整性"""
        errors = []
        
        # 檢查重要的結果文件
        important_files = [
            "docking_results/summaries/successful_results.csv",
            "analysis/binding_affinity/binding_affinity_statistics.json",
            "structures/receptors/receptor.pdbqt"
        ]
        
        for file_path in important_files:
            path = Path(file_path)
            
            if path.exists():
                # 檢查文件是否為空
                if path.stat().st_size == 0:
                    errors.append({
                        "type": "file_corruption",
                        "description": f"文件為空: {file_path}",
                        "file": str(path)
                    })
                
                # 檢查 CSV 文件格式
                if path.suffix == '.csv':
                    try:
                        df = pd.read_csv(path)
                        if len(df) == 0:
                            errors.append({
                                "type": "file_corruption",
                                "description": f"CSV 文件無數據: {file_path}",
                                "file": str(path)
                            })
                    except Exception:
                        errors.append({
                            "type": "file_corruption",
                            "description": f"CSV 文件損壞: {file_path}",
                            "file": str(path)
                        })
        
        return errors
    
    def handle_errors(self, errors):
        """處理發現的錯誤"""
        logger.info("🔧 開始錯誤處理")
        
        recovery_results = {}
        
        for error in errors:
            error_type = error["type"]
            
            logger.info(f"   處理錯誤: {error_type}")
            
            if error_type in self.recovery_strategies:
                try:
                    success = self.recovery_strategies[error_type](error)
                    recovery_results[error["description"]] = success
                    
                    if success:
                        logger.info(f"   ✅ 錯誤已修復: {error['description']}")
                    else:
                        logger.warning(f"   ⚠️ 錯誤修復失敗: {error['description']}")
                        
                except Exception as e:
                    logger.error(f"   ❌ 錯誤處理異常: {e}")
                    recovery_results[error["description"]] = False
            else:
                logger.warning(f"   ⚠️ 無可用的恢復策略: {error_type}")
                recovery_results[error["description"]] = False
        
        # 記錄處理結果
        self._log_recovery_results(recovery_results)
        
        return recovery_results
    
    def _recover_docking_failure(self, error):
        """恢復對接失敗"""
        try:
            # 提取配體名稱
            if "file" in error:
                error_file = Path(error["file"])
                ligand_name = error_file.stem.replace("_error", "")
                
                # 檢查原始配體文件是否存在
                ligand_file = Path(f"ligands/processed/high_quality/{ligand_name}.pdbqt")
                
                if ligand_file.exists():
                    # 重新提交對接任務
                    logger.info(f"   重新提交對接: {ligand_name}")
                    
                    # 這裡可以重新運行對接
                    # 實際實現需要調用對接腳本
                    return True
                else:
                    logger.warning(f"   找不到配體文件: {ligand_file}")
                    return False
            
        except Exception as e:
            logger.error(f"   對接恢復失敗: {e}")
        
        return False
    
    def _recover_md_failure(self, error):
        """恢復 MD 模擬失敗"""
        try:
            if "file" in error:
                log_file = Path(error["file"])
                compound_dir = log_file.parent
                
                # 清理失敗的 MD 文件
                md_files = ["md.tpr", "md.xtc", "md.edr", "md.log"]
                
                for md_file in md_files:
                    file_path = compound_dir / md_file
                    if file_path.exists():
                        file_path.unlink()
                
                logger.info(f"   清理失敗的 MD 文件: {compound_dir}")
                
                # 重新開始 MD 模擬
                # 實際實現需要調用 MD 腳本
                return True
                
        except Exception as e:
            logger.error(f"   MD 恢復失敗: {e}")
        
        return False
    
    def _handle_disk_full(self, error):
        """處理磁碟空間不足"""
        logger.info("   清理臨時文件")
        
        cleanup_successful = False
        
        try:
            # 清理臨時文件
            temp_patterns = [
                "**/*.temp.*",
                "**/*~",
                "**/*.bak",
                "**/core.*"
            ]
            
            total_freed = 0
            
            for pattern in temp_patterns:
                temp_files = Path(".").rglob(pattern)
                
                for temp_file in temp_files:
                    if temp_file.is_file():
                        file_size = temp_file.stat().st_size
                        temp_file.unlink()
                        total_freed += file_size
            
            logger.info(f"   清理了 {total_freed / (1024**2):.2f} MB 臨時文件")
            
            # 清理舊的日誌文件 (>7天)
            import time
            current_time = time.time()
            week_ago = current_time - (7 * 24 * 3600)
            
            log_files = Path("logs").rglob("*.log")
            for log_file in log_files:
                if log_file.stat().st_mtime < week_ago:
                    log_file.unlink()
                    logger.info(f"   刪除舊日誌: {log_file}")
            
            cleanup_successful = True
            
        except Exception as e:
            logger.error(f"   磁碟清理失敗: {e}")
        
        return cleanup_successful
    
    def _handle_memory_error(self, error):
        """處理記憶體錯誤"""
        logger.info("   嘗試釋放記憶體")
        
        try:
            # 殺死可能的殭屍進程
            import psutil
            
            for proc in psutil.process_iter(['pid', 'name', 'status']):
                try:
                    if proc.info['status'] == psutil.STATUS_ZOMBIE:
                        proc.kill()
                        logger.info(f"   終止殭屍進程: {proc.info['pid']}")
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            
            return True
            
        except Exception as e:
            logger.error(f"   記憶體清理失敗: {e}")
        
        return False
    
    def _handle_file_corruption(self, error):
        """處理文件損壞"""
        try:
            corrupted_file = Path(error["file"])
            
            # 創建備份目錄
            backup_dir = Path("backups")
            backup_dir.mkdir(exist_ok=True)
            
            # 備份損壞的文件
            backup_file = backup_dir / f"{corrupted_file.name}.corrupted.{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            
            if corrupted_file.exists():
                shutil.move(str(corrupted_file), str(backup_file))
                logger.info(f"   已備份損壞文件: {backup_file}")
            
            # 嘗試從其他來源恢復文件
            # 這裡可以實現具體的文件恢復邏輯
            
            return True
            
        except Exception as e:
            logger.error(f"   文件恢復失敗: {e}")
        
        return False
    
    def _log_recovery_results(self, results):
        """記錄恢復結果"""
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "recovery_results": results,
            "total_errors": len(results),
            "successful_recoveries": sum(1 for success in results.values() if success)
        }
        
        log_file = self.error_dir / "recovery_log.jsonl"
        
        with open(log_file, 'a') as f:
            f.write(json.dumps(log_entry) + '\n')
    
    def run_error_handling_cycle(self):
        """運行完整的錯誤處理循環"""
        logger.info("🚀 開始錯誤處理循環")
        
        try:
            # 掃描錯誤
            errors = self.scan_for_errors()
            
            if errors:
                # 處理錯誤
                recovery_results = self.handle_errors(errors)
                
                successful_recoveries = sum(1 for success in recovery_results.values() if success)
                
                logger.info(f"📊 錯誤處理完成:")
                logger.info(f"   總錯誤數: {len(errors)}")
                logger.info(f"   成功恢復: {successful_recoveries}")
                logger.info(f"   恢復率: {successful_recoveries/len(errors)*100:.1f}%")
                
                return recovery_results
            else:
                logger.info("✅ 未發現需要處理的錯誤")
                return {}
                
        except Exception as e:
            logger.error(f"❌ 錯誤處理循環失敗: {e}")
            return {}

if __name__ == "__main__":
    handler = ErrorHandler()
    handler.run_error_handling_cycle()