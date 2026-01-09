# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
虛擬篩選進度監控儀表板
"""

import os
import time
import json
import psutil
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime, timedelta
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ProgressMonitor:
    def __init__(self):
        self.monitoring_dir = Path("analysis/monitoring")
        self.monitoring_dir.mkdir(parents=True, exist_ok=True)
        
        self.start_time = time.time()
        
    def check_system_resources(self):
        """檢查系統資源使用情況"""
        logger.info("🖥️ 檢查系統資源")
        
        # CPU 使用率
        cpu_percent = psutil.cpu_percent(interval=1)
        
        # 記憶體使用
        memory = psutil.virtual_memory()
        
        # 磁碟使用
        disk = psutil.disk_usage('/')
        
        resources = {
            "timestamp": datetime.now().isoformat(),
            "cpu_percent": cpu_percent,
            "memory_percent": memory.percent,
            "memory_available_gb": memory.available / (1024**3),
            "disk_used_percent": (disk.used / disk.total) * 100,
            "disk_free_gb": disk.free / (1024**3)
        }
        
        logger.info(f"   CPU: {cpu_percent:.1f}%")
        logger.info(f"   記憶體: {memory.percent:.1f}% (可用: {resources['memory_available_gb']:.1f} GB)")
        logger.info(f"   磁碟: {resources['disk_used_percent']:.1f}% (可用: {resources['disk_free_gb']:.1f} GB)")
        
        return resources
    
    def check_docking_progress(self):
        """檢查對接進度"""
        logger.info("🔗 檢查對接進度")
        
        # 計算總配體數
        ligand_dir = Path("ligands/processed/high_quality")
        if ligand_dir.exists():
            total_ligands = len(list(ligand_dir.glob("*.pdbqt")))
        else:
            total_ligands = 0
        
        # 計算完成的對接數
        docking_dir = Path("docking_results/individual")
        if docking_dir.exists():
            completed_dockings = len(list(docking_dir.glob("*_result.pdbqt")))
            failed_dockings = len(list(Path("docking_results/failed").glob("*.pdbqt")))
        else:
            completed_dockings = 0
            failed_dockings = 0
        
        # 計算進度
        if total_ligands > 0:
            progress_percent = (completed_dockings / total_ligands) * 100
            success_rate = completed_dockings / (completed_dockings + failed_dockings) if (completed_dockings + failed_dockings) > 0 else 0
        else:
            progress_percent = 0
            success_rate = 0
        
        docking_status = {
            "total_ligands": total_ligands,
            "completed_dockings": completed_dockings,
            "failed_dockings": failed_dockings,
            "progress_percent": progress_percent,
            "success_rate": success_rate * 100
        }
        
        logger.info(f"   總配體數: {total_ligands}")
        logger.info(f"   完成對接: {completed_dockings} ({progress_percent:.1f}%)")
        logger.info(f"   失敗對接: {failed_dockings}")
        logger.info(f"   成功率: {success_rate*100:.1f}%")
        
        return docking_status
    
    def check_md_progress(self):
        """檢查 MD 模擬進度"""
        logger.info("🧬 檢查 MD 模擬進度")
        
        md_dir = Path("md_simulations")
        if not md_dir.exists():
            return {"md_simulations": 0, "completed_md": 0}
        
        # 計算 MD 目錄數量
        md_dirs = [d for d in md_dir.iterdir() if d.is_dir()]
        
        # 計算完成的 MD 模擬
        completed_md = 0
        for md_subdir in md_dirs:
            if (md_subdir / "md.xtc").exists():
                completed_md += 1
        
        md_status = {
            "md_simulations": len(md_dirs),
            "completed_md": completed_md,
            "md_progress_percent": (completed_md / len(md_dirs) * 100) if md_dirs else 0
        }
        
        logger.info(f"   MD 模擬數: {len(md_dirs)}")
        logger.info(f"   完成 MD: {completed_md}")
        logger.info(f"   MD 進度: {md_status['md_progress_percent']:.1f}%")
        
        return md_status
    
    def check_file_sizes(self):
        """檢查重要目錄的文件大小"""
        logger.info("📁 檢查文件大小")
        
        important_dirs = [
            "ligands",
            "docking_results", 
            "md_simulations",
            "analysis"
        ]
        
        dir_sizes = {}
        
        for dir_name in important_dirs:
            dir_path = Path(dir_name)
            if dir_path.exists():
                total_size = sum(f.stat().st_size for f in dir_path.rglob('*') if f.is_file())
                dir_sizes[dir_name] = total_size / (1024**3)  # GB
                logger.info(f"   {dir_name}: {dir_sizes[dir_name]:.2f} GB")
            else:
                dir_sizes[dir_name] = 0
        
        return dir_sizes
    
    def check_running_processes(self):
        """檢查正在運行的相關進程"""
        logger.info("⚙️ 檢查運行中的進程")
        
        target_processes = ["vina", "gmx", "python"]
        running_processes = {}
        
        for proc in psutil.process_iter(['pid', 'name', 'cpu_percent', 'memory_percent']):
            try:
                proc_name = proc.info['name'].lower()
                for target in target_processes:
                    if target in proc_name:
                        if target not in running_processes:
                            running_processes[target] = []
                        
                        running_processes[target].append({
                            "pid": proc.info['pid'],
                            "name": proc.info['name'],
                            "cpu_percent": proc.info['cpu_percent'],
                            "memory_percent": proc.info['memory_percent']
                        })
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        
        for process_type, processes in running_processes.items():
            logger.info(f"   {process_type}: {len(processes)} 個進程運行中")
        
        return running_processes
    
    def estimate_remaining_time(self, docking_status):
        """估算剩餘時間"""
        if docking_status["progress_percent"] <= 0:
            return {"estimated_hours": "未知", "estimated_completion": "未知"}
        
        elapsed_time = time.time() - self.start_time
        
        # 基於已完成的百分比估算總時間
        if docking_status["progress_percent"] > 0:
            estimated_total_time = elapsed_time / (docking_status["progress_percent"] / 100)
            remaining_time = estimated_total_time - elapsed_time
            
            estimated_hours = remaining_time / 3600
            estimated_completion = datetime.now() + timedelta(seconds=remaining_time)
            
            return {
                "estimated_hours": estimated_hours,
                "estimated_completion": estimated_completion.strftime("%Y-%m-%d %H:%M:%S")
            }
        
        return {"estimated_hours": "未知", "estimated_completion": "未知"}
    
    def generate_progress_report(self):
        """生成進度報告"""
        logger.info("📊 生成進度報告")
        
        # 收集所有狀態信息
        system_resources = self.check_system_resources()
        docking_status = self.check_docking_progress()
        md_status = self.check_md_progress()
        file_sizes = self.check_file_sizes()
        running_processes = self.check_running_processes()
        time_estimates = self.estimate_remaining_time(docking_status)
        
        # 組合完整報告
        full_report = {
            "timestamp": datetime.now().isoformat(),
            "system_resources": system_resources,
            "docking_status": docking_status,
            "md_status": md_status,
            "file_sizes": file_sizes,
            "running_processes": running_processes,
            "time_estimates": time_estimates,
            "overall_status": self._determine_overall_status(docking_status, md_status)
        }
        
        # 保存報告
        report_file = self.monitoring_dir / f"progress_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, 'w') as f:
            json.dump(full_report, f, indent=2, default=str)
        
        # 保存最新報告
        latest_report_file = self.monitoring_dir / "latest_progress_report.json"
        with open(latest_report_file, 'w') as f:
            json.dump(full_report, f, indent=2, default=str)
        
        logger.info(f"📄 進度報告已保存: {report_file}")
        
        return full_report
    
    def _determine_overall_status(self, docking_status, md_status):
        """判斷整體狀態"""
        if docking_status["progress_percent"] < 10:
            return "對接初始階段"
        elif docking_status["progress_percent"] < 50:
            return "對接進行中"
        elif docking_status["progress_percent"] < 95:
            return "對接即將完成"
        elif docking_status["progress_percent"] >= 95 and md_status["md_simulations"] == 0:
            return "對接完成，等待分析"
        elif md_status["md_simulations"] > 0 and md_status["completed_md"] < md_status["md_simulations"]:
            return "MD 模擬進行中"
        else:
            return "所有計算完成"
    
    def generate_dashboard_plot(self, report):
        """生成監控儀表板圖表"""
        logger.info("📈 生成監控儀表板")
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle(f"CPV 虛擬篩選進度監控 - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 
                    fontsize=16)
        
        # 1. 對接進度
        docking_data = report["docking_status"]
        axes[0, 0].pie([docking_data["completed_dockings"], 
                       docking_data["total_ligands"] - docking_data["completed_dockings"]], 
                      labels=["已完成", "未完成"], 
                      autopct='%1.1f%%',
                      colors=['green', 'lightgray'])
        axes[0, 0].set_title("對接進度")
        
        # 2. 系統資源
        resources = report["system_resources"]
        resource_names = ["CPU", "記憶體", "磁碟"]
        resource_values = [resources["cpu_percent"], 
                          resources["memory_percent"], 
                          resources["disk_used_percent"]]
        
        bars = axes[0, 1].bar(resource_names, resource_values)
        axes[0, 1].set_ylabel("使用率 (%)")
        axes[0, 1].set_title("系統資源使用")
        axes[0, 1].set_ylim(0, 100)
        
        # 為不同使用率著色
        for bar, value in zip(bars, resource_values):
            if value > 80:
                bar.set_color('red')
            elif value > 60:
                bar.set_color('orange')
            else:
                bar.set_color('green')
        
        # 3. 文件大小分佈
        file_sizes = report["file_sizes"]
        if any(file_sizes.values()):
            axes[0, 2].pie(file_sizes.values(), 
                          labels=file_sizes.keys(), 
                          autopct='%1.1f%%')
            axes[0, 2].set_title("存儲空間分佈")
        
        # 4. MD 模擬進度
        md_data = report["md_status"]
        if md_data["md_simulations"] > 0:
            axes[1, 0].pie([md_data["completed_md"], 
                           md_data["md_simulations"] - md_data["completed_md"]], 
                          labels=["已完成", "進行中"], 
                          autopct='%1.1f%%',
                          colors=['blue', 'lightblue'])
        else:
            axes[1, 0].text(0.5, 0.5, "尚未開始\nMD 模擬", 
                           ha='center', va='center', transform=axes[1, 0].transAxes)
        axes[1, 0].set_title("MD 模擬進度")
        
        # 5. 成功率統計
        success_rate = docking_data.get("success_rate", 0)
        axes[1, 1].bar(["對接成功率"], [success_rate])
        axes[1, 1].set_ylabel("成功率 (%)")
        axes[1, 1].set_title("對接成功率")
        axes[1, 1].set_ylim(0, 100)
        
        # 6. 狀態摘要
        axes[1, 2].text(0.1, 0.8, f"整體狀態: {report['overall_status']}", 
                       transform=axes[1, 2].transAxes, fontsize=12, weight='bold')
        
        if "estimated_hours" in report["time_estimates"] and isinstance(report["time_estimates"]["estimated_hours"], (int, float)):
            axes[1, 2].text(0.1, 0.6, f"預計剩餘: {report['time_estimates']['estimated_hours']:.1f} 小時", 
                           transform=axes[1, 2].transAxes, fontsize=10)
            axes[1, 2].text(0.1, 0.4, f"預計完成: {report['time_estimates']['estimated_completion']}", 
                           transform=axes[1, 2].transAxes, fontsize=10)
        
        axes[1, 2].text(0.1, 0.2, f"總配體: {docking_data['total_ligands']:,}", 
                       transform=axes[1, 2].transAxes, fontsize=10)
        axes[1, 2].set_title("執行摘要")
        axes[1, 2].axis('off')
        
        plt.tight_layout()
        
        # 保存圖表
        dashboard_file = self.monitoring_dir / "monitoring_dashboard.png"
        plt.savefig(dashboard_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"📊 監控儀表板已保存: {dashboard_file}")
    
    def run_monitoring_cycle(self):
        """運行一次完整的監控循環"""
        logger.info("🔄 開始監控循環")
        
        try:
            # 生成進度報告
            report = self.generate_progress_report()
            
            # 生成儀表板圖表
            self.generate_dashboard_plot(report)
            
            # 檢查警告條件
            self._check_warnings(report)
            
            logger.info("✅ 監控循環完成")
            return True
            
        except Exception as e:
            logger.error(f"❌ 監控循環失敗: {e}")
            return False
    
    def _check_warnings(self, report):
        """檢查警告條件"""
        warnings = []
        
        # 檢查磁碟空間
        if report["system_resources"]["disk_used_percent"] > 85:
            warnings.append("磁碟空間不足 (>85%)")
        
        # 檢查記憶體使用
        if report["system_resources"]["memory_percent"] > 90:
            warnings.append("記憶體使用過高 (>90%)")
        
        # 檢查對接成功率
        if report["docking_status"]["success_rate"] < 80 and report["docking_status"]["completed_dockings"] > 100:
            warnings.append("對接成功率過低 (<80%)")
        
        if warnings:
            logger.warning("⚠️ 檢測到警告:")
            for warning in warnings:
                logger.warning(f"   - {warning}")
            
            # 保存警告記錄
            warning_file = self.monitoring_dir / "warnings.log"
            with open(warning_file, 'a') as f:
                f.write(f"{datetime.now().isoformat()}: {', '.join(warnings)}\n")

if __name__ == "__main__":
    monitor = ProgressMonitor()
    monitor.run_monitoring_cycle()