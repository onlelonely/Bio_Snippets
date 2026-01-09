# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
對接結果分析和排序
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import logging
from scipy import stats

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DockingAnalyzer:
    def __init__(self, results_dir="docking_results/summaries"):
        self.results_dir = Path(results_dir)
        self.analysis_dir = Path("analysis/binding_affinity")
        self.analysis_dir.mkdir(parents=True, exist_ok=True)
        
        self.results_df = None
        self.statistics = {}
        
    def load_results(self):
        """載入對接結果"""
        results_file = self.results_dir / "successful_results.csv"
        
        if not results_file.exists():
            logger.error(f"❌ 結果文件不存在: {results_file}")
            return False
        
        try:
            self.results_df = pd.read_csv(results_file)
            logger.info(f"📊 載入 {len(self.results_df)} 個成功對接結果")
            
            # 基本數據檢查
            if "binding_affinity" not in self.results_df.columns:
                logger.error("❌ 結果文件缺少 binding_affinity 欄位")
                return False
            
            # 移除無效數據
            original_count = len(self.results_df)
            self.results_df = self.results_df.dropna(subset=["binding_affinity"])
            self.results_df = self.results_df[self.results_df["binding_affinity"] != 0]
            
            if len(self.results_df) < original_count:
                logger.warning(f"⚠️ 移除了 {original_count - len(self.results_df)} 個無效結果")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ 載入結果文件失敗: {e}")
            return False
    
    def calculate_statistics(self):
        """計算統計信息"""
        logger.info("📈 計算統計信息")
        
        if self.results_df is None or self.results_df.empty:
            logger.error("❌ 沒有可用的結果數據")
            return
        
        binding_affinities = self.results_df["binding_affinity"]
        
        self.statistics = {
            "total_compounds": len(self.results_df),
            "mean_binding_affinity": float(binding_affinities.mean()),
            "std_binding_affinity": float(binding_affinities.std()),
            "min_binding_affinity": float(binding_affinities.min()),
            "max_binding_affinity": float(binding_affinities.max()),
            "median_binding_affinity": float(binding_affinities.median()),
            "q25_binding_affinity": float(binding_affinities.quantile(0.25)),
            "q75_binding_affinity": float(binding_affinities.quantile(0.75)),
            "percentiles": {
                "p1": float(binding_affinities.quantile(0.01)),
                "p5": float(binding_affinities.quantile(0.05)),
                "p10": float(binding_affinities.quantile(0.10)),
                "p90": float(binding_affinities.quantile(0.90)),
                "p95": float(binding_affinities.quantile(0.95)),
                "p99": float(binding_affinities.quantile(0.99))
            }
        }
        
        # 計算 Z-scores
        self.results_df["z_score"] = stats.zscore(binding_affinities)
        
        # 定義顯著結合者 (< -7 kcal/mol 或 Z-score < -2)
        strong_binders = self.results_df[
            (self.results_df["binding_affinity"] < -7.0) | 
            (self.results_df["z_score"] < -2.0)
        ]
        
        self.statistics["strong_binders"] = len(strong_binders)
        self.statistics["strong_binder_rate"] = len(strong_binders) / len(self.results_df)
        
        logger.info(f"   總化合物數: {self.statistics['total_compounds']}")
        logger.info(f"   平均結合能: {self.statistics['mean_binding_affinity']:.2f} kcal/mol")
        logger.info(f"   標準差: {self.statistics['std_binding_affinity']:.2f}")
        logger.info(f"   最佳結合能: {self.statistics['min_binding_affinity']:.2f} kcal/mol")
        logger.info(f"   強結合者: {self.statistics['strong_binders']} ({self.statistics['strong_binder_rate']:.1%})")
    
    def check_positive_controls(self):
        """檢查陽性對照排名"""
        logger.info("🎯 檢查陽性對照排名")
        
        # 載入已知抑制劑信息
        controls_file = Path("ligands/controls/known_inhibitors.json")
        if not controls_file.exists():
            logger.warning("⚠️ 未找到陽性對照信息")
            return {}
        
        with open(controls_file, 'r') as f:
            inhibitor_data = json.load(f)
        
        positive_controls = [comp["name"].lower() for comp in 
                           inhibitor_data.get("positive_controls", [])]
        
        control_rankings = {}
        
        for i, row in self.results_df.iterrows():
            ligand_name = row["ligand_name"].lower()
            
            for control_name in positive_controls:
                if control_name in ligand_name:
                    rank = i + 1  # 1-based ranking
                    percentile = (rank / len(self.results_df)) * 100
                    
                    control_rankings[control_name] = {
                        "rank": rank,
                        "percentile": percentile,
                        "binding_affinity": row["binding_affinity"],
                        "z_score": row["z_score"]
                    }
                    
                    logger.info(f"   {control_name}: 排名 {rank} ({percentile:.1f}%), "
                              f"結合能 {row['binding_affinity']:.2f} kcal/mol")
        
        if not control_rankings:
            logger.warning("⚠️ 未在結果中找到陽性對照")
        
        return control_rankings
    
    def generate_visualizations(self):
        """生成視覺化圖表"""
        logger.info("📊 生成視覺化圖表")
        
        if self.results_df is None or self.results_df.empty:
            logger.error("❌ 沒有數據可視覺化")
            return
        
        # 設置繪圖風格
        plt.style.use('seaborn-v0_8')
        fig_dir = Path("figures")
        fig_dir.mkdir(exist_ok=True)
        
        # 1. 結合能分佈直方圖
        plt.figure(figsize=(10, 6))
        plt.hist(self.results_df["binding_affinity"], bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel("Binding Affinity (kcal/mol)")
        plt.ylabel("Frequency")
        plt.title("Distribution of Binding Affinities")
        plt.axvline(self.statistics["mean_binding_affinity"], color='red', 
                   linestyle='--', label=f'Mean: {self.statistics["mean_binding_affinity"]:.2f}')
        plt.axvline(-7.0, color='green', linestyle='--', label='Strong Binder Threshold')
        plt.legend()
        plt.tight_layout()
        plt.savefig(fig_dir / "binding_affinity_distribution.png", dpi=300)
        plt.close()
        
        # 2. 累積分佈圖
        plt.figure(figsize=(10, 6))
        sorted_affinities = np.sort(self.results_df["binding_affinity"])
        y_vals = np.arange(1, len(sorted_affinities) + 1) / len(sorted_affinities)
        plt.plot(sorted_affinities, y_vals)
        plt.xlabel("Binding Affinity (kcal/mol)")
        plt.ylabel("Cumulative Probability")
        plt.title("Cumulative Distribution of Binding Affinities")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fig_dir / "binding_affinity_cumulative.png", dpi=300)
        plt.close()
        
        # 3. Z-score 分佈
        plt.figure(figsize=(10, 6))
        plt.hist(self.results_df["z_score"], bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel("Z-Score")
        plt.ylabel("Frequency")
        plt.title("Z-Score Distribution")
        plt.axvline(-2, color='red', linestyle='--', label='Significant Threshold (Z < -2)')
        plt.axvline(0, color='black', linestyle='-', alpha=0.5, label='Mean')
        plt.legend()
        plt.tight_layout()
        plt.savefig(fig_dir / "z_score_distribution.png", dpi=300)
        plt.close()
        
        # 4. 排名 vs 結合能散點圖
        plt.figure(figsize=(12, 6))
        ranks = range(1, min(1001, len(self.results_df) + 1))  # 前1000名
        top_affinities = self.results_df["binding_affinity"].head(1000)
        
        plt.scatter(ranks, top_affinities, alpha=0.6, s=20)
        plt.xlabel("Rank")
        plt.ylabel("Binding Affinity (kcal/mol)")
        plt.title("Binding Affinity vs Rank (Top 1000)")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fig_dir / "rank_vs_affinity.png", dpi=300)
        plt.close()
        
        logger.info(f"📁 圖表已保存至 {fig_dir}")
    
    def identify_outliers(self, method="iqr", threshold=1.5):
        """識別異常值"""
        logger.info(f"🔍 使用 {method} 方法識別異常值")
        
        binding_affinities = self.results_df["binding_affinity"]
        
        if method == "iqr":
            Q1 = binding_affinities.quantile(0.25)
            Q3 = binding_affinities.quantile(0.75)
            IQR = Q3 - Q1
            
            lower_bound = Q1 - threshold * IQR
            upper_bound = Q3 + threshold * IQR
            
            outliers = self.results_df[
                (self.results_df["binding_affinity"] < lower_bound) |
                (self.results_df["binding_affinity"] > upper_bound)
            ]
            
        elif method == "zscore":
            outliers = self.results_df[abs(self.results_df["z_score"]) > threshold]
        
        else:
            logger.error(f"❌ 不支持的異常值檢測方法: {method}")
            return pd.DataFrame()
        
        logger.info(f"   發現 {len(outliers)} 個異常值")
        
        if len(outliers) > 0:
            outliers_file = self.analysis_dir / f"outliers_{method}.csv"
            outliers.to_csv(outliers_file, index=False)
            logger.info(f"   異常值已保存至: {outliers_file}")
        
        return outliers
    
    def save_analysis_results(self):
        """保存分析結果"""
        logger.info("💾 保存分析結果")
        
        # 保存統計信息
        stats_file = self.analysis_dir / "binding_affinity_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(self.statistics, f, indent=2)
        
        # 保存陽性對照排名
        control_rankings = self.check_positive_controls()
        if control_rankings:
            controls_file = self.analysis_dir / "positive_control_rankings.json"
            with open(controls_file, 'w') as f:
                json.dump(control_rankings, f, indent=2)
        
        # 保存帶 Z-score 的完整結果
        if self.results_df is not None:
            enhanced_results_file = self.analysis_dir / "results_with_statistics.csv"
            self.results_df.to_csv(enhanced_results_file, index=False)
        
        logger.info(f"📁 分析結果已保存至 {self.analysis_dir}")
    
    def run_full_analysis(self):
        """執行完整的結果分析"""
        logger.info("🚀 開始結果分析")
        
        try:
            # 載入結果
            if not self.load_results():
                return False
            
            # 計算統計信息
            self.calculate_statistics()
            
            # 檢查陽性對照
            self.check_positive_controls()
            
            # 生成視覺化
            self.generate_visualizations()
            
            # 識別異常值
            self.identify_outliers(method="iqr")
            self.identify_outliers(method="zscore", threshold=2)
            
            # 保存結果
            self.save_analysis_results()
            
            logger.info("🎉 結果分析完成！")
            return True
            
        except Exception as e:
            logger.error(f"❌ 分析過程發生錯誤: {e}")
            return False

if __name__ == "__main__":
    analyzer = DockingAnalyzer()
    success = analyzer.run_full_analysis()
    
    if success:
        print("✅ 結果分析完成")
    else:
        print("❌ 結果分析失敗")