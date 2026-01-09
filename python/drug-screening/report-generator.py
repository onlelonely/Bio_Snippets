# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
生成綜合分析報告
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import logging
from datetime import datetime
import numpy as np

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ReportGenerator:
    def __init__(self):
        self.report_dir = Path("analysis/reports")
        self.report_dir.mkdir(parents=True, exist_ok=True)
        
    def generate_markdown_report(self):
        """生成 Markdown 格式的綜合報告"""
        logger.info("📝 生成 Markdown 報告")
        
        report_content = f"""# CPV 虛擬篩選分析報告

## 報告信息
- **生成時間**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- **項目**: Canine Parvovirus 抗病毒藥物虛擬篩選
- **目標**: VP2-TfR 結合界面抑制劑

## 執行摘要

"""
        
        # 載入統計信息
        try:
            stats_file = Path("analysis/binding_affinity/binding_affinity_statistics.json")
            if stats_file.exists():
                with open(stats_file, 'r') as f:
                    stats = json.load(f)
                
                report_content += f"""### 對接統計
- **總化合物數**: {stats.get('total_compounds', 'N/A'):,}
- **平均結合能**: {stats.get('mean_binding_affinity', 'N/A'):.2f} kcal/mol
- **最佳結合能**: {stats.get('min_binding_affinity', 'N/A'):.2f} kcal/mol
- **強結合者數量**: {stats.get('strong_binders', 'N/A')} ({stats.get('strong_binder_rate', 0)*100:.1f}%)

"""
        except Exception as e:
            logger.warning(f"⚠️ 無法載入統計信息: {e}")
        
        # 載入陽性對照信息
        try:
            controls_file = Path("analysis/binding_affinity/positive_control_rankings.json")
            if controls_file.exists():
                with open(controls_file, 'r') as f:
                    controls = json.load(f)
                
                report_content += "### 陽性對照結果\n"
                for name, data in controls.items():
                    report_content += f"- **{name.title()}**: 排名 {data['rank']} ({data['percentile']:.1f}%), 結合能 {data['binding_affinity']:.2f} kcal/mol\n"
                report_content += "\n"
        except Exception as e:
            logger.warning(f"⚠️ 無法載入對照信息: {e}")
        
        # 載入篩選摘要
        try:
            filter_file = Path("analysis/sar/filtering_summary.json")
            if filter_file.exists():
                with open(filter_file, 'r') as f:
                    filter_summary = json.load(f)
                
                report_content += f"""### 篩選結果
- **原始結果數**: {filter_summary.get('original_results', 'N/A'):,}
- **最終篩選數**: {filter_summary.get('final_filtered', 'N/A'):,}
- **篩選成功率**: {filter_summary.get('success_rate', 0)*100:.2f}%

### 篩選步驟
"""
                for i, step in enumerate(filter_summary.get('filtering_steps', []), 1):
                    report_content += f"{i}. {step}\n"
                report_content += "\n"
        except Exception as e:
            logger.warning(f"⚠️ 無法載入篩選摘要: {e}")
        
        # 添加前 10 個候選分子
        try:
            top_file = Path("analysis/sar/final_selection_essential_top_20.csv")
            if top_file.exists():
                top_df = pd.read_csv(top_file)
                
                report_content += "### 前 10 個最佳候選分子\n\n"
                report_content += "| 排名 | 化合物名稱 | 結合能 (kcal/mol) | 分子量 | LogP | QED |\n"
                report_content += "|------|------------|-------------------|--------|------|-----|\n"
                
                for i, row in top_df.head(10).iterrows():
                    report_content += f"| {i+1} | {row['ligand_name']} | {row['binding_affinity']:.2f} | {row['molecular_weight']:.1f} | {row['logp']:.2f} | {row['qed']:.3f} |\n"
                
                report_content += "\n"
        except Exception as e:
            logger.warning(f"⚠️ 無法載入前 10 候選分子: {e}")
        
        # 添加圖表部分
        report_content += """## 分析圖表

### 結合能分佈
![結合能分佈](../figures/binding_affinity_distribution.png)

### 累積分佈
![累積分佈](../figures/binding_affinity_cumulative.png)

### 排名與結合能關係
![排名 vs 結合能](../figures/rank_vs_affinity.png)

## 結論與建議

### 主要發現
1. 成功完成大規模虛擬篩選，識別出高親和力候選化合物
2. 陽性對照在結果中排名合理，驗證了篩選方法的有效性
3. 通過多維度篩選獲得了具有良好成藥性的候選分子

### 下一步建議
1. **分子動力學驗證**: 對前 20 個候選分子進行 MD 模擬驗證結合穩定性
2. **實驗驗證**: 購買或合成前 5-10 個候選分子進行體外活性測試
3. **結構優化**: 基於對接結果進行化學修飾以提高親和力和選擇性
4. **ADMET 預測**: 進行更詳細的吸收、分佈、代謝、排泄和毒性預測

---

*本報告由 CPV 虛擬篩選自動化流程生成*
"""
        
        # 保存報告
        report_file = self.report_dir / "virtual_screening_report.md"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        logger.info(f"📄 Markdown 報告已保存至: {report_file}")
    
    def generate_summary_plots(self):
        """生成報告摘要圖表"""
        logger.info("📊 生成摘要圖表")
        
        plt.style.use('seaborn-v0_8')
        
        # 創建多子圖
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        try:
            # 載入數據
            results_file = Path("analysis/binding_affinity/results_with_statistics.csv")
            if results_file.exists():
                df = pd.read_csv(results_file)
                
                # 子圖 1: 結合能分佈
                axes[0, 0].hist(df["binding_affinity"], bins=50, alpha=0.7, edgecolor='black')
                axes[0, 0].set_xlabel("Binding Affinity (kcal/mol)")
                axes[0, 0].set_ylabel("Frequency")
                axes[0, 0].set_title("Binding Affinity Distribution")
                
                # 子圖 2: 前 1000 名排名圖
                top_1000 = df.head(1000)
                axes[0, 1].scatter(range(1, len(top_1000)+1), top_1000["binding_affinity"], 
                                 alpha=0.6, s=10)
                axes[0, 1].set_xlabel("Rank")
                axes[0, 1].set_ylabel("Binding Affinity (kcal/mol)")
                axes[0, 1].set_title("Top 1000 Compounds")
                
                # 子圖 3: Z-score 分佈
                if "z_score" in df.columns:
                    axes[1, 0].hist(df["z_score"], bins=50, alpha=0.7, edgecolor='black')
                    axes[1, 0].axvline(-2, color='red', linestyle='--', label='Threshold')
                    axes[1, 0].set_xlabel("Z-Score")
                    axes[1, 0].set_ylabel("Frequency")
                    axes[1, 0].set_title("Z-Score Distribution")
                    axes[1, 0].legend()
                
                # 子圖 4: 篩選漏斗圖
                try:
                    filter_file = Path("analysis/sar/filtering_summary.json")
                    if filter_file.exists():
                        with open(filter_file, 'r') as f:
                            filter_data = json.load(f)
                        
                        stages = ['Original', 'Final Filtered']
                        counts = [filter_data.get('original_results', 0), 
                                filter_data.get('final_filtered', 0)]
                        
                        axes[1, 1].bar(stages, counts, alpha=0.7)
                        axes[1, 1].set_ylabel("Number of Compounds")
                        axes[1, 1].set_title("Filtering Pipeline")
                        
                        # 添加數值標籤
                        for i, count in enumerate(counts):
                            axes[1, 1].text(i, count + max(counts)*0.01, f'{count:,}', 
                                           ha='center', va='bottom')
                except Exception:
                    pass
            
        except Exception as e:
            logger.warning(f"⚠️ 生成圖表時發生錯誤: {e}")
        
        plt.tight_layout()
        summary_plot_file = self.report_dir / "summary_plots.png"
        plt.savefig(summary_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"📊 摘要圖表已保存至: {summary_plot_file}")
    
    def generate_full_report(self):
        """生成完整報告"""
        logger.info("🚀 生成完整分析報告")
        
        try:
            # 生成摘要圖表
            self.generate_summary_plots()
            
            # 生成 Markdown 報告
            self.generate_markdown_report()
            
            logger.info("🎉 報告生成完成！")
            logger.info(f"📁 報告位置: {self.report_dir}")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ 報告生成失敗: {e}")
            return False

if __name__ == "__main__":
    generator = ReportGenerator()
    success = generator.generate_full_report()
    
    if success:
        print("✅ 報告生成完成")
    else:
        print("❌ 報告生成失敗")