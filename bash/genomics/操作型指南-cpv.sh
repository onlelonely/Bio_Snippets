#!/bin/bash
# ---------------------------------------------
# Title: 操作型指南_CPV
# Description: From: Source/3. Efforts/Res_Drug_&_Transcriptome/待整合/操作型指南_CPV.md (40 blocks)
# ---------------------------------------------

# --- Part 1 ---
#!/bin/bash
# 創建標準化項目目錄結構

echo "🏗️ 創建 CPV 虛擬篩選項目目錄結構..."

# 主要目錄
mkdir -p {structures,ligands,docking_results,md_simulations,analysis,figures,documentation,scripts,logs}

# 子目錄
mkdir -p structures/{original,processed,receptors}
mkdir -p ligands/{raw,processed,controls}
mkdir -p docking_results/{individual,summaries,failed}
mkdir -p md_simulations/{inputs,trajectories,analysis}
mkdir -p analysis/{binding_affinity,sar,statistics,visualizations}
mkdir -p scripts/{preparation,docking,analysis,validation}
mkdir -p logs/{docking,md,errors}

# 創建 README 文件
cat > directory_structure.txt << 'EOF'
📁 CPV 虛擬篩選項目目錄結構
├── structures/           # 蛋白質結構文件
│   ├── original/        # 原始 PDB 文件
│   ├── processed/       # 預處理後的結構
│   └── receptors/       # 對接用受體文件
├── ligands/             # 配體文件
│   ├── raw/            # 原始化合物庫
│   ├── processed/      # 轉換後的 PDBQT 文件
│   └── controls/       # 陽性和陰性對照
├── docking_results/     # 對接結果
│   ├── individual/     # 個別對接結果
│   ├── summaries/      # 結果匯總
│   └── failed/         # 失敗記錄
├── md_simulations/      # 分子動力學模擬
│   ├── inputs/         # 輸入文件
│   ├── trajectories/   # 軌跡文件
│   └── analysis/       # 分析結果
├── analysis/           # 數據分析
│   ├── binding_affinity/  # 結合親和力分析
│   ├── sar/              # 結構活性關係
│   ├── statistics/       # 統計分析
│   └── visualizations/   # 視覺化圖表
├── figures/            # 發表用圖表
├── documentation/      # 文檔記錄
├── scripts/           # 自動化腳本
│   ├── preparation/   # 準備階段腳本
│   ├── docking/      # 對接腳本
│   ├── analysis/     # 分析腳本
│   └── validation/   # 驗證腳本
└── logs/             # 日誌文件
    ├── docking/      # 對接日誌
    ├── md/          # MD 模擬日誌
    └── errors/      # 錯誤記錄
EOF

echo "✅ 項目目錄結構創建完成！"
echo "📖 詳細說明請查看 directory_structure.txt"

# 設置權限
chmod +x scripts/*.sh 2>/dev/null || true
chmod +x scripts/*.py 2>/dev/null || true

echo "🎯 請執行 'python scripts/check_environment.py' 檢查環境"

# --- Part 2 ---
python scripts/check_environment.py

# --- Part 3 ---
bash scripts/setup_project.sh

# --- Part 4 ---
ls -la  # 確認目錄結構
   cat directory_structure.txt  # 查看說明

# --- Part 5 ---
python scripts/analysis/analyze_binding_site.py

# --- Part 6 ---
python scripts/analysis/compile_inhibitors.py

# --- Part 7 ---
cat analysis/binding_site_analysis.json
   cat ligands/controls/known_inhibitors.csv

# --- Part 8 ---
# 啟動 PyMOL 並執行上述命令
   pymol structures/original/6OAS.pdb

# --- Part 9 ---
python scripts/preparation/prepare_receptor.py

# --- Part 10 ---
python scripts/preparation/identify_binding_site.py

# --- Part 11 ---
ls structures/receptors/
   cat docking_results/vina_config.txt
   pymol analysis/visualizations/binding_site.pml

# --- Part 12 ---
# 創建原始文件目錄
   mkdir -p ligands/raw
   
   # 手動下載 ZINC20 數據至 ligands/raw/
   # 下載對照化合物至 ligands/controls/

# --- Part 13 ---
python scripts/preparation/convert_ligands.py

# --- Part 14 ---
python scripts/preparation/validate_ligands.py

# --- Part 15 ---
ls ligands/processed/
   cat ligands/processed/conversion_stats.json
   cat ligands/processed/validation_statistics.json

# --- Part 16 ---
python scripts/docking/validate_docking.py

# --- Part 17 ---
# 正常執行
   python scripts/docking/run_docking.py
   
   # 恢復中斷的對接
   python scripts/docking/run_docking.py --resume

# --- Part 18 ---
# 檢查成功率
   tail -f logs/docking.log
   
   # 查看統計信息
   cat docking_results/summaries/docking_statistics.json

# --- Part 19 ---
ls docking_results/summaries/
   head -20 docking_results/summaries/top_1000_results.csv

# --- Part 20 ---
# 搜索 Nitazoxanide 在結果中的排名
   grep -n "Nitazoxanide\|nitazoxanide" docking_results/summaries/successful_results.csv

# --- Part 21 ---
# 查看前 20 個最佳結果
   head -20 docking_results/summaries/top_1000_results.csv

# --- Part 22 ---
python scripts/analysis/analyze_results.py

# --- Part 23 ---
python scripts/analysis/filter_hits.py

# --- Part 24 ---
python scripts/analysis/generate_reports.py

# --- Part 25 ---
ls analysis/
   cat analysis/sar/final_selection_essential_top_20.csv
   cat analysis/reports/virtual_screening_report.md

# --- Part 26 ---
# 使用 PyMOL 手動合併最佳對接姿態與受體
   pymol structures/receptors/receptor.pdbqt docking_results/individual/top_compound_result.pdbqt
   # 另存為 complex.pdb

# --- Part 27 ---
# 蛋白質拓撲
   gmx pdb2gmx -f complex.pdb -o processed.gro -p topol.top
   
   # 配體拓撲 (需要手動創建或使用 CGenFF)

# --- Part 28 ---
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
   gmx solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top
   gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
   gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral

# --- Part 29 ---
# 對前10個候選分子進行MD模擬
   python scripts/validation/run_md_simulations.py

# --- Part 30 ---
python scripts/validation/cross_validate.py

# --- Part 31 ---
ls md_simulations/
   cat md_simulations/simulation_summary.json

# --- Part 32 ---
cat analysis/validation/cross_validation_results.json

# --- Part 33 ---
df -h
   du -sh md_simulations/ docking_results/ ligands/

# --- Part 34 ---
ps aux | grep -E "(vina|gmx|python)"
   htop  # 查看 CPU 和記憶體使用

# --- Part 35 ---
find . -name "*.log" -mtime -1 -exec grep -l "ERROR\|FAILED" {} \;
   tail -50 logs/errors/*

# --- Part 36 ---
find docking_results/individual -name "*_result.pdbqt" | wc -l
   find md_simulations -name "md.xtc" | wc -l

# --- Part 37 ---
# 運行一次監控
   python scripts/monitoring/monitor_progress.py
   
   # 設置定期監控 (每小時)
   crontab -e
   # 添加: 0 * * * * cd /path/to/cpv && python scripts/monitoring/monitor_progress.py

# --- Part 38 ---
python scripts/monitoring/error_handler.py

# --- Part 39 ---
ls analysis/monitoring/
   cat analysis/monitoring/latest_progress_report.json

# --- Part 40 ---
# 在支援GUI的環境中檢視
   xdg-open analysis/monitoring/monitoring_dashboard.png