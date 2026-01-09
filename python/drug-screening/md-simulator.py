# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
自動化分子動力學模擬管道
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

class MDSimulator:
    def __init__(self, 
                 top_compounds_file="analysis/sar/final_selection_essential_top_20.csv",
                 receptor_file="structures/receptors/receptor.pdbqt",
                 docking_dir="docking_results/individual"):
        
        self.top_compounds_file = Path(top_compounds_file)
        self.receptor_file = Path(receptor_file)
        self.docking_dir = Path(docking_dir)
        self.md_dir = Path("md_simulations")
        self.md_dir.mkdir(parents=True, exist_ok=True)
        
        # MD 參數
        self.simulation_time = 50  # ns
        self.temperature = 300     # K
        self.pressure = 1.0        # bar
        
    def load_top_compounds(self):
        """載入前N個化合物"""
        if not self.top_compounds_file.exists():
            logger.error(f"❌ 化合物文件不存在: {self.top_compounds_file}")
            return []
        
        import pandas as pd
        df = pd.read_csv(self.top_compounds_file)
        
        # 取前10個進行MD模擬
        top_compounds = df.head(10)["ligand_name"].tolist()
        
        logger.info(f"📋 載入 {len(top_compounds)} 個化合物進行 MD 模擬")
        return top_compounds
    
    def prepare_complex_structure(self, ligand_name):
        """準備蛋白質-配體複合體結構"""
        logger.info(f"🔧 準備複合體結構: {ligand_name}")
        
        compound_dir = self.md_dir / ligand_name
        compound_dir.mkdir(exist_ok=True)
        
        # 查找對接結果
        docked_ligand = self.docking_dir / f"{ligand_name}_result.pdbqt"
        
        if not docked_ligand.exists():
            logger.error(f"❌ 未找到對接結果: {docked_ligand}")
            return None
        
        # 轉換格式
        receptor_pdb = compound_dir / "receptor.pdb"
        ligand_pdb = compound_dir / "ligand.pdb"
        
        # 轉換受體
        cmd_receptor = [
            "obabel", str(self.receptor_file), 
            "-O", str(receptor_pdb)
        ]
        
        # 轉換配體 (取第一個構象)
        cmd_ligand = [
            "obabel", str(docked_ligand),
            "-O", str(ligand_pdb),
            "-m1"  # 只取第一個構象
        ]
        
        try:
            subprocess.run(cmd_receptor, check=True, capture_output=True)
            subprocess.run(cmd_ligand, check=True, capture_output=True)
            
            # 合併為複合體
            complex_pdb = compound_dir / "complex.pdb"
            
            with open(complex_pdb, 'w') as outfile:
                # 寫入受體
                with open(receptor_pdb, 'r') as f:
                    for line in f:
                        if line.startswith(('ATOM', 'HETATM')):
                            outfile.write(line)
                
                # 寫入配體
                with open(ligand_pdb, 'r') as f:
                    for line in f:
                        if line.startswith(('ATOM', 'HETATM')):
                            # 修改殘基名稱為 LIG
                            if len(line) >= 21:
                                line = line[:17] + "LIG" + line[20:]
                            outfile.write(line)
                
                outfile.write("END\n")
            
            logger.info(f"✅ 複合體結構已準備: {complex_pdb}")
            return complex_pdb
            
        except subprocess.CalledProcessError as e:
            logger.error(f"❌ 格式轉換失敗: {e}")
            return None
    
    def generate_topology(self, complex_pdb, compound_dir):
        """生成 GROMACS 拓撲"""
        logger.info("📐 生成拓撲文件")
        
        # 蛋白質拓撲
        try:
            cmd = [
                "gmx", "pdb2gmx",
                "-f", str(complex_pdb),
                "-o", str(compound_dir / "processed.gro"),
                "-p", str(compound_dir / "topol.top"),
                "-ff", "oplsaa",  # 力場選擇
                "-water", "spce",
                "-ignh"           # 忽略氫原子
            ]
            
            result = subprocess.run(cmd, input="1\n1\n", text=True, 
                                  capture_output=True, cwd=compound_dir)
            
            if result.returncode != 0:
                logger.error(f"❌ 拓撲生成失敗: {result.stderr}")
                return False
            
            logger.info("✅ 拓撲生成成功")
            return True
            
        except Exception as e:
            logger.error(f"❌ 拓撲生成異常: {e}")
            return False
    
    def create_mdp_files(self, compound_dir):
        """創建 MDP 參數文件"""
        
        # 能量最小化
        em_mdp = f"""
; Energy minimization
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
"""
        
        # NVT 平衡
        nvt_mdp = f"""
; NVT equilibration
define          = -DPOSRES
integrator      = md
nsteps          = 50000
dt              = 0.002
nstxout         = 500
nstvout         = 500
nstenergy       = 500
nstlog          = 500
continuation    = no
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4
ns_type         = grid
nstlist         = 5
cutoff-scheme   = Verlet
rcoulomb        = 1.0
rvdw            = 1.0
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1     0.1
ref_t           = {self.temperature}     {self.temperature}
pbc             = xyz
gen_vel         = yes
gen_temp        = {self.temperature}
gen_seed        = -1
"""
        
        # NPT 平衡
        npt_mdp = f"""
; NPT equilibration
define          = -DPOSRES
integrator      = md
nsteps          = 50000
dt              = 0.002
nstxout         = 500
nstvout         = 500
nstenergy       = 500
nstlog          = 500
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4
ns_type         = grid
nstlist         = 5
cutoff-scheme   = Verlet
rcoulomb        = 1.0
rvdw            = 1.0
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1     0.1
ref_t           = {self.temperature}     {self.temperature}
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = {self.pressure}
compressibility = 4.5e-5
pbc             = xyz
gen_vel         = no
"""
        
        # 生產運行
        md_mdp = f"""
; Production MD
integrator      = md
nsteps          = {int(self.simulation_time * 1000000 / 2)}  ; {self.simulation_time} ns
dt              = 0.002
nstxout         = 5000
nstvout         = 5000
nstxtcout       = 1000
nstenergy       = 1000
nstlog          = 1000
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4
ns_type         = grid
nstlist         = 5
cutoff-scheme   = Verlet
rcoulomb        = 1.0
rvdw            = 1.0
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1     0.1
ref_t           = {self.temperature}     {self.temperature}
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = {self.pressure}
compressibility = 4.5e-5
pbc             = xyz
gen_vel         = no
"""
        
        # 寫入文件
        mdp_files = {
            "em.mdp": em_mdp,
            "nvt.mdp": nvt_mdp,
            "npt.mdp": npt_mdp,
            "md.mdp": md_mdp
        }
        
        for filename, content in mdp_files.items():
            with open(compound_dir / filename, 'w') as f:
                f.write(content)
        
        logger.info("📝 MDP 文件已創建")
    
    def run_simulation_step(self, step_name, mdp_file, input_gro, output_prefix, compound_dir):
        """運行單個模擬步驟"""
        logger.info(f"⚡ 運行 {step_name}")
        
        tpr_file = compound_dir / f"{output_prefix}.tpr"
        
        # 創建 TPR 文件
        grompp_cmd = [
            "gmx", "grompp",
            "-f", str(mdp_file),
            "-c", str(input_gro),
            "-p", "topol.top",
            "-o", str(tpr_file),
            "-maxwarn", "2"
        ]
        
        try:
            result = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=compound_dir)
            
            if result.returncode != 0:
                logger.error(f"❌ {step_name} GROMPP 失敗: {result.stderr}")
                return False
            
            # 運行模擬
            if step_name == "Energy Minimization":
                mdrun_cmd = ["gmx", "mdrun", "-deffnm", output_prefix, "-v"]
            else:
                mdrun_cmd = ["gmx", "mdrun", "-deffnm", output_prefix, "-v", "-nb", "GPU"]
            
            result = subprocess.run(mdrun_cmd, capture_output=True, text=True, cwd=compound_dir)
            
            if result.returncode != 0:
                logger.error(f"❌ {step_name} MDRUN 失敗: {result.stderr}")
                return False
            
            logger.info(f"✅ {step_name} 完成")
            return True
            
        except Exception as e:
            logger.error(f"❌ {step_name} 異常: {e}")
            return False
    
    def run_full_simulation(self, ligand_name):
        """運行完整的 MD 模擬流程"""
        logger.info(f"🚀 開始 {ligand_name} 的 MD 模擬")
        
        compound_dir = self.md_dir / ligand_name
        
        try:
            # 1. 準備結構
            complex_pdb = self.prepare_complex_structure(ligand_name)
            if not complex_pdb:
                return False
            
            # 2. 生成拓撲
            if not self.generate_topology(complex_pdb, compound_dir):
                return False
            
            # 3. 添加溶劑盒子
            logger.info("💧 添加溶劑")
            
            # 定義盒子
            editconf_cmd = [
                "gmx", "editconf",
                "-f", "processed.gro",
                "-o", "newbox.gro",
                "-c", "-d", "1.0", "-bt", "cubic"
            ]
            
            result = subprocess.run(editconf_cmd, capture_output=True, cwd=compound_dir)
            if result.returncode != 0:
                logger.error("❌ 盒子定義失敗")
                return False
            
            # 添加水
            solvate_cmd = [
                "gmx", "solvate",
                "-cp", "newbox.gro",
                "-cs", "spc216.gro",
                "-o", "solvated.gro",
                "-p", "topol.top"
            ]
            
            result = subprocess.run(solvate_cmd, capture_output=True, cwd=compound_dir)
            if result.returncode != 0:
                logger.error("❌ 溶劑化失敗")
                return False
            
            # 4. 創建 MDP 文件
            self.create_mdp_files(compound_dir)
            
            # 5. 運行模擬步驟
            simulation_steps = [
                ("Energy Minimization", "em.mdp", "solvated.gro", "em"),
                ("NVT Equilibration", "nvt.mdp", "em.gro", "nvt"),
                ("NPT Equilibration", "npt.mdp", "nvt.gro", "npt"),
                ("Production MD", "md.mdp", "npt.gro", "md")
            ]
            
            for step_name, mdp_file, input_gro, output_prefix in simulation_steps:
                if not self.run_simulation_step(step_name, mdp_file, input_gro, output_prefix, compound_dir):
                    logger.error(f"❌ {ligand_name} 模擬在 {step_name} 階段失敗")
                    return False
            
            logger.info(f"🎉 {ligand_name} MD 模擬完成！")
            return True
            
        except Exception as e:
            logger.error(f"❌ {ligand_name} 模擬過程發生錯誤: {e}")
            return False
    
    def analyze_md_results(self, ligand_name):
        """分析 MD 模擬結果"""
        logger.info(f"📊 分析 {ligand_name} MD 結果")
        
        compound_dir = self.md_dir / ligand_name
        
        if not (compound_dir / "md.xtc").exists():
            logger.error(f"❌ 未找到軌跡文件: {compound_dir / 'md.xtc'}")
            return {}
        
        analysis_results = {}
        
        try:
            # RMSD 分析
            rmsd_cmd = [
                "gmx", "rms",
                "-s", "md.tpr",
                "-f", "md.xtc",
                "-o", "rmsd.xvg",
                "-tu", "ns"
            ]
            
            result = subprocess.run(rmsd_cmd, input="4\n4\n", text=True, 
                                  capture_output=True, cwd=compound_dir)
            
            if result.returncode == 0:
                # 解析 RMSD 數據
                rmsd_file = compound_dir / "rmsd.xvg"
                if rmsd_file.exists():
                    rmsd_values = []
                    with open(rmsd_file, 'r') as f:
                        for line in f:
                            if not line.startswith(('@', '#')):
                                parts = line.strip().split()
                                if len(parts) >= 2:
                                    rmsd_values.append(float(parts[1]))
                    
                    if rmsd_values:
                        analysis_results["rmsd_mean"] = sum(rmsd_values) / len(rmsd_values)
                        analysis_results["rmsd_max"] = max(rmsd_values)
                        analysis_results["rmsd_final"] = rmsd_values[-1]
            
            # 回轉半徑分析
            rg_cmd = [
                "gmx", "gyrate",
                "-s", "md.tpr",
                "-f", "md.xtc",
                "-o", "gyrate.xvg"
            ]
            
            result = subprocess.run(rg_cmd, input="1\n", text=True,
                                  capture_output=True, cwd=compound_dir)
            
            # 能量分析
            energy_cmd = [
                "gmx", "energy",
                "-s", "md.tpr",
                "-f", "md.edr",
                "-o", "energy.xvg"
            ]
            
            result = subprocess.run(energy_cmd, input="10\n0\n", text=True,
                                  capture_output=True, cwd=compound_dir)
            
            logger.info(f"✅ {ligand_name} 分析完成")
            
        except Exception as e:
            logger.warning(f"⚠️ {ligand_name} 分析部分失敗: {e}")
        
        return analysis_results
    
    def run_parallel_simulations(self, max_workers=2):
        """並行運行多個 MD 模擬"""
        top_compounds = self.load_top_compounds()
        
        if not top_compounds:
            logger.error("❌ 沒有化合物可供模擬")
            return False
        
        logger.info(f"🚀 開始並行 MD 模擬，使用 {max_workers} 個進程")
        
        results = {}
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交模擬任務
            future_to_compound = {
                executor.submit(self.run_full_simulation, compound): compound
                for compound in top_compounds
            }
            
            # 收集結果
            for future in as_completed(future_to_compound):
                compound = future_to_compound[future]
                
                try:
                    success = future.result()
                    results[compound] = success
                    
                    if success:
                        # 分析結果
                        analysis = self.analyze_md_results(compound)
                        results[f"{compound}_analysis"] = analysis
                        
                        logger.info(f"✅ {compound} 模擬和分析完成")
                    else:
                        logger.error(f"❌ {compound} 模擬失敗")
                        
                except Exception as e:
                    logger.error(f"❌ {compound} 處理異常: {e}")
                    results[compound] = False
        
        # 保存結果摘要
        summary_file = self.md_dir / "simulation_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        successful_simulations = sum(1 for k, v in results.items() 
                                   if not k.endswith('_analysis') and v)
        
        logger.info(f"📊 MD 模擬摘要:")
        logger.info(f"   總化合物數: {len(top_compounds)}")
        logger.info(f"   成功模擬: {successful_simulations}")
        logger.info(f"   成功率: {successful_simulations/len(top_compounds)*100:.1f}%")
        logger.info(f"💾 結果摘要: {summary_file}")
        
        return successful_simulations > 0

if __name__ == "__main__":
    simulator = MDSimulator()
    success = simulator.run_parallel_simulations(max_workers=2)
    
    if success:
        print("✅ MD 模擬完成")
    else:
        print("❌ MD 模擬失敗")