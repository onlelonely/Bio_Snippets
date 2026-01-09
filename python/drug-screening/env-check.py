# ---------------------------------------------
# Source: Drug Screening Pipeline (CPV)
# ---------------------------------------------

"""
環境檢查腳本 - 驗證所有必需軟體已正確安裝
"""

import subprocess
import sys
import importlib
from pathlib import Path

def check_command(cmd, version_flag="--version"):
    """檢查命令列工具是否可用"""
    try:
        result = subprocess.run([cmd, version_flag], 
                              capture_output=True, text=True, timeout=10)
        return True, result.stdout.strip()
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        return False, "Not found"

def check_python_package(package_name):
    """檢查 Python 套件是否已安裝"""
    try:
        importlib.import_module(package_name)
        return True, "Available"
    except ImportError:
        return False, "Not installed"

def main():
    print("🔍 CPV 虛擬篩選環境檢查")
    print("=" * 50)
    
    # 檢查命令列工具
    tools = [
        ("vina", "--help"),
        ("obabel", "--version"),
        ("gmx", "--version"),
        ("pymol", "-h")
    ]
    
    print("\n📦 命令列工具:")
    for tool, flag in tools:
        available, info = check_command(tool, flag)
        status = "✅" if available else "❌"
        print(f"{status} {tool:<15} {info[:50]}")
    
    # 檢查 Python 套件
    packages = [
        "pandas", "numpy", "matplotlib", "seaborn",
        "rdkit", "Bio", "MDAnalysis", "scipy"
    ]
    
    print("\n🐍 Python 套件:")
    for package in packages:
        available, info = check_python_package(package)
        status = "✅" if available else "❌"
        print(f"{status} {package:<15} {info}")
    
    print("\n" + "=" * 50)
    print("環境檢查完成！")

if __name__ == "__main__":
    main()