# ---------------------------------------------
# Title: AlphaFold
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/AI & Structural Biology/AlphaFold.md
# ---------------------------------------------

import requests

# 從 AlphaFold Database 下載結構
def download_alphafold_structure(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    with open(f"{uniprot_id}_alphafold.pdb", "wb") as f:
        f.write(response.content)
    return f"{uniprot_id}_alphafold.pdb"

# 使用範例
structure_file = download_alphafold_structure("P53_HUMAN")