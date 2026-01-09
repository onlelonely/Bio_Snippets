#!/bin/bash
# ---------------------------------------------
# Title: Development Environment Setup for H2A-Pipe Project
# Description: From: Source/3. Efforts/Res_Ab_developability _ML/Archives/Development Environment Setup for H2A-Pipe Project.md (9 blocks)
# ---------------------------------------------

# --- Part 1 ---
# Save environment for reproducibility
conda env export > environment.yml

# Platform-independent export
conda env export --from-history > environment-minimal.yml

# --- Part 2 ---
# Setup ColabFold for local use
pip install "colabfold[alphafold]"

# Download model parameters
python -m colabfold.download

# Create prediction script
cat > scripts/run_alphafold.py << 'EOF'
import sys
from colabfold import batch

def predict_structure(sequence, output_dir):
    batch.run(
        queries=[("antibody", sequence)],
        output_dir=output_dir,
        num_models=5,
        use_templates=True,
        amber_relax=True
    )

if __name__ == "__main__":
    predict_structure(sys.argv[1], sys.argv[2])
EOF

# --- Part 3 ---
# Create local BLAST databases
mkdir -p data/blast_db

# Download and format antibody databases
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz
for file in nr.*.tar.gz; do tar -xzf $file; done

# Create custom database from project sequences
makeblastdb \
  -in data/processed/all_antibodies.fasta \
  -dbtype prot \
  -out data/blast_db/h2a_antibodies \
  -title "H2A-Pipe Antibody Database"

# --- Part 4 ---
# Check CUDA availability
nvidia-smi

# Configure PyTorch for GPU
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
python -c "import torch; print(f'CUDA device: {torch.cuda.get_device_name(0)}')"

# Set memory growth for TensorFlow (if used)
export TF_FORCE_GPU_ALLOW_GROWTH=true

# Monitor GPU usage
gpustat -i 1  # Install with: pip install gpustat

# --- Part 5 ---
# Initialize git with proper ignores
git init
cat > .gitignore << 'EOF'
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
env/
venv/
ENV/

# Jupyter
.ipynb_checkpoints
*.ipynb_checkpoints

# Data (use DVC instead)
data/raw/*
data/processed/*
!data/raw/.gitkeep
!data/processed/.gitkeep

# Models (use DVC)
models/*.pkl
models/*.pt
models/*.h5

# Logs
logs/
*.log

# Cache
.cache/
cache/

# IDE
.vscode/
.idea/
*.swp
*.swo

# Environment
.env
.env.local

# OS
.DS_Store
Thumbs.db
EOF

# Configure DVC for data versioning
dvc init
dvc add data/raw/gdpa1.csv
git add data/raw/gdpa1.csv.dvc .gitignore
git commit -m "Initialize DVC with GDPa1 dataset"

# --- Part 6 ---
pip install pre-commit
pre-commit install
pre-commit run --all-files  # Test on all files

# --- Part 7 ---
# Install gcloud CLI
curl https://sdk.cloud.google.com | bash

# Initialize and authenticate
gcloud init
gcloud auth login

# Create project
gcloud projects create h2a-pipe-dev --name="H2A-Pipe Development"
gcloud config set project h2a-pipe-dev

# Enable necessary APIs
gcloud services enable compute.googleapis.com
gcloud services enable storage.googleapis.com
gcloud services enable aiplatform.googleapis.com

# Create VM with GPU
gcloud compute instances create h2a-dev \
  --machine-type=n1-standard-8 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=pytorch-latest-gpu \
  --image-project=deeplearning-platform-release \
  --boot-disk-size=100GB \
  --zone=us-central1-a

# --- Part 8 ---
# .env file (never commit!)
ANTHROPIC_API_KEY=sk-ant-...
OPENAI_API_KEY=sk-...
HUGGINGFACE_TOKEN=hf_...
DATABASE_URL=postgresql://...
REDIS_URL=redis://...
AWS_ACCESS_KEY_ID=...
AWS_SECRET_ACCESS_KEY=...
GITHUB_TOKEN=ghp_...
NCBI_API_KEY=...

# --- Part 9 ---
# 1. Set up Claude Code and environment
pip install claude-code
claude-code init h2a-pipe --template bioinformatics

# 2. Download Ginkgo dataset
claude-code task "Download Ginkgo GDPa1 dataset from HuggingFace and create initial analysis notebook"

# 3. Generate baseline model
claude-code implement "Create XGBoost baseline model for antibody thermal stability prediction using GDPa1 data"

# 4. Start documentation
claude-code document --init "Create README with project overview, installation, and usage instructions"