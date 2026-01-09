#!/bin/bash
# ---------------------------------------------
# Title: 設定指南
# Description: From: Source/3. Efforts/Res_Ab_developability _ML/Methods/設定指南.md (30 blocks)
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
env/
virtualenv/
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
# System requirements
- Python 3.11+ 
- Node.js 18+ (for MCP servers)
- Git
- 8GB+ RAM recommended
- CUDA-capable GPU (optional, for deep learning)

# --- Part 10 ---
# Install via pip (recommended)
pip install claude-code

# Or install from source
git clone https://github.com/anthropic/claude-code
cd claude-code
pip install -e .

# Verify installation
claude-code --version

# --- Part 11 ---
# Set up API key
export ANTHROPIC_API_KEY="your-api-key-here"

# Or create config file
cat > ~/.claude-code/config.yaml << EOF
api_key: "your-api-key-here"
model: "claude-opus-4-1-20250805"
max_tokens: 4096
temperature: 0.3  # Lower for more deterministic code generation
EOF

# --- Part 12 ---
# Create H2A-Pipe project structure
claude-code init h2a-pipe --template bioinformatics

# This creates:
h2a-pipe/
├── .claude-code/
│   ├── config.yaml         # Project-specific Claude settings
│   ├── prompts/            # Custom prompt templates
│   └── context/            # Project context files
├── src/
│   ├── antibody/           # Antibody analysis modules
│   ├── ml_models/          # Machine learning pipelines
│   └── utils/              # Utility functions
├── data/
│   ├── raw/                # Raw sequence data
│   ├── processed/          # Processed features
│   └── models/             # Trained models
├── notebooks/              # Jupyter notebooks
├── tests/                  # Unit tests
└── requirements.txt        # Dependencies

# --- Part 13 ---
# Generate antibody feature extractor
claude-code generate "Create a comprehensive antibody feature extractor class that includes:
- Amino acid composition
- CDR region identification and characterization  
- Hydrophobicity and charge calculations
- Predicted glycosylation sites
- Use Biopython and include type hints"

# Implement cross-species transfer learning
claude-code implement "Build a transfer learning pipeline that:
- Takes human antibody model (trained on GDPa1)
- Adapts to canine antibodies
- Uses domain adaptation techniques
- Includes confidence estimation
- Saves adapted models"

# Debug existing code
claude-code debug src/ml_models/transfer_learning.py \
  --error "KeyError in feature alignment" \
  --context "Processing canine antibody sequences"

# Optimize performance
claude-code optimize src/antibody/feature_extractor.py \
  --metric "execution_time" \
  --constraint "memory < 4GB" \
  --suggestion "Consider batch processing and caching"

# --- Part 14 ---
# Sequence analysis pipeline
claude-code create-pipeline "Antibody sequence analysis" \
  --steps "
    1. Parse FASTA files
    2. Quality control (length, valid amino acids)
    3. CDR identification using IMGT numbering
    4. Extract 350-dimensional feature vectors
    5. Save as HDF5 for efficient loading" \
  --input "data/raw/canine_antibodies.fasta" \
  --output "data/processed/canine_features.h5"

# ML model development
claude-code ml-model "XGBoost for antibody stability prediction" \
  --features "sequence, structure, physicochemical" \
  --target "thermal_stability" \
  --validation "5-fold cross-validation" \
  --metrics "R2, RMSE, Spearman correlation" \
  --save "models/stability_predictor_v1.pkl"

# Data visualization
claude-code visualize "Cross-species prediction results" \
  --plot-types "heatmap, scatter, violin" \
  --comparisons "human vs canine vs feline" \
  --save-format "publication-ready PDF"

# --- Part 15 ---
# Literature integration
claude-code research "Find and summarize papers on:
- Canine antibody therapeutics
- Cross-species antibody engineering
- Transfer learning in bioinformatics
Output as structured JSON with key findings"

# Hypothesis generation
claude-code hypothesize \
  --data "results/cross_species_analysis.csv" \
  --question "Why do canine antibodies show lower prediction accuracy?" \
  --suggest "Testable experiments using computational methods only"

# Paper writing assistance
claude-code write-methods \
  --pipeline "src/pipelines/h2a_main.py" \
  --style "Nature Methods" \
  --include "Equations, parameter justification, validation approach"

# --- Part 16 ---
# Install
npm install -g @anthropic/mcp-server-filesystem

# Configure for H2A-Pipe
cat > ~/.claude-code/mcp-config.json << EOF
{
  "mcpServers": {
    "filesystem": {
      "command": "mcp-server-filesystem",
      "args": ["--root", "/path/to/h2a-pipe"],
      "permissions": {
        "read": ["data/", "results/", "src/"],
        "write": ["results/", "logs/", "models/"]
      }
    }
  }
}
EOF

# --- Part 17 ---
# Install PostgreSQL MCP server
npm install -g @anthropic/mcp-server-postgres

# Configure
{
  "postgres": {
    "command": "mcp-server-postgres",
    "args": [
      "--connection-string", "postgresql://user:pass@localhost/antibody_db",
      "--schema", "public"
    ],
    "description": "Antibody sequence and annotation database"
  }
}

# --- Part 18 ---
# Custom MCP server for Python execution
npm install -g @anthropic/mcp-server-python

# Configure with bioinformatics environment
{
  "python-bio": {
    "command": "mcp-server-python",
    "args": [
      "--conda-env", "h2a-pipe",
      "--packages", "biopython,pandas,numpy,scikit-learn,torch"
    ]
  }
}

# --- Part 19 ---
# Download and process GDPa1 data
claude-code task "Download Ginkgo GDPa1 dataset from HuggingFace and:
1. Parse antibody sequences and biophysical properties
2. Create train/validation/test splits (70/15/15)
3. Generate summary statistics
4. Save processed data as Parquet files
Include error handling and progress bars"

# Expected output structure
data/processed/
├── gdpa1_train.parquet
├── gdpa1_val.parquet
├── gdpa1_test.parquet
└── gdpa1_stats.json

# --- Part 20 ---
# Create transfer learning pipeline
claude-code implement --interactive "Build complete transfer learning pipeline:

Phase 1: Train on human data (GDPa1)
- Use ensemble of RF, XGBoost, and NN
- Predict 10 biophysical properties
- Save feature importance

Phase 2: Adapt to canine
- Implement DANN (Domain Adversarial Neural Network)
- Use MMD for distribution matching
- Fine-tune on limited canine data (if available)

Phase 3: Evaluate
- Zero-shot performance
- Few-shot learning curves
- Confidence calibration

Use modular design with clear interfaces"

# --- Part 21 ---
# Create FastAPI server for predictions
claude-code api "Create FastAPI server for H2A-Pipe with:
- POST /predict endpoint accepting FASTA sequences
- Batch prediction support
- Species selection (human/canine/feline)
- Confidence intervals
- Result caching with Redis
- OpenAPI documentation
- Rate limiting and authentication"

# Deploy configuration
claude-code deploy --platform "google-cloud-run" \
  --dockerfile \
  --requirements \
  --github-actions

# --- Part 22 ---
# DVC (Data Version Control) setup
dvc init
dvc add data/raw/gdpa1.csv
dvc add models/baseline_rf.pkl
git add data/raw/gdpa1.csv.dvc models/baseline_rf.pkl.dvc
git commit -m "Add initial data and model versions"

# Track experiments
dvc exp run -n "canine_transfer_v1" \
  python src/train.py --species canine --epochs 100

# --- Part 23 ---
# Generate comprehensive tests
claude-code test "Create pytest tests for:
- Feature extraction functions
- Model prediction pipeline
- Data validation
- API endpoints
Include fixtures for sample sequences and expected outputs"

# Run tests with coverage
pytest tests/ --cov=src --cov-report=html

# Continuous testing during development
claude-code watch src/ --test --on-change

# --- Part 24 ---
claude-code troubleshoot "Memory error when processing 10,000 sequences" \
  --suggest "
  - Implement batch processing with configurable batch size
  - Use memory-mapped files for large arrays
  - Profile memory usage with memory_profiler
  - Consider Dask for distributed processing"

# --- Part 25 ---
claude-code optimize-performance \
  --profile src/features/extractor.py \
  --bottleneck "Sequence alignment taking 90% of time" \
  --suggest "
  - Parallelize with multiprocessing
  - Use Numba JIT compilation
  - Cache alignment results
  - Consider approximate methods"

# --- Part 26 ---
claude-code improve-model \
  --current-performance "R2=0.35 for canine" \
  --target "R2>0.50" \
  --strategies "
  - Increase training data augmentation
  - Implement species-specific encoders
  - Use evolutionary information
  - Ensemble multiple approaches"

# --- Part 27 ---
# Start simple, iterate complexity
claude-code iterate "Simple antibody classifier" \
  --versions "
  v1: Binary classification (functional/non-functional)
  v2: Add multi-class (stability levels)
  v3: Add regression (exact values)
  v4: Add uncertainty quantification"

# --- Part 28 ---
# Auto-generate documentation
claude-code document src/ \
  --format "sphinx" \
  --include "API reference, tutorials, theory" \
  --examples "From notebooks/*.ipynb"

# --- Part 29 ---
# Create shareable analysis
claude-code notebook "Cross-species analysis" \
  --include "
  - Motivation and background
  - Data loading and preprocessing
  - Model training with explanations
  - Results visualization
  - Interactive widgets for exploration" \
  --export "HTML for sharing with collaborators"

# --- Part 30 ---
# Project initialization
claude-code init [project-name] --template bioinformatics

# Code generation
claude-code generate [description]
claude-code implement [specification]
claude-code api [endpoints]

# Analysis and optimization  
claude-code analyze [file/directory]
claude-code optimize [target] --metric [metric]
claude-code profile [script]

# Testing and validation
claude-code test [module] --coverage
claude-code validate [predictions] --biological-constraints

# Documentation
claude-code document [source] --format [sphinx|markdown]
claude-code notebook [analysis] --interactive

# Deployment
claude-code deploy --platform [gcp|aws|azure]
claude-code containerize --dockerfile

# Research tasks
claude-code research [query] --papers --code-examples
claude-code hypothesize --data [file] --question [question]
claude-code visualize [results] --publication-ready