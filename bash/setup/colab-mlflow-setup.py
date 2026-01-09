# ---------------------------------------------
# Title: 設定指南
# Description: From: Source/3. Efforts/Res_Ab_developability _ML/Methods/設定指南.md (6 blocks)
# ---------------------------------------------

# --- Part 1 ---
# colab_setup.ipynb - First cell
!pip install claude-code
!git clone https://github.com/your-username/h2a-pipe
%cd h2a-pipe
!pip install -r requirements.txt

# Mount Google Drive for persistent storage
from google.colab import drive
drive.mount('/content/drive')

# Setup Claude Code
import os
os.environ['ANTHROPIC_API_KEY'] = 'your-key-here'

# GPU check
import torch
print(f"GPU Available: {torch.cuda.is_available()}")
print(f"GPU Name: {torch.cuda.get_device_name(0)}")

# --- Part 2 ---
# src/utils/mlflow_config.py
import mlflow
import mlflow.pytorch
from mlflow.tracking import MlflowClient

def setup_mlflow():
    mlflow.set_tracking_uri("http://localhost:5000")
    mlflow.set_experiment("h2a-pipe-experiments")
    
def log_model_performance(model, metrics, params):
    with mlflow.start_run():
        # Log parameters
        for key, value in params.items():
            mlflow.log_param(key, value)
        
        # Log metrics
        for key, value in metrics.items():
            mlflow.log_metric(key, value)
        
        # Log model
        mlflow.pytorch.log_model(model, "model")

# --- Part 3 ---
# Alternative to MLflow
import wandb

wandb.init(
    project="h2a-pipe",
    config={
        "learning_rate": 0.001,
        "architecture": "transformer",
        "dataset": "gdpa1",
        "species": "canine"
    }
)

# Log metrics during training
wandb.log({"loss": loss, "accuracy": accuracy})

# --- Part 4 ---
# src/utils/secrets.py
import os
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables
env_path = Path(__file__).parent.parent.parent / '.env'
load_dotenv(env_path)

class Config:
    ANTHROPIC_API_KEY = os.getenv('ANTHROPIC_API_KEY')
    DATABASE_URL = os.getenv('DATABASE_URL')
    
    @classmethod
    def validate(cls):
        """Ensure all required secrets are present"""
        required = ['ANTHROPIC_API_KEY', 'DATABASE_URL']
        missing = [k for k in required if not getattr(cls, k)]
        if missing:
            raise ValueError(f"Missing required secrets: {missing}")

Config.validate()

# --- Part 5 ---
# .claude-code/prompts/data_handling.txt
When working with biological sequence data:
1. Use generators for large FASTA files
2. Implement chunked processing for memory efficiency
3. Cache computed features with joblib
4. Use HDF5 for large numerical arrays
5. Implement checkpointing for long-running tasks

Example pattern:

# --- Part 6 ---
# .claude-code/prompts/validation.txt
For bioinformatics predictions, implement:

1. Biological validity checks:
   - CDR lengths within expected ranges
   - Hydrophobicity patterns consistent with antibodies
   - No stop codons in sequences

2. Statistical validation:
   - Permutation tests for feature importance
   - Bootstrap confidence intervals
   - Cross-species correlation analysis

3. Computational reproducibility:
   - Fixed random seeds
   - Environment versioning
   - Detailed logging of all parameters