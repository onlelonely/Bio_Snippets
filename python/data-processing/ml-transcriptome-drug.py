# ---------------------------------------------
# Title: ml_transcriptome_drug
# Description: From: Source/3. Efforts/Res_Drug_&_Transcriptome/討論記錄/ml_transcriptome_drug.md (4 blocks)
# ---------------------------------------------

# --- Part 1 ---
class DrugFeatureEngineering:
    """藥物分子特徵提取"""
    
    def extract_molecular_features(self, smiles):
        """從SMILES提取分子特徵"""
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski
        
        mol = Chem.MolFromSmiles(smiles)
        
        features = {
            # Lipinski's Rule of Five
            'MW': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'HBA': Descriptors.NumHAcceptors(mol),
            'TPSA': Descriptors.TPSA(mol),
            
            # 其他重要描述符
            'RotatableBonds': Descriptors.NumRotatableBonds(mol),
            'AromaticRings': Descriptors.NumAromaticRings(mol),
            'QED': Descriptors.qed(mol),  # 藥物相似性評分
            
            # 複雜度
            'BertzCT': Descriptors.BertzCT(mol),
            'MolecularComplexity': len(mol.GetBonds())
        }
        
        return features
    
    def extract_fingerprints(self, smiles):
        """分子指紋特徵"""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        
        fingerprints = {
            'morgan_fp': AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
            'maccs_keys': AllChem.GetMACCSKeysFingerprint(mol),
            'topological_fp': Chem.RDKFingerprint(mol)
        }
        
        return fingerprints

# --- Part 2 ---
class DeepTranscriptomeModel(nn.Module):
    """專門設計用於轉錄體數據的深度網路"""
    
    def __init__(self, input_dim, num_drugs, dropout=0.3):
        super().__init__()
        
        # 基因表現編碼器
        self.gene_encoder = nn.Sequential(
            nn.Linear(input_dim, 2048),
            nn.BatchNorm1d(2048),
            nn.ReLU(),
            nn.Dropout(dropout),
            
            nn.Linear(2048, 1024),
            nn.BatchNorm1d(1024),
            nn.ReLU(),
            nn.Dropout(dropout),
            
            nn.Linear(1024, 512),
            nn.BatchNorm1d(512),
            nn.ReLU()
        )
        
        # 藥物特徵編碼器
        self.drug_encoder = nn.Sequential(
            nn.Linear(2048, 512),  # Morgan fingerprint dimension
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(dropout),
            
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.ReLU()
        )
        
        # 交互作用預測器
        self.interaction_predictor = nn.Sequential(
            nn.Linear(512 + 256, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(dropout),
            
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(dropout),
            
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            
            nn.Linear(128, 1),
            nn.Sigmoid()
        )
        
    def forward(self, gene_expr, drug_features):
        gene_embedding = self.gene_encoder(gene_expr)
        drug_embedding = self.drug_encoder(drug_features)
        
        # 結合兩種特徵
        combined = torch.cat([gene_embedding, drug_embedding], dim=1)
        
        # 預測藥物反應
        response = self.interaction_predictor(combined)
        
        return response

# --- Part 3 ---
import torch_geometric
from torch_geometric.nn import GCNConv, global_mean_pool

class DrugTargetGNN(nn.Module):
    """圖神經網路用於藥物-標靶預測"""
    
    def __init__(self, num_features, hidden_dim=256):
        super().__init__()
        
        # 圖卷積層
        self.conv1 = GCNConv(num_features, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, hidden_dim)
        self.conv3 = GCNConv(hidden_dim, 128)
        
        # 預測層
        self.predictor = nn.Sequential(
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )
        
    def forward(self, x, edge_index, batch):
        # 圖卷積
        x = F.relu(self.conv1(x, edge_index))
        x = F.dropout(x, p=0.2, training=self.training)
        
        x = F.relu(self.conv2(x, edge_index))
        x = F.dropout(x, p=0.2, training=self.training)
        
        x = F.relu(self.conv3(x, edge_index))
        
        # 圖池化
        x = global_mean_pool(x, batch)
        
        # 預測
        return self.predictor(x)

# --- Part 4 ---
def analyze_feature_importance(model, feature_names):
    """分析哪些特徵最重要"""
    
    if hasattr(model, 'feature_importances_'):
        # Tree-based models
        importances = model.feature_importances_
    elif hasattr(model, 'coef_'):
        # Linear models
        importances = np.abs(model.coef_[0])
    else:
        # Use SHAP for any model
        import shap
        explainer = shap.Explainer(model)
        shap_values = explainer(X_test)
        importances = np.abs(shap_values.values).mean(0)
    
    # 排序並視覺化
    feature_importance_df = pd.DataFrame({
        'feature': feature_names,
        'importance': importances
    }).sort_values('importance', ascending=False)
    
    # Top 20特徵
    top_features = feature_importance_df.head(20)
    
    # 生物學解釋
    biological_interpretation = {
        'TFRC_expression': '鐵受體表現量 - CPV進入細胞的關鍵',
        'interferon_pathway_score': '干擾素通路活性 - 抗病毒免疫',
        'cell_cycle_score': '細胞週期活性 - CPV複製需求',
        'iron_metabolism_coherence': '鐵代謝基因協同性'
    }
    
    return top_features, biological_interpretation