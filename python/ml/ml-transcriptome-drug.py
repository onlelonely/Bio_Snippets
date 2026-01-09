# ---------------------------------------------
# Title: ml_transcriptome_drug
# Description: From: Source/3. Efforts/Res_Drug_&_Transcriptome/討論記錄/ml_transcriptome_drug.md (7 blocks)
# ---------------------------------------------

# --- Part 1 ---
class TranscriptomeFeatureEngineering:
    """轉錄體數據的特徵工程"""
    
    def __init__(self, expression_data):
        self.expr_data = expression_data  # genes x samples
        
    def extract_basic_features(self):
        """基礎統計特徵"""
        features = {
            # 表現量特徵
            'mean_expression': np.mean(self.expr_data, axis=1),
            'variance': np.var(self.expr_data, axis=1),
            'cv': np.std(self.expr_data, axis=1) / np.mean(self.expr_data, axis=1),
            
            # 動態範圍
            'dynamic_range': np.max(self.expr_data, axis=1) - np.min(self.expr_data, axis=1),
            'fold_change': np.max(self.expr_data, axis=1) / (np.min(self.expr_data, axis=1) + 1e-6),
            
            # 分位數特徵
            'q25': np.percentile(self.expr_data, 25, axis=1),
            'q75': np.percentile(self.expr_data, 75, axis=1),
            'iqr': np.percentile(self.expr_data, 75, axis=1) - np.percentile(self.expr_data, 25, axis=1)
        }
        return features
    
    def extract_pathway_features(self):
        """通路層級特徵（更有生物意義）"""
        pathway_features = {}
        
        # 使用GSEA或ssGSEA計算通路活性分數
        pathways = {
            'immune_response': ['IFNB1', 'IFNA1', 'IL6', 'TNF'],
            'iron_metabolism': ['TFRC', 'TF', 'FTH1', 'FTL'],
            'cell_cycle': ['CDKN1A', 'CDKN2A', 'PCNA', 'MKI67'],
            'apoptosis': ['CASP3', 'CASP9', 'BCL2', 'BAX']
        }
        
        for pathway_name, genes in pathways.items():
            # 計算通路活性分數（簡化版ssGSEA）
            pathway_genes_expr = self.expr_data[self.expr_data.index.isin(genes)]
            pathway_features[f'{pathway_name}_score'] = np.mean(pathway_genes_expr, axis=0)
            pathway_features[f'{pathway_name}_coherence'] = np.corrcoef(pathway_genes_expr)[0,1:].mean()
        
        return pathway_features
    
    def extract_network_features(self):
        """網路拓撲特徵"""
        import networkx as nx
        
        # 建立共表現網路
        corr_matrix = np.corrcoef(self.expr_data)
        G = nx.from_numpy_array(np.abs(corr_matrix) > 0.7)
        
        network_features = {
            'degree_centrality': nx.degree_centrality(G),
            'betweenness_centrality': nx.betweenness_centrality(G),
            'closeness_centrality': nx.closeness_centrality(G),
            'clustering_coefficient': nx.clustering(G),
            'pagerank': nx.pagerank(G)
        }
        
        # 模組特徵
        from sklearn.cluster import SpectralClustering
        clustering = SpectralClustering(n_clusters=10, affinity='precomputed')
        modules = clustering.fit_predict(corr_matrix)
        
        network_features['module_membership'] = modules
        network_features['module_connectivity'] = self.calculate_module_connectivity(modules, corr_matrix)
        
        return network_features

# --- Part 2 ---
# PCA - 最快但可能丟失非線性關係
from sklearn.decomposition import PCA
pca = PCA(n_components=50)
pca_features = pca.fit_transform(expression_data)

# Truncated SVD - 適合稀疏數據
from sklearn.decomposition import TruncatedSVD
svd = TruncatedSVD(n_components=50)
svd_features = svd.fit_transform(expression_data)

# NMF - 保證非負，適合基因表現數據
from sklearn.decomposition import NMF
nmf = NMF(n_components=50)
nmf_features = nmf.fit_transform(expression_data)

# --- Part 3 ---
# UMAP - 保持局部和全局結構
import umap
reducer = umap.UMAP(n_components=50, n_neighbors=15, min_dist=0.1)
umap_features = reducer.fit_transform(expression_data)

# t-SNE - 主要用於視覺化
from sklearn.manifold import TSNE
tsne = TSNE(n_components=2, perplexity=30)
tsne_features = tsne.fit_transform(expression_data)

# AutoEncoder - 深度學習降維
import torch
import torch.nn as nn

class TranscriptomeAutoEncoder(nn.Module):
    def __init__(self, input_dim=20000, latent_dim=128):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 2048),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(2048, 512),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(512, latent_dim)
        )
        
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 512),
            nn.ReLU(),
            nn.Linear(512, 2048),
            nn.ReLU(),
            nn.Linear(2048, input_dim)
        )
    
    def forward(self, x):
        latent = self.encoder(x)
        reconstructed = self.decoder(latent)
        return reconstructed, latent

# --- Part 4 ---
class TraditionalMLModels:
    """傳統ML模型比較"""
    
    def train_all_models(self, X_train, y_train, X_test, y_test):
        from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
        from sklearn.svm import SVC
        from xgboost import XGBClassifier
        from lightgbm import LGBMClassifier
        from catboost import CatBoostClassifier
        
        models = {
            'RandomForest': RandomForestClassifier(
                n_estimators=200,
                max_depth=20,
                min_samples_split=5
            ),
            
            'XGBoost': XGBClassifier(
                n_estimators=200,
                max_depth=6,
                learning_rate=0.1,
                subsample=0.8
            ),
            
            'LightGBM': LGBMClassifier(
                n_estimators=200,
                num_leaves=31,
                learning_rate=0.1
            ),
            
            'CatBoost': CatBoostClassifier(
                iterations=200,
                depth=6,
                learning_rate=0.1,
                verbose=False
            ),
            
            'SVM': SVC(
                kernel='rbf',
                C=1.0,
                probability=True
            )
        }
        
        results = {}
        for name, model in models.items():
            model.fit(X_train, y_train)
            predictions = model.predict_proba(X_test)[:, 1]
            results[name] = self.evaluate_model(y_test, predictions)
            
        return results

# --- Part 5 ---
class CPVDrugPredictionPipeline:
    """CPV藥物預測的完整流程"""
    
    def __init__(self):
        self.feature_extractor = TranscriptomeFeatureEngineering()
        self.drug_encoder = DrugFeatureEngineering()
        self.models = {}
        
    def prepare_features(self, gdpx_data, lopac_drugs):
        """準備訓練數據"""
        
        # 1. 提取已知抗病毒藥物的特徵
        known_antivirals = ['Ribavirin', 'Nitazoxanide', 'Remdesivir']
        positive_samples = []
        
        for drug in known_antivirals:
            drug_response = gdpx_data[drug]
            
            # 轉錄體特徵
            trans_features = self.feature_extractor.extract_basic_features(drug_response)
            pathway_features = self.feature_extractor.extract_pathway_features(drug_response)
            
            # 藥物特徵
            drug_features = self.drug_encoder.extract_molecular_features(drug.smiles)
            
            # 組合特徵
            combined = {**trans_features, **pathway_features, **drug_features}
            positive_samples.append(combined)
        
        # 2. 負樣本（非抗病毒藥物）
        negative_samples = self.prepare_negative_samples(gdpx_data)
        
        return positive_samples, negative_samples
    
    def train_ensemble(self, X_train, y_train):
        """訓練集成模型"""
        
        # 1. XGBoost（通常表現最好）
        self.models['xgboost'] = XGBClassifier(
            n_estimators=300,
            max_depth=8,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8
        )
        
        # 2. Random Forest（穩定性好）
        self.models['rf'] = RandomForestClassifier(
            n_estimators=500,
            max_depth=None,
            min_samples_split=5
        )
        
        # 3. Neural Network（捕捉複雜關係）
        self.models['nn'] = self.build_neural_network()
        
        # 訓練所有模型
        for name, model in self.models.items():
            print(f"Training {name}...")
            model.fit(X_train, y_train)
        
        return self.models
    
    def predict_cpv_drugs(self, candidate_drugs):
        """預測CPV藥物效果"""
        
        predictions = {}
        
        for drug in candidate_drugs:
            # 提取特徵
            features = self.extract_drug_features(drug)
            
            # 集成預測
            pred_scores = []
            for model in self.models.values():
                pred_scores.append(model.predict_proba(features)[0, 1])
            
            # 加權平均（可以根據驗證集調整權重）
            weights = [0.4, 0.3, 0.3]  # XGBoost, RF, NN
            final_score = np.average(pred_scores, weights=weights)
            
            predictions[drug] = {
                'score': final_score,
                'individual_scores': pred_scores,
                'confidence': np.std(pred_scores)  # 用標準差評估預測一致性
            }
        
        return predictions

# --- Part 6 ---
cpv_recommended_pipeline = {
    "Phase 1": {
        "方法": "PCA + XGBoost",
        "原因": "快速、可解釋、不需要大量數據",
        "預期準確率": "75-85%"
    },
    
    "Phase 2": {
        "方法": "AutoEncoder + 集成學習",
        "原因": "更好的特徵學習、更穩定",
        "預期準確率": "85-90%"
    },
    
    "Phase 3": {
        "方法": "GNN + Attention機制",
        "原因": "最大化利用所有資訊",
        "預期準確率": "90%+"
    }
}

# --- Part 7 ---
# 立即可以開始的簡單實作
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from sklearn.metrics import roc_auc_score

# 載入GDPx數據
gdpx_data = pd.read_csv('gdpx_processed.csv')

# 簡單特徵工程
features = []
for drug in gdpx_data.columns:
    drug_features = {
        'mean_expr': gdpx_data[drug].mean(),
        'std_expr': gdpx_data[drug].std(),
        'top10_genes': gdpx_data[drug].nlargest(10).mean(),
        'bottom10_genes': gdpx_data[drug].nsmallest(10).mean()
    }
    features.append(drug_features)

# 標記已知抗病毒藥物
labels = [1 if drug in known_antivirals else 0 for drug in gdpx_data.columns]

# 訓練XGBoost
X_train, X_test, y_train, y_test = train_test_split(features, labels)
model = XGBClassifier(n_estimators=100)
model.fit(X_train, y_train)

# 評估
predictions = model.predict_proba(X_test)[:, 1]
auc = roc_auc_score(y_test, predictions)
print(f"AUC: {auc:.3f}")

# 預測新藥物
new_drug_predictions = model.predict_proba(new_drug_features)