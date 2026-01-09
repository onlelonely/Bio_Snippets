# ---------------------------------------------
# Title: Res_Ab_developability_ML
# Description: From: Source/3. Efforts/Res_Ab_developability _ML/Res_Ab_developability_ML.md (6 blocks)
# ---------------------------------------------

# --- Part 1 ---
class SpeciesAwareAttention(nn.Module):
    """
    物種感知注意力機制
    動態調整人類與小鼠抗體的特徵權重
    """
    def __init__(self, d_model=768, n_heads=12):
        super().__init__()
        self.d_model = d_model
        self.n_heads = n_heads
        
        # 區域識別器 - 識別保守vs變異區域
        self.region_classifier = nn.Sequential(
            nn.Linear(d_model, 256),
            nn.ReLU(),
            nn.Linear(256, 4)  # VH, VL, CH, CL
        )
        
        # 物種特異性注意力頭
        self.human_attention = nn.MultiheadAttention(d_model, n_heads//2)
        self.mouse_attention = nn.MultiheadAttention(d_model, n_heads//2)
        self.conserved_attention = nn.MultiheadAttention(d_model, n_heads)
        
        # 動態融合門控
        self.fusion_gate = nn.Sequential(
            nn.Linear(d_model * 3, d_model),
            nn.Sigmoid()
        )
    
    def forward(self, x, species_embedding, sequence_mask=None):
        """
        x: [seq_len, batch, d_model] - 序列特徵
        species_embedding: [batch, d_model] - 物種編碼
        sequence_mask: [batch, seq_len] - 區域標記(VH/VL/CH/CL)
        """
        # Step 1: 識別序列區域
        region_scores = self.region_classifier(x)  # [seq_len, batch, 4]
        region_weights = F.softmax(region_scores, dim=-1)
        
        # Step 2: 計算區域特異性注意力
        # 保守區域（VH/VL）使用共享注意力
        conserved_mask = (sequence_mask == 0) | (sequence_mask == 1)  # VH or VL
        conserved_attn, _ = self.conserved_attention(
            x * conserved_mask.unsqueeze(-1), x, x
        )
        
        # 物種特異區域（CH/CL）使用分離注意力
        variable_mask = ~conserved_mask
        human_attn, _ = self.human_attention(
            x * variable_mask.unsqueeze(-1), x, x
        )
        mouse_attn, _ = self.mouse_attention(
            x * variable_mask.unsqueeze(-1), x, x
        )
        
        # Step 3: 動態融合
        species_weight = torch.sigmoid(species_embedding).unsqueeze(0)  # [1, batch, d_model]
        
        # 根據物種混合human/mouse注意力
        species_specific_attn = (species_weight * human_attn + 
                                 (1 - species_weight) * mouse_attn)
        
        # Step 4: 區域加權組合
        output = (region_weights[..., :2].sum(-1, keepdim=True) * conserved_attn +
                 region_weights[..., 2:].sum(-1, keepdim=True) * species_specific_attn)
        
        return output, region_weights

# --- Part 2 ---
class RegionAwareMasking:
    """
    根據AHo numbering為不同區域創建mask
    """
    def __init__(self, aho_annotations):
        self.aho_annotations = aho_annotations
        self.region_boundaries = self._parse_aho_regions()
    
    def _parse_aho_regions(self):
        """解析AHo編號，識別CDR和框架區"""
        boundaries = {
            'VH': {
                'FR1': (1, 25),
                'CDR1': (26, 35),
                'FR2': (36, 49),
                'CDR2': (50, 65),
                'FR3': (66, 104),
                'CDR3': (105, 117),
                'FR4': (118, 128)
            },
            'VL': {
                'FR1': (1, 23),
                'CDR1': (24, 34),
                'FR2': (35, 48),
                'CDR2': (49, 56),
                'FR3': (57, 93),
                'CDR3': (94, 102),
                'FR4': (103, 113)
            }
        }
        return boundaries
    
    def create_species_differential_mask(self, sequence, species='human'):
        """
        創建物種差異mask，標記需要特別注意的位置
        """
        mask = torch.zeros(len(sequence))
        
        # CDR區域 - 高變異但功能保守
        cdr_positions = self._get_cdr_positions(sequence)
        mask[cdr_positions] = 0.5  # 中等權重
        
        # Fc區域 - 物種差異大
        if species == 'mouse':
            # 小鼠特異性位點
            mouse_specific = [221, 224, 235, 268]  # 示例位置
            mask[mouse_specific] = 1.0  # 高權重
        
        # Vernier區 - 影響VH-VL配對
        vernier_positions = [2, 48, 67, 69, 71, 93]  # VH
        mask[vernier_positions] = 0.7
        
        return mask

# --- Part 3 ---
class HierarchicalSpeciesEmbedding(nn.Module):
    """
    層次化物種編碼器
    捕捉物種間的演化關係和功能相似性
    """
    def __init__(self, d_model=768):
        super().__init__()
        
        # Level 1: 物種家族 (靈長類 vs 囓齒類)
        self.family_embedding = nn.Embedding(2, d_model//4)
        
        # Level 2: 具體物種
        self.species_embedding = nn.Embedding(10, d_model//4)  # 支援10個物種
        
        # Level 3: 抗體亞型
        self.isotype_embedding = nn.Embedding(8, d_model//4)  # IgG1-4 for human/mouse
        
        # Level 4: 實驗條件（可選）
        self.condition_embedding = nn.Linear(5, d_model//4)  # pH, temp, etc.
        
        # 融合層
        self.fusion = nn.Sequential(
            nn.Linear(d_model, d_model),
            nn.LayerNorm(d_model),
            nn.ReLU(),
            nn.Linear(d_model, d_model)
        )
        
        # 物種特異性參數
        self.species_specific_params = nn.ModuleDict({
            'human': nn.Linear(d_model, 14),  # 14個開發性指標
            'mouse': nn.Linear(d_model, 14),
            'canine': nn.Linear(d_model, 14)
        })
    
    def forward(self, family_id, species_id, isotype_id, conditions=None):
        # 組合多層次編碼
        family_emb = self.family_embedding(family_id)
        species_emb = self.species_embedding(species_id)
        isotype_emb = self.isotype_embedding(isotype_id)
        
        if conditions is not None:
            cond_emb = self.condition_embedding(conditions)
            combined = torch.cat([family_emb, species_emb, isotype_emb, cond_emb], dim=-1)
        else:
            combined = torch.cat([family_emb, species_emb, isotype_emb, 
                                 torch.zeros_like(family_emb)], dim=-1)
        
        # 生成最終物種編碼
        species_encoding = self.fusion(combined)
        
        return species_encoding
    
    def get_species_specific_projection(self, species_encoding, target_species='mouse'):
        """獲取物種特異性投影矩陣"""
        if target_species in self.species_specific_params:
            return self.species_specific_params[target_species](species_encoding)
        else:
            # 使用最近鄰物種的參數
            return self.species_specific_params['human'](species_encoding)

# --- Part 4 ---
class CrossSpeciesAntibodyModel(nn.Module):
    """
    跨物種抗體預測模型
    保留人類VH/VL知識，適應物種特異性
    """
    def __init__(self, pretrained_human_model, target_species='mouse'):
        super().__init__()
        
        # 凍結人類模型的VH/VL編碼器
        self.frozen_vh_encoder = copy.deepcopy(pretrained_human_model.vh_encoder)
        self.frozen_vl_encoder = copy.deepcopy(pretrained_human_model.vl_encoder)
        for param in self.frozen_vh_encoder.parameters():
            param.requires_grad = False
        for param in self.frozen_vl_encoder.parameters():
            param.requires_grad = False
        
        # 可訓練的適配器層
        self.vh_adapter = AdapterLayer(d_model=768, reduction_factor=16)
        self.vl_adapter = AdapterLayer(d_model=768, reduction_factor=16)
        
        # 物種特異性Fc編碼器（完全可訓練）
        self.fc_encoder = SpeciesSpecificFcEncoder(target_species=target_species)
        
        # 物種感知注意力
        self.species_attention = SpeciesAwareAttention()
        
        # 預測頭 - 不同指標可能需要不同程度的適應
        self.prediction_heads = nn.ModuleDict({
            # 保守指標（主要依賴VH/VL）
            'Tm1': nn.Linear(768, 1),  # 共享權重
            'Titer': nn.Linear(768, 1),
            
            # 物種特異指標（需要更多適應）
            'SEC_Monomer': SpeciesAdaptiveHead(768, 1, target_species),
            'AC_SINS_pH6': SpeciesAdaptiveHead(768, 1, target_species),
            'HIC': SpeciesAdaptiveHead(768, 1, target_species)
        })
    
    def forward(self, vh_seq, vl_seq, hc_seq, lc_seq, species_info):
        # Step 1: 提取保守的V區特徵
        vh_features_frozen = self.frozen_vh_encoder(vh_seq)
        vl_features_frozen = self.frozen_vl_encoder(vl_seq)
        
        # Step 2: 通過適配器微調
        vh_features = vh_features_frozen + self.vh_adapter(vh_features_frozen)
        vl_features = vl_features_frozen + self.vl_adapter(vl_features_frozen)
        
        # Step 3: 物種特異性Fc編碼
        fc_features = self.fc_encoder(hc_seq, lc_seq, species_info)
        
        # Step 4: 整合所有特徵
        combined_features = torch.cat([vh_features, vl_features, fc_features], dim=1)
        
        # Step 5: 物種感知注意力
        attended_features, attention_weights = self.species_attention(
            combined_features, species_info['embedding']
        )
        
        # Step 6: 生成預測
        predictions = {}
        for metric_name, head in self.prediction_heads.items():
            predictions[metric_name] = head(attended_features)
        
        return predictions, attention_weights

# --- Part 5 ---
class AdapterLayer(nn.Module):
    """
    輕量級適配器，最小化參數同時保持適應能力
    """
    def __init__(self, d_model, reduction_factor=16):
        super().__init__()
        d_bottleneck = d_model // reduction_factor
        
        self.down_project = nn.Linear(d_model, d_bottleneck)
        self.activation = nn.ReLU()
        self.up_project = nn.Linear(d_bottleneck, d_model)
        
        # 殘差縮放因子（可學習）
        self.scale = nn.Parameter(torch.ones(1) * 0.1)
        
    def forward(self, x):
        residual = x
        x = self.down_project(x)
        x = self.activation(x)
        x = self.up_project(x)
        return residual + self.scale * x

# --- Part 6 ---
def progressive_validation_pipeline():
    """
    三階段驗證策略
    """
    # Stage 1: 驗證V區特徵的跨物種保守性
    vh_vl_conservation_test = {
        'metric': 'cosine_similarity',
        'threshold': 0.85,  # 預期VH/VL特徵高度相似
        'datasets': ['human_test', 'mouse_bashour']
    }
    
    # Stage 2: 驗證物種特異性適應
    species_specific_test = {
        'metrics': ['HIC', 'AC_SINS'],  # 物種差異較大的指標
        'evaluation': 'correlation_with_computational_DPs',
        'expected_improvement': '>15%'  # 相比直接遷移
    }
    
    # Stage 3: 端到端預測性能
    end_to_end_test = {
        'test_set': 'mouse_holdout_20%',
        'metrics': 'all_14_GDPa1_assays',
        'success_criteria': {
            'average_R2': '>0.6',
            'critical_metrics_R2': '>0.7'  # Tm, Titer
        }
    }
    
    return [vh_vl_conservation_test, species_specific_test, end_to_end_test]