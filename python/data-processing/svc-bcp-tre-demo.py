# ---------------------------------------------
# Title: Svc_BCP TRE Demo
# Description: From: Source/3. Efforts/_Archives/Svc_BCP TRE Demo.md (26 blocks)
# ---------------------------------------------

# --- Part 1 ---
# 建立肺癌-口腔癌欄位對應字典
field_mapping = {
    # 人口學資料
    'lung_age': 'oral_age',
    'lung_gender': 'oral_sex',
    'lung_education': 'oral_edu_level',
    
    # 生活習慣
    'lung_smoking_years': 'oral_smoke_duration',
    'lung_alcohol': 'oral_drinking_habit',
    
    # 新增口腔癌特有風險因子
    'betel_nut_chewing': 'oral_betel_years',
    'oral_hygiene': 'oral_dental_visits',
    'HPV_status': 'oral_hpv_result'
}

# --- Part 2 ---
combined_cohort = {
    'cohort_id': 'unique_identifier',
    'cancer_type': ['lung', 'oral', 'control'],
    'diagnosis_status': [0, 1],  # 0: control, 1: case
    'common_features': [...],     # 共同風險因子
    'specific_features': {
        'lung': [...],            # 肺癌特有
        'oral': [...]             # 口腔癌特有
    }
}

# --- Part 3 ---
lung_cancer_snps = {
    # CHRNA3-CHRNA5-CHRNB4 cluster (強相關)
    'rs16969968': {'chr': 15, 'pos': 78898723, 'RAF': 0.35, 'OR': 1.32},
    'rs1051730': {'chr': 15, 'pos': 78894339, 'RAF': 0.33, 'OR': 1.30},
    
    # TERT region
    'rs2736100': {'chr': 5, 'pos': 1286516, 'RAF': 0.50, 'OR': 1.12},
    
    # CYP1A1
    'rs4646903': {'chr': 15, 'pos': 75041917, 'RAF': 0.35, 'OR': 1.45},
    
    # TP53
    'rs1042522': {'chr': 17, 'pos': 7676154, 'RAF': 0.27, 'OR': 1.09},
    
    # 加入10-15個已驗證的風險位點
}

# --- Part 4 ---
oral_cancer_snps = {
    # ADH1B
    'rs1229984': {'chr': 4, 'pos': 100239319, 'RAF': 0.70, 'OR': 0.59},
    
    # ALDH2
    'rs671': {'chr': 12, 'pos': 112241766, 'RAF': 0.30, 'OR': 1.67},
    
    # CYP1A1
    'rs4646903': {'chr': 15, 'pos': 75041917, 'RAF': 0.35, 'OR': 1.89},
    
    # PLCE1
    'rs2274223': {'chr': 10, 'pos': 96087873, 'RAF': 0.27, 'OR': 1.34},
    
    # 加入8-10個口腔癌相關位點
}

# --- Part 5 ---
shared_risk_snps = {
    'rs4646903': {'lung_OR': 1.45, 'oral_OR': 1.89},  # CYP1A1
    'rs1048943': {'lung_OR': 1.20, 'oral_OR': 1.35},  # CYP1A1
    'rs2031920': {'lung_OR': 1.15, 'oral_OR': 1.28},  # CYP2E1
}

# --- Part 6 ---
def simulate_genotypes(n_samples, snp_info, case_control_status):
    """
    根據疾病狀態生成真實的基因型分布
    Cases: 增加風險等位基因頻率
    Controls: 使用群體等位基因頻率
    """
    genotypes = []
    for sample in range(n_samples):
        if case_control_status[sample] == 1:  # Case
            # 根據OR調整等位基因頻率
            adjusted_freq = calculate_case_frequency(
                snp_info['RAF'], 
                snp_info['OR']
            )
        else:  # Control
            adjusted_freq = snp_info['RAF']
        
        # 生成基因型 (0/0, 0/1, 1/1)
        genotype = np.random.choice(
            [0, 1, 2], 
            p=hardy_weinberg_proportions(adjusted_freq)
        )
        genotypes.append(genotype)
    
    return genotypes

# --- Part 7 ---
# 單一癌症GWAS
gwas_lung = GWAS(
    phenotype='lung_cancer',
    covariates=['age', 'sex', 'smoking_status', 'PC1', 'PC2'],
    model='logistic'
)

gwas_oral = GWAS(
    phenotype='oral_cancer',
    covariates=['age', 'sex', 'betel_nut', 'alcohol', 'PC1', 'PC2'],
    model='logistic'
)

# Meta-analysis找共同風險位點
meta_analysis = MetaGWAS([gwas_lung, gwas_oral])

# --- Part 8 ---
# 方法1: 簡單加權法
PRS_simple = Σ(βi × dosei)

# 方法2: LDpred2
PRS_LDpred = LDpred2_auto(
    sumstats=gwas_results,
    ld_reference='EAS_LD_panel',
    h2_init=0.01,
    sparse=True
)

# 方法3: PRS-CS
PRS_CS = PRScs(
    summary_stats=gwas_results,
    ld_reference='1000G_EAS',
    phi=1e-4
)

# --- Part 9 ---
# 肺癌PRS (基於15個SNPs)
PRS_lung = calculate_PRS(
    snps=lung_cancer_snps,
    weights=effect_sizes_lung,
    genotypes=individual_genotypes
)

# 口腔癌PRS (基於12個SNPs)
PRS_oral = calculate_PRS(
    snps=oral_cancer_snps,
    weights=effect_sizes_oral,
    genotypes=individual_genotypes
)

# 整合PRS
PRS_combined = α × PRS_lung + β × PRS_oral + γ × PRS_shared

# --- Part 10 ---
prs_categories = {
    'Low': 'PRS < 20th percentile',
    'Intermediate': '20th ≤ PRS < 80th percentile',
    'High': 'PRS ≥ 80th percentile',
    'Very High': 'PRS ≥ 95th percentile'
}

# --- Part 11 ---
integrated_features = {
    # Level 1: 基因體層次
    'genomic': {
        'PRS_lung': continuous,
        'PRS_oral': continuous,
        'risk_snp_count': discrete,
        'protective_snp_count': discrete
    },
    
    # Level 2: 環境暴露層次
    'environmental': {
        'smoking_pack_years': continuous,
        'alcohol_frequency': ordinal,
        'betel_nut_years': continuous,
        'air_pollution_index': continuous,
        'occupational_hazard': ordinal
    },
    
    # Level 3: 生活型態層次
    'lifestyle': {
        'diet_score': continuous,
        'exercise_frequency': ordinal,
        'sleep_quality': ordinal,
        'stress_level': ordinal
    },
    
    # Level 4: 臨床指標層次
    'clinical': {
        'BMI': continuous,
        'FEV1_FVC_ratio': continuous,
        'inflammatory_markers': continuous,
        'oral_lesions': binary
    }
}

# --- Part 12 ---
# G×E交互作用項
interaction_features = {
    'PRS_smoking': PRS_lung × smoking_status,
    'PRS_alcohol': PRS_oral × alcohol_consumption,
    'CHRNA5_smoking': rs16969968 × pack_years,
    'ALDH2_alcohol': rs671 × drinking_frequency,
    'CYP1A1_pollution': rs4646903 × air_quality_index
}

# --- Part 13 ---
# Step 1: 分層特徵選擇
selected_genomic = select_features(genomic_features, method='LASSO')
selected_environmental = select_features(env_features, method='RFE')
selected_clinical = select_features(clinical_features, method='mRMR')

# Step 2: 整合特徵選擇
integrated_selection = ensemble_feature_selection([
    selected_genomic,
    selected_environmental,
    selected_clinical
], voting_threshold=0.6)

# Step 3: 交互作用項篩選
interaction_selection = interaction_filter(
    features=interaction_features,
    method='mutual_information',
    threshold=0.1
)

# --- Part 14 ---
def analyze_risk_deterioration(model, prs_score, environmental_factors):
    """
    分析相同PRS下，不良環境因子如何影響風險
    """
    # 基準風險（良好生活習慣）
    baseline_env = {
        'smoking': 0,
        'alcohol': 0,
        'exercise': 1,
        'diet_quality': 1
    }
    baseline_risk = model.predict_proba([prs_score, baseline_env])
    
    # 逐項惡化分析
    deterioration_analysis = {}
    for factor in environmental_factors:
        worst_case_env = baseline_env.copy()
        worst_case_env[factor] = get_worst_value(factor)
        worst_risk = model.predict_proba([prs_score, worst_case_env])
        
        deterioration_analysis[factor] = {
            'baseline_risk': baseline_risk,
            'deteriorated_risk': worst_risk,
            'risk_increase': worst_risk - baseline_risk,
            'relative_increase': worst_risk / baseline_risk
        }
    
    return deterioration_analysis

# --- Part 15 ---
# 3D風險表面圖
def plot_risk_surface(model, prs_range, env_factor_range):
    """
    繪製PRS與環境因子的風險交互作用3D表面
    """
    PRS, ENV = np.meshgrid(prs_range, env_factor_range)
    RISK = model.predict_proba(np.c_[PRS.ravel(), ENV.ravel()])
    
    fig = go.Figure(data=[go.Surface(
        x=PRS,
        y=ENV,
        z=RISK.reshape(PRS.shape),
        colorscale='RdYlBu_r'
    )])
    
    fig.update_layout(
        title='PRS-環境因子風險交互作用',
        scene=dict(
            xaxis_title='PRS Score',
            yaxis_title='Environmental Risk',
            zaxis_title='Cancer Risk'
        )
    )
    return fig

# --- Part 16 ---
class TemporalRiskModel:
    def __init__(self):
        self.trajectory_model = XGBRegressor(
            objective='reg:squarederror',
            n_estimators=1000,
            early_stopping_rounds=50
        )
    
    def predict_risk_trajectory(self, individual_profile, time_points):
        """預測個體在不同時間點的風險變化"""
        trajectories = []
        for t in time_points:
            # 考慮年齡增長
            age_adjusted_profile = adjust_for_age(individual_profile, t)
            # 考慮累積暴露
            exposure_adjusted = adjust_for_cumulative_exposure(
                age_adjusted_profile, t
            )
            risk = self.trajectory_model.predict(exposure_adjusted)
            trajectories.append(risk)
        
        return trajectories

# --- Part 17 ---
"""
挑戰：在沒有基因檢測的情況下，從詳細的人口學和問卷資料推斷PRS風險分層
理論基礎：某些表型特徵和生活習慣可能與遺傳風險相關
"""

class PRSImputationModel:
    def __init__(self):
        # 使用深度學習捕捉複雜的非線性關係
        self.imputation_network = Sequential([
            Dense(512, activation='relu', input_dim=n_features),
            Dropout(0.3),
            Dense(256, activation='relu'),
            BatchNormalization(),
            Dense(128, activation='relu'),
            Dropout(0.2),
            Dense(64, activation='relu'),
            Dense(3, activation='softmax')  # Low/Medium/High PRS
        ])
        
        # 輔助的梯度提升模型
        self.gb_model = HistGradientBoostingClassifier(
            max_iter=1000,
            learning_rate=0.05,
            max_depth=8
        )

# --- Part 18 ---
prs_predictive_features = {
    # 家族史特徵（最強預測因子）
    'family_history': {
        'n_affected_relatives': count,
        'closest_affected_relative': ordinal,  # 1st/2nd/3rd degree
        'age_at_diagnosis_relatives': continuous,
        'multiple_cancer_types': binary
    },
    
    # 族裔背景（與等位基因頻率相關）
    'ethnicity_features': {
        'paternal_origin': categorical,
        'maternal_origin': categorical,
        'ethnicity_admixture': continuous  # 如果有資料
    },
    
    # 早發性疾病指標
    'early_onset_indicators': {
        'age_first_symptoms': continuous,
        'precancerous_lesions': binary,
        'benign_tumors_history': count
    },
    
    # 藥物反應（可能反映基因變異）
    'pharmacogenomic_hints': {
        'alcohol_flush_reaction': binary,  # ALDH2相關
        'smoking_addiction_level': ordinal,  # CHRNA5相關
        'caffeine_sensitivity': ordinal
    },
    
    # 生理特徵
    'physiological_markers': {
        'baseline_inflammation': continuous,
        'metabolic_rate': continuous,
        'immune_response_pattern': categorical
    }
}

# --- Part 19 ---
def quantify_prediction_uncertainty(model, features):
    """
    提供PRS預測的信心區間
    """
    # Monte Carlo Dropout
    predictions = []
    for _ in range(100):
        pred = model.predict(features, training=True)  # Keep dropout on
        predictions.append(pred)
    
    mean_pred = np.mean(predictions, axis=0)
    std_pred = np.std(predictions, axis=0)
    
    confidence_levels = {
        'High_confidence': std_pred < 0.1,
        'Medium_confidence': 0.1 <= std_pred < 0.2,
        'Low_confidence': std_pred >= 0.2
    }
    
    return mean_pred, confidence_levels

# --- Part 20 ---
def interactive_risk_calculator():
    """
    使用者輸入資料，即時計算並視覺化風險
    """
    # 輸入介面
    st.sidebar.header("個人資料輸入")
    
    # 基本資料
    age = st.sidebar.slider("年齡", 20, 80, 50)
    sex = st.sidebar.selectbox("性別", ["男", "女"])
    
    # 家族史
    family_history = st.sidebar.selectbox(
        "肺癌家族史", 
        ["無", "二級親屬", "一級親屬", "多位親屬"]
    )
    
    # 生活習慣
    smoking = st.sidebar.selectbox(
        "吸菸史",
        ["從未", "<10包年", "10-20包年", ">20包年"]
    )
    
    # 環境暴露
    occupation = st.sidebar.selectbox(
        "職業暴露",
        ["低風險", "中風險", "高風險"]
    )
    
    # 計算風險
    if st.button("計算風險"):
        # 推斷PRS
        imputed_prs = model.impute_prs(user_data)
        
        # 計算綜合風險
        risk_score = model.calculate_risk(user_data, imputed_prs)
        
        # 視覺化結果
        display_risk_gauge(risk_score)
        show_risk_factors_breakdown()
        provide_recommendations()

# --- Part 21 ---
def vcf_interactive_viewer():
    """展示VCF檔案結構和內容"""
    # 顯示VCF header
    st.code(vcf_header, language='text')
    
    # 互動式變異位點展示
    selected_snp = st.selectbox("選擇SNP", snp_list)
    
    # 顯示該SNP的詳細資訊
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("位置", f"chr{snp_info['chr']}:{snp_info['pos']}")
    with col2:
        st.metric("風險等位基因頻率", f"{snp_info['RAF']:.2%}")
    with col3:
        st.metric("Odds Ratio", f"{snp_info['OR']:.2f}")
    
    # 顯示基因型分布
    plot_genotype_distribution(selected_snp)

# --- Part 22 ---
def risk_deterioration_simulator():
    """模擬環境因子對PRS風險的影響"""
    
    # 設定基準PRS
    prs_level = st.select_slider(
        "PRS風險等級",
        options=["低", "中", "高"]
    )
    
    # 環境因子調整
    st.subheader("調整環境因子")
    smoking = st.slider("吸菸(包年)", 0, 50, 0)
    alcohol = st.slider("飲酒頻率", 0, 7, 0)
    exercise = st.slider("運動(小時/週)", 0, 20, 5)
    
    # 即時計算並顯示風險變化
    baseline_risk = calculate_baseline_risk(prs_level)
    adjusted_risk = calculate_adjusted_risk(
        prs_level, smoking, alcohol, exercise
    )
    
    # 動態風險條
    risk_change = adjusted_risk - baseline_risk
    st.metric(
        "風險變化",
        f"{adjusted_risk:.1%}",
        delta=f"{risk_change:+.1%}"
    )
    
    # 風險軌跡圖
    plot_risk_trajectory(baseline_risk, adjusted_risk)

# --- Part 23 ---
def explainability_dashboard():
    """整合式模型解釋介面"""
    
    tab1, tab2, tab3 = st.tabs(["SHAP", "LIME", "特徵重要性"])
    
    with tab1:
        # SHAP瀑布圖
        individual_id = st.selectbox("選擇個案", case_list)
        plot_shap_waterfall(individual_id)
        
    with tab2:
        # LIME解釋
        show_lime_explanation(individual_id)
        
    with tab3:
        # 特徵重要性比較
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("肺癌模型")
            plot_feature_importance('lung')
        with col2:
            st.subheader("口腔癌模型")
            plot_feature_importance('oral')

# --- Part 24 ---
# 測試模型在不同癌症類型的泛化能力
cross_cancer_validation = {
    'lung_to_oral': test_lung_model_on_oral_data(),
    'oral_to_lung': test_oral_model_on_lung_data(),
    'combined_model': test_unified_model()
}

# --- Part 25 ---
# 評估從人口學資料推斷PRS的準確度
imputation_metrics = {
    'accuracy': classification_accuracy,
    'auc': roc_auc_score,
    'correlation': pearson_correlation,
    'calibration': calibration_error
}

# --- Part 26 ---
screening_recommendations = {
    'Very_High_Risk': {
        'PRS': 'Top 5%',
        'Environmental': 'High',
        'Screening': 'Annual LDCT + Oral examination',
        'Intervention': 'Intensive lifestyle modification'
    },
    'High_Risk': {
        'PRS': 'Top 20%',
        'Environmental': 'Moderate-High',
        'Screening': 'Biennial LDCT',
        'Intervention': 'Lifestyle counseling'
    },
    'Moderate_Risk': {
        'PRS': 'Middle 60%',
        'Environmental': 'Moderate',
        'Screening': 'Standard guidelines',
        'Intervention': 'General health promotion'
    }
}