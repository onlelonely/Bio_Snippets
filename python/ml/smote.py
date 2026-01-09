# ---------------------------------------------
# Title: SMOTE
# Description: From: Source/1. Atlas/📊 Methods & Statistics/Data Processing/SMOTE.md
# ---------------------------------------------

import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.linear_model import LassoCV
from sklearn.linear_model import Lasso
from sklearn.feature_selection import RFE
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from imblearn.over_sampling import SMOTE
from sklearn.utils import shuffle
import numpy as np

def preprocess_data(file_path, info_path, Vthreshold=0.95):
    data = pd.read_csv(file_path)
    info = pd.read_csv(info_path)

    labels = data.iloc[:, 0].str.replace('.CEL', '')
    feature_data = data.iloc[:, 1:]

    columns_to_drop = [col for col in feature_data.columns if not col.startswith('T')]
    feature_data = feature_data.drop(columns=columns_to_drop)

    selector = VarianceThreshold(threshold=Vthreshold)
    filtered_data = selector.fit_transform(feature_data)
    mask = selector.get_support()
    original_columns = feature_data.columns
    filtered_columns = original_columns[mask]
    filtered_data_df = pd.DataFrame(filtered_data, columns=filtered_columns)

    scaler = StandardScaler()
    filtered_data_scaled = scaler.fit_transform(filtered_data_df)
    
    info['Drug_Response'] = info['Drug'] + '_' + info['Response']
    mapping = {'ASI_Good': 1, 'ARB_Good': 2, 'ASI_Poor': 3, 'ARB_Poor': 4}
    info['Numerical_Group'] = info['Drug_Group'].map(mapping)
    groups = info['Numerical_Group']

    return filtered_data_scaled, groups, filtered_columns


def run_smote(iter_time, features, groups, fold_1, fold_2):
    dataset_stats = []
    results = []
    
    pca = PCA(n_components=2)
    features_pca = pca.fit_transform(features)
    
    for iteration in range(iter_time):
        smote_2_fold = SMOTE(sampling_strategy=fold_1, k_neighbors=2)
        features_resampled_2_fold, groups_resampled_2_fold = smote_2_fold.fit_resample(features, groups)

        smote_10_fold = SMOTE(sampling_strategy=fold_2, k_neighbors=4)
        features_resampled_10_fold, groups_resampled_10_fold = smote_10_fold.fit_resample(features_resampled_2_fold, groups_resampled_2_fold)

        features_resampled, groups_resampled = shuffle(features_resampled_10_fold, groups_resampled_10_fold)

        features_resampled_pca = pca.transform(features_resampled)

        X_train_resampled, X_test_resampled, y_train_resampled, y_test_resampled = train_test_split(features_resampled_pca, groups_resampled, test_size=0.2, random_state=iteration)

        clf_post = RandomForestClassifier()
        clf_post.fit(X_train_resampled, y_train_resampled)

        accuracy = clf_post.score(X_test_resampled, y_test_resampled)

        results.append((accuracy, features_resampled, groups_resampled))

        mean_features_resampled = np.mean(features_resampled, axis=0)
        std_features_resampled = np.std(features_resampled, axis=0)

        dataset_stats.append((mean_features_resampled, std_features_resampled, accuracy, features_resampled, groups_resampled))
        
    return dataset_stats, results


def find_most_dissimilar_datasets(dataset_stats):
    def dissimilarity_cosine(dataset1, dataset2):
        mean1, std1, _, _, _ = dataset1
        mean2, std2, _, _, _ = dataset2
        cosine_sim_mean = cosine_similarity(mean1.reshape(1, -1), mean2.reshape(1, -1))
        cosine_sim_std = cosine_similarity(std1.reshape(1, -1), std2.reshape(1, -1))
        return (1 - cosine_sim_mean) + (1 - cosine_sim_std)

    max_distance = 0
    most_dissimilar_datasets = (None, None)

    for i in range(len(dataset_stats)):
        for j in range(i + 1, len(dataset_stats)):
            distance = dissimilarity_cosine(dataset_stats[i], dataset_stats[j])
            if distance > max_distance:
                max_distance = distance
                most_dissimilar_datasets = (dataset_stats[i], dataset_stats[j])

    return most_dissimilar_datasets

# t-SNE visualization
def plot_tsne(features, labels, title, map, save_filename=None):
    tsne = TSNE(n_components=2, random_state=0)
    features_tsne = tsne.fit_transform(features)
    
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(features_tsne[:, 0], features_tsne[:, 1], c=labels, cmap='viridis', alpha=0.6)
    
    # Reverse the mapping dictionary
    reverse_mapping = {v: k for k, v in map.items()}
    
    # Create legend labels list
    legend_labels = [reverse_mapping[i] for i in np.sort(np.unique(labels))]
    
    # Add legend
    legend1 = plt.legend(handles=scatter.legend_elements()[0], title="Classes", labels=legend_labels)
    plt.title(title)
    
    if save_filename:
        plt.savefig(save_filename, format='png')  # Save the plot as a PNG file
    plt.show()

# Main Execution
file_path = 'expr_sample.csv'
info_path = 'traits.csv'
iter_time = 100
fold_1 = {1: 18, 2: 18, 3: 8, 4: 30}
fold_2 = {1: 250, 2: 250, 3: 250, 4: 250}
mapping = {'ASI_Good': 1, 'ARB_Good': 2, 'ASI_Poor': 3, 'ARB_Poor': 4}

features, groups, filtered_columns = preprocess_data(file_path, info_path, Vthreshold=0.9)
dataset_stats, results = run_smote(iter_time, features, groups, fold_1, fold_2)
most_dissimilar_datasets = find_most_dissimilar_datasets(dataset_stats)

# Extract the selected results
selected_results = [most_dissimilar_datasets[0][-3:], most_dissimilar_datasets[1][-3:]]

# Apply t-SNE and plot for the selected 2 most dissimilar datasets
for i, result in enumerate(selected_results):
    _, features_resampled, groups_resampled = result
    plot_tsne(features_resampled, groups_resampled, title=f"t-SNE for Selected Dataset {i+1}", map=mapping, save_filename=f"tsne_plot_dataset_{i+1}.png")

# RFE, and save the results before and after RFE
for i, result in enumerate(selected_results):
    accuracy, features_resampled, groups_resampled = result
    resampled_data = pd.DataFrame(features_resampled, columns=filtered_columns)  # Assuming features_resampled is an array-like structure
    resampled_data.insert(0, 'Labels', groups_resampled)
    resampled_data.to_csv(f'selected_{i+1}_resampled_data.csv', index=False)

    lasso_cv = LassoCV(cv=5, max_iter=10000) 
    lasso_cv.fit(features_resampled, groups_resampled)
    # Get the best alpha value
    best_alpha = lasso_cv.alpha_
    print("Optimal alpha:", best_alpha)

    # Create Lasso regressor
    lasso = Lasso(alpha=best_alpha, max_iter=1000)
    # Create the RFE model and select top N features
    rfe = RFE(lasso, n_features_to_select=50) 
    fit = rfe.fit(features_resampled, groups_resampled)

    # Transform data using RFE
    filtered_data_rfe = fit.transform(features_resampled)

    # Plot t-SNE after RFE
    plot_tsne(filtered_data_rfe, groups_resampled, title=f"t-SNE after RFEfor Selected Dataset {i+1}", map=mapping)
    plt.show()

    # Get the column names of the features that were selected by RFE
    selected_columns = features_resampled.columns[selector.get_support()][fit.support_]

    # Create a DataFrame with the filtered data from RFE and the labels
    filtered_data_rfe_df = pd.DataFrame(filtered_data_rfe, columns=selected_columns)
    filtered_data_rfe_df.insert(0, 'Labels', groups_resampled)  # Insert 'Labels' column at the beginning

		# Save the DataFrame to a CSV file
    filtered_data_rfe_df.to_csv(f"filtered_data for Selected Dataset_{i+1}.csv", index=False)