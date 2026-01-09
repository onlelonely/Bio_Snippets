# ---------------------------------------------
# Title: Target Encoding
# Description: Target encoding for categorical variables with smoothing and cross-validation
# Input: DataFrame with categorical features
# Output: Encoded DataFrame
# ---------------------------------------------

import pandas as pd
from category_encoders import TargetEncoder
from sklearn.model_selection import KFold

# Basic target encoding
def basic_target_encoding(df, feature_col, target_col, smoothing=10.0, min_samples_leaf=2):
    """
    Apply basic target encoding to a categorical feature.

    Args:
        df: DataFrame containing the feature and target
        feature_col: Name of the categorical feature column
        target_col: Name of the target column
        smoothing: Smoothing parameter (higher = more regularization)
        min_samples_leaf: Minimum samples to trigger smoothing

    Returns:
        DataFrame with encoded feature
    """
    encoder = TargetEncoder(smoothing=smoothing, min_samples_leaf=min_samples_leaf)
    df[f'{feature_col}_encoded'] = encoder.fit_transform(df[feature_col], df[target_col])
    return df, encoder


# Cross-validated target encoding (prevents target leakage)
def cv_target_encoding(df, feature_col, target_col, n_splits=5, random_state=42):
    """
    Apply target encoding with KFold cross-validation to prevent target leakage.

    Args:
        df: DataFrame containing the feature and target
        feature_col: Name of the categorical feature column
        target_col: Name of the target column
        n_splits: Number of cross-validation folds
        random_state: Random seed for reproducibility

    Returns:
        DataFrame with encoded feature
    """
    df = df.copy()
    df[f'{feature_col}_encoded'] = 0.0
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    for train_idx, val_idx in kf.split(df):
        encoder = TargetEncoder()
        encoder.fit(df.loc[train_idx, feature_col], df.loc[train_idx, target_col])
        df.loc[val_idx, f'{feature_col}_encoded'] = encoder.transform(df.loc[val_idx, feature_col])

    return df


# Example usage
if __name__ == "__main__":
    # Sample data
    df = pd.DataFrame({
        'city': ['A', 'B', 'A', 'C', 'A', 'B'],
        'target': [1, 0, 1, 0, 0, 1]
    })

    # Basic encoding
    df_basic, encoder = basic_target_encoding(df.copy(), 'city', 'target')
    print("Basic encoding:")
    print(df_basic)

    # CV encoding (for larger datasets)
    # df_cv = cv_target_encoding(df, 'city', 'target')
