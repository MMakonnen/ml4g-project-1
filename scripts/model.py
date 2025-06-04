import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import lightgbm as lgb
from scipy.stats import spearmanr

X_1_train_path = "data/X1-train/features.tsv"
y_1_train_path = "data/X1-train/y.tsv"

X_1_val_path = "data/X1-val/features.tsv"
y_1_val_path = "data/X1-val/y.tsv"

X_2_train_path = "data/X2-train/features.tsv"
y_2_train_path = "data/X2-train/y.tsv"

X_2_val_path = "data/X2-val/features.tsv"
y_2_val_path = "data/X2-val/y.tsv"

X_3_test_path_orig = "data/CAGE-train/X3_test_info.tsv"
X_3_test_path = "data/X3-test/features.tsv"

cols_to_keep = [
    "DNase_avg_int",
    "DNase_avg_peaks",
    "H3K4me1_avg_int",
    "H3K4me1_avg_peaks",
    "H3K4me3_avg_int",
    "H3K4me3_avg_peaks",
    "H3K27ac_avg_int",
    "H3K27ac_avg_peaks",
]

# Parameters for LightGBM with Huber loss
lightgbm_params = {
    "boosting_type": "gbdt",
    "objective": "huber",
    "alpha": 0.9,
    "learning_rate": 0.40,
    "num_leaves": 73,
    "n_estimators": 105,
}


def extend_df(df, cols_to_keep):
    # Add a binary column for strand
    df["strand_binary"] = df["strand"].map({"+": 1, "-": 0})

    # gene length
    df["gene_length"] = df["gene_end"] - df["gene_start"]

    # transcription site length
    df["trans_site_len"] = df["TSS_end"] - df["TSS_start"]

    # ratio transcription site length & land gene length
    df["trans_gene_ratio"] = df["trans_site_len"] / df["gene_length"]

    new_cols = ["strand_binary", "gene_length", "trans_site_len", "trans_gene_ratio"]

    return df[cols_to_keep + new_cols]


def train_model():
    X_1_train = pd.read_csv(X_1_train_path, sep="\t")
    y_1_train = pd.read_csv(y_1_train_path, sep="\t")

    X_1_val = pd.read_csv(X_1_val_path, sep="\t")
    y_1_val = pd.read_csv(y_1_val_path, sep="\t")

    X_2_train = pd.read_csv(X_2_train_path, sep="\t")
    y_2_train = pd.read_csv(y_2_train_path, sep="\t")

    X_2_val = pd.read_csv(X_2_val_path, sep="\t")
    y_2_val = pd.read_csv(y_2_val_path, sep="\t")

    X_1_train = extend_df(X_1_train, cols_to_keep)
    X_1_val = extend_df(X_1_val, cols_to_keep)
    X_2_train = extend_df(X_2_train, cols_to_keep)
    X_2_val = extend_df(X_2_val, cols_to_keep)

    # standardize data

    # Standardize X_1 data
    scaler_1 = StandardScaler()
    X_1_train = scaler_1.fit_transform(X_1_train)
    X_1_val = scaler_1.transform(X_1_val)

    # Standardize X_2 data
    scaler_2 = StandardScaler()
    X_2_train = scaler_2.fit_transform(X_2_train)
    X_2_val = scaler_2.transform(X_2_val)

    # stack data

    # Stack the training data and validation data for X and y
    X_train = np.vstack([X_1_train, X_2_train])
    X_val = np.vstack([X_1_val, X_2_val])

    # Combine y data for training and validation
    y_train = pd.concat([y_1_train, y_2_train], ignore_index=True)
    y_val = pd.concat([y_1_val, y_2_val], ignore_index=True)

    # Prepare the LightGBM dataset
    train_data = lgb.Dataset(X_train, label=y_train)
    val_data = lgb.Dataset(X_val, label=y_val, reference=train_data)

    # Train the model
    model = lgb.train(
        lightgbm_params,
        train_data,
        valid_sets=[val_data],
    )

    # Predict on validation set
    y_pred_val = model.predict(X_val)

    # Calculate Spearman correlation for the predictions on validation set
    spearman_corr, _ = spearmanr(y_val, y_pred_val)
    print(f"Spearman correlation on validation set: {spearman_corr}")

    return model


def eval_model(model):
    # load initial test data for initial gene name order (before feature engineering)
    original_test = pd.read_csv(X_3_test_path_orig, sep="\t")
    gene_names_orig_order = original_test["gene_name"]

    # load feature engineered test data (with different order of rows) and extend it
    X_test = pd.read_csv(X_3_test_path, sep="\t")
    gene_names_eng_order = X_test["gene_name"]
    X_test = extend_df(X_test, cols_to_keep)

    # standardize data
    scaler_test = StandardScaler()
    X_test = scaler_test.fit_transform(X_test)

    raw_preds = model.predict(X_test)

    preds_unsorted = pd.concat(
        [gene_names_eng_order, pd.Series(raw_preds)], ignore_index=True, axis=1
    )

    # Rename the columns
    preds_unsorted.columns = ["gene_name", "gex_predicted"]

    # sort the predicted genes according to intially order for test data
    final_preds = (
        preds_unsorted.set_index("gene_name").loc[gene_names_orig_order].reset_index()
    )

    return final_preds


def export(preds, name):
    save_dir = "."
    file_name = "gex_predicted.csv"
    zip_name = f"{name}_Project1.zip"
    save_path = f"{save_dir}/{zip_name}"
    compression_options = dict(method="zip", archive_name=file_name)

    preds[["gene_name", "gex_predicted"]].to_csv(
        save_path, compression=compression_options
    )


def pipeline(argv):
    model = train_model()
    preds = eval_model(model)
    if len(argv) > 2:
        export(preds, argv[2] + "_" + argv[1])
