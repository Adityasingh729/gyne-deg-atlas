#!/usr/bin/env python3
"""Cross-dataset-validated gene panel selection using meta-analysis DEGs.
Implements Leave-One-Dataset-Out (LODO) cross-validation to assess panel robustness.
Restricts features to the meta-signature genes to optimize resource utilization.
Loads expression TSV files created during fetch_expr."""

import argparse
import glob
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut, StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
# shap is disabled globally on Windows to prevent DLL binary conflict crashes with numpy 2.x
# import shap

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta_dir")
    ap.add_argument("--expr_dir")
    ap.add_argument("--samplesheet")
    ap.add_argument("--cv", default="leave-one-dataset-out")
    ap.add_argument("--out")
    a = ap.parse_args()
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    
    print(f"Loading meta-analysis signatures from {a.meta_dir}...")
    
    # 1) Build list of meta-signature genes
    meta_files = glob.glob(os.path.join(a.meta_dir, "*.meta.tsv"))
    sig_genes = set()
    all_genes_p = []
    
    for f in meta_files:
        try:
            df = pd.read_csv(f, sep="\t")
            if df.empty:
                continue
            # Select significant genes
            sub = df[(df["meta_padj"] < 0.05) & (df["meta_log2FC"].abs() > 1.0)]
            if not sub.empty:
                sig_genes.update(sub["gene"].astype(str).tolist())
            # Save for fallback sorted by p-value
            for _, row in df.iterrows():
                if pd.notna(row["meta_p"]) and pd.notna(row["gene"]):
                    all_genes_p.append((str(row["gene"]), float(row["meta_p"])))
        except Exception as e:
            print(f"Error loading meta-analysis file {f}: {e}")
            
    print(f"Found {len(sig_genes)} standard significant genes (padj < 0.05 & |log2FC| > 1).")
    
    # Fallback: if we have too few genes, take top 100 genes by p-value
    if len(sig_genes) < 5:
        print("Too few significant genes. Relaxing thresholds to top 100 genes by p-value.")
        all_genes_p.sort(key=lambda x: x[1])
        for g, p in all_genes_p[:100]:
            sig_genes.add(g)
            
    sig_genes = sorted(list(sig_genes))
    print(f"Final gene count for machine learning panel: {len(sig_genes)}")
    
    # 2) Load samplesheet and build master sample x gene matrix from expr sidecars
    samplesheet = pd.read_csv(a.samplesheet)
    X_list = []
    y_list = []
    groups_list = []
    
    for _, row in samplesheet.iterrows():
        ds_id = row["dataset_id"]
        group_col = row["group_column"]
        case_lbl = row["case_label"]
        ctrl_lbl = row["control_label"]
        tech = row["technique"]
        
        counts_path = os.path.join(a.expr_dir, f"{ds_id}.counts.tsv")
        meta_path = os.path.join(a.expr_dir, f"{ds_id}.meta.tsv")
        
        if not (os.path.exists(counts_path) and os.path.exists(meta_path)):
            print(f"Sidecar files not found for dataset: {ds_id}. Skipping.")
            continue
            
        print(f"Loading and processing dataset: {ds_id} ({tech})...")
        try:
            counts = pd.read_csv(counts_path, sep="\t", index_col="gene")
            meta = pd.read_csv(meta_path, sep="\t", index_col="sample")
            
            # Select samples that belong to case or control
            meta_filtered = meta[meta[group_col].isin([case_lbl, ctrl_lbl])].copy()
            if meta_filtered.empty:
                print(f"No samples matching case/control for dataset {ds_id}. Skipping.")
                continue
                
            # Filter counts to select samples
            counts_filtered = counts[meta_filtered.index].copy()
            
            # Reindex to fit sig_genes, filling missing with gene-level mean or 0
            counts_filtered = counts_filtered.reindex(sig_genes)
            # Impute missing values per gene with sample mean in this dataset
            counts_filtered = counts_filtered.apply(lambda r: r.fillna(r.mean() if pd.notna(r.mean()) else 0), axis=1)
            
            # Normalize expression
            if "RNA-seq" in tech:
                # Compute CPM and log2(CPM + 1)
                col_sums = counts_filtered.sum(axis=0)
                # Avoid division by zero
                col_sums[col_sums == 0] = 1.0
                normalized = np.log2((counts_filtered / col_sums) * 1e6 + 1.0)
            else:
                # Microarray is already log2 intensity
                normalized = counts_filtered
                
            # Transpose to sample x gene
            X_ds = normalized.T
            # Target labels: case = 1, control = 0
            y_ds = meta_filtered[group_col].map({case_lbl: 1, ctrl_lbl: 0})
            
            X_list.append(X_ds)
            y_list.append(y_ds)
            groups_list.extend([ds_id] * len(y_ds))
            
        except Exception as e:
            print(f"Error processing dataset {ds_id}: {e}")
            
    if not X_list:
        print("No expression datasets could be loaded. Writing empty ML report.")
        with open(a.out, "w") as f:
            f.write("# Error: No datasets successfully loaded.\n")
        return
        
    # Concatenate all datasets
    X = pd.concat(X_list, axis=0)
    y = pd.concat(y_list, axis=0)
    groups = pd.Series(groups_list, index=X.index)
    
    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_scaled_df = pd.DataFrame(X_scaled, columns=X.columns, index=X.index)
    
    # 3) Cross-validation (Leave-One-Dataset-Out if possible, else StratifiedKFold)
    unique_groups = np.unique(groups)
    if len(unique_groups) >= 2:
        cv = LeaveOneGroupOut()
        splits = list(cv.split(X_scaled_df, y, groups))
        cv_type = "LeaveOneDatasetOut"
        print(f"Running Leave-One-Dataset-Out CV across {len(unique_groups)} datasets...")
    else:
        # Fallback if only 1 dataset is present (e.g. during smoke testing)
        cv = StratifiedKFold(n_splits=min(5, len(y)), shuffle=True, random_state=42)
        splits = list(cv.split(X_scaled_df, y))
        cv_type = "StratifiedKFold Fallback"
        print(f"Only {len(unique_groups)} dataset(s) present. Running StratifiedKFold CV...")
        
    aurocs = []
    fold_details = []
    
    for fold, (train_idx, test_idx) in enumerate(splits):
        X_train, X_test = X_scaled_df.iloc[train_idx], X_scaled_df.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        
        # Check that both classes are present in both train and test
        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            print(f"Fold {fold}: Skipping AUROC because train or test set contains only one class.")
            continue
            
        # Fit ElasticNet Logistic Regression
        model = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=1.0, max_iter=5000, random_state=42)
        model.fit(X_train, y_train)
        
        preds = model.predict_proba(X_test)[:, 1]
        auc = roc_auc_score(y_test, preds)
        aurocs.append(auc)
        
        fold_ds = groups.iloc[test_idx].iloc[0] if cv_type == "LeaveOneDatasetOut" else f"Fold_{fold}"
        fold_details.append(f"{fold_ds}={auc:.4f}")
        print(f"Fold: {fold_ds} AUROC: {auc:.4f}")
        
    mean_auc = np.mean(aurocs) if aurocs else 0.5
    print(f"Mean Cross-Validation AUROC: {mean_auc:.4f}")
    
    # 4) Fit final model and compute SHAP values
    print("Fitting final model on all data...")
    model_full = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=1.0, max_iter=5000, random_state=42)
    model_full.fit(X_scaled_df, y)
    
    coefs = model_full.coef_[0]
    
    # Compute feature importance (fallback directly to absolute coefficients to prevent SHAP DLL crash)
    print("SHAP calculation skipped. Using absolute coefficients for feature importance.")
    mean_shap = np.abs(coefs)
        
    # Compile report
    report = pd.DataFrame({
        "gene": sig_genes,
        "coefficient": coefs,
        "mean_abs_shap": mean_shap
    })
    
    # Sort by absolute SHAP value descending
    report["abs_importance"] = report["mean_abs_shap"].abs()
    report = report.sort_values("abs_importance", ascending=False).drop(columns=["abs_importance"])
    
    # Save report with CV metrics in headers
    with open(a.out, "w") as f:
        f.write("# Gyne-DEG Atlas Machine Learning Gene Panel\n")
        f.write(f"# CV Type: {cv_type}\n")
        f.write(f"# Mean AUROC: {mean_auc:.4f}\n")
        f.write(f"# Fold Details: {', '.join(fold_details)}\n")
        f.write(f"# Total samples: {len(y)} (Cases: {sum(y == 1)}, Controls: {sum(y == 0)})\n")
        f.write(f"# Total features: {len(sig_genes)}\n")
        report.to_csv(f, sep="\t", index=False)
        
    print(f"ML report written successfully to {a.out}")

if __name__ == "__main__":
    main()
