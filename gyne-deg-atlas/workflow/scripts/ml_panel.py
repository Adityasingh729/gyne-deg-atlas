#!/usr/bin/env python3
"""Cross-DATASET-validated gene panel on meta-DEGs.
Leave-one-dataset-out CV (NOT leave-one-sample-out) is what makes this honest.
8 GB-safe: restrict features to the meta-signature genes before fitting."""
import argparse, glob, os
# import pandas as pd
# from sklearn.linear_model import LogisticRegressionCV  # elastic net
# from sklearn.model_selection import LeaveOneGroupOut
# import shap

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--meta_dir"); ap.add_argument("--expr_dir")
    ap.add_argument("--samplesheet"); ap.add_argument("--cv", default="leave-one-dataset-out")
    ap.add_argument("--out")
    a=ap.parse_args()
    # TODO:
    #  1) build sample x gene matrix from expr (restricted to meta-signature genes)
    #  2) group = dataset_id ; LeaveOneGroupOut() -> AUROC per held-out dataset
    #  3) report parsimonious panel + SHAP importances
    raise SystemExit("ml_panel.py is a stub - implement leave-one-dataset-out elastic net.")

if __name__ == "__main__":
    main()
