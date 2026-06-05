#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import celltypist
from celltypist import models
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_h5ad", required=True)
    parser.add_argument("--out_prefix", required=True, default="sample")
    args = parser.parse_args()

    print("Loading AnnData object...")
    adata = sc.read_h5ad(args.input_h5ad)

    print("Loading CellTypist model...")
    models.download_models(force_update=False, model=['Immune_All_Low.pkl'])
    model = models.Model.load(model='Immune_All_Low.pkl')

    print("Predicting cell types...")
    if adata.raw is not None:
        adata_pred = adata.raw.to_adata()
    else:
        adata_pred = adata.copy()

    if adata_pred.X.max() > 100:
        sc.pp.normalize_total(adata_pred, target_sum=1e4)
        sc.pp.log1p(adata_pred)

    predictions = celltypist.annotate(adata_pred, model=model, majority_voting=True)
    
    adata.obs['majority_voting'] = predictions.predicted_labels['majority_voting']

    print("Saving outputs...")
    counts = adata.obs['majority_voting'].value_counts().reset_index()
    counts.columns = ['Cell_Type', 'Count']
    counts.to_csv(f"{args.out_prefix}_celltypes.csv", index=False)

    sc.pl.umap(adata, color=['majority_voting'], show=False, title='CellTypist Annotation')
    plt.savefig(f"{args.out_prefix}_annotation_umap.png", bbox_inches='tight', dpi=300)

    adata.write(f"{args.out_prefix}_annotated.h5ad")

if __name__ == "__main__":
    main()
