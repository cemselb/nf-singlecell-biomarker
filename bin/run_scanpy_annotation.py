#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import celltypist
from celltypist import models
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="ML Cell Annotation using CellTypist.")
    parser.add_argument("--input_h5ad", required=True, help="Path to input H5AD file")
    parser.add_argument("--out_prefix", required=True, default="sample", help="Prefix for outputs")
    args = parser.parse_args()

    print("Loading AnnData object...")
    adata = sc.read_h5ad(args.input_h5ad)

    print("Downloading and loading CellTypist model...")
    models.download_models(force_update=False, model=['Immune_All_Low.pkl'])
    model = models.Model.load(model='Immune_All_Low.pkl')

    print("Predicting cell types...")
    ###Celltypist requires raw counts normalized to 10k
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    predictions = celltypist.annotate(adata, model=model, majority_voting=True)
    adata = predictions.to_adata()

    print("Saving cell type counts...")
    counts = adata.obs['majority_voting'].value_counts().reset_index()
    counts.columns = ['Cell_Type', 'Count']
    counts.to_csv(f"{args.out_prefix}_celltypes.csv", index=False)

    print("Plotting annotated UMAP...")
    sc.pl.umap(adata, color=['majority_voting'], show=False, title='CellTypist Annotation')
    plt.savefig(f"{args.out_prefix}_annotation_umap.png", bbox_inches='tight', dpi=300)

    print("Saving annotated object...")
    adata.write(f"{args.out_prefix}_annotated.h5ad")
    print("Annotation complete.")

if __name__ == "__main__":
    main()
