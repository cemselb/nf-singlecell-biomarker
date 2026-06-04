#!/usr/bin/env python3

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import sys

def main():
    parser = argparse.ArgumentParser(description="Run PCA and UMAP using Scanpy.")
    parser.add_argument("--input_h5ad", required=True, help="Path to the filtered .h5ad file.")
    parser.add_argument("--out_prefix", default="sample", help="Prefix for output files.")
    args = parser.parse_args()

    try:
        adata = sc.read_h5ad(args.input_h5ad)
    except Exception as e:
        print(f"[ERROR] Failed to load data: {e}")
        sys.exit(1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata) 

    sc.pl.umap(adata, color=['leiden'], show=False)
    plt.savefig(f"{args.out_prefix}_umap.png", bbox_inches='tight')

    output_file = f"{args.out_prefix}_dimred.h5ad"
    adata.write(output_file)

if __name__ == "__main__":
    main()
