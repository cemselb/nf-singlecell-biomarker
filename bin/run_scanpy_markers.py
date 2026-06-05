#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Find cluster biomarkers using Scanpy.")
    parser.add_argument("--input_h5ad", required=True, help="Path to input anndata H5AD file")
    parser.add_argument("--out_prefix", required=True, default="sample", help="Prefix for output files")
    args = parser.parse_args()

    print("Loading AnnData object...")
    adata = sc.read_h5ad(args.input_h5ad)

    if 'leiden' not in adata.obs:
        print("Running neighborhood graph and Leiden clustering...")
        # Assuming PCA was run in the dimred step
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        sc.tl.leiden(adata, resolution=0.5)

    print("Calculating DGE...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    print("Extracting marker table...")
    # Extract the results into a flat Pandas DataFrame
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    df_list = []
    for group in groups:
        df_temp = pd.DataFrame({
            'cluster': group,
            'gene': result['names'][group],
            'scores': result['scores'][group],
            'logfoldchanges': result['logfoldchanges'][group],
            'pvals': result['pvals'][group],
            'pvals_adj': result['pvals_adj'][group]
        })
        df_list.append(df_temp)
        
    marker_df = pd.concat(df_list, ignore_index=True)
    
    sig_markers = marker_df[marker_df['pvals_adj'] < 0.05].sort_values(by=['cluster', 'logfoldchanges'], ascending=[True, False])
    
    print("Saving marker table...")
    sig_markers.to_csv(f"{args.out_prefix}_top_markers.csv", index=False)

    print("Generating marker visualisation...")
    sc.pl.rank_genes_groups_matrixplot(
        adata, 
        n_genes=3, 
        values_to_plot='logfoldchanges', 
        cmap='bwr', 
        vmin=-2, 
        vmax=2, 
        min_logfoldchange=0.5,
        show=False
    )
    
    plt.savefig(f"{args.out_prefix}_marker_heatmap.png", bbox_inches='tight', dpi=300)
    print("Scanpy marker discovery complete.")

if __name__ == "__main__":
    main()
