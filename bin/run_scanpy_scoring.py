#!/usr/bin/env python3

import argparse
import scanpy as sc
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Calculate module scores for curated signatures.")
    parser.add_argument("--input_h5ad", required=True, help="Path to input H5AD file")
    parser.add_argument("--out_prefix", required=True, default="sample", help="Prefix for outputs")
    args = parser.parse_args()

    print("Loading AnnData object...")
    adata = sc.read_h5ad(args.input_h5ad)

    inflammatory_genes = ['IL32', 'CCL5', 'GNLY', 'NKG7', 'CD14', 'LYZ', 'HLA-DRA']
    
    valid_genes = [gene for gene in inflammatory_genes if gene in adata.var_names]

    print(f"Scoring Inflammatory Signature using {len(valid_genes)} genes...")
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name='Inflammation_Score')

    print("Generating UMAP visualisation...")
    sc.pl.umap(adata, color='Inflammation_Score', cmap='magma', show=False, title="Inflammatory Signature Score")
    plt.savefig(f"{args.out_prefix}_signature_umap.png", bbox_inches='tight', dpi=300)
    plt.close()

    print("Generating violinpPlot...")
    if 'leiden' in adata.obs:
        sc.pl.violin(adata, keys='Inflammation_Score', groupby='leiden', show=False, rotation=45)
        plt.savefig(f"{args.out_prefix}_signature_violin.png", bbox_inches='tight', dpi=300)
    
    print("Saving scored obj...")
    adata.write(f"{args.out_prefix}_scored.h5ad")
    print("Scoring complete.")

if __name__ == "__main__":
    main()
