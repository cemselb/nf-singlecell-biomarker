#!/usr/bin/env python3

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Run standard scRNA-seq QC using Scanpy.")
    parser.add_argument("--matrix_dir", required=True, help="Path to the 10x Genomics filtered_gene_bc_matrices dir.")
    parser.add_argument("--min_genes", type=int, default=200, help="Minimum number of genes expressed per cell.")
    parser.add_argument("--min_cells", type=int, default=3, help="Minimum number of cells a gene must be expressed in.")
    parser.add_argument("--max_mito", type=float, default=5.0, help="Maximum allowed percentage of mitochondrial reads.")
    parser.add_argument("--out_prefix", default="scanpy_qc", help="Prefix for output files.")
    
    args = parser.parse_args()

    print(f"[INFO] Loading 10x matrix from {args.matrix_dir}...")
    try:
        adata = sc.read_10x_mtx(args.matrix_dir, var_names='gene_symbols', cache=True)
    except Exception as e:
        print(f"[ERROR] Failed to load data: {e}")
        sys.exit(1)

    #basic filtering
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    #mito qc
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    #plotting pre-filter metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                 jitter=0.4, multi_panel=True, show=False)
    plt.savefig(f"{args.out_prefix}_pre_filter_violins.png", bbox_inches='tight')

    #strict filtering based on mitochondrial threshold
    print(f"[INFO] Filtering cells with > {args.max_mito}% mitochondrial reads...")
    adata = adata[adata.obs.pct_counts_mt < args.max_mito, :]

    #save
    output_file = f"{args.out_prefix}_filtered.h5ad"
    adata.write(output_file)
    print(f"[INFO] QC complete. Filtered object saved to {output_file}.")

if __name__ == "__main__":
    main()
