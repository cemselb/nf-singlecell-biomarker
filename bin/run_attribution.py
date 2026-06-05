#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import torch
import torch.nn as nn
from captum.attr import IntegratedGradients
import matplotlib.pyplot as plt
import seaborn as sns

### a lightweight proxy model to simulate
class MockFM(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(input_dim, 64), nn.ReLU(), nn.Linear(64, 1))
    def forward(self, x):
        return self.net(x)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_h5ad", required=True)
    parser.add_argument("--out_prefix", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_h5ad)
    
    X = torch.tensor(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X, dtype=torch.float32)
    
    model = MockFM(X.shape[1])
    ig = IntegratedGradients(model)
    
    attributions = ig.attribute(X, target=0, n_steps=50)
    attr_np = attributions.detach().numpy()
    
    mean_attr = pd.Series(attr_np.mean(axis=0), index=adata.var_names).sort_values(ascending=False)
    mean_attr.to_csv(f"{args.out_prefix}_attributions.csv")

    ### plot
    plt.figure(figsize=(8, 6))
    sns.barplot(x=mean_attr.head(20).values, y=mean_attr.head(20).index, palette='viridis')
    plt.title("Top 20 Genes by Attribution Importance")
    plt.savefig(f"{args.out_prefix}_attr_heatmap.png", bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    main()
