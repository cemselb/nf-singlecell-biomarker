#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list = list(
  make_option(c("-i", "--input_rds"), type="character", default=NULL,
              help="Path to input RDS file"),
  make_option(c("-o", "--out_prefix"), type="character", default="sample",
              help="Prefix for outputs")
)
opt = parse_args(OptionParser(option_list=option_list))

cat("Loading Seurat obj...\n")
seurat_obj <- readRDS(opt$input_rds)

inflammatory_genes <- c('IL32', 'CCL5', 'GNLY', 'NKG7', 'CD14', 'LYZ', 'HLA-DRA')

valid_genes <- intersect(inflammatory_genes, rownames(seurat_obj))

cat(sprintf("Scoring Inflammatory Signature using %d genes...\n", length(valid_genes)))
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(valid_genes),
  name = 'Inflammation_Score'
)

score_col <- "Inflammation_Score1"

cat("Generating UMAP visualisation...\n")
p1 <- FeaturePlot(seurat_obj, features = score_col, cols = c("lightgrey", "darkred")) + 
  ggtitle("Inflammatory Signature Score")
ggsave(paste0(opt$out_prefix, "_signature_umap.png"), plot = p1, width = 8, height = 6, dpi = 300)

cat("Generating violin plot...\n")
p2 <- VlnPlot(seurat_obj, features = score_col, pt.size = 0) + 
  ggtitle("Inflammation Score by Cluster") +
  theme(legend.position = "none")
ggsave(paste0(opt$out_prefix, "_signature_violin.png"), plot = p2, width = 8, height = 6, dpi = 300)

cat("Saving scored obj...\n")
saveRDS(seurat_obj, file = paste0(opt$out_prefix, "_scored.rds"))
cat("Scoring complete.\n")
