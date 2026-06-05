#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

option_list = list(
  make_option(c("-i", "--input_rds"), type="character", default=NULL,
              help="Path to input Seurat RDS file"),
  make_option(c("-o", "--out_prefix"), type="character", default="sample",
              help="Prefix for output files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input_rds)) {
  print_help(opt_parser)
  stop("Input RDS file must be supplied.", call.=FALSE)
}

cat("Loading Seurat obj...\n")
seurat_obj <- readRDS(opt$input_rds)

if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  cat("Running clustering...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
}

cat("Calculating DE for all clusters...\n")
markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25,
                          verbose = FALSE)

top_markers <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

cat("Saving marker table...\n")
write.csv(markers, file = paste0(opt$out_prefix, "_top_markers.csv"), row.names = FALSE)

cat("Generating marker visualisation...\n")
top_3_genes <- top_markers %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC) %>% pull(gene) %>% unique()

p <- DotPlot(seurat_obj, features = top_3_genes) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Cluster Biomarkers")

ggsave(paste0(opt$out_prefix, "_marker_heatmap.png"), plot = p, width = 10, height = 6, dpi = 300)

cat("Seurat marker discovery complete.\n")
