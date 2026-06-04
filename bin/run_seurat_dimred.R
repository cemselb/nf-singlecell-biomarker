#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list = list(
  make_option(c("-i", "--input_rds"), type="character", default=NULL, 
              help="Path to the filtered Seurat .rds file"),
  make_option(c("-o", "--out_prefix"), type="character", default="sample", 
              help="Prefix for output files")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input_rds)){
  print_help(opt_parser)
  stop("Input RDS file must be supplied.", call.=FALSE)
}

seurat_obj <- readRDS(opt$input_rds)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

umap_plot <- DimPlot(seurat_obj, reduction = "umap")
ggsave(filename = paste0(opt$out_prefix, "_umap.png"), plot = umap_plot, width = 8, height = 6)

output_file <- paste0(opt$out_prefix, "_dimred.rds")
saveRDS(seurat_obj, file = output_file)
