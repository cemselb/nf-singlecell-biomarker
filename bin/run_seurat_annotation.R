#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
})

option_list = list(
  make_option(c("-i", "--input_rds"), type="character", default=NULL,
              help="Path to input RDS file"),
  make_option(c("-o", "--out_prefix"), type="character", default="sample",
              help="Prefix for outputs")
)
opt = parse_args(OptionParser(option_list=option_list))

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(opt$input_rds)
ref <- HumanPrimaryCellAtlasData()

cat("Predicting cell types with SingleR...\n")
##Extract normalised counts
counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
predictions <- SingleR(test = counts, ref = ref, labels = ref$label.main)

seurat_obj <- AddMetaData(seurat_obj, metadata = predictions$labels, col.name = "SingleR_Labels")

cat("Saving cell type counts...\n")
counts_df <- as.data.frame(table(seurat_obj$SingleR_Labels))
colnames(counts_df) <- c("Cell_Type", "Count")
write.csv(counts_df, file = paste0(opt$out_prefix, "_celltypes.csv"), row.names = FALSE)

cat("Plotting annotated UMAP...\n")
p <- DimPlot(seurat_obj, reduction = "umap", group.by = "SingleR_Labels", label = TRUE, repel = TRUE) + 
  ggtitle("SingleR Automated Annotation")
ggsave(paste0(opt$out_prefix, "_annotation_umap.png"), plot = p, width = 8, height = 6, dpi = 300)

cat("Saving annotated object...\n")
saveRDS(seurat_obj, file = paste0(opt$out_prefix, "_annotated.rds"))
cat("Annotation complete.\n")
