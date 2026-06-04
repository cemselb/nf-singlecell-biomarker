#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

#args

option_list = list(
  make_option(c("-m", "--matrix_dir"), type="character", default=NULL, 
              help="Path to the 10x Genomics filtered_gene_bc_matrices dir", metavar="character"),
  make_option(c("-g", "--min_genes"), type="integer", default=200, 
              help="Minimum number of genes expressed per cell [default= %default]", metavar="integer"),
  make_option(c("-c", "--min_cells"), type="integer", default=3, 
              help="Minimum number of cells a gene must be expressed in [default= %default]", metavar="integer"),
  make_option(c("-t", "--max_mito"), type="numeric", default=5.0, 
              help="Maximum allowed percentage of mitochondrial reads [default= %default]", metavar="numeric"),
  make_option(c("-o", "--out_prefix"), type="character", default="seurat_qc", 
              help="Prefix for output files [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$matrix_dir)){
  print_help(opt_parser)
  stop("Matrix directory argument must be supplied.", call.=FALSE)
}

message("[INFO] Loading 10x matrix...")
#seurat expects a dir containing matrix.mtx, features.tsv, barcodes.tsv
counts <- Read10X(data.dir = opt$matrix_dir)

#Seurat obj
message("[INFO] Creating Seurat object...")
seurat_obj <- CreateSeuratObject(counts = counts, 
                                 min.cells = opt$min_cells, 
                                 min.features = opt$min_genes, 
                                 project = "scRNA_Pipeline")

#mito qc
message("[INFO] Calculating mitochondrial percentages...")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

#plotting pre-filter metrics
message("[INFO] Generating QC plots...")
vln_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = paste0(opt$out_prefix, "_pre_filter_violins.png"), plot = vln_plot, width = 12, height = 5)

#strict filtering
message(sprintf("[INFO] Filtering cells with > %.1f%% mitochondrial reads...", opt$max_mito))
seurat_obj <- subset(seurat_obj, subset = percent.mt < opt$max_mito)

#save
output_file <- paste0(opt$out_prefix, "_filtered.rds")
saveRDS(seurat_obj, file = output_file)
message("[INFO] QC complete. Filtered object saved to ", output_file)
