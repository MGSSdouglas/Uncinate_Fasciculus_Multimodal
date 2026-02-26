library(scCustomize)
library(Seurat)
library(Matrix)
library(DropletUtils)



# Step 1 : Set workdir
setwd("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib12B/")


# Step 1: Read CellBender file
cb_data <- Read_CellBender_h5_Mat("lib12B_cellbender_filtered.h5")

# Step 2: Create Seurat object (optional, just to inspect)
obj <- CreateSeuratObject(cb_data)

# Step 3: Extract matrix, features, and barcodes
sparse_mat <- obj@assays$RNA$counts
barcodes <- colnames(sparse_mat)
features <- rownames(sparse_mat)

# Step 4: Prepare a features data.frame
features_df <- data.frame(
  gene_ids = features,
  gene_symbols = features,
  feature_type = "Gene Expression"
)

# Step 5: Write to Cell Ranger-compatible format
write10xCounts(
  path = "converted_cellranger_format",
  x = sparse_mat,
  gene.id = features_df$gene_ids,
  gene.symbol = features_df$gene_symbols,
  gene.type = "Gene Expression",
  type = "sparse", #output type
  barcodes = barcodes,
  version = "3"
)



#Step 6 - read 10x folder

data_dir <- 'converted_cellranger_format'
expression_matrix <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = expression_matrix)
