library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(anndata)
library(scCustomize)
library(reticulate)


libs_UF_all_harmony_integrated <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")

#add metadata
covariates_df <- read.csv("~/scratch/20250714_libs1-12B_analysis_v4/20250807_subjectinfo_snRNAseq.csv")

# Impute missing pH
#pH of S517 was NA so imputed as mean of all pH of other subjects
#pH_mean <- mean(na.omit(covariates_df$pH))
#covariates_df[which(covariates_df$`Brain ID` == "517"),]$pH <-  pH_mean


# === Add Sequencing Batch Information ===

# Define lookup table: mapping library IDs (orig.ident) to sequencing batches
seq_batch_assignments <- c(
  "lib1_UF" = "SeqBatch1", "lib2_UF" = "SeqBatch1", "lib3_UF" = "SeqBatch1",
  "lib4_UF" = "SeqBatch1", "lib5_UF" = "SeqBatch1", "lib6b_UF" = "SeqBatch1",
  "lib7A_UF" = "SeqBatch2", "lib8A_UF" = "SeqBatch2", "lib9A_UF" = "SeqBatch2",
  "lib10B_UF" = "SeqBatch2", "lib11A_UF" = "SeqBatch2", "lib12B_UF" = "SeqBatch2"
)

seq_batch_vector <- seq_batch_assignments[libs_UF_all_harmony_integrated$orig.ident]
names(seq_batch_vector) <- colnames(libs_UF_all_harmony_integrated)


# Add SeqBatch metadata to Seurat object
libs_UF_all_harmony_integrated <- AddMetaData(
  object = libs_UF_all_harmony_integrated,
  metadata = factor(seq_batch_vector),
  col.name = "SeqBatch"
)

#rename orig ident to batch for consistency (batch refers to individual pools, seqbatch refers to the sequencing runs)
libs_UF_all_harmony_integrated$Batch <- libs_UF_all_harmony_integrated$orig.ident


#######################################################################################################################

# === Rename Subject IDs in Metadata ===
covariates_df <- covariates_df %>%
  dplyr::rename(SubjectID = Genotype.ID) 

  
# === Prepare Seurat Metadata for Merge ===
# Create SubjectID column in Seurat metadata
libs_UF_all_harmony_integrated$SubjectID <- libs_UF_all_harmony_integrated$AnyDoublet_Individual_Assignment

# Extract current Seurat metadata and prepare for join
meta <- libs_UF_all_harmony_integrated@meta.data %>%
  mutate(cell_barcode = rownames(.)) %>%          # Add barcodes as column
  left_join(covariates_df, by = "SubjectID")   # Join on SubjectID

# === Add Merged Metadata Back to Seurat Object ===
libs_UF_all_harmony_integrated <- AddMetaData(libs_UF_all_harmony_integrated, metadata = meta)

# === Final Check ===
stopifnot(all(rownames(libs_UF_all_harmony_integrated@meta.data) == colnames(libs_UF_all_harmony_integrated)))  # should be TRUE

saveRDS(libs_UF_all_harmony_integrated, "~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")

########################################################################################################################################

#convert to sce
sce <- as.SingleCellExperiment(libs_UF_all_harmony_integrated)
saveRDS(sce, "~/scratch/20250714_libs1-12B_analysis_v4/sce_20250815.rds")


#done locally
library(zellkonverter)
library(reticulate)
library(anndata)

# Create a new virtual environment
reticulate::virtualenv_create("zell_env")

# Install required Python packages
reticulate::virtualenv_install("zell_env", packages = c("anndata", "h5py", "numpy", "scipy"))

# Use it
reticulate::use_virtualenv("zell_env", required = TRUE)

#load sce
sce <- readRDS("~/Desktop/MGSS/snRNAseq/sce_20250815.rds")

#convert to anndata
ace_ann <- SCE2AnnData(
  sce)

# Write to compressed h5ad file
write_h5ad(ace_ann,'sce_anndata.h5ad',compression='gzip')


########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
######## #load in allen institute "map my cells" #########################################################


## load updated object

MTG_annot <- read.csv("~/scratch/20250714_libs1-12B_analysis_v4/sce_anndata_10xHumanMTGSEA-AD(CCN20230505)_DeepGenerativeMapping_UTC_1754076889297/sce_anndata_10xHumanMTGSEA-AD(CCN20230505)_DeepGenerativeMapping_UTC_1754076889297.csv",comment.char="#")
MTG_annot_meta <- as.data.frame(cbind(MTG_annot$subclass_name, MTG_annot$supertype_name))
rownames(MTG_annot_meta) <- MTG_annot$cell_id
colnames(MTG_annot_meta) <- c("Allen_subclass_name", "Allen_superclass_name")

libs_UF_all_harmony_integrated <- AddMetaData(libs_UF_all_harmony_integrated,MTG_annot_meta)

DimPlot(libs_UF_all_harmony_integrated, group.by = "Allen_subclass_name")
DimPlot(libs_UF_all_harmony_integrated, group.by = "Allen_superclass_name")

whole_brain_annot <- read.csv("~/scratch/20250714_libs1-12B_analysis_v4/sce_anndata_10xWholeHumanBrain(CCN202210140)_HierarchicalMapping_UTC_1754075731883/sce_anndata_10xWholeHumanBrain(CCN202210140)_HierarchicalMapping_UTC_1754075731883.csv",comment.char="#")
whole_brain_annot_meta <- as.data.frame(cbind(whole_brain_annot$cluster_name, whole_brain_annot$supercluster_name))
rownames(whole_brain_annot_meta) <- whole_brain_annot$cell_id
colnames(whole_brain_annot_meta) <- c("WholeBrainAllen_class_name", "WholeBrainAllen_superclass_name")

libs_UF_all_harmony_integrated <- AddMetaData(libs_UF_all_harmony_integrated,whole_brain_annot_meta)


DimPlot(libs_UF_all_harmony_integrated, group.by = "WholeBrainAllen_class_name")
DimPlot(libs_UF_all_harmony_integrated, group.by = "WholeBrainAllen_superclass_name")

saveRDS(libs_UF_all_harmony_integrated, "~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")

