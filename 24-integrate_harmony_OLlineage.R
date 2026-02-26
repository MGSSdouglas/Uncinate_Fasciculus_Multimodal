library(Seurat)
library(dplyr)
library(harmony)
library(RhpcBLASctl)
library(ggplot2)

# Limit threading for BLAS and OpenMP to prevent nested threads
omp_set_num_threads(1)
blas_set_num_threads(1)

library(future)
options(future.globals.maxSize = 15 * 20 * 1024^3)


libs_UF_OL_lineage <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_OL_lineage_sctransformed.rds")


# Step 3: Run PCA and Harmony
libs_UF_OL_lineage <- RunPCA(libs_UF_OL_lineage)


libs_UF_OL_lineage <- IntegrateLayers(
  object = libs_UF_OL_lineage,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "harmony", 
  verbose = TRUE)


elbow <- ElbowPlot(libs_UF_OL_lineage, ndims = 100)
ggsave("~/scratch/20250714_libs1-12B_analysis_v4/UF_merged_SCT_OL_lineage_Elbow.pdf", elbow, width = 12, height = 8)


saveRDS(libs_UF_OL_lineage,"~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_OL_lineage_harmony_integrated.rds")

