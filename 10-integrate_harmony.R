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


libs_UF_all <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_all_sctransformed.rds")

# Remove mitochondrial genes from variable features
#hvgs <- VariableFeatures(libs_UF_all)
#hvgs_clean <- hvgs[!grepl("^MT-", hvgs)]
#VariableFeatures(libs_UF_all) <- hvgs_clean

# Step 3: Run PCA and Harmony
#libs_UF_all <- RunPCA(libs_UF_all, features = hvgs_clean)

libs_UF_all <- RunPCA(libs_UF_all, assay = "SCT")


libs_UF_all <- IntegrateLayers(
  object = libs_UF_all,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "harmony", 
  verbose = TRUE)


elbow <- ElbowPlot(libs_UF_all, ndims = 50)
ggsave("~/scratch/20250714_libs1-12B_analysis_v4/UF_merged_SCT_Elbow.pdf", elbow, width = 12, height = 8)


saveRDS(libs_UF_all,"~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_all_harmony_integrated.rds")
