library(MetaNeighbor)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(anndata)
library(SeuratDisk)
library(Matrix)
library(scclusteval)
library(mclust)

#run on jupyterhub
#code from hbhl 2023 workshop


#########################################################################################################################################################
#Allen MTG

allenMTG <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/Reference_MTG_RNAseq_all-nuclei.2022-06-07.rds")
allenMTG <- NormalizeData(object = allenMTG)
allenMTG <- FindVariableFeatures(object = allenMTG)
allenMTG <- ScaleData(object = allenMTG)
vargenes <- VariableFeatures(allenMTG)

allenMTG_sce <- as.SingleCellExperiment(allenMTG)
rm(allenMTG)

#Look at the major cell types and clusters, also distribution of genes and UMIs
allenMTG$subclass_label %>% table



pretrained_model_major_cluster = MetaNeighbor::trainModel(
  var_genes = vargenes,
  dat = allenMTG_sce,
  study_id = rep("allenMTG", dim(allenMTG)[2]),
  cell_type = allenMTG$subclass_label)

# save model
saveRDS(pretrained_model_major_cluster, "~/scratch/20250714_libs1-12B_analysis_v4/pretainedModelAllenMTG.rds")



libs_UF_all_harmony_integrated <- readRDS('/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds')
pretrained_model_major_cluster <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/pretainedModelAllenMTG.rds")

#try subsetting out the mixed cluster first

sce <- as.SingleCellExperiment(libs_UF_all_harmony_integrated)
sce <- sce[, sce$clusters_broad != "Mixed"]

rm(libs_UF_all_harmony_integrated)


#One versus all
aurocs = MetaNeighborUS(
  trained_model = pretrained_model_major_cluster, dat = sce,
  study_id = rep("UF", dim(sce)[2]), 
  cell_type = sce$celltype,
  fast_version = TRUE
)
tryCatch({plotHeatmapPretrained(aurocs, margins = c(11,5))}, error = function(error_condition) {})



#One versus best
best_hits = MetaNeighborUS(
  trained_model = pretrained_model_major_cluster, dat = sce,
  study_id = rep("UF", dim(sce)[2]), 
  cell_type = sce$celltype,
  one_vs_best = TRUE,
  fast_version = TRUE
)

tryCatch({plotHeatmapPretrained(best_hits, margins = c(11,7))}, 
         error = function(error_condition) {})


celltype_vec <- sce$celltype
celltype_vec <- factor(celltype_vec, levels = c("Astro1", "Astro2", "COP", "Endo", "Excit1", "Excit2", "Excit3", "Excit4", "Excit5", "Excit6",
                                                   "Excit7", "Excit8", "Excit9", "Inhib1", "Inhib2", "Inhib3", "Inhib4", "Inhib5", "Inhib6", 
                                                   "Micro1", "Micro2", "OL1", "OL2", "OL3", "OPC", "Pericytes"))
names(celltype_vec) <- rownames(colData(sce))

allenMTG_vec <- sce$Allen_subclass_name
names(allenMTG_vec) <- rownames(colData(sce))  

clusters_broad_vec <- sce$clusters_broad
clusters_broad_vec <- factor(clusters_broad_vec, levels = c("Astro", "COP", "Endo", "Excit", "Inhib", "Micro", "OL", "OPC", "Pericytes"))
names(clusters_broad_vec) <- rownames(colData(sce))

  
plot1 <- PairWiseJaccardSetsHeatmap(celltype_vec, allenMTG_vec)
adjustedRandIndex(celltype_vec, allenMTG_vec) %>% print

plot2 <- PairWiseJaccardSetsHeatmap(clusters_broad_vec, allenMTG_vec)
adjustedRandIndex(clusters_broad_vec, allenMTG_vec) %>% print


#########################################################################################################################################################


#Jakel2018
exp_Jakel <- read.table("~/scratch/20250714_libs1-12B_analysis_v4/GSE118257_MSCtr_snRNA_ExpressionMatrix_R 2.txt")
metadata_jakel <- read.table("~/scratch/20250714_libs1-12B_analysis_v4/GSE118257_MSCtr_snRNA_FinalAnnotationTable (1).txt", header = TRUE)
rownames(metadata_jakel) <- metadata_jakel$Detected

Seurat_Jakel <- CreateSeuratObject(exp_Jakel, project = "Jakel")
Seurat_Jakel <- AddMetaData(Seurat_Jakel, metadata = metadata_jakel)

Seurat_Jakel <- NormalizeData(object = Seurat_Jakel)
Seurat_Jakel <- FindVariableFeatures(object = Seurat_Jakel)
Seurat_Jakel <- ScaleData(object = Seurat_Jakel)
vargenes <- VariableFeatures(Seurat_Jakel)

Jakel_sce <- as.SingleCellExperiment(Seurat_Jakel)
rm(Seurat_Jakel)

Jakel_sce$Celltypes %>% table


pretrained_model_Jakel = MetaNeighbor::trainModel(
  var_genes = vargenes,
  dat = Jakel_sce,
  study_id = rep("Jakel_sce", dim(Jakel_sce)[2]),
  cell_type = Jakel_sce$Celltypes)

# save model
saveRDS(pretrained_model_Jakel, "~/scratch/20250714_libs1-12B_analysis_v4/pretainedModelJakel.rds")



libs_UF_all_harmony_integrated <- readRDS('/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds')
pretrained_model_Jakel <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/pretainedModelJakel.rds")

sce <- as.SingleCellExperiment(libs_UF_all_harmony_integrated)
sce <- sce[, sce$clusters_broad != "Mixed"]

rm(libs_UF_all_harmony_integrated)

celltype_vec <- sce$celltype
celltype_vec <- factor(celltype_vec, levels = c("Astro1", "Astro2", "COP", "Endo", "Excit1", "Excit2", "Excit3", "Excit4", "Excit5", "Excit6",
                                                "Excit7", "Excit8", "Excit9", "Inhib1", "Inhib2", "Inhib3", "Inhib4", "Inhib5", "Inhib6", 
                                                "Micro1", "Micro2", "OL1", "OL2", "OL3", "OPC", "Pericytes"))
names(celltype_vec) <- rownames(colData(sce))




#One versus all
aurocs = MetaNeighborUS(
  trained_model = pretrained_model_Jakel, dat = sce,
  study_id = rep("UF", dim(sce)[2]), 
  cell_type = sce$celltype,
  fast_version = TRUE
)
tryCatch({plotHeatmapPretrained(aurocs, margins = c(15,6))}, error = function(error_condition) {})



#One versus best
best_hits = MetaNeighborUS(
  trained_model = pretrained_model_Jakel, dat = sce,
  study_id = rep("UF", dim(sce)[2]), 
  cell_type = sce$celltype,
  one_vs_best = TRUE,
  fast_version = TRUE
)

tryCatch({plotHeatmapPretrained(best_hits, margins = c(15,6))}, 
         error = function(error_condition) {})

