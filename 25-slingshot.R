# === Load Libraries ===
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(viridis)
library(ggplot2)
library(tradeSeq)
library(lmerTest)
library(dplyr)
library(broom.mixed)
library(purrr)
library(tidyr)
library(rstatix)
library(multcomp)
#library(glmGamPoi)

#Parallelization 
#library(future)
#plan("multicore", workers = 10)

#options(future.globals.maxSize = 100 * 1024^3)

# === Load Seurat object and metadata ===

libs_UF_harmony_integrated_OL_lineage <-  readRDS("~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_OL_lineage_harmony_integrated.rds")

libs_UF_harmony_integrated_OL_lineage <- JoinLayers(libs_UF_harmony_integrated_OL_lineage, assay = "RNA")

metadata_full_seurat <- read.csv("~/scratch/20250714_libs1-12B_analysis_v4/metadata_20250815.csv")
metadata_OLseurat <- filter(metadata_full_seurat, clusters_broad %in% c("OPC", "COP", "OL"))
rownames(metadata_OLseurat) <- metadata_OLseurat$X.1

libs_UF_harmony_integrated_OL_lineage <- AddMetaData(libs_UF_harmony_integrated_OL_lineage, metadata_OLseurat)

# Run PCA if not already done
libs_UF_harmony_integrated_OL_lineage <- RunPCA(libs_UF_harmony_integrated_OL_lineage, verbose = TRUE, reduction.name = "pca")

#try reclustering

libs_UF_harmony_integrated_OL_lineage <- FindNeighbors(libs_UF_harmony_integrated_OL_lineage, reduction = "harmony", dims = 1:25)
#libs_UF_harmony_integrated_OL_lineage <- FindClusters(libs_UF_harmony_integrated_OL_lineage, resolution = 0.1)
libs_UF_harmony_integrated_OL_lineage <- RunUMAP(libs_UF_harmony_integrated_OL_lineage, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")

p <- DimPlot(libs_UF_harmony_integrated_OL_lineage, reduction = "umap.harmony")


#Idents(libs_UF_harmony_integrated_OL_lineage) <- libs_UF_harmony_integrated_OL_lineage$celltype


#https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html#basic_processing_with_seurat_pipeline
#https://nbisweden.github.io/workshop-scRNAseq/labs/seurat/seurat_07_trajectory.html

#--------------------------------------------------
# Step 2: Convert to SingleCellExperiment and attach embeddings
#--------------------------------------------------

# Convert to SCE
sce <- as.SingleCellExperiment(libs_UF_harmony_integrated_OL_lineage)

# Add PCA and UMAP embeddings from Seurat
reducedDims(sce)$PCA <- Embeddings(libs_UF_harmony_integrated_OL_lineage, "pca")
reducedDims(sce)$UMAP <- Embeddings(libs_UF_harmony_integrated_OL_lineage, "umap.harmony")

# Add cluster labels from the 'celltype' column
colData(sce)$clusters_broad <- libs_UF_harmony_integrated_OL_lineage$clusters_broad
colData(sce)$celltype <- libs_UF_harmony_integrated_OL_lineage$celltype

#--------------------------------------------------
# Step 3: Run Slingshot using PCA (accurate fitting)
#--------------------------------------------------

sce <- slingshot(
  sce,
  clusterLabels = sce$clusters_broad,          # masked labels
  reducedDim    = "PCA",        
  start.clus    = "OPC", 
  allow.breaks  = FALSE,        
  stretch       = 0.1,           
  shrink        = 0.7,
  reweight=TRUE)

slingLineages(sce)  # Lists lineage(s) as cluster sequences
slingCurves(sce)    # Returns the actual curve data

saveRDS(sce,"~/scratch/20250714_libs1-12B_analysis_v4/OL_lineage_sce_slingshot_20250815.rds")

#--------------------------------------------------
# Step 4: Extract pseudotime and cell weights
#--------------------------------------------------

pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)

#  store pseudotime in Seurat object
libs_UF_harmony_integrated_OL_lineage$pseudotime <- pseudotime

saveRDS(libs_UF_harmony_integrated_OL_lineage, "~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_harmony_integrated_OL_lineage.rds")
#--------------------------------------------------
# Step 5: Visualize pseudotime on UMAP
#--------------------------------------------------

colors <- viridis(100)[cut(pseudotime, breaks = 100)]

plot(reducedDims(sce)$UMAP,
     col = colors,
     pch = 16,
     asp = 1,
     main = "Slingshot Pseudotime (UMAP)")
lines(SlingshotDataSet(sce), lwd = 2, col = "black", reducedDim = "UMAP")



saveRDS(sce,"~/scratch/20250714_libs1-12B_analysis_v4/OL_lineage_sce_slingshot_20250815.rds")


#
# #--------------------------------------------------
# # Step 8: Plotting and group analysis
# #--------------------------------------------------
#
 libs_UF_harmony_integrated_OL_lineage <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_harmony_integrated_OL_lineage.rds")


DimPlot(
  libs_UF_harmony_integrated_OL_lineage,
  reduction = "umap.harmony",
  group.by = "Group",  # replace with your actual group column name
  pt.size = 0.5
) + ggtitle("UMAP by Group")


#umap
FeaturePlot(
  libs_UF_harmony_integrated_OL_lineage,
  features = "pseudotime",
  reduction = "umap.harmony",
  pt.size = 0.5
) + scale_color_viridis_c() + ggtitle("UMAP by Pseudotime")

#pca
FeaturePlot(libs_UF_harmony_integrated_OL_lineage, reduction = "pca", features = "pseudotime", label = TRUE, repel = TRUE)


FeaturePlot(
  libs_UF_harmony_integrated_OL_lineage,
  features = "Age",
  reduction = "umap.harmony",
  pt.size = 0.5
) + scale_color_viridis_c() + ggtitle("UMAP by Age")




# Extract UMAP coordinates and metadata
umap_df <- as.data.frame(Embeddings(libs_UF_harmony_integrated_OL_lineage, "umap.harmony"))
umap_df$pseudotime <- libs_UF_harmony_integrated_OL_lineage$pseudotime
umap_df$group <- libs_UF_harmony_integrated_OL_lineage$Group  # change this if needed
umap_df$age <- as.numeric(libs_UF_harmony_integrated_OL_lineage$Age)  # change this if needed

ggplot(umap_df, aes(x = umapharmony_1, y = umapharmony_2, color = pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  facet_wrap(~ group) +
  theme_minimal() +
  labs(title = "UMAP by Pseudotime, Faceted by Group") +   theme(plot.title = element_text(hjust = 0.5))



df <- data.frame(
  pseudotime = libs_UF_harmony_integrated_OL_lineage$pseudotime,
  group = libs_UF_harmony_integrated_OL_lineage$Group  # adjust if different
)

###########################################################################################################################################

ggplot(df, aes(x = group, y = pseudotime, fill = group)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Pseudotime by Group", y = "Pseudotime", x = "")

meta <- libs_UF_harmony_integrated_OL_lineage@meta.data
summary(lmer(pseudotime ~ Group + Age + Sex + pH  + (1|SubjectID), data = meta))  # adjust as needed



results <- meta %>%
  group_by(celltype) %>%
  group_modify(~{
    fit <- lmer(pseudotime ~ Group + Age + Sex + pH + (1|SubjectID),
                data = .x, REML = FALSE)
    broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
  }) %>%
  ungroup()

# View
results


ggplot(df, aes(x = group, y = pseudotime, color = group)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_summary(fun = median, geom = "point", shape = 18, size = 2, color = "black") +
  theme_minimal() +
  labs(title = "Pseudotime by Group with Median Overlay",
       x = "Group",
       y = "Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5)) +  theme(text = element_text(size = 14)) + theme(legend.position = "none")



ggplot(df, aes(x = pseudotime, fill = group, color = group)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Pseudotime Distributions by Group",
       x = "Pseudotime",
       y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))


#########################################################################################################################################

# Construct a data frame with relevant columns
df <- data.frame(
  pseudotime = meta$pseudotime,
  Sample = meta$Brain.ID,         # subject ID
  Condition = meta$Group,   # e.g., "Case", "Control"
  celltype = meta$celltype      
)

# Summarize median pseudotime per subject and cluster
pt_summary <- df %>%
  group_by(celltype, Sample, Condition) %>%
  summarise(median_PT = median(pseudotime), .groups = "drop")


# Run Wilcoxon test per cluster (Case vs Control)
wilcox_results <- pt_summary %>%
  group_by(celltype) %>%
  wilcox_test(median_PT ~ Condition) %>%
  adjust_pvalue(method = "fdr") %>%
 add_significance("p.adj")


# Run Wilcoxon test all together  (Case vs Control)
wilcox_results_all <- pt_summary %>%
  wilcox_test(median_PT ~ Condition) 

pt_summary$celltype <- factor(pt_summary$celltype, levels = c("OPC", "COP", "OL1", "OL2", "OL3"))



ggplot(pt_summary, aes(x = Condition, y = median_PT, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, aes(color = Condition)) +
  facet_wrap(~ celltype) +
  theme_minimal() +
  labs(title = "Median Pseudotime by Condition and Cluster",
       y = "Median Pseudotime", x = "Condition") +
 theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") + theme(text = element_text(size = 14))




#########################################################################################################################################

#compute pseudotime by clus

PseudoByClus_median <- meta %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(median_PT_clus = median(pseudotime, na.rm = TRUE), .groups = "drop")

PseudoByClus_mean <- meta %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(mean_PT_clus = mean(pseudotime, na.rm = TRUE), .groups = "drop")

lme_celltype = lmer(pseudotime ~ celltype +  (1|SubjectID), data=meta)
anova(lme_celltype)


meta$celltype <- factor(meta$celltype, levels = c("OPC", "COP", "OL1", "OL2", "OL3"))

#########################################################################################################################################

summary(glht(lme_celltype, linfct=mcp(celltype = "Tukey")), test = adjusted(type = "bonferroni"))

Idents(libs_UF_harmony_integrated_OL_lineage) <- libs_UF_harmony_integrated_OL_lineage@meta.data$celltype

#plot pseudo by cluster with colors
meds <- meta %>%
  group_by(celltype) %>%
  summarise(med = median(pseudotime), .groups = "drop")

ggplot(meta, aes(x = celltype, y = pseudotime, fill = celltype, color = celltype)) +
  geom_boxplot(outlier.size = 0.7, outlier.alpha = 0.6, outlier.stroke = 0,
               alpha = 0.8, color = "black", outlier.color = NA, whisker.colour = "grey50" ) +
  geom_point(data = meds, aes(y = med), size = 2, color = "black") +
  geom_text(data = meds, aes(y = med, label = sprintf("%.2f", med)),
            vjust = -3, color = "black") +
  theme_minimal() +
  labs(title = "Pseudotime by Cluster",
       x = "Cluster", y = "Pseudotime") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")


libs_UF_harmony_integrated_OL_lineage <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_harmony_integrated_OL_lineage.rds")
