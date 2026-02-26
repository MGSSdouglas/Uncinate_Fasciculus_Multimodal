library(Seurat)
library(sctransform)
library(dplyr)
library(ggplot2)

libs_UF_all_harmony_integrated <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")

libs_UF_all_harmony_integrated <- PrepSCTFindMarkers(
  libs_UF_all_harmony_integrated,
  assay = "SCT",
  verbose = TRUE
)


# Wilcoxon test with SCTransform
Wilcox_markers <- FindAllMarkers(
  libs_UF_all_harmony_integrated,
  assay = "SCT",
  test.use = "wilcox",        # No MAST dependency
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = TRUE
)


# Get top 10 markers per scType cluster
top_markers_harmony <- Wilcox_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  top_n(n = 10, wt = avg_log2FC)

# Heatmap of top scType markers
p <- DoHeatmap(libs_UF_all_harmony_integrated,
               features = top_markers_harmony$gene,
               group.by = "celltype",
               assay = "SCT")


write.csv(Wilcox_markers, "~/scratch/20250714_libs1-12B_analysis_v4/20250815_Wilcox_markers_res0point5_withpericytes.csv", row.names = FALSE)
saveRDS(libs_UF_all_harmony_integrated, "~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")
ggsave("~/scratch/20250714_libs1-12B_analysis_v4/UF_libs_all_singlets_merge_SCT_harmony_markers_heatmap.pdf", p, width = 40, height = 40)
