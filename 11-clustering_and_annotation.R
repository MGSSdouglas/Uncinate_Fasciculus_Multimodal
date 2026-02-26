library(Seurat)
library(clustree)
library(dplyr)
library(HGNChelper)
library(openxlsx)
library(braincellann)
library(ggplot2)
library(patchwork)
library(scales)  # for comma labels
library(scCustomize)

#library(RhpcBLASctl)

#run on jupyter hub

libs_UF_all_harmony_integrated <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_all_harmony_integrated.rds")

DefaultAssay(libs_UF_all_harmony_integrated) <- "RNA"
libs_UF_all_harmony_integrated <- JoinLayers(libs_UF_all_harmony_integrated)

DefaultAssay(libs_UF_all_harmony_integrated) <- "SCT"


libs_UF_all_harmony_integrated <- FindNeighbors(libs_UF_all_harmony_integrated, reduction = "harmony", dims = 1:25)
libs_UF_all_harmony_integrated <- FindClusters(libs_UF_all_harmony_integrated, resolution = c(0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
libs_UF_all_harmony_integrated <- RunUMAP(libs_UF_all_harmony_integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

p <- DimPlot(libs_UF_all_harmony_integrated, reduction = "umap.harmony")


######## CLUSTREE ########
p <- clustree(libs_UF_all_harmony_integrated@meta.data, prefix = "SCT_snn_res.")
Idents(libs_UF_all_harmony_integrated) <- libs_UF_all_harmony_integrated@meta.data$SCT_snn_res.0.5

######## SCTYPE ########
source("~/scratch/20250327_libs1-12B_analysis_v1/sctype_wrapper.R"); 

libs_UF_all_harmony_integrated <- run_sctype(libs_UF_all_harmony_integrated, assay = "SCT", scaled = TRUE, known_tissue_type="Brain", name="sctype_celltypes")

p <- DimPlot(libs_UF_all_harmony_integrated, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = "sctype_celltypes")

#####################################################################################################################################################

p1 <- ggplot(libs_UF_all_harmony_integrated@meta.data, aes(x = celltype, fill = AnyDoublet_Individual_Assignment)) + geom_bar(position = "fill") + theme_minimal() + theme(axis.text.x = element_text(angle = 45))
#p2 <- ggplot(libs_UF_all_harmony_integrated@meta.data, aes(x = celltype, fill = Group)) + geom_bar(position = "fill")+ theme_minimal()
p3 <- ggplot(libs_UF_all_harmony_integrated@meta.data, aes(x = celltype, fill = orig.ident)) + geom_bar(position = "fill")+ theme_minimal() + theme(axis.text.x = element_text(angle = 45))

p1 + p3


ggplot(libs_UF_all_harmony_integrated@meta.data, aes(x = AnyDoublet_Individual_Assignment, fill = celltype)) +
  geom_bar(position = "fill") + theme_minimal() 


######## MARKERS ########

main_markers <- c("ALDH1L1", "AQP4", "MOG", "PDGFRA", "CSPG4", "CX3CR1", "CSF1R", "CLDN5", "PDGFRB", "RBFOX3", "SLC17A7", "GAD1")

excitatory_layer_markers <- c("SLC17A7", "SATB2", "FXYD6", "GSG1L", "RASGRF2", "CUX2", "RELN", "TLE1","HTR2C", "FOXP2")

inhibitory_subtype_markers <- c("GAD1", "GAD2", "SST", "PVALB", "VIP", "LAMP5", "LHX6", "ADARB2", "CCK", "NPY", "CALB1", "CALB2")

glial_markers <- c("GFAP", "GJA1", "ALDH1L1", "AQP4", "MAG", "MOG", "OLIG1", "OLIG2", "SOX10", "MYT1", "PDGFRA", "ZFPM2", "ITPR2", "TCF7L2", "MRC1", "CX3CR1", "SPI1", "VIM", "CLDN5")

OL_markers <- c("MBP", "MOG", "MOBP", "PLP1", "MAL", "OPALIN", "CNP", "SOX10", "OLIG2", "PDGFRA", "CSPG4", "PCDH15", "VCAN", "BCAS1", "TNR", "BCAN")


#DotPlot(libs_UF_all_harmony_integrated, features = main_markers) + RotatedAxis()


#ol markers

main <- VlnPlot(libs_UF_all_harmony_integrated, features = main_markers)
ggsave("my_violin_plot.png", plot = main, width = 18, height = 10, units = "in", dpi = 300)

ol <- VlnPlot(libs_UF_all_harmony_integrated, features = OL_markers)
ggsave("my_violin_plot_OL.png", plot = ol, width = 18, height = 10, units = "in", dpi = 300)

glia <- VlnPlot(libs_UF_all_harmony_integrated, features = glial_markers)
ggsave("my_violin_plot_glia.png", plot = glia, width = 18, height = 10, units = "in", dpi = 300)

inhib <- VlnPlot(libs_UF_all_harmony_integrated, features = inhibitory_subtype_markers)
ggsave("my_violin_plot_inhib.png", plot = inhib, width = 18, height = 10, units = "in", dpi = 300)

excit <- VlnPlot(libs_UF_all_harmony_integrated, features = excitatory_layer_markers)
ggsave("my_violin_plot_excit.png", plot = excit, width = 18, height = 10, units = "in", dpi = 300)



#ggsave("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/plot.png")

FeaturePlot(libs_UF_all_harmony_integrated, features = "percent.mt")
#to look for immune OL
FeaturePlot(libs_UF_all_harmony_integrated, c("CD74", "APOE"))

#pericytes
FeaturePlot(libs_UF_all_harmony_integrated, "PDGFRB")


VlnPlot(libs_UF_all_harmony_integrated, features = "percent.mt")
p1 <- VlnPlot(libs_UF_all_harmony_integrated, features = "nCount_RNA")
p2 <- VlnPlot(libs_UF_all_harmony_integrated, features = "nFeature_RNA")
p3 <- VlnPlot(libs_UF_all_harmony_integrated, features = "nCount_SCT")
p4 <- VlnPlot(libs_UF_all_harmony_integrated, features = "nFeature_SCT")

p1+p2
p3+p4

#manually select perictyes
plot <- DimPlot(libs_UF_all_harmony_integrated, reduction = "umap.harmony")
select.cells <- CellSelector(plot = plot)


Idents(libs_UF_all_harmony_integrated) <- libs_UF_all_harmony_integrated@meta.data$SCT_snn_res.0.5


clusters_broad <- c("0"="OPC", "1"="OL", "2" = "OL", "3"="Micro", "4"="Astro", "5" = "OL", 
                    "6" = "Astro", "7" = "Excit", "8" = "Excit", "9"= "Excit", "10" = "Inhib", 
                    "11" = "Inhib", "12" = "Excit", "13"= "Inhib", "14" = "Inhib", "15" = "Excit", 
                    "16" = "Micro", "17" = "Excit", "18" = "Inhib", "19" = "Excit", "20" = "Mixed",
                    "21"="Endo", "22" ="Excit", "23"="Inhib", "24"="COP","25"="Excit", "26"="Mixed")


libs_UF_all_harmony_integrated <- RenameIdents(libs_UF_all_harmony_integrated, clusters_broad)

#add in pericytes
# 5) Add the new level, then assign it
Idents(libs_UF_all_harmony_integrated, cells = select.cells) <- "Pericytes"

libs_UF_all_harmony_integrated$clusters_broad <- Idents(libs_UF_all_harmony_integrated)
DimPlot(libs_UF_all_harmony_integrated)

##########################################################################################################################################

Idents(libs_UF_all_harmony_integrated) <- libs_UF_all_harmony_integrated@meta.data$SCT_snn_res.0.5

name_map <- c("0"="OPC", "1"="OL1", "2" = "OL2", "3"="Micro1", "4"="Astro1", "5" = "OL3", 
              "6" = "Astro2", "7" = "Excit1", "8" = "Excit2", "9"= "Excit3", "10" = "Inhib1", 
              "11" = "Inhib2", "12" = "Excit4", "13"= "Inhib3", "14" = "Inhib4", "15" = "Excit5", 
              "16" = "Micro2", "17" = "Excit6", "18" = "Inhib5", "19" = "Excit7", "20" = "Mixed1",
              "21"="Endo", "22" ="Excit8", "23"="Inhib6", "24"="COP","25"="Excit9", "26"="Mixed2")

libs_UF_all_harmony_integrated <- RenameIdents(libs_UF_all_harmony_integrated, name_map)

#add in pericytes based on marker genes and UMAP coordinates
Idents(libs_UF_all_harmony_integrated, cells = select.cells) <- "Pericytes"

libs_UF_all_harmony_integrated$celltype <- Idents(libs_UF_all_harmony_integrated)
DimPlot(libs_UF_all_harmony_integrated)


saveRDS(libs_UF_all_harmony_integrated,"~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")


###################################################################################################################################################################################

#more plotting 
#trying different dotplot options
libs_UF_all_harmony_integrated <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")


topmarkers <- c("GFAP", "AQP4", "PDGFRA", "GPR17", "MBP", "P2RY12", "CX3CR1", "CLDN5", "PDGFRB", "SNAP25", "SLC17A7", "GAD1")

clusterorder <- c("Astro1", "Astro2", "OPC", "COP", "OL1", "OL2", "OL3", "Micro1", "Micro2", "Endo", "Pericytes", "Excit1", "Excit2", "Excit3", 
                  "Excit4", "Excit5", "Excit6", "Excit7", "Excit8", "Excit9", "Inhib1", "Inhib2", "Inhib3", "Inhib4", "Inhib5", "Inhib6", "Mixed1", "Mixed2")

libs_UF_all_harmony_integrated$celltype <- factor(libs_UF_all_harmony_integrated@meta.data$celltype, levels = clusterorder)
Idents(libs_UF_all_harmony_integrated) <- libs_UF_all_harmony_integrated$celltype


DotPlot(libs_UF_all_harmony_integrated, topmarkers, idents = c("Astro1", "Astro2", "OPC", "COP", "OL1", "OL2", "OL3", "Micro1", "Micro2", "Endo", "Pericytes", 
                                                                   "Excit1", "Excit2", "Excit3", "Excit4", "Excit5", "Excit6", "Excit7", "Excit8", "Excit9", 
                                                                   "Inhib1", "Inhib2", "Inhib3", "Inhib4", "Inhib5", "Inhib6")) + FontSize(x.text = 10, y.text = 11, x.title = 14, y.title = 14)                                                               

DotPlot(libs_UF_all_harmony_integrated, topmarkers) + FontSize(x.text = 10, y.text = 11, x.title = 14, y.title = 14)                                                               


Clustered_DotPlot(seurat_object = libs_UF_all_harmony_integrated, features = c("GFAP", "AQP4", "PDGFRA", "GPR17", "MBP", "CX3CR1", "P2RY12", "CLDN5", "PDGFRB", "SNAP25", "SLC17A7", "GAD1"), flip = TRUE)


topmarkers <- c("GFAP","AQP4","PDGFRA","GPR17","MBP","P2RY12",
                "CX3CR1","CLDN5","PDGFRB","SNAP25","SLC17A7","GAD1")

#with mixed clusters
clusterorder <- c("Astro1","Astro2","OPC","COP","OL1","OL2","OL3",
                  "Micro1","Micro2","Endo","Pericytes",
                  "Excit1","Excit2","Excit3","Excit4","Excit5","Excit6","Excit7","Excit8","Excit9",
                  "Inhib1","Inhib2","Inhib3","Inhib4","Inhib5","Inhib6","Mixed1","Mixed2")

#want reverse of order above, hence rev()
# ensure factor order & identities
libs_UF_all_harmony_integrated$celltype <- factor(
  libs_UF_all_harmony_integrated@meta.data$celltype,
  levels = rev(clusterorder)
)
Idents(libs_UF_all_harmony_integrated) <- libs_UF_all_harmony_integrated$celltype

#without mixed clusters
idents_vec <- c("Astro1","Astro2","OPC","COP","OL1","OL2","OL3","Micro1","Micro2","Endo","Pericytes",
                "Excit1","Excit2","Excit3","Excit4","Excit5","Excit6","Excit7","Excit8","Excit9",
                "Inhib1","Inhib2","Inhib3","Inhib4","Inhib5","Inhib6")

# ---- DotPlot ----
p_dot <- DotPlot(
  object   = libs_UF_all_harmony_integrated,
  features = topmarkers,
  idents   = idents_vec,
  dot.scale = 4.5
) +
  Seurat::FontSize(x.text = 12.5, y.text = 12, x.title = 14, y.title = 14) +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
  labs(x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

y_levels <- levels(p_dot$data$id)

# ---- cell counts (same ordering) ----
counts_df <- libs_UF_all_harmony_integrated@meta.data %>%
  filter(celltype %in% y_levels) %>%
  count(celltype, name = "n") %>%
  mutate(
    celltype = factor(celltype, levels = y_levels),
    frac = n / sum(n)
  )

# ---- bar plot (on the LEFT) ----
p_bar <- ggplot(counts_df, aes(x = n, y = celltype)) +
  geom_col(width = 0.7, fill = "grey") +
  # put the labels inside or just to the left of the bar
  geom_text(aes(label = scales::comma(n)), hjust = 1.25, size = 2.75) +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0))) +
  labs(x = "Cells", y = "") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + scale_y_discrete(position = "right") +  theme(axis.text.y = element_text(color = "black"))

# show y labels only once (on the bars)
p_dot <- p_dot + theme(axis.text.y = element_blank())

# ---- combine ----
p_combined <- p_bar + p_dot + plot_layout(widths = c(4.5, 3.2))
p_combined

