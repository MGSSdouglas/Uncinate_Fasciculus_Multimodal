library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)

########################################################################################################################################

lib1_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/converted_cellranger_format/")
lib1_UF_obj <- CreateSeuratObject(counts = lib1_UF, project = "lib1_UF", min.cells = 3, min.features = 200)


demultiplex_metadata_lib1 <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib1.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib1) <- demultiplex_metadata_lib1$Barcode
demultiplex_metadata_lib1 <- demultiplex_metadata_lib1[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib1) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib1_UF_obj <- AddMetaData(lib1_UF_obj, metadata = demultiplex_metadata_lib1)
lib1_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib1_UF_obj <- subset(lib1_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib1_UF_obj <- subset(lib1_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib1_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib1_UF_obj, pattern = "^MT-")
View(lib1_UF_obj@meta.data)
VlnPlot(lib1_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib1_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib1_UF_obj)


saveRDS(lib1_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib1_UF_obj.rds")


########################################################################################################################################
lib2_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib2/converted_cellranger_format/")
lib2_UF_obj <- CreateSeuratObject(counts = lib2_UF, project = "lib2_UF", min.cells = 3, min.features = 200)
demultiplex_metadata_lib2 <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib2.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib2) <- demultiplex_metadata_lib2$Barcode
demultiplex_metadata_lib2 <- demultiplex_metadata_lib2[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib2) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib2_UF_obj <- AddMetaData(lib2_UF_obj, metadata = demultiplex_metadata_lib2)
lib2_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib2_UF_obj <- subset(lib2_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib2_UF_obj <- subset(lib2_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib2_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib2_UF_obj, pattern = "^MT-")
View(lib2_UF_obj@meta.data)
VlnPlot(lib2_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib2_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib2_UF_obj)
saveRDS(lib2_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib2_UF_obj.rds")
########################################################################################################################################

lib3_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib3/converted_cellranger_format/")
lib3_UF_obj <- CreateSeuratObject(counts = lib3_UF, project = "lib3_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib3 <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib3.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib3) <- demultiplex_metadata_lib3$Barcode
demultiplex_metadata_lib3 <- demultiplex_metadata_lib3[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib3) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib3_UF_obj <- AddMetaData(lib3_UF_obj, metadata = demultiplex_metadata_lib3)
lib3_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib3_UF_obj <- subset(lib3_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib3_UF_obj <- subset(lib3_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib3_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib3_UF_obj, pattern = "^MT-")
View(lib3_UF_obj@meta.data)
VlnPlot(lib3_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib3_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib3_UF_obj)

saveRDS(lib3_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib3_UF_obj.rds")
########################################################################################################################################

lib4_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib4/converted_cellranger_format/")
lib4_UF_obj <- CreateSeuratObject(counts = lib4_UF, project = "lib4_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib4 <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib4.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib4) <- demultiplex_metadata_lib4$Barcode
demultiplex_metadata_lib4 <- demultiplex_metadata_lib4[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib4) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib4_UF_obj <- AddMetaData(lib4_UF_obj, metadata = demultiplex_metadata_lib4)
lib4_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib4_UF_obj <- subset(lib4_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib4_UF_obj <- subset(lib4_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib4_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib4_UF_obj, pattern = "^MT-")
View(lib4_UF_obj@meta.data)
VlnPlot(lib4_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib4_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib4_UF_obj)

saveRDS(lib4_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib4_UF_obj.rds")
########################################################################################################################################

lib5_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib5/converted_cellranger_format/")
lib5_UF_obj <- CreateSeuratObject(counts = lib5_UF, project = "lib5_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib5 <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib5.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib5) <- demultiplex_metadata_lib5$Barcode
demultiplex_metadata_lib5 <- demultiplex_metadata_lib5[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib5) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib5_UF_obj <- AddMetaData(lib5_UF_obj, metadata = demultiplex_metadata_lib5)
lib5_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib5_UF_obj <- subset(lib5_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib5_UF_obj <- subset(lib5_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib5_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib5_UF_obj, pattern = "^MT-")
View(lib5_UF_obj@meta.data)
VlnPlot(lib5_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib5_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib5_UF_obj)

saveRDS(lib5_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib5_UF_obj.rds")
########################################################################################################################################

lib6b_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib6b/converted_cellranger_format/")
lib6b_UF_obj <- CreateSeuratObject(counts = lib6b_UF, project = "lib6b_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib6b <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib6b.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib6b) <- demultiplex_metadata_lib6b$Barcode
demultiplex_metadata_lib6b <- demultiplex_metadata_lib6b[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib6b) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib6b_UF_obj <- AddMetaData(lib6b_UF_obj, metadata = demultiplex_metadata_lib6b)
lib6b_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib6b_UF_obj <- subset(lib6b_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib6b_UF_obj <- subset(lib6b_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib6b_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib6b_UF_obj, pattern = "^MT-")
View(lib6b_UF_obj@meta.data)
VlnPlot(lib6b_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib6b_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib6b_UF_obj)

saveRDS(lib6b_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib6b_UF_obj.rds")
########################################################################################################################################

lib7A_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib7A/converted_cellranger_format/")
lib7A_UF_obj <- CreateSeuratObject(counts = lib7A_UF, project = "lib7A_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib7A <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib7A.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib7A) <- demultiplex_metadata_lib7A$Barcode
demultiplex_metadata_lib7A <- demultiplex_metadata_lib7A[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib7A) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib7A_UF_obj <- AddMetaData(lib7A_UF_obj, metadata = demultiplex_metadata_lib7A)
lib7A_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib7A_UF_obj <- subset(lib7A_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib7A_UF_obj <- subset(lib7A_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib7A_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib7A_UF_obj, pattern = "^MT-")
View(lib7A_UF_obj@meta.data)
VlnPlot(lib7A_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib7A_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib7A_UF_obj)

saveRDS(lib7A_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib7A_UF_obj.rds")
########################################################################################################################################

lib8A_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib8A/converted_cellranger_format/")
lib8A_UF_obj <- CreateSeuratObject(counts = lib8A_UF, project = "lib8A_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib8A <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib8A.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib8A) <- demultiplex_metadata_lib8A$Barcode
demultiplex_metadata_lib8A <- demultiplex_metadata_lib8A[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib8A) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib8A_UF_obj <- AddMetaData(lib8A_UF_obj, metadata = demultiplex_metadata_lib8A)
lib8A_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib8A_UF_obj <- subset(lib8A_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib8A_UF_obj <- subset(lib8A_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib8A_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib8A_UF_obj, pattern = "^MT-")
View(lib8A_UF_obj@meta.data)
VlnPlot(lib8A_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib8A_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib8A_UF_obj)

saveRDS(lib8A_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib8A_UF_obj.rds")
########################################################################################################################################

lib9A_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib9A/converted_cellranger_format/")
lib9A_UF_obj <- CreateSeuratObject(counts = lib9A_UF, project = "lib9A_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib9A <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib9A.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib9A) <- demultiplex_metadata_lib9A$Barcode
demultiplex_metadata_lib9A <- demultiplex_metadata_lib9A[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib9A) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib9A_UF_obj <- AddMetaData(lib9A_UF_obj, metadata = demultiplex_metadata_lib9A)
lib9A_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib9A_UF_obj <- subset(lib9A_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib9A_UF_obj <- subset(lib9A_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib9A_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib9A_UF_obj, pattern = "^MT-")
View(lib9A_UF_obj@meta.data)
VlnPlot(lib9A_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib9A_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib9A_UF_obj)

saveRDS(lib9A_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib9A_UF_obj.rds")
########################################################################################################################################

lib10B_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib10B/converted_cellranger_format/")
lib10B_UF_obj <- CreateSeuratObject(counts = lib10B_UF, project = "lib10B_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib10B <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib10B.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib10B) <- demultiplex_metadata_lib10B$Barcode
demultiplex_metadata_lib10B <- demultiplex_metadata_lib10B[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib10B) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib10B_UF_obj <- AddMetaData(lib10B_UF_obj, metadata = demultiplex_metadata_lib10B)
lib10B_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib10B_UF_obj <- subset(lib10B_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib10B_UF_obj <- subset(lib10B_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib10B_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib10B_UF_obj, pattern = "^MT-")
View(lib10B_UF_obj@meta.data)
VlnPlot(lib10B_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib10B_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib10B_UF_obj)

saveRDS(lib10B_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib10B_UF_obj.rds")
########################################################################################################################################

lib11A_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib11A/converted_cellranger_format/")
lib11A_UF_obj <- CreateSeuratObject(counts = lib11A_UF, project = "lib11A_UF", min.cells = 3, min.features = 200)

#to rerun once have genoype
demultiplex_metadata_lib11A <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib11A.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib11A) <- demultiplex_metadata_lib11A$Barcode
demultiplex_metadata_lib11A <- demultiplex_metadata_lib11A[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib11A) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib11A_UF_obj <- AddMetaData(lib11A_UF_obj, metadata = demultiplex_metadata_lib11A)
lib11A_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib11A_UF_obj <- subset(lib11A_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib11A_UF_obj <- subset(lib11A_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib11A_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib11A_UF_obj, pattern = "^MT-")
View(lib11A_UF_obj@meta.data)
VlnPlot(lib11A_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib11A_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib11A_UF_obj)

saveRDS(lib11A_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib11A_UF_obj.rds")
########################################################################################################################################

lib12B_UF <- Read10X("/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib12B/converted_cellranger_format/")
lib12B_UF_obj <- CreateSeuratObject(counts = lib12B_UF, project = "lib12B_UF", min.cells = 3, min.features = 200)

demultiplex_metadata_lib12B <- as.data.frame(read_delim(file = '~/scratch/20250714_libs1-12B_analysis_v4/combined_results_w_combined_assignments_lib12B.tsv', delim = "\t" ))
rownames(demultiplex_metadata_lib12B) <- demultiplex_metadata_lib12B$Barcode
demultiplex_metadata_lib12B <- demultiplex_metadata_lib12B[,c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")]
colnames(demultiplex_metadata_lib12B) <- c("AnyDoublet_DropletType", "AnyDoublet_Individual_Assignment")

#add meta data
lib12B_UF_obj <- AddMetaData(lib12B_UF_obj, metadata = demultiplex_metadata_lib12B)
lib12B_UF_obj@meta.data %>% head

#filter out doublets and unassigned cells
lib12B_UF_obj <- subset(lib12B_UF_obj, subset = AnyDoublet_Individual_Assignment != "doublet")
lib12B_UF_obj <- subset(lib12B_UF_obj, subset = AnyDoublet_Individual_Assignment != "unassigned")

# % MT reads
lib12B_UF_obj[["percent.mt"]] <- PercentageFeatureSet(lib12B_UF_obj, pattern = "^MT-")
View(lib12B_UF_obj@meta.data)
VlnPlot(lib12B_UF_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lib12B_UF_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
View(lib12B_UF_obj)

saveRDS(lib12B_UF_obj, "/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib12B_UF_obj.rds")

########################################################################################################################################
################################################################# MERGE ################################################################
########################################################################################################################################

library(RhpcBLASctl)
# Control threading for memory efficiency
blas_set_num_threads(1)
omp_set_num_threads(1)


#load in libraries
lib1_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib1_UF_obj.rds")
lib2_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib2_UF_obj.rds")
lib3_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib3_UF_obj.rds")
lib4_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib4_UF_obj.rds")
lib5_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib5_UF_obj.rds")
lib6b_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib6b_UF_obj.rds")
lib7A_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib7A_UF_obj.rds")
lib8A_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib8A_UF_obj.rds")
lib9A_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib9A_UF_obj.rds")
lib10B_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib10B_UF_obj.rds")
lib11A_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib11A_UF_obj.rds")
lib12B_UF_obj <- readRDS("/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/lib12B_UF_obj.rds")




libs_UF_all <- merge(x = lib1_UF_obj, y = list(lib2_UF_obj, lib3_UF_obj, lib4_UF_obj, lib5_UF_obj, lib6b_UF_obj, lib7A_UF_obj, lib8A_UF_obj, lib9A_UF_obj, lib10B_UF_obj, lib11A_UF_obj, lib12B_UF_obj))

saveRDS(libs_UF_all, '~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_all.rds')




