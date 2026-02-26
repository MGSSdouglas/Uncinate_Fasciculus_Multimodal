library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

#Goal is to check for any obvious outliers and look at overall data structure 
#run in jupyter hub

libs_UF_all_harmony_integrated <- readRDS('/home/kelperl/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds')

#note: AnyDoublet_Individual_Assignment == subjectID

#pseudobulk by cell type x subject x batch
########################################################################################################################################################################


# Step 1: Make sure Seurat object has:
# cluster ids 
# sample ids
# group and other covariates of interest


# Step 2: Pseudobulk using AverageExpression
avg_exp <- AggregateExpression(
  libs_UF_all_harmony_integrated,assays = "RNA",
  group.by = c("celltype", "AnyDoublet_Individual_Assignment", "Batch"),return.seurat = TRUE # pseudobulk by cluster & subject & batch
)

# Step 3: Extract matrix for a specific cluster (e.g., "OPC")
# The result is a list per assay — e.g., avg_exp$RNA is a matrix
avg_exp <- NormalizeData(avg_exp)%>% FindVariableFeatures() %>% ScaleData()%>% RunPCA(npcs = 10)


DimPlot(avg_exp, reduction = "pca", group.by = "celltype", label = TRUE)



#pseudobulk by subject x batch
#####################################################################################################################################################################
# Step 1: Make sure Seurat object has:
# cluster ids 
# sample ids
# group and other covariates of interest

# If needed, set cluster identities
# Idents(seurat_obj) <- "cluster_id"

# Step 2: Pseudobulk using AverageExpression
avg_exp <- AggregateExpression(
  libs_UF_all_harmony_integrated,assays = "RNA",
  group.by = c("AnyDoublet_Individual_Assignment", "Batch"),return.seurat = TRUE # pseudobulk by subject & batch
)

# Step 3: Extract matrix for a specific cluster (e.g., "OPC")
# The result is a list per assay — e.g., avg_exp$RNA is a matrix
avg_exp <- NormalizeData(avg_exp)%>% FindVariableFeatures() %>% ScaleData()%>% RunPCA(npcs = 10)


# draw points only
p <- DimPlot(avg_exp, reduction = "pca", label = FALSE, pt.size = 2)

# then add labels with ggrepel params (Seurat passes ... to ggrepel)
p <- LabelClusters(
  plot         = p,
  id           = "ident",
  repel        = TRUE,
  max.overlaps = Inf,     # or a larger number like 100
  box.padding  = 0.2,
  point.padding= 0.1
)
p


#pseudobulk by subject only
#####################################################################################################################################################################

# Step 1: Make sure Seurat object has:
# cluster ids 
# sample ids
# group and other covariates of interest


# If needed, set cluster identities
# Idents(seurat_obj) <- "cluster_id"

# Step 2: Pseudobulk using AverageExpression
avg_exp <- AggregateExpression(
  libs_UF_all_harmony_integrated,assays = "RNA",
  group.by = "AnyDoublet_Individual_Assignment", return.seurat = TRUE #  pseudobulk by subject 
)

# Step 3: Extract matrix for a specific cluster (e.g., "OPC")
# The result is a list per assay — e.g., avg_exp$RNA is a matrix
avg_exp <- NormalizeData(avg_exp)%>% FindVariableFeatures() %>% ScaleData()%>% RunPCA(npcs = 10)



# draw points only
p <- DimPlot(avg_exp, reduction = "pca", label = FALSE, pt.size = 2)

# then add labels with ggrepel params (Seurat passes ... to ggrepel)
p <- LabelClusters(
  plot         = p,
  id           = "ident",
  repel        = TRUE,
  max.overlaps = Inf,     # or a larger number like 100
  box.padding  = 0.2,
  point.padding= 0.1
)
p


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Perform statistical tests to see if significant differences between groups 

#get subject info
subjects <- unique(libs_UF_all_harmony_integrated$SubjectID)

#read in covariates
covar_csv <- read.csv(covar_csv)

covariates_df <-covar_csv %>% filter(Genotype.ID %in% subjects)  %>%
  group_by(Group) %>%
  summarise(
    mean_age = mean(Age),
    se_age = std.error(Age),
    mean_pmi = mean(PMI),
    se_pmi = std.error(PMI),
    mean_ph = mean(as.numeric((pH))),
    se_ph = std.error(as.numeric((pH))))


tmp <- covar_csv %>% filter(Genotype.ID %in% subjects)

table(tmp$Group, tmp$Sex)
table(tmp$Group, tmp$Known.recent.antidepressant)


kw_PMI <- kruskal.test(PMI ~ Group, data = tmp)
kw_age <- kruskal.test(Age ~ Group, data = tmp)
kw_pH <- kruskal.test(pH ~ Group, data = tmp)


t_PMI <- t.test(PMI ~ Group, data = tmp)
t_age <- t.test(Age ~ Group, data = tmp)
t_pH <- t.test(pH ~ Group, data = tmp)


