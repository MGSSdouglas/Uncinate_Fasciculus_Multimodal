#############################################
## Age DEGs downstream analysis: Individual + Broad  
#############################################

library(dplyr)
library(stringr)
library(ggplot2)       
library(AnnotationDbi) 
library(patchwork)     
library(clusterProfiler)
library(enrichplot)
library("org.Hs.eg.db")
library(DOSE)

##inputs
df_summary_ind   <- read.csv("~/AgeDE_celltype_subjectBatch_ge5cells_minSubj27_QWTRUE_GlobalAge_3pass_20260214_ALL_CELLTYPES.csv")
df_summary_broad <- read.csv("~/AgeDE_BROAD_subjectBatch_ge10cells_minSubj27_QWTRUE_GlobalAge_3pass_20260214_ALL_CLUSTERS.csv") 

## Myelin heatmaps (individual + broad) 
df_summary_ind_myelin <- df_summary_ind[df_summary_ind$gene %in% c("MBP","PLP1","CNP","PLLP","MOG","MAG","MOBP"), ] %>%
  dplyr::filter(cluster_id %in% c("OL1","OL2","OL3"))

ggplot(df_summary_ind_myelin, aes(x = cluster_id, y = gene, fill = logFC)) +
  geom_tile() +
  theme_bw() +
  labs(x = "Cluster", y = "Gene", fill = "LogFC", title = "Myelin-constituent genes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_gradient2(low = "#075AFF", high = "white") + theme(axis.title = element_text(size = 14)) +  theme(text = element_text(size = 14))

df_summary_broad_myelin <- df_summary_broad[df_summary_broad$gene %in% c("MBP","PLP1","CNP","PLLP","MOG","MAG","MOBP"), ] %>%
  dplyr::filter(cluster_id %in% "OL")

ggplot(df_summary_broad_myelin, aes(x = gene, y = cluster_id, fill = logFC)) +
  geom_tile() + theme_bw() +
  labs(x= "Gene", fill = "LogFC", title =  "Myelin-constituent genes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_gradient2(low = "#075AFF", high = "white") +   theme(axis.title = element_text(size = 14))+  theme(text = element_text(size = 14))


#######################################################################################################################
# re-filter ind to padj < 0.1
df_summary_ind <- df_summary_ind %>%
  arrange(adj.P.Val) %>% filter(adj.P.Val < 0.1)

df_summary_broad <- df_summary_broad %>%
  arrange(adj.P.Val) %>% filter(adj.P.Val < 0.1)

#######################################################################################################################
# GSEA: individual clusters
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

ensembl <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = df_summary_ind$gene,
                                 keytype = "SYMBOL", column = "ENSEMBL")
df_summary_ind$ensembl <- ensembl

genelistFC <- df_summary_ind$logFC
names(genelistFC) <- df_summary_ind$ensembl
genelistFC <- genelistFC[!is.na(names(genelistFC))]
genelistFC <- genelistFC[unique(names(genelistFC))]
genelistFC <- sort(genelistFC, decreasing = TRUE)

gse <- gseGO(geneList = genelistFC,
             ont = "ALL", keyType = "ENSEMBL",
             nPerm = 10000, minGSSize = 10, maxGSSize = 500,
             pvalueCutoff = 0.05, verbose = TRUE,
             OrgDb = organism, pAdjustMethod = "none")

if (".sign" %in% colnames(as.data.frame(gse@result))) {
  dotplot(gse, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign)
} else {
  dotplot(gse, showCategory = 20)
}

edox2 <- pairwise_termsim(gse)
p1 <- treeplot(edox2, fontsize = 3, hclust_method = "centroid")

if (".sign" %in% colnames(as.data.frame(gse@result))) {
  dotplot(gse, showCategory = 30, split = ".sign") + facet_grid(. ~ .sign) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text  = element_text(size = 10),
      axis.title  = element_text(size = 10),
      plot.margin = margin(5, 5, 5, 40)
    )
} else {
  dotplot(gse, showCategory = 30) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text  = element_text(size = 10),
      axis.title  = element_text(size = 10),
      plot.margin = margin(5, 5, 5, 40)
    )
}

cnetplot(gse, categorySize = "pvalue", foldChange = genelistFC, showCategory = 10)
ridgeplot(gse, showCategory = 20) + guides(fill = guide_colorbar(title = "p-value")) +  
  labs(x = "enrichment distribution",title = "Biological Processes Enriched in Age DEGs") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face= "bold"))



#######################################################################################################################
# GSEA: broad clusters
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

ensembl <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = df_summary_broad$gene,
                                 keytype = "SYMBOL", column = "ENSEMBL")
df_summary_broad$ensembl <- ensembl

genelistFC <- df_summary_broad$logFC
names(genelistFC) <- df_summary_broad$ensembl
genelistFC <- genelistFC[!is.na(names(genelistFC))]
genelistFC <- genelistFC[unique(names(genelistFC))]
genelistFC <- sort(genelistFC, decreasing = TRUE)

gse <- gseGO(geneList = genelistFC,
             ont = "ALL", keyType = "ENSEMBL",
             nPerm = 10000, minGSSize = 10, maxGSSize = 500,
             pvalueCutoff = 0.05, verbose = TRUE,
             OrgDb = organism, pAdjustMethod = "none")

edox2 <- pairwise_termsim(gse)
p2 <- treeplot(edox2, fontsize = 3, hclust_method = "centroid", nCluster = 4)

if (".sign" %in% colnames(as.data.frame(gse@result))) {
  dotplot(gse, showCategory = 30, split = ".sign") + facet_grid(. ~ .sign) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text  = element_text(size = 10),
      axis.title  = element_text(size = 10),
      plot.margin = margin(5, 5, 5, 40)
    )
} else {
  dotplot(gse, showCategory = 30) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text  = element_text(size = 10),
      axis.title  = element_text(size = 10),
      plot.margin = margin(5, 5, 5, 40)
    )
}

ridgeplot(gse, showCategory = 30) + xlab("Enrichment distribution") +
  guides(fill = guide_colorbar(title = "p-value"))

#######################################################################################################################
# Count DEGs up/down per cluster (broad + individual) and heatmaps

# Broad
deg_counts_broad <- df_summary_broad %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down")) %>%
  group_by(cluster_id, direction) %>%
  summarise(n_DEGs_broad = n(), .groups = "drop")

deg_counts_complete_broad <- deg_counts_broad %>%
  tidyr::complete(cluster_id, direction, fill = list(n_DEGs_broad = 0))

p2 <- ggplot(deg_counts_complete_broad, aes(x = direction, y = cluster_id, fill = n_DEGs_broad)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(label = n_DEGs_broad), color = "black", size = 4) +
  scale_fill_gradient2(low = "white", high = "firebrick",
                       midpoint = median(deg_counts_complete_broad$n_DEGs_broad, na.rm = TRUE)) +
  labs(x = NULL, y = "Cluster", fill = "# DEGs", title = "Age DEGs - broad clusters") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16))

# Individual
deg_counts_ind <- df_summary_ind %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down")) %>%
  group_by(cluster_id, direction) %>%
  summarise(n_DEGs_ind = n(), .groups = "drop")

deg_counts_complete_ind <- deg_counts_ind %>%
  tidyr::complete(cluster_id, direction, fill = list(n_DEGs_ind = 0))

p1 <- ggplot(deg_counts_complete_ind, aes(x = direction, y = cluster_id, fill = n_DEGs_ind)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(label = n_DEGs_ind), color = "black", size = 4) +
  scale_fill_gradient2(low = "white", high = "firebrick",
                       midpoint = median(deg_counts_complete_ind$n_DEGs_ind, na.rm = TRUE)) +
  labs(x = NULL, y = "Cluster", title = "Age DEGs - individual clusters") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  theme(legend.position = "none")

p1 + p2
