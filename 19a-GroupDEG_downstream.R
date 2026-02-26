#############################################
#Downstream processing of group DEG
#############################################
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(disgenet2r)
library(clusterProfiler)
library(enrichplot)
library("org.Hs.eg.db")
library(DOSE)
library(AnnotationDbi)

#Inputs
res       <- read.csv("~/GroupDE_celltype_subjectBatch_ge5cells_20260214_ZEROCOUNTS.csv")
res_broad <- read.csv("~/GroupDE_BROAD_subjectBatch_ge10cells_20260214_ZEROCOUNTS.csv")

##cFilter individual-level DEGs from 'res'
df_summary_individual <- res %>%
  arrange(adj.P.Val) %>%
  filter(adj.P.Val < 0.1)

df_summary_individual$n_ctrl_rem_0nan <- df_summary_individual$n_CTRL_subjects - df_summary_individual$n_zero_CTRL_subjects
df_summary_individual$n_dsca_rem_0nan <- df_summary_individual$n_DSCA_subjects  - df_summary_individual$n_zero_DSCA_subjects

df_summary_individual <- dplyr::filter(df_summary_individual, n_ctrl_rem_0nan >= 12)
df_summary_individual <- dplyr::filter(df_summary_individual, n_dsca_rem_0nan >= 15)

## broad-level DEGs from 'res_broad'
df_summary_broad <- res_broad %>%
  arrange(adj.P.Val) %>%
  filter(adj.P.Val < 0.1)

## NOTE: Adjust column names below to match res_broad if they differ from res
df_summary_broad$n_ctrl_rem_0nan <- df_summary_broad$n_CTRL_subjects - df_summary_broad$n_zero_CTRL_subjects
df_summary_broad$n_dsca_rem_0nan <- df_summary_broad$n_DSCA_subjects  - df_summary_broad$n_zero_DSCA_subjects

df_summary_broad <- dplyr::filter(df_summary_broad, n_ctrl_rem_0nan >= 12)
df_summary_broad <- dplyr::filter(df_summary_broad, n_dsca_rem_0nan >= 15)

##= Defined DEG_all — union of broad and individual DEGs.
DEG_all <- bind_rows(df_summary_individual, df_summary_broad) %>%
  distinct()

## ---- DisGeNET API key ----
#api_key <- input API key derived from disgenet
Sys.setenv(DISGENET_API_KEY = api_key)

#############################################
## DisGeNET: Broad
#############################################
results <- gene2disease(
  gene     = df_summary_broad$gene,    # FIX: was undefined; now uses df_summary_broad
  database = "PSYGENET",
  score    = c(0, 1),
  verbose  = TRUE
)

qres <- results@qresult
  tab <- qres[, c("gene_symbol","disease_name","score","diseaseUMLSCUI")] %>%
    unique() %>% arrange(desc(score))
  
  plot(results, type = "Network", prop = 10, verbose = TRUE)

  ggplot(tab, aes(x = gene_symbol, y = disease_name, fill = score)) +
    geom_tile() + theme_bw() +
    labs(x = "Gene", y = "Psychiatric Disease", fill = "Score", title = "Gene-Disease Heatmap (Broad)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient(low = "skyblue", high = "skyblue4")

## Fisher test (Broad)
de_genes       <- unique(df_summary_broad$gene)
universe_genes <- unique(res_broad$gene)
a <- 3
M <- 1549  #  # total PsyGeNET genes 

U    <- length(universe_genes)
n_DE <- length(de_genes)

if (M > U)    stop("Invalid input: PsyGeNET size M cannot exceed universe size U.")
if (a > M)    stop("Invalid input: overlap a cannot exceed M.")
if (a > n_DE) stop("Invalid input: overlap a cannot exceed number of DE genes.")
if (n_DE > U) stop("Invalid input: number of DE genes cannot exceed universe size U.")

b <- M - a
c <- n_DE - a
d <- U - n_DE - b

tab_fisher <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                     dimnames = list(DE = c("yes","no"), PsyGeNET = c("yes","no")))
tab_fisher

ft <- fisher.test(tab_fisher, alternative = "greater")
ft
ft$estimate
ft$p.value

#############################################
## DisGeNET: Individual
#############################################
results <- gene2disease(
  gene     = df_summary_individual$gene,  # FIX: was df_summary$gene; now df_summary_individual
  database = "PSYGENET",
  score    = c(0, 1),
  verbose  = TRUE
)

qres <- results@qresult
  tab <- qres[, c("gene_symbol","disease_name","score","diseaseUMLSCUI")] %>%
    unique() %>% arrange(desc(score))
  
  plot(results, type = "Network", prop = 10, verbose = TRUE)
  plot(results, type = "Heatmap", limit = 100, nchars = 60, interactive = TRUE, verbose = TRUE)

## Fisher test (Individual)
de_genes       <- unique(df_summary_individual$gene)  # FIX: was df_summary$gene
universe_genes <- unique(res$gene)

a <- 4
M <- 1549  # total PsyGeNET genes 


U    <- length(universe_genes)
n_DE <- length(de_genes)

if (M > U)    stop("Invalid input: PsyGeNET size M cannot exceed universe size U.")
if (a > M)    stop("Invalid input: overlap a cannot exceed M.")
if (a > n_DE) stop("Invalid input: overlap a cannot exceed number of DE genes.")
if (n_DE > U) stop("Invalid input: number of DE genes cannot exceed universe size U.")

b <- M - a
c <- n_DE - a
d <- U - n_DE - b

tab_fisher <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                     dimnames = list(DE = c("yes","no"), PsyGeNET = c("yes","no")))
tab_fisher

ft <- fisher.test(tab_fisher, alternative = "greater")
ft
ft$estimate
ft$p.value

#############################################
## DisGeNET: All genes in DEG_all
#############################################
results <- gene2disease(
  gene     = DEG_all$gene,  # FIX: DEG_all is now defined above
  database = "PSYGENET",
  score    = c(0, 1),
  verbose  = TRUE
)

qres <- results@qresult
  tab <- qres[, c("gene_symbol","disease_name","score","diseaseUMLSCUI")] %>%
    unique() %>% arrange(desc(score))
  
  plot(results, type = "Network", prop = 10, verbose = TRUE)
  plot(results, type = "Heatmap", limit = 100, nchars = 60, interactive = TRUE, verbose = TRUE)
  
  ggplot(tab, aes(x = gene_symbol, y = stringr::str_wrap(disease_name, width = 40), fill = score)) +
    geom_tile() +
    theme_bw() +
    labs(x = "Gene", y = "Psychiatric Disease", fill = "Score", title = "Gene-Disease Heatmap") +
    scale_fill_gradient(low = "skyblue", high = "skyblue4") +
    theme(
      plot.title  = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title  = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10)
    )

## Fisher test (DEG_all)
de_genes       <- unique(DEG_all$gene)
universe_genes <- union(unique(res$gene), unique(res_broad$gene))
a <- 7   
M <- 1549  # total PsyGeNET genes 

U    <- length(universe_genes)
n_DE <- length(de_genes)

if (M > U)    stop("Invalid input: PsyGeNET size M cannot exceed universe size U.")
if (a > M)    stop("Invalid input: overlap a cannot exceed M.")
if (a > n_DE) stop("Invalid input: overlap a cannot exceed number of DE genes.")
if (n_DE > U) stop("Invalid input: number of DE genes cannot exceed universe size U.")

b <- M - a
c <- n_DE - a
d <- U - n_DE - b

tab_fisher <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                     dimnames = list(DE = c("yes","no"), PsyGeNET = c("yes","no")))
tab_fisher

ft <- fisher.test(tab_fisher, alternative = "greater")
ft
ft$estimate
ft$p.value

#############################################
## Plots: Broad vs Individual
#############################################

theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", linewidth = 0.6),
      axis.text.x        = element_text(size = base_size,     color = "black"),
      axis.text.y        = element_text(size = base_size,     color = "black"),
      axis.title.x       = element_text(size = base_size + 1, color = "black", face = "bold"),
      axis.title.y       = element_text(size = base_size + 1, color = "black", face = "bold"),
      plot.title         = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      legend.title       = element_text(size = base_size - 1, face = "bold"),
      legend.text        = element_text(size = base_size - 1),
      legend.key.size    = unit(0.4, "cm")
    )
}

## ---- Genes to label ----
genes_to_label <- c("NECTIN3", "MIR646HG")

## ---- Prepare data ----
df_broad_for_plot <- df_summary_broad %>%
  mutate(neglog10p = -log10(adj.P.Val),
         label     = ifelse(gene %in% genes_to_label, gene, NA))

df_ind_for_plot <- df_summary_individual %>%
  mutate(neglog10p = -log10(adj.P.Val),
         label     = ifelse(gene %in% genes_to_label, gene, NA))

## ---- Broad clusters ----
pbroad <- ggplot(df_broad_for_plot, aes(x = logFC, y = cluster_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(aes(fill = neglog10p),
             shape = 21, color = "grey20", size = 2.25, stroke = 0.4, alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3, fontface = "italic", color = "black",
    box.padding = 1, point.padding = 0.8,
    nudge_x = 0.4, nudge_y = 0.1,
    segment.color = "grey60", segment.size = 0.3,
    na.rm = TRUE, max.overlaps = 20
  ) +
  scale_fill_gradient(low = "skyblue2", high = "skyblue4",
                      name = expression(-log[10]~italic(p)[adj])) +
  xlim(-2.5, 2.5) +
  theme_publication() +
  labs(x = expression(log[2]~"Fold Change"), y = "Cell Type (Broad)")

## ---- Individual clusters ----
pind <- ggplot(df_ind_for_plot, aes(x = logFC, y = cluster_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(aes(fill = neglog10p),
             shape = 21, color = "grey20", size = 2.25, stroke = 0.4, alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3, fontface = "italic", color = "black",
    box.padding = 1, point.padding = 0.8,
    nudge_x = 0.4, nudge_y = 0.1,
    segment.color = "grey60", segment.size = 0.3,
    na.rm = TRUE, max.overlaps = 20
  ) +
  scale_fill_gradient(low = "skyblue2", high = "skyblue4",
                      name = expression(-log[10]~italic(p)[adj])) +
  xlim(-2.5, 2.5) +
  theme_publication() +
  labs(x = expression(log[2]~"Fold Change"), y = "Cell Type (Individual)")

## ---- Combine ----
(pind | pbroad) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Differentially Expressed Genes by Cell Type",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    )
  )
#############################################
## Myelin genes
#############################################
res_myelin <- res[res$gene %in% c("MBP","PLP1","CNP","PLLP","MOG","MAG","MOBP"), ] %>%
  dplyr::filter(cluster_id %in% c("OL1","OL2","OL3"))

ggplot(res_myelin, aes(x = gene, y = cluster_id, fill = logFC)) +
  geom_tile() + theme_bw() +
  labs(x = "Gene", fill = "LogFC", title = "Myelin-constituent genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")

#############################################
## GSEA GO on DEG_all
#############################################
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

ensembl <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = DEG_all$gene, keytype = "SYMBOL", column = "ENSEMBL")
DEG_all$ensembl <- ensembl

genelistFC <- DEG_all$logFC
names(genelistFC) <- DEG_all$ensembl

genelistFC <- genelistFC[!is.na(names(genelistFC))]
genelistFC <- genelistFC[unique(names(genelistFC))]
genelistFC <- sort(genelistFC, decreasing = TRUE)

## NOTE: 'nPerm' is deprecated in newer versions of clusterProfiler; use 'nPermSimple' instead.
gse <- gseGO(geneList      = genelistFC,
             ont           = "BP",
             keyType       = "ENSEMBL",
             nPermSimple   = 10000,   # FIX: nPerm deprecated; replaced with nPermSimple
             minGSSize     = 5,
             maxGSSize     = 800,
             pvalueCutoff  = 0.05,
             verbose       = TRUE,
             OrgDb         = organism,
             pAdjustMethod = "none")

ridgeplot(gse, showCategory = 15) +
  labs(x = "enrichment distribution",
       title = "Biological Processes Enriched in Group DEGs") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))

if (".sign" %in% colnames(as.data.frame(gse@result))) {
  dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)
} else {
  dotplot(gse, showCategory = 10)
}