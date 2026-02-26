############################
library(muscat)
library(edgeR)
library(limma)
library(tidyverse)
library(WGCNA)
library(tidyr)
library(DESeq2)
library(cluster)
library(reshape2)
library(RColorBrewer)
library(flashClust)
library(SingleCellExperiment)
library(scCustomize)
library(stringr)
# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
#enableWGCNAThreads(nThreads = 8)

#used https://github.com/MgssGroup/snRNASeq_public/blob/main/10_WGCNA/wgcna_InN.R as template

#run on jupyter hub


seurat_obj <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")
#set RNA active assay
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- subset(seurat_obj, subset = clusters_broad == "OL")

seurat_obj_pb <- AggregateExpression(seurat_obj, assays = "RNA", group.by = "SubjectID")

###################################################################################################
#read in subjectinfo
subinfo <- read.csv("~/scratch/20250714_libs1-12B_analysis_v4/20250807_subjectinfo_snRNAseq.csv")

#subinfo<- merge(subinfo, seqbatches, by.x = "Genotype.ID", by.y = "SubjectID")
subinfo$Genotype.ID <- gsub("_", "-", subinfo$Genotype.ID)
subinfo$Sex_bin <- gsub("Male", 0, subinfo$Sex)
subinfo$Sex_bin <- gsub("Female", 1, subinfo$Sex_bin)

counts <- seurat_obj_pb[["RNA"]]

#filter for low counts 
filter <- rowSums(counts) > 5
counts_select <- counts[filter,]
nrow(counts_select)
#set covariates in meta file
subinfo$Group <-factor(subinfo$Group, levels = c("DS-CA", "CTRL"))
subinfo$Sex <- as.factor(subinfo$Sex)
subinfo$Sex_bin <- as.factor(subinfo$Sex_bin)
subinfo$Age <- as.numeric(subinfo$Age)
subinfo$pH <- as.numeric(subinfo$pH)
subinfo$PMI <- as.numeric(subinfo$PMI)
datTraits <- as.data.frame(subinfo)
rownames(datTraits) <- datTraits$Genotype.ID

#check format for meta and counts
datTraits <- datTraits[order(rownames(datTraits)),]  
datTraits <- datTraits[colnames(counts_select),]  
counts_select <- counts_select[,order(colnames(counts_select))]  
all(rownames(datTraits) == colnames(counts_select))

rm(seurat_obj)
#normalize and transform DESeqDataSetFromMatrix without defining a specific model / variance stabilizing transformation
counts_select <- round(counts_select) %>% as.data.frame()
dds <- DESeqDataSetFromMatrix(countData=counts_select, colData=datTraits, design = ~ 1)
rld <- rlog(dds, blind = TRUE)

#remove batch effect using limma::removeBatcheffect
cov.matrix <- cbind(datTraits$pH, datTraits$PMI)
designK <- model.matrix(~ Group + Age + Sex, data = datTraits) 
normalized <- removeBatchEffect(assay(rld), covariates=cov.matrix, design=designK)
datExpr <- t(normalized)


######################################################################################################################
############################################ Power analysis check ####################################################
powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr, powerVector=powers, verbose=5, blockSize= 5000, networkType = "signed")

save.image("~/scratch/20250714_libs1-12B_analysis_v4/myWGCNA_OL.Rdata")

#load("~/scratch/20250714_libs1-12B_analysis_v4/myWGCNA_OL.Rdata")

par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

######################################################################################################################
############################################ MODULE construction #####################################################

PWR=16
minModule=30
cor <- WGCNA::cor

#check scale free
adjacency = adjacency(datExpr, power = PWR, type="signed")
TOM = TOMsimilarity(adjacency,TOMType = "signed")
dissTOM = 1-TOM
k=as.vector(apply(adjacency, 2, sum, na.rm=T))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

#build modules
net = blockwiseModules(datExpr, corType="pearson", maxBlockSize = 40000, networkType="signed",power=PWR, minModuleSize= minModule, nThreads=30, TOMType = "signed", TOMDenom = "min", deepSplit=2, verbose=5, mergeCutHeight=0.15, reassignThreshold = 1e-6, detectCutHeight = 0.995, numericLabels=TRUE, saveTOMs=FALSE, pamStage=TRUE, pamRespectsDendro=FALSE)
moduleLabelsAutomatic=net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
MEsAutomatic=net$MEs
unique(moduleColorsAutomatic)
table(moduleColorsAutomatic)
write.table(moduleColorsAutomatic, "DG_colors.txt",sep="\t",quote=F)
save.image("~/scratch/20250714_libs1-12B_analysis_v4/myWGCNA_OL.Rdata")

#plot dendogram
svg("~/scratch/20250714_libs1-12B_analysis_v4/OL_Dendrogram.Module.svg")
plotDendroAndColors(net$dendrograms[[1]], moduleColorsAutomatic, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

#list hub gene for each module
chooseTopHubInEachModule(datExpr, moduleColorsAutomatic, omitColors = "grey", power = PWR, type = "signed")

#KMEs
KMEs<-signedKME(datExpr, net$MEs,corFnc = "cor", corOptions = "use = 'p'")
kme=data.frame(rownames(counts_select), moduleColorsAutomatic, KMEs)
colnames(kme)[1]="Symbol"
rownames(kme)=NULL
write.table(kme,"~/scratch/20250714_libs1-12B_analysis_v4/OL_KME_DG.txt",sep="\t",quote=F)

#correlate module and phenotype
moduleColorsIEGG=moduleColorsAutomatic
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsIEGG)$eigengenes
MEsIEGG = MEs0
MEsIEGG$MEgrey=NULL

#set phenotype of interest
datTraits$Group_bin <- gsub("CTRL", 0, datTraits$Group)
datTraits$Group_bin <- gsub("DS-CA", 1, datTraits$Group_bin)
datTraits$Sex_bin <- gsub("Male", 0, datTraits$Sex)
datTraits$Sex_bin <- gsub("Female", 1, datTraits$Sex_bin)
datCovariate <- datTraits %>% dplyr::select(Age,Group_bin, Sex_bin, pH, PMI)

modTraitCor= cor(MEsIEGG,datCovariate,method="pearson")
write.table(modTraitCor,"~/scratch/20250714_libs1-12B_analysis_v4/OL_modTraitCor_DG.txt",sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitP,"~/scratch/20250714_libs1-12B_analysis_v4/OL_modTraitP_DG.txt",sep="\t",quote=F)
textMatrix = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
svg("~/scratch/20250714_libs1-12B_analysis_v4/OL_Heatmap_DatTraits.svg", height=2000, width=2000)
par(mar = c(3,9, 2, 1))
#labeledHeatmap(Matrix = modTraitCor, xLabels = c("Age", "Group", "Sex", "pH", "PMI"), yLabels = names(MEsIEGG), ySymbols = names(MEsIEGG), colorLabels =FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("OL Module Association"))
labeledHeatmap(Matrix = modTraitCor[,1:3], xLabels = c("Age", "Group", "Sex"), yLabels = names(MEsIEGG), ySymbols = names(MEsIEGG), colorLabels =FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix[,1:3], setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("OL Module Association"))

dev.off()

#eigenGeneNet
MEList=moduleEigengenes(datExpr,colors=moduleColorsAutomatic,softPower = PWR,impute = TRUE)
MEs = MEList$eigengenes
#MEs <- dplyr::select(MEs, -"MEgrey")  
MET=orderMEs(MEs)

#GGplotInput
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsAutomatic,softPower = PWR,impute = TRUE)$eigengenes
MEs0$Rows=colnames(counts_select)
write.table(MEs0, "~/scratch/20250714_libs1-12B_analysis_v4/OL_DG_Matrix_module_correlation.txt",sep="\t",quote=F)

#adjacency matrix
Adj = adjacency(datExpr, power = PWR,type="signed",corFnc = "cor", corOptions = "use = 'p'")
moduleOutput <- data.frame(rownames(counts_select))
moduleOutput[,2]<- moduleColorsAutomatic
intraCon <- intramodularConnectivity(Adj, moduleColorsAutomatic)
moduleOutput[,3]<-intraCon$kWithin
colnames(moduleOutput) <- c("Gene", "ModuleColor", "kWithin")
write.table(moduleOutput, "~/scratch/20250714_libs1-12B_analysis_v4/OL_ModuleOutput_DG.txt", sep="\t", quote=F)

#display the corelation values with a heatmap plot
#INCLUE THE NEXT LINE TO SAVE TO FILE
pdf(file="OLheatmap.pdf")

save.image("~/scratch/20250714_libs1-12B_analysis_v4/myWGCNA_OL.Rdata")
load("~/scratch/20250714_libs1-12B_analysis_v4/myWGCNA_OL.Rdata")

# =========================
# FGSEA
# =========================
suppressPackageStartupMessages({
  library(fgsea)
  library(WGCNA)
  library(stringr)
  library(ggplot2)
})

# ---- 1) Map numeric kME columns (kME1..kME0) to color names; drop grey ----
stopifnot(exists("KMEs"))
kme_nums <- as.integer(sub("^kME", "", colnames(KMEs)))
kme_cols <- labels2colors(kme_nums)

keep_idx <- kme_nums != 0 & kme_cols != "grey"
KMEs_color <- KMEs[, keep_idx, drop = FALSE]
colnames(KMEs_color) <- paste0("kME_", kme_cols[keep_idx])

# quick sanity: list available colors
avail_colors <- sub("^kME_", "", colnames(KMEs_color))

# ---- 2) Load pathways (GO BP; HGNC symbols) and tidy names ----
pathways <- fgsea::gmtPathways('~/scratch/20250714_libs1-12B_analysis_v4/GO_Biological_Process_2025.txt')
names(pathways) <- stringr::str_replace(names(pathways), "\\s*\\([^\\)]+\\)$", "")

# ---- 3) Build ranks from kME for module across ALL genes ----
# select module of interest by color
target_color <- "orange"
target_col   <- paste0("kME_", target_color)

if (!target_col %in% colnames(KMEs_color)) {
  stop(sprintf(
    "Module '%s' not found. Available modules: %s",
    target_color, paste(avail_colors, collapse = ", ")
  ))
}

stats <- KMEs_color[, target_col]
stopifnot(!is.null(names(stats)) || !is.null(rownames(KMEs_color)))
names(stats) <- rownames(KMEs_color)

# clean and sort (descending order)
stats <- stats[!is.na(stats)]
stats <- stats[!duplicated(names(stats))]
stats <- sort(stats, decreasing = TRUE)

# ---- 4) Filter pathways by overlap and size ----
minSize <- 10; maxSize <- 500
paths_f <- lapply(pathways, function(gs) intersect(gs, names(stats)))
paths_f <- paths_f[lengths(paths_f) >= minSize & lengths(paths_f) <= maxSize]

if (length(paths_f) == 0) {
  overlaps <- vapply(pathways, function(gs) sum(gs %in% names(stats)), integer(1))
  stop(sprintf(
    "No GO terms passed size filters (minSize=%d, maxSize=%d).\nStats genes: %d | median overlap: %d.\nTry minSize=5 or check ID matching.",
    minSize, maxSize, length(stats), median(overlaps)
  ))
}

# ---- 5) Run FGSEA (multilevel) ----
res <- fgsea::fgseaMultilevel(pathways = paths_f, stats = stats, eps = 0.0)
res <- res[order(res$padj, -abs(res$NES)), ]
res$module_color <- target_color
# after you create `res`
res <- res[order(res$padj, -abs(res$NES)), ]
res$module_color <- target_color

# convert to data.frame so column subsetting works as expected
res_out <- as.data.frame(res)

# collapse list-columns (e.g., leadingEdge) to semicolon-joined strings
list_cols <- vapply(res_out, is.list, logical(1))
if (any(list_cols)) {
  res_out[list_cols] <- lapply(res_out[list_cols], function(x)
    vapply(x, function(item) paste0(as.character(item), collapse = ";"),
           FUN.VALUE = character(1))
  )
}

# write
outdir <- "~/scratch/20250714_libs1-12B_analysis_v4/fgsea"
dir.create(outdir, showWarnings = FALSE)
outfile <- file.path(outdir, paste0("fgsea_OL_black", target_color, ".tsv"))
write.table(res_out, outfile, sep = "\t", row.names = FALSE, quote = FALSE)


library(dplyr)
library(ggplot2)
library(stringr)
library(fgsea)



# --- pick top-N pathways (by FDR then |NES|) ---
topN <- 25
res_top <- res_out %>%
  arrange(padj, -abs(NES)) %>%filter(padj < 0.05) %>%
  slice_head(n = topN) %>%
  mutate(
    pathway_wrapped = str_wrap(pathway, width = 45),
    pathway_wrapped = factor(pathway_wrapped, levels = rev(pathway_wrapped)),
    neglog10padj = -log10(padj),
    dir = ifelse(NES >= 0, "Up", "Down")
  )

# 1) Dot plot: NES on x, pathways on y, point size = -log10(FDR)
p2 <- ggplot(res_top, aes(x = NES, y = pathway_wrapped)) +
  geom_point(aes(size = neglog10padj,
                 fill = dir),
             shape = 21, color = "black") +
  scale_size_continuous(name = "-log10(FDR)") +
  scale_fill_manual(values = c(Down = "steelblue", Up = "tomato")) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL,
       title = sprintf("Top %d pathways for OL orange module", topN)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())+  theme(plot.title = element_text(hjust = 0.5))


# 3) Classic FGSEA summary table plot for the top pathways
#    Uses your original 'paths_f' (list of gene sets) and 'stats' (named stat vector)
top_names <- res_top$pathway
fgsea::plotGseaTable(pathways = paths_f[top_names],
                     stats = stats,
                     fgseaRes = res_out[match(top_names, res_out$pathway), ])

