# Variance partitioning for pseudobulk SingleCellExperiment (sce)
# - Loads a pseudobulk sce 
# - Cleans metadata, coerces factor/continuous variables
# - Computes log2(FPM+1) for modeling
# - Checks multicollinearity with canCorPairs()
# - Runs variance partitioning and saves plots + summaries
# ==============================

  library(SingleCellExperiment)
  library(DESeq2)
  library(variancePartition)
  library(dplyr)
  library(Seurat)

#run in jupyter hub

seurat_obj <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")

#pseudobulk
subject_batch_id <- paste0(seurat_obj$SubjectID, "_", seurat_obj$Batch)
names(subject_batch_id) <- colnames(seurat_obj)
seurat_obj <- AddMetaData(seurat_obj, subject_batch_id, col.name = "subject_batch_id")


subinfo <- read.csv("~/scratch/20250714_libs1-12B_analysis_v4/20250807_subjectinfo_snRNAseq.csv")


subinfo <- subinfo %>%
  mutate(
    CompletedHighSchool = case_when(
      Education.level == "High school diploma obtained (or professional training complete)" ~ "Yes",
      Education.level == "Graduate level diploma obtained (Master/Doctorate)" ~ "Yes",
      Education.level == "High school diploma incomplete" ~ "No",
      Education.level ==  "Bachelor obtained" ~ "Yes",
      Education.level == "Less than 7 years of education"  ~ "No",
      Education.level ==  "Jr. High school diploma incomplete" ~ "No",
      Education.level ==  "College diploma obtained" ~ "Yes",
      Education.level == "DNK" ~NA,
      Education.level ==  "N/A" ~ NA))
    
    
    

# contingency table
tab <- table(seurat_obj$SubjectID, seurat_obj$SeqBatch)

# turn into a data frame
df <- as.data.frame.matrix(tab)

# assign the batch label (whichever column is non-zero)
df$SeqBatch <- apply(df, 1, function(x) names(x)[which.max(x)])

# keep only the label
seqbatches <- data.frame(
  SubjectID = rownames(df),
  SeqBatch  = df$SeqBatch
)

#remove mixed clusters  
seurat_obj <- subset(seurat_obj, idents = c("Mixed1", "Mixed2"), invert = TRUE)


#merge subject info

subinfo<- merge(subinfo, seqbatches, by.x = "Genotype.ID", by.y = "SubjectID")
subinfo$SubjectID <- gsub("_", "-", subinfo$Genotype.ID)


counts <- seurat_obj[["RNA"]]$counts
genes_to_keep <- rowSums(counts > 0) >= 10


seurat_obj_pb <- AggregateExpression(seurat_obj, assays = "RNA", group.by = c("SubjectID", "Batch"), return.seurat = TRUE, features =  rownames(counts)[genes_to_keep])

meta <- seurat_obj_pb@meta.data

info <- merge(meta, subinfo, by= "SubjectID")
rownames(info) <- info$orig.ident


# ------------------------------
# User inputs
# ------------------------------
outDir       <- "~/scratch/20250714_libs1-12B_analysis_v4/varpartition"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

counts_pb <- seurat_obj_pb@assays[["RNA"]]$counts


# Finalize rownames of colData to match counts columns exactly
stopifnot(identical(rownames(info), colnames(counts_pb)))  # critical for DESeq2

# ------------------------------
# Outputs + modeling choices
# ------------------------------
outDir <- "~/scratch/20250714_libs1-12B_analysis_v4/varpartition"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Choose covariates (tune as needed)
# Ensure these columns exist in `info` after the join
# Choose which variables to include in the variance partitioning formula

form <- ~ Group + Sex + Age + Batch + SeqBatch + Refr.Delay..h. + Known.substance.abuse.dependence + Known.recent.antidepressant + PMI + pH 

# ------------------------------
# Factors:
fac_vars <- c("Group", "Sex", "Batch", "SeqBatch", "Known.substance.abuse.dependence", "Known.recent.antidepressant")
for (v in intersect(fac_vars, colnames(info))) {
  info[[v]] <- as.character(info[[v]])
  info[[v]][is.na(info[[v]])] <- "Unknown"
  info[[v]] <- factor(info[[v]])
}


# Numeric:
num_vars <- c("Age", "PMI", "pH", "Refr.Delay..h.")
for (v in intersect(num_vars, colnames(info))) {
  info[[v]] <- suppressWarnings(as.numeric(info[[v]]))
  if (anyNA(info[[v]])) info[[v]][is.na(info[[v]])] <- mean(info[[v]], na.rm = TRUE)
}

# ------------------------------
# Build DESeq2 object to get size factors → FPM → log2(FPM+1)
# ------------------------------
if (any(counts_pb < 0, na.rm = TRUE)) stop("Counts contain negative values.")

dds <- DESeqDataSetFromMatrix(
  countData = counts_pb,
  colData   = info,
  design    = ~ 1
)
dds <- estimateSizeFactors(dds)

fpm_mat <- fpm(dds)
keep    <- rowSums(fpm_mat > 1) >= 0.5 * ncol(fpm_mat)
quantLog <- log2(fpm_mat[keep, , drop = FALSE] + 1)

message("Genes kept after FPM>1 filter: ", sum(keep), " / ", nrow(fpm_mat))

# ------------------------------
# Multicollinearity diagnostics (canonical correlations)
# ------------------------------
C <- canCorPairs(form, info)
saveRDS(C, file = file.path(outDir, "UF_pb_cca_result.rds"))
print(C)


   # shrink row/col label fonts and adjust margins
   plotCorrMatrix(
     C,
     cexRow   = 1,              # y-axis label size
     cexCol   = 1,              # x-axis label size
     margins  = c(15, 15),        # bottom/left margins (lines); increase if labels clip
     key.xlab = "correlation",    # legend title (optional)
   )
   
   
# ------------------------------
# Variance partitioning
# ------------------------------
do_varpart <- TRUE
if (do_varpart) {
  varPart <- fitExtractVarPartModel(quantLog, form, info)
  vp_sorted <- sortCols(varPart)
  medians <- apply(vp_sorted, 2, median)
  
  saveRDS(varPart, file = file.path(outDir, "UF_pb_varPart.rds"))
  write.csv(
    data.frame(covariate = names(medians), median_variance = medians),
    file = file.path(outDir, "UF_pb_median_variance_by_covariate.csv"),
    row.names = FALSE
  )
  pdf(file.path(outDir, "UF_pb_varPart_plot.pdf"), width = 9, height = 6)
  print(plotVarPart(vp_sorted))
  dev.off()
}
