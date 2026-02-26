################################################################################
# GROUP DEG per BROAD CLUSTER (pseudobulk limma-voom)
# Design: ~ 0 + group_id + Age + Sex + pH + Batch
# Contrast: DS-CA vs CTRL
# Method: 3-pass voomWithQualityWeights + duplicateCorrelation pipeline
################################################################################

  library(muscat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(scuttle)
  library(openxlsx)
  library(rlang)

#run in jupyter hub

# ========================
#  topTable wrapper
# ========================
safe_topTable_1coef <- function(fit2, number = Inf, sort.by = "P") {
  tt <- tryCatch(
    limma::topTable(fit2, number = number, sort.by = sort.by),
    error = function(e) NULL
  )
  if (!is.null(tt) && nrow(tt) > 0) return(tt)
  
  stopifnot(!is.null(fit2$coefficients))
  stopifnot(ncol(fit2$coefficients) == 1)
  
  rn <- rownames(fit2$coefficients)
  pv <- fit2$p.value[, 1]
  
  tt <- data.frame(
    logFC     = fit2$coefficients[, 1],
    AveExpr   = fit2$Amean,
    t         = fit2$t[, 1],
    P.Value   = pv,
    adj.P.Val = p.adjust(pv, method = "BH"),
    B         = fit2$lods[, 1],
    row.names = rn,
    check.names = FALSE
  )
  
  if (!is.null(sort.by) && toupper(sort.by) == "P") {
    finite_mask <- is.finite(tt$P.Value)
    if (any(finite_mask)) {
      ord <- order(tt$P.Value[finite_mask], decreasing = FALSE)
      tt  <- rbind(tt[finite_mask, , drop = FALSE][ord, , drop = FALSE],
                   tt[!finite_mask, , drop = FALSE])
    }
  }
  if (is.finite(number)) tt <- utils::head(tt, number)
  tt
}

# ========================
# Parameters
# ========================
alpha               <- 0.1 
min_cells_per_key   <- 10 #min number of cells to be included
min_CTRL_subjects   <- 12 #min CTRL
min_DSCA_subjects   <- 15 #min DS-CA
use_quality_weights <- TRUE

date_tag <- format(Sys.Date(), "%Y%m%d")

# ========================
# Load data
# ========================
sce <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/sce_20250815.rds")

# Set cluster variable
stopifnot("clusters_broad" %in% colnames(colData(sce)))
colData(sce)$cluster_id <- as.character(colData(sce)$clusters_broad)

# Exclude Mixed cluster 
keep_cell <- !is.na(colData(sce)$cluster_id) & 
  !(colData(sce)$cluster_id %in% c("Mixed"))
sce <- sce[, keep_cell]

# Gene filtering
sce <- sce[rowSums(counts(sce) > 0) > 0, ]

# Cell QC: remove low detected-gene cells
qc <- perCellQCMetrics(sce)
low_detected <- isOutlier(metric = qc$detected, nmads = 3, log = TRUE, type = "lower")
sce <- sce[, !low_detected]

# Gene QC
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# Compute size factors and log-normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

# ========================
# Setup IDs
# ========================
stopifnot(all(c("SubjectID","Group","Batch") %in% colnames(colData(sce))))
colData(sce)$sample_id <- as.character(colData(sce)$SubjectID)
sce <- prepSCE(sce, kid = "cluster_id", gid = "Group", sid = "sample_id", drop = FALSE)
colData(sce)$subject_batch_id <- paste0(colData(sce)$sample_id, "_", colData(sce)$Batch)

# ========================
# Helper: collapse metadata by subject×Batch
# ========================
collapse_md_by_key <- function(sce_sub, key_col = "subject_batch_id") {
  md0 <- as.data.frame(SummarizedExperiment::colData(sce_sub))
  need_cols <- intersect(
    c("Age","Sex","pH","Batch","sample_id","subject_batch_id","group_id"),
    names(md0)
  )
  md0 <- md0[, need_cols, drop = FALSE]
  key_sym <- rlang::sym(key_col)
  
  md0 %>%
    dplyr::group_by(!!key_sym) %>%
    dplyr::summarise(
      Age       = mean(as.numeric(Age), na.rm = TRUE),
      pH        = mean(as.numeric(pH),  na.rm = TRUE),
      Sex       = dplyr::first(Sex[!is.na(Sex)]),
      Batch     = dplyr::first(Batch[!is.na(Batch)]),
      sample_id = dplyr::first(sample_id[!is.na(sample_id)]),
      group_id  = dplyr::first(group_id[!is.na(group_id)]),
      .groups   = "drop"
    ) %>%
    as.data.frame()
}

# ========================
# Main loop
# ========================
all_results <- list()
zeros_list  <- list()

broad_clusters <- sort(unique(as.character(colData(sce)$cluster_id)))

for (cl in broad_clusters) {
  message("Analyzing broad cluster: ", cl)
  
  sce_sub <- sce[, sce$cluster_id == cl]
  
  # Filter by minimum cells per subject×batch
  keepers <- names(which(table(colData(sce_sub)$subject_batch_id) >= min_cells_per_key))
  if (length(keepers) < 2L) {
    message("  Skip: insufficient subject_batch_id with sufficient cells")
    next
  }
  sce_sub <- sce_sub[, colData(sce_sub)$subject_batch_id %in% keepers]
  
  # Check group balance
  tab <- table(sce_sub$sample_id, sce_sub$group_id)
  n_CTRL  <- if ("CTRL"  %in% colnames(tab)) sum(tab[, "CTRL"]  > 0) else 0
  n_DS_CA <- if ("DS-CA" %in% colnames(tab)) sum(tab[, "DS-CA"] > 0) else 0
  
  if (n_CTRL < min_CTRL_subjects || n_DS_CA < min_DSCA_subjects) {
    message("  Skip: insufficient subjects per group")
    next
  }
  
  # Pseudobulk aggregation
  pb <- muscat::aggregateData(
    x = sce_sub, assay = "counts", fun = "sum",
    by = c("cluster_id", "subject_batch_id")
  )
  counts <- assay(pb, cl)
  
  # Match metadata to pseudobulk columns
  md_key  <- collapse_md_by_key(sce_sub, "subject_batch_id")
  key_vals <- colnames(counts)
  row_idx <- match(key_vals, md_key$subject_batch_id)
  if (anyNA(row_idx)) {
    message("  Skip: metadata mismatch")
    next
  }
  md <- md_key[row_idx, , drop = FALSE]
  rownames(md) <- key_vals
  
  # Convert to factors (make.names converts DS-CA to DS.CA)
  md$group_id <- factor(make.names(as.character(md$group_id)))
  if ("CTRL" %in% levels(md$group_id)) md$group_id <- stats::relevel(md$group_id, ref = "CTRL")
  md$sample_id <- factor(md$sample_id)
  md$Batch     <- droplevels(factor(md$Batch))
  md$Sex       <- droplevels(factor(md$Sex))
  md$Age       <- suppressWarnings(as.numeric(md$Age))
  md$pH        <- suppressWarnings(as.numeric(md$pH))
  
  # Remove rows with missing covariates
  covars <- c("group_id", "Age", "Sex", "pH", "Batch")
  ok <- stats::complete.cases(md[, covars, drop = FALSE])
  if (!all(ok)) {
    md     <- droplevels(md[ok, , drop = FALSE])
    counts <- counts[, ok, drop = FALSE]
  }
  
  if (nlevels(md$group_id) < 2) {
    message("  Skip: insufficient groups")
    next
  }
  if (nlevels(md$Batch) < 2) {
    message("  Skip: insufficient batch levels")
    next
  }
  
  stopifnot(ncol(counts) == nrow(md), identical(colnames(counts), rownames(md)))
  
  # Design matrix
  design <- model.matrix(~ 0 + group_id + Age + Sex + pH + Batch, data = md)
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  
  # Check estimability
  ne <- limma::nonEstimable(design)
  if (!is.null(ne)) {
    message("  Skip: design not estimable")
    next
  }
  
  # edgeR normalization and filtering
  y <- DGEList(counts)
  keep <- filterByExpr(y, design)
  if (!any(keep)) {
    message("  Skip: no genes pass filtering")
    next
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  # Check residual degrees of freedom
  resid_df <- ncol(y$counts) - qr(design)$rank
  if (resid_df <= 0) {
    message("  Skip: no residual df")
    next
  }
  
  # 3-pass voomWithQualityWeights + duplicateCorrelation
  reps_per_subj <- table(md$sample_id)
  use_dup_cor <- any(reps_per_subj >= 2)
  dupCor_val <- NA_real_
  
  if (!use_dup_cor) {
    v <- if (use_quality_weights) {
      limma::voomWithQualityWeights(y, design, plot = FALSE)
    } else {
      limma::voom(y, design, plot = FALSE)
    }
    fit <- limma::lmFit(v, design)
    
  } else {
    # Pass 1
    v1 <- if (use_quality_weights) {
      limma::voomWithQualityWeights(y, design, plot = FALSE)
    } else {
      limma::voom(y, design, plot = FALSE)
    }
    
    corfit1 <- limma::duplicateCorrelation(v1, design, block = md$sample_id)
    rho1 <- corfit1$consensus
    
    if (!is.finite(rho1) || abs(rho1) >= 0.99) {
      message("  Invalid correlation, fitting without block")
      use_dup_cor <- FALSE
      v <- v1
      fit <- limma::lmFit(v, design)
      
    } else {
      # Pass 2
      v2 <- if (use_quality_weights) {
        limma::voomWithQualityWeights(y, design, plot = FALSE,
                                      block = md$sample_id, correlation = rho1)
      } else {
        limma::voom(y, design, plot = FALSE,
                    block = md$sample_id, correlation = rho1)
      }
      
      corfit2 <- limma::duplicateCorrelation(v2, design, block = md$sample_id)
      rho2 <- corfit2$consensus
      if (!is.finite(rho2) || abs(rho2) >= 0.99) {
        rho2 <- rho1
      }
      dupCor_val <- rho2
      
      # Pass 3
      v <- if (use_quality_weights) {
        limma::voomWithQualityWeights(y, design, plot = FALSE,
                                      block = md$sample_id, correlation = dupCor_val)
      } else {
        limma::voom(y, design, plot = FALSE,
                    block = md$sample_id, correlation = dupCor_val)
      }
      
      fit <- limma::lmFit(v, design, block = md$sample_id, correlation = dupCor_val)
    }
  }
  
  # Contrast: DS.CA - CTRL
  lhs <- "group_idDS.CA"
  rhs <- "group_idCTRL"
  if (!all(c(lhs, rhs) %in% colnames(design))) {
    message("  Skip: missing group columns in design")
    next
  }
  
  cmat <- limma::makeContrasts(contrasts = paste(lhs, "-", rhs), levels = design)
  fitc <- limma::contrasts.fit(fit, cmat)
  fit2 <- limma::eBayes(fitc, trend = TRUE, robust = TRUE)
  
  # Extract results
  tt <- safe_topTable_1coef(fit2, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  
  # Annotate results
  tt$cluster_id       <- cl
  tt$n_cols_kept      <- ncol(y)
  tt$n_CTRL_subjects  <- n_CTRL
  tt$n_DSCA_subjects  <- n_DS_CA
  tt$dupCor           <- dupCor_val
  tt$used_dupCor      <- use_dup_cor
  tt$used_qw          <- use_quality_weights
  tt$significant_FDR  <- tt$adj.P.Val < alpha
  tt$model            <- "GroupDE"
  
  first <- intersect(c("gene", "cluster_id", "model"), colnames(tt))
  tt <- tt[, c(first, setdiff(colnames(tt), first))]
  all_results[[cl]] <- tt
  
  # Zero-count diagnostics
  counts_fit <- y$counts
  gvec       <- as.character(md$group_id)
  subj_vec   <- as.character(md$sample_id)
  
  ctrl_mask_cols <- gvec == "CTRL"
  dsca_mask_cols <- gvec == "DS.CA"
  
  n_zero_CTRL_cols <- rowSums(counts_fit[, ctrl_mask_cols, drop = FALSE] == 0)
  n_zero_DSCA_cols <- rowSums(counts_fit[, dsca_mask_cols, drop = FALSE] == 0)
  n_CTRL_cols      <- sum(ctrl_mask_cols)
  n_DSCA_cols      <- sum(dsca_mask_cols)
  
  # Subject-level collapse
  subj_factor <- factor(subj_vec, levels = unique(subj_vec))
  counts_by_subject <- t(rowsum(t(counts_fit), group = subj_factor, reorder = FALSE))
  
  subj_levels <- levels(subj_factor)
  subj_group  <- factor(
    sapply(subj_levels, function(s) unique(na.omit(gvec[subj_vec == s]))[1]),
    levels = levels(md$group_id)
  )
  
  ctrl_mask_subj <- subj_group == "CTRL"
  dsca_mask_subj <- subj_group == "DS.CA"
  
  n_zero_CTRL_subjects <- if (any(ctrl_mask_subj)) rowSums(counts_by_subject[, ctrl_mask_subj, drop = FALSE] == 0) else 0L
  n_zero_DSCA_subjects <- if (any(dsca_mask_subj)) rowSums(counts_by_subject[, dsca_mask_subj, drop = FALSE] == 0) else 0L
  n_CTRL_subjects_post <- sum(ctrl_mask_subj)
  n_DSCA_subjects_post <- sum(dsca_mask_subj)
  
  zeros_df <- data.frame(
    gene                   = rownames(counts_fit),
    cluster_id             = cl,
    model                  = "GroupDE",
    n_zero_CTRL_cols       = as.integer(n_zero_CTRL_cols),
    n_zero_DSCA_cols       = as.integer(n_zero_DSCA_cols),
    n_CTRL_cols            = n_CTRL_cols,
    n_DSCA_cols            = n_DSCA_cols,
    pct_zero_CTRL_cols     = if (n_CTRL_cols > 0) n_zero_CTRL_cols / n_CTRL_cols else NA_real_,
    pct_zero_DSCA_cols     = if (n_DSCA_cols > 0) n_zero_DSCA_cols / n_DSCA_cols else NA_real_,
    n_zero_CTRL_subjects   = as.integer(n_zero_CTRL_subjects),
    n_zero_DSCA_subjects   = as.integer(n_zero_DSCA_subjects),
    n_CTRL_subjects        = n_CTRL_subjects_post,
    n_DSCA_subjects        = n_DSCA_subjects_post,
    pct_zero_CTRL_subjects = if (n_CTRL_subjects_post > 0) n_zero_CTRL_subjects / n_CTRL_subjects_post else NA_real_,
    pct_zero_DSCA_subjects = if (n_DSCA_subjects_post > 0) n_zero_DSCA_subjects / n_DSCA_subjects_post else NA_real_,
    stringsAsFactors       = FALSE
  )
  
  tt_key <- tt[, c("gene", "logFC", "P.Value", "adj.P.Val", "model")]
  zeros_df <- dplyr::left_join(zeros_df, tt_key, by = c("gene", "model"))
  zeros_list[[cl]] <- zeros_df
}

stopifnot(length(all_results) > 0)

# ========================
# Export results
# ========================
res <- dplyr::bind_rows(all_results) %>%
  dplyr::mutate(
    FDR_global        = p.adjust(P.Value, method = "BH"),
    significant_global = FDR_global < alpha
  )

tag <- paste0("GroupDE_BROAD_subjectBatch_ge", min_cells_per_key, "cells_", date_tag)

# Combined CSV
write.csv(res, paste0(tag, "_ALL_CLUSTERS.csv"), row.names = FALSE)

# Excel workbook (one sheet per cluster)
wb <- openxlsx::createWorkbook()
by_cl <- split(res, res$cluster_id)

for (cl in names(by_cl)) {
  sheet <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", cl), 1, 31)
  if (sheet %in% names(wb)) sheet <- paste0(substr(sheet, 1, 28), "_", sample(1:9, 1))
  openxlsx::addWorksheet(wb, sheet)
  out <- by_cl[[cl]]
  first <- intersect(c("gene", "cluster_id", "model"), colnames(out))
  out <- out[, c(first, setdiff(colnames(out), first))]
  openxlsx::writeData(wb, sheet, out)
}
openxlsx::saveWorkbook(wb, paste0(tag, "_BY_CLUSTER.xlsx"), overwrite = TRUE)

# Summary table
deg_summary <- res %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(
    n_genes_tested   = dplyr::n(),
    n_sig_adjP       = sum(adj.P.Val < alpha, na.rm = TRUE),
    n_up_adjP        = sum(adj.P.Val < alpha & logFC > 0, na.rm = TRUE),
    n_down_adjP      = sum(adj.P.Val < alpha & logFC < 0, na.rm = TRUE),
    n_sig_globalFDR  = sum(FDR_global < alpha, na.rm = TRUE),
    n_up_globalFDR   = sum(FDR_global < alpha & logFC > 0, na.rm = TRUE),
    n_down_globalFDR = sum(FDR_global < alpha & logFC < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_sig_adjP))

write.csv(deg_summary, paste0(tag, "_SUMMARY.csv"), row.names = FALSE)
print(deg_summary)

# Zero-counts CSV
zeros_all <- dplyr::bind_rows(zeros_list) %>%
  dplyr::left_join(res[, c("cluster_id", "gene", "FDR_global", "model")],
                   by = c("cluster_id", "gene", "model")) %>%
  dplyr::distinct(cluster_id, gene, model, .keep_all = TRUE)

write.csv(zeros_all, paste0(tag, "_ZEROCOUNTS.csv"), row.names = FALSE)

message("\nAnalysis complete. Files written:")
message("  - ", tag, "_ALL_CLUSTERS.csv")
message("  - ", tag, "_BY_CLUSTER.xlsx")
message("  - ", tag, "_SUMMARY.csv")
message("  - ", tag, "_ZEROCOUNTS.csv")