################################################################################
# AGE EFFECT per BROAD CLUSTER (pseudobulk limma-voom)
# Design: ~ Age_z + group_id + Sex + pH + Batch
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
  library(ggplot2)

#run in jupyter hub

# ========================
# Safe topTable wrapper
# ========================
safe_topTable <- function(fit2, coef_name, number = Inf, sort.by = "P") {
  tt <- tryCatch(
    limma::topTable(fit2, coef = coef_name, number = number, sort.by = sort.by),
    error = function(e) NULL
  )
  if (!is.null(tt) && nrow(tt) > 0) return(tt)
  
  if (!coef_name %in% colnames(fit2$coefficients)) {
    stop("safe_topTable: coef '", coef_name, "' not found in fit2$coefficients.")
  }
  
  rn <- rownames(fit2$coefficients)
  pv <- fit2$p.value[, coef_name]
  
  tt <- data.frame(
    logFC     = fit2$coefficients[, coef_name],
    AveExpr   = fit2$Amean,
    t         = fit2$t[, coef_name],
    P.Value   = pv,
    adj.P.Val = p.adjust(pv, method = "BH"),
    B         = fit2$lods[, coef_name],
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
alpha               <- 0.10
min_cells_per_sbid  <- 10 #min number of cells to be included
min_subjects_needed <- 27
use_quality_weights <- TRUE

# Plotting options
do_plots               <- TRUE
plot_top_n_per_cluster <- Inf
plot_output_dir        <- "AgeEffect_PLOTS_AgeScatter_BROAD"

if (isTRUE(do_plots)) {
  if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)
}

date_tag <- format(Sys.Date(), "%Y%m%d")

# ========================
# Load data
# ========================
sce <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/sce_20250815.rds")

# Set cluster variable
stopifnot("clusters_broad" %in% colnames(colData(sce)))
colData(sce)$cluster_id <- as.character(colData(sce)$clusters_broad)

# Exclude Mixed clusters 
drop_mask <- !is.na(colData(sce)$cluster_id) & 
  !(colData(sce)$cluster_id %in% c("Mixed"))
sce <- sce[, drop_mask]

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
# collapse metadata by subject×Batch
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
all_results_linear <- list()
broad_clusters <- sort(unique(as.character(colData(sce)$cluster_id)))

for (cl in broad_clusters) {
  message("Analyzing cluster: ", cl)
  
  sce_sub <- sce[, sce$cluster_id == cl]
  
  # Check minimum subjects
  n_subjects_pre <- length(unique(as.character(colData(sce_sub)$sample_id)))
  if (n_subjects_pre < min_subjects_needed) {
    message("  Skip: insufficient subjects")
    next
  }
  
  # Filter by minimum cells per subject×batch
  keep_keys <- names(which(table(colData(sce_sub)$subject_batch_id) >= min_cells_per_sbid))
  if (length(keep_keys) < 2L) {
    message("  Skip: <2 subject_batch_id with sufficient cells")
    next
  }
  sce_sub <- sce_sub[, colData(sce_sub)$subject_batch_id %in% keep_keys]
  
  n_subjects_after <- length(unique(as.character(colData(sce_sub)$sample_id)))
  if (n_subjects_after < min_subjects_needed) {
    message("  Skip: insufficient subjects after filtering")
    next
  }
  
  # Pseudobulk aggregation
  pb <- muscat::aggregateData(
    x = sce_sub, assay = "counts", fun = "sum",
    by = c("cluster_id", "subject_batch_id")
  )
  counts_mat <- assay(pb, cl)
  
  # Match metadata to pseudobulk columns
  md_key <- collapse_md_by_key(sce_sub, "subject_batch_id")
  key_vals <- colnames(counts_mat)
  row_idx <- match(key_vals, md_key$subject_batch_id)
  if (anyNA(row_idx)) {
    message("  Skip: metadata mismatch")
    next
  }
  md <- md_key[row_idx, , drop = FALSE]
  rownames(md) <- key_vals
  
  # Convert to factors
  md$sample_id <- factor(md$sample_id)
  md$Batch     <- droplevels(factor(md$Batch))
  md$Sex       <- droplevels(factor(md$Sex))
  md$group_id  <- droplevels(factor(md$group_id))
  if ("CTRL" %in% levels(md$group_id)) md$group_id <- stats::relevel(md$group_id, ref = "CTRL")
  md$Age <- suppressWarnings(as.numeric(md$Age))
  md$pH  <- suppressWarnings(as.numeric(md$pH))
  
  # Remove rows with missing covariates
  covars_pre <- c("Age","group_id","Sex","pH","Batch","sample_id")
  ok_pre <- stats::complete.cases(md[, covars_pre, drop = FALSE])
  if (!all(ok_pre)) {
    md         <- droplevels(md[ok_pre, , drop = FALSE])
    counts_mat <- counts_mat[, ok_pre, drop = FALSE]
  }
  
  if (nlevels(md$Batch) < 2) {
    message("  Skip: insufficient batch levels")
    next
  }
  
  # Global age z-scoring
  Age_z <- as.numeric(scale(md$Age))
  if (!all(is.finite(Age_z)) || sd(Age_z, na.rm = TRUE) == 0) {
    message("  Skip: invalid Age_z")
    next
  }
  md$Age_z <- Age_z
  
  # Final completeness check
  covars <- c("Age_z","group_id","Sex","pH","Batch","sample_id")
  ok <- stats::complete.cases(md[, covars, drop = FALSE])
  if (!all(ok)) {
    md         <- droplevels(md[ok, , drop = FALSE])
    counts_mat <- counts_mat[, ok, drop = FALSE]
  }
  
  n_subjects_final <- length(unique(as.character(md$sample_id)))
  if (n_subjects_final < min_subjects_needed) {
    message("  Skip: insufficient subjects after Age_z filtering")
    next
  }
  
  stopifnot(ncol(counts_mat) == nrow(md), identical(colnames(counts_mat), rownames(md)))
  
  # Design matrix
  design <- model.matrix(~ Age_z + group_id + Sex + pH + Batch, data = md)
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  
  # Check estimability
  ne <- limma::nonEstimable(design)
  if (!is.null(ne)) {
    message("  Skip: design not estimable")
    next
  }
  
  # edgeR normalization and filtering
  y <- DGEList(counts_mat)
  keep_genes <- filterByExpr(y, design = design)
  if (!any(keep_genes)) {
    message("  Skip: no genes pass filtering")
    next
  }
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
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
  dupCor_val  <- NA_real_
  
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
  
  # eBayes and extract results
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  
  coef_name <- "Age_z"
  if (!coef_name %in% colnames(fit2$coefficients)) {
    message("  Skip: Age_z not in coefficients")
    next
  }
  
  tt <- safe_topTable(fit2, coef_name = coef_name, number = Inf, sort.by = "P")
  
  # Annotate results
  tt$gene            <- rownames(tt)
  tt$cluster_id      <- cl
  tt$n_cols_kept     <- ncol(y$counts)
  tt$n_genes_tested  <- nrow(y$counts)
  tt$n_subjects      <- n_subjects_final
  tt$dupCor          <- dupCor_val
  tt$used_dupCor     <- use_dup_cor
  tt$used_qw         <- use_quality_weights
  tt$significant     <- tt$adj.P.Val < alpha
  
  first <- intersect(c("gene","cluster_id"), colnames(tt))
  tt <- tt[, c(first, setdiff(colnames(tt), first))]
  all_results_linear[[cl]] <- tt
  
  # Optional plotting
  if (isTRUE(do_plots)) {
    sig_genes <- tt %>%
      dplyr::filter(adj.P.Val < alpha) %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::pull(gene)
    
    if (length(sig_genes) > 0) {
      sel_genes <- if (is.infinite(plot_top_n_per_cluster)) sig_genes else utils::head(sig_genes, plot_top_n_per_cluster)
      
      facets_per_page <- 24
      gene_pages <- split(sel_genes, ceiling(seq_along(sel_genes) / facets_per_page))
      n_pages <- length(gene_pages)
      
      page_idx <- 1L
      for (genes_page in gene_pages) {
        df_plot <- lapply(genes_page, function(g) {
          data.frame(
            gene = g,
            expr = as.numeric(v$E[g, , drop = TRUE]),
            Age  = md$Age,
            stringsAsFactors = FALSE
          )
        }) %>% dplyr::bind_rows()
        
        df_plot$gene <- factor(df_plot$gene, levels = genes_page)
        
        p <- ggplot(df_plot, aes(x = Age, y = expr)) +
          geom_point(alpha = 0.8, size = 1.5) +
          geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "black") +
          facet_wrap(~ gene, scales = "free_y", ncol = 4) +
          labs(
            title = paste0("Age vs expression — ", cl),
            subtitle = paste0("FDR < ", alpha, " | page ", page_idx, "/", n_pages),
            x = "Age (years)", y = "Expression (voom logCPM)"
          ) +
          theme_bw(base_size = 11) +
          theme(panel.grid.minor = element_blank())
        
        height_in <- max(5, 2.2 * 6)
        outfile <- if (n_pages > 1) {
          file.path(plot_output_dir, sprintf("AgeScatter_%s_p%02d_%s.png", cl, page_idx, date_tag))
        } else {
          file.path(plot_output_dir, sprintf("AgeScatter_%s_%s.png", cl, date_tag))
        }
        
        ggsave(outfile, plot = p, width = 12, height = height_in, dpi = 300)
        page_idx <- page_idx + 1L
      }
    }
  }
}

# ========================
# Export results
# ========================
stopifnot(length(all_results_linear) > 0)

res_lin <- dplyr::bind_rows(all_results_linear) %>%
  dplyr::mutate(
    FDR_global         = p.adjust(P.Value, method = "BH"),
    significant_global = FDR_global < alpha
  )

tag <- paste0("AgeDE_BROAD_subjectBatch_ge", min_cells_per_sbid,
              "cells_minSubj", min_subjects_needed,
              "_QW", use_quality_weights,
              "_GlobalAge_3pass_", date_tag)

# Combined CSV
write.csv(res_lin, paste0(tag, "_ALL_CLUSTERS.csv"), row.names = FALSE)

# Excel workbook (one sheet per cluster)
wb_lin <- openxlsx::createWorkbook()
for (cl in names(all_results_linear)) {
  sheet <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", cl), 1, 31)
  if (sheet %in% names(wb_lin)) sheet <- paste0(substr(sheet, 1, 28), "_", sample(1:9, 1))
  openxlsx::addWorksheet(wb_lin, sheet)
  out <- all_results_linear[[cl]]
  first <- intersect(c("gene","cluster_id"), colnames(out))
  out <- out[, c(first, setdiff(colnames(out), first))]
  openxlsx::writeData(wb_lin, sheet, out)
}
openxlsx::saveWorkbook(wb_lin, paste0(tag, "_BY_CLUSTER.xlsx"), overwrite = TRUE)

# Summary table
deg_summary_lin <- res_lin %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(
    n_genes_tested   = max(n_genes_tested, na.rm = TRUE),
    n_sig_adjP       = sum(adj.P.Val < alpha, na.rm = TRUE),
    n_up_adjP        = sum(adj.P.Val < alpha & logFC > 0, na.rm = TRUE),
    n_down_adjP      = sum(adj.P.Val < alpha & logFC < 0, na.rm = TRUE),
    n_sig_globalFDR  = sum(FDR_global < alpha, na.rm = TRUE),
    n_up_globalFDR   = sum(FDR_global < alpha & logFC > 0, na.rm = TRUE),
    n_down_globalFDR = sum(FDR_global < alpha & logFC < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>% dplyr::arrange(dplyr::desc(n_sig_adjP))

write.csv(deg_summary_lin, paste0(tag, "_SUMMARY.csv"), row.names = FALSE)
print(deg_summary_lin)

message("\nAnalysis complete. Files written:")
message("  - ", tag, "_ALL_CLUSTERS.csv")
message("  - ", tag, "_BY_CLUSTER.xlsx")
message("  - ", tag, "_SUMMARY.csv")
if (isTRUE(do_plots)) {
  message("  - Plots: ", plot_output_dir)
}