################################################################################
# VIOLIN PLOT VIZ FOR DEG RESULTS
# Subject-level log2(CPM) per cluster — matches the logic of DEG pseudobulk code
# One dot per subject, violin plot (outline only) over subject-level distribution
################################################################################

  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(edgeR)
  library(ggplot2)
  library(dplyr)

# ========================
sce_path     <- "~/scratch/20250714_libs1-12B_analysis_v4/sce_20250815.rds"
deg_csv_path <- "GroupDE_celltype_subjectBatch_ge5cells_20260214_ALL_CELLTYPES.csv"

genes_of_interest    <- c("NECTIN3", "MIR646HG", "TRIM36")
#clusters_of_interest <- c("OL1", "OL3", "Inhib3")
clusters_of_interest <- c("OPC", "OL")

# minimum cells a subject must have in this cluster. 
#5 for individual clusters, 10 for broad clusters
min_cells_per_subject <- 10    
alpha                 <- 0.1

out_dir <- "./violin_plots"
dir.create(out_dir, showWarnings = FALSE)

# Plot aesthetics
group_colors       <- c("CTRL" = "#F8766D", "DS-CA" = "#00BFC4")
point_size         <- 2.0
point_alpha        <- 0.9
point_jitter_width <- 0.08
violin_linewidth   <- 0.7
violin_trim        <- TRUE
prior_count        <- 0.5     # pseudocount for log2(CPM), same as edgeR default

# ========================
# Load SCE — raw counts only, no normalization needed here
# (we compute log2CPM per subject from raw counts below)
# ========================
message("Loading SCE...")
sce <- readRDS(sce_path)

#use if individual clusters by uncommenting
#stopifnot("celltype" %in% colnames(colData(sce)))
#colData(sce)$cluster_id <- as.character(colData(sce)$celltype)

#use if broad clusters by uncommenting
stopifnot("clusters_broad" %in% colnames(colData(sce)))
colData(sce)$cluster_id <- as.character(colData(sce)$clusters_broad)


keep_cell <- !is.na(colData(sce)$cluster_id) &
  !(colData(sce)$cluster_id %in% c("Mixed1", "Mixed2")) #or Mixed for broad clusters
sce <- sce[, keep_cell]
sce <- sce[rowSums(counts(sce) > 0) > 0, ]

# Cell-level outlier removal (same as DEG script)
qc          <- scuttle::perCellQCMetrics(sce)
low_detected <- scuttle::isOutlier(metric = qc$detected, nmads = 3, log = TRUE, type = "lower")
sce <- sce[, !low_detected]

stopifnot(all(c("SubjectID", "Group") %in% colnames(colData(sce))))

# ========================
# Load DEG results
# ========================
message("Loading DEG results...")
deg <- read.csv(deg_csv_path, stringsAsFactors = FALSE)

# ========================
# significance label
# ========================
sig_label <- function(deg_row, alpha) {
  if (nrow(deg_row) == 0) return("not tested")
  fdr_val <- deg_row$adj.P.Val
  lfc     <- deg_row$logFC
  pval    <- deg_row$P.Value
  stars <- dplyr::case_when(
    fdr_val < 0.001 ~ "***",
    fdr_val < 0.01  ~ "**",
    fdr_val < alpha ~ "*",
    TRUE            ~ "ns"
  )
  sprintf("logFC=%.2f  p=%.2g  FDR=%.2g  %s", lfc, pval, fdr_val, stars)
}

# ========================
# main plot function
#   - Subset SCE to cluster
#   - Filter subjects with < min_cells_per_subject cells
#   - Aggregate raw counts per subject (gene counts + cluster-wide lib size)
#   - Compute log2(CPM) with prior_count via edgeR::cpm()
#   - Plot: violin (outline) + jittered subject-level dots
# ========================
make_violin <- function(sce, gene, cluster, deg,
                        min_cells_per_subject, alpha,
                        group_colors, prior_count,
                        point_size, point_alpha, point_jitter_width,
                        violin_linewidth, violin_trim) {
  
  # ---- subset to cluster ----
  sce_cl <- sce[, sce$cluster_id == cluster]
  if (ncol(sce_cl) == 0) {
    warning("No cells in cluster '", cluster, "'. Skipping.")
    return(NULL)
  }
  if (!gene %in% rownames(sce_cl)) {
    warning("Gene '", gene, "' not found in cluster '", cluster, "'. Skipping.")
    return(NULL)
  }
  
  # ---- filter subjects with too few cells in this cluster ----
  subj_vec     <- as.character(colData(sce_cl)$SubjectID)
  keep_subjects <- names(which(table(subj_vec) >= min_cells_per_subject))
  sce_cl <- sce_cl[, subj_vec %in% keep_subjects]
  if (ncol(sce_cl) == 0) {
    warning("No subjects with >= ", min_cells_per_subject,
            " cells in cluster '", cluster, "'. Skipping.")
    return(NULL)
  }
  subj_vec <- as.character(colData(sce_cl)$SubjectID)  # refresh
  grp_vec  <- as.character(colData(sce_cl)$Group)
  
  # ---- aggregate raw counts per subject ----
  C <- counts(sce_cl)  # genes x cells sparse matrix
  
  # Gene-level counts summed per subject
  gene_counts_cell <- as.numeric(C[gene, , drop = TRUE])
  gene_counts_subj <- tapply(gene_counts_cell, subj_vec, sum, na.rm = TRUE)
  
  # Cluster-wide library sizes per subject (sum across ALL genes in this cluster)
  libsize_cell <- Matrix::colSums(C)
  libsize_subj <- tapply(libsize_cell, subj_vec, sum, na.rm = TRUE)
  
  # Subject -> group mapping (one row per subject)
  map_df <- data.frame(
    sample_id = subj_vec,
    group_id  = grp_vec,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(group_id)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  # Align across subjects
  subjects <- intersect(names(gene_counts_subj), names(libsize_subj))
  subjects <- intersect(subjects, map_df$sample_id)
  if (!length(subjects)) {
    warning("No subjects remain after alignment for cluster '", cluster, "'. Skipping.")
    return(NULL)
  }
  gene_counts_subj <- gene_counts_subj[subjects]
  libsize_subj     <- libsize_subj[subjects]
  
  # ---- compute log2(CPM) per subject using cluster-wide lib size ----
  # counts for this gene / total cluster counts for this subject, log-scaled.
  y <- edgeR::DGEList(
    matrix(gene_counts_subj, nrow = 1, dimnames = list(gene, subjects))
  )
  y$samples$lib.size     <- as.numeric(libsize_subj[subjects])
  y$samples$norm.factors <- 1   # no TMM — keep consistent with pseudobulk DEG model
  logcpm <- as.numeric(edgeR::cpm(y, log = TRUE, prior.count = prior_count))
  
  # ---- build plot data frame ----
  df <- data.frame(
    sample_id = subjects,
    expr      = logcpm,
    group_id  = map_df$group_id[match(subjects, map_df$sample_id)],
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(group_id), is.finite(expr))
  
  if (!nrow(df)) {
    warning("No plottable rows for gene '", gene, "' in cluster '", cluster, "'. Skipping.")
    return(NULL)
  }
  
  # ---- restrict to groups in color palette ----
  present_groups <- intersect(names(group_colors), unique(df$group_id))
  df <- dplyr::filter(df, group_id %in% present_groups)
  df$group_id <- factor(df$group_id, levels = present_groups)
  
  # ---- guard: violin needs >= 2 subjects with >= 2 distinct values ----
  dens_ok <- df |>
    dplyr::group_by(group_id) |>
    dplyr::summarise(n = dplyr::n(), nuniq = dplyr::n_distinct(expr), .groups = "drop") |>
    dplyr::filter(n >= 2, nuniq >= 2) |>
    dplyr::pull(group_id) |>
    as.character()
  df_violin <- dplyr::filter(df, group_id %in% dens_ok)
  
  colors_present <- group_colors[present_groups]
  
  # ---- plot ----
  ggplot() +
    
    # Violin: colored outline, transparent fill
    { if (nrow(df_violin) > 0)
      geom_violin(
        data    = df_violin,
        mapping = aes(x = group_id, y = expr, color = group_id),
        trim    = violin_trim, scale = "width", na.rm = TRUE,
        fill    = NA, linewidth = violin_linewidth
      )
    } +
    
    # One jittered dot per subject
    geom_jitter(
      data    = df,
      mapping = aes(x = group_id, y = expr, color = group_id),
      width   = point_jitter_width, height = 0,
      size    = point_size, alpha = point_alpha
    ) +
    
    scale_color_manual(values = colors_present, name = NULL) +
    
    labs(
      title = paste0(gene, " | ", cluster),
      x     = NULL,
      y     = "log2(CPM)"
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      axis.text.x     = element_text(size = 11, face = "bold")
    )
}

# ========================
# select clusters to iterate
# ========================
if (is.null(clusters_of_interest) || length(clusters_of_interest) == 0) {
  clusters_to_plot <- sort(unique(deg$cluster_id))
} else {
  clusters_to_plot <- clusters_of_interest
}

# ========================
# make and save plots
# ========================
for (cl in clusters_to_plot) {
  for (gn in genes_of_interest) {
    message("Plotting: ", gn, " | ", cl)
    
    p <- make_violin(
      sce                   = sce,
      gene                  = gn,
      cluster               = cl,
      deg                   = deg,
      min_cells_per_subject = min_cells_per_subject,
      alpha                 = alpha,
      group_colors          = group_colors,
      prior_count           = prior_count,
      point_size            = point_size,
      point_alpha           = point_alpha,
      point_jitter_width    = point_jitter_width,
      violin_linewidth      = violin_linewidth,
      violin_trim           = violin_trim
    )
    
    if (is.null(p)) next
    
    fname <- file.path(
      out_dir,
      paste0("violin_", gsub("[^A-Za-z0-9_]", "_", gn),
             "_", gsub("[^A-Za-z0-9_]", "_", cl), ".pdf")
    )
    
    ggsave(fname, plot = p, width = 4, height = 5, device = "pdf")
    message("  Saved: ", fname)
  }
}

message("\nDone. All plots saved to: ", out_dir)