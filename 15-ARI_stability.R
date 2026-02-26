library(Seurat)
#calculate clustering stability by randomly subsetting seurat object with defined cell types to 75% of the data, 
#then reclustering it with the same parameters 100 times and and calculating the average adjusted rand index (ARI)
#to see how similar they are

#load seurat object
libs_UF_all_harmony_integrated <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")


compute_subset_ari <- function(
    seurat_obj,
    orig_clusters = NULL,          # defaults to Idents(seurat_obj)
    reduction = "harmony",         # e.g., "harmony" or "pca"
    dims = 1:25,
    prop = 0.75,
    resolution = 0.5,
    n_iter = 100,
    seed = 123
) {
  # --- deps ---
  require(Seurat)
  require(mclust)
  
  # --- original labels (from full object) ---
  if (is.null(orig_clusters)) {
    orig_clusters <- Idents(seurat_obj)
  }
  if (length(orig_clusters) != ncol(seurat_obj)) {
    stop("`orig_clusters` must be a vector with one label per cell in `seurat_obj` (same length and cell names).")
  }
  # Make sure names are cell barcodes
  orig_clusters <- setNames(as.character(orig_clusters), colnames(seurat_obj))
  
  # --- storage ---
  ari_vals <- numeric(n_iter)
  iters <- seq_len(n_iter)
  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  
  # --- iterate ---
  for (i in iters) {
    set.seed(seed + i - 1L)
    keep <- sample(colnames(seurat_obj), size = floor(prop * ncol(seurat_obj)))
    
    # Subset preserves reductions for kept cells
    sobj_sub <- subset(seurat_obj, cells = keep)
    
    # Rebuild graph & clusters on the subset using the specified reduction
    sobj_sub <- FindNeighbors(sobj_sub, reduction = reduction, dims = dims, verbose = FALSE)
    sobj_sub <- FindClusters(sobj_sub, resolution = resolution, verbose = FALSE)
    
    # Compare: original labels vs new labels on the *same* subset cells
    new_lab  <- as.character(Idents(sobj_sub))
    names(new_lab) <- colnames(sobj_sub)
    old_lab  <- orig_clusters[colnames(sobj_sub)]
    
    # ARI (robust to label permutations)
    ari_vals[i] <- mclust::adjustedRandIndex(old_lab, new_lab)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  mean_ari <- mean(ari_vals)
  sd_ari   <- stats::sd(ari_vals)
  se_ari   <- sd_ari / sqrt(n_iter)
  ci95     <- c(lower = mean_ari - 1.96 * se_ari, upper = mean_ari + 1.96 * se_ari)
  
  out <- list(
    mean_ari = mean_ari,
    sd_ari   = sd_ari,
    n_iter   = n_iter,
    prop     = prop,
    resolution = resolution,
    reduction  = reduction,
    dims       = dims,
    ci95     = ci95,
    ari_per_iter = ari_vals
  )
  class(out) <- c("subsetARI", class(out))
  return(out)
}


res <- compute_subset_ari(
  seurat_obj   = libs_UF_all_harmony_integrated,
  reduction    = "harmony",
  dims         = 1:25,
  prop         = 0.75,
  resolution   = 0.5,
  n_iter       = 100,
  seed         = 123
)

saveRDS(res, "~/scratch/20250714_libs1-12B_analysis_v4/ARI_Stability_analysis.rds")

