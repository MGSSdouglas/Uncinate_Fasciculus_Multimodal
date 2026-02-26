library(dplyr)
library(tidyr)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(rstatix)
library(boot)
library(parallel)
library(ggpubr)


#based on: https://github.com/MgssGroup/snRNASeq_public/blob/main/5_Finalized_downstream_analysis/Finalized_scripts/1.3_celltype_props_case_control.R

#run in jupyter hub

# === Load Seurat obj and covariates ===
libs_UF_all_harmony_integrated <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/obj_20250815_libs_UF_all_harmony_integrated.rds")

meta_all <- libs_UF_all_harmony_integrated@meta.data

# 2. Count number of cells per sample_id and cluster
cell_counts <- meta_all %>%
  group_by(SubjectID, celltype) %>%
  summarise(n_cells = n(), .groups = "drop")

# 3. Calculate total number of cells per sample
total_counts <- meta_all %>%
  group_by(SubjectID) %>%
  summarise(total_cells = n(), .groups = "drop")

# 4. Merge and compute proportions
cluster_proportions <- cell_counts %>%
  left_join(total_counts, by = "SubjectID") %>%
  dplyr::mutate(Proportion = n_cells / total_cells) %>%
  dplyr::select(SubjectID, celltype, Proportion)

# 5. Pivot wider: 1 row per subject, 1 column per cluster
subject_proportions <- cluster_proportions %>%
  pivot_wider(names_from = celltype, values_from = Proportion, values_fill = 0)

# 6. Add Group
subject_proportions <- subject_proportions %>%
  left_join(meta_all %>% select(SubjectID, Group) %>% distinct(), by = "SubjectID") 

######################################################## GROUP ########################################################

# Set parameters
R <- 10000        # Number of bootstrap resamples
parallel <- "multicore" # Bootstrapping parallelization ("multicore" or "no")
ncpus <- 1       # Number of CPUs for bootstrapping

# Wilcoxon test per cluster
wilcoxon_proportions <- function(cluster_proportions, conditions) {
  myData <- data.frame(Proportion = cluster_proportions, Condition = as.factor(conditions))
  res_test <- wilcox_test(formula = Proportion ~ Condition, data = myData, alternative = "two.sided", paired = FALSE, detailed = TRUE)
  return(res_test)
}

# Bootstrapped Wilcoxon test
w_test <- function(d, i) {
  d <- d[i, ]
  res <- wilcox_test(formula = Proportion ~ Condition, data = d, paired = FALSE, alternative = "two.sided", detailed = TRUE)
  return(res$estimate)
}

# Input: subject_proportions

# Run Wilcoxon tests
print("Running Wilcoxon tests for cell-type proportions...")

results_wilcoxon <- lapply(select(subject_proportions, where(is.numeric)), function(x) {
  wilcoxon_proportions(x, subject_proportions$Group)
}) %>% bind_rows(.id = "Cluster")

# Run bootstrapped Wilcoxon tests
print("Running bootstrapped Wilcoxon tests for cell-type proportions...")

results_wilcoxon_booted <- lapply(select(subject_proportions, where(is.numeric)), function(x) {
  myData <- data.frame(Proportion = x, Condition = subject_proportions$Group)
  booted_res <- boot(myData, statistic = w_test, R = R, parallel = parallel, ncpus = ncpus)
  
  p1 <- 2 * sum(booted_res$t[,1] > 0) / R
  p2 <- 2 * sum(booted_res$t[,1] < 0) / R
  
  p_value <- ifelse(p1 < 1, p1, p2)
  
  tibble(estimate_booted = booted_res$t0, booted_p = p_value)
}) %>% bind_rows(.id = "Cluster")

# Merge Wilcoxon and bootstrapped results
results_wilcoxonsummary <- full_join(results_wilcoxon, results_wilcoxon_booted, by = "Cluster") %>%
  mutate(
    p_adjust = p.adjust(p, method = "BH"),
    booted_p_adjust = p.adjust(booted_p, method = "BH")
  )

# Save results
write.csv(results_wilcoxonsummary, "~/scratch/20250714_libs1-12B_analysis_v4/20250815_results_wilcoxonsummary_celltypeprops.csv", row.names = FALSE)
print("Results saved")

# Plot boxplots per cluster
plots_list <- list()

for (name in names(select(subject_proportions, where(is.numeric)))) {
  print(name)

  myData <- data.frame(
    Proportion = subject_proportions[[name]],
    Condition = as.factor(subject_proportions$Group),
    Subject = subject_proportions$SubjectID
  )

  p <- ggplot(data = myData, aes(x = Proportion, y = Condition, color = Condition)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.6) +
    coord_flip() +
    ggtitle(paste0(name)) +
    theme_classic(base_size = 16)

  plots_list[[name]] <- p
}

#remove mixed clusters
plots_list <- plots_list[-c(21,25)]

# Optionally, save plots
# Example:
# for (n in names(plots_list)) {
#   ggsave(filename = paste0("plot_", n, ".pdf"), plot = plots_list[[n]], width = 8, height = 6)
# }

library(ggpubr)

combined_plot <- ggarrange(
  plotlist = plots_list, ncol = 7,nrow = 4,
  common.legend = TRUE,
  legend = "none")

combined_plot


library(ggpubr)
library(ggplot2)

# Extract the titles from each plot
titles <- sapply(plots_list, function(p) p$labels$title)

# Order the plots by the titles alphabetically
plots_list_ordered <- plots_list[order(titles)]

# Combine in alphabetical order
combined_plot <- ggarrange(
  plotlist = plots_list_ordered,
  ncol = 7, nrow = 4,
  common.legend = TRUE,
  legend = "none"
)

combined_plot


# Save
ggsave("combined_proportion.pdf", combined_plot, width = 20, height = 20)


