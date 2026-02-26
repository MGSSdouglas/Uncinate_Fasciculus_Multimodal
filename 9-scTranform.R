library(Seurat)
library(sctransform)
library(dplyr)
library(glmGamPoi)

#Parallelization for sctransform - doesn't work for harmony - results in errors
library(future)
plan("multicore", workers = 10)

options(future.globals.maxSize = 100 * 1024^3)


libs_UF_all <- readRDS("~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_all.rds")

#filter out S251 as not enough cells (n=15 cells)
libs_UF_all <- subset(x = libs_UF_all, subset = AnyDoublet_Individual_Assignment != "R2_251")


#mitochondrial
libs_UF_all[["percent.mt"]] <- PercentageFeatureSet(libs_UF_all, pattern = "^MT-")
libs_UF_all <- subset(libs_UF_all, subset = percent.mt < 5)

# Filter out unwanted genes from the RNA assay (raw counts)
keep_genes <- setdiff(rownames(libs_UF_all), "MALAT1")
libs_UF_all <- subset(libs_UF_all, features = keep_genes)


libs_UF_all <- libs_UF_all %>%
  SCTransform(vars.to.regress = "percent.mt", return.only.var.genes = FALSE, verbose = TRUE)


saveRDS(libs_UF_all, "~/scratch/20250714_libs1-12B_analysis_v4/libs_UF_all_sctransformed.rds")
