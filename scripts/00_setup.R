# Project setup
# Run this file once from the repository root, then run the numbered R scripts
# sequentially in the same R session.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

# Optional local configuration. This file is ignored by git.
if (file.exists(file.path(project_root, "config", "local_paths.R"))) {
  source(file.path(project_root, "config", "local_paths.R"))
}

if (!exists("data_file")) {
  data_file <- Sys.getenv(
    "PANKBASE_RDS",
    unset = file.path(project_root, "data", "060425_scRNA_v3.3.rds")
  )
}

if (!exists("result_dir")) {
  result_dir <- Sys.getenv(
    "PANKBASE_RESULTS",
    unset = file.path(project_root, "results")
  )
}

if (!exists("motif_dir")) {
  motif_dir <- Sys.getenv("SCENIC_MOTIF_DIR", unset = "")
}

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
for (d in c(
  "SCN_donor_barplots_all_celltypes_diabetes_status",
  "SCN_donor_barplots_all_celltypes_HbA1c_group",
  "SCN_pseudobulk_DESeq2_diabetes_group",
  "SCN_targeted_pseudobulk",
  "SCN_UMAP_diabetes_status",
  "SCN_UMAP_HbA1c_group",
  "SCN9A_continuous_GSEA",
  "CytoTRACE2",
  "SCENIC",
  "SCN9A_regulon_continuous_scatter_HPAP",
  "SCN9A_HPAP_top10_regulon_UMAP_violin",
  "Cytoscape_Top10_Regulon_Networks",
  "Beta_cell_Monocle3",
  "monocle2"
)) {
  dir.create(file.path(result_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# Core packages used across the workflow
library(Seurat)
library(data.table)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(openxlsx)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(enrichplot)
library(ComplexHeatmap)
library(igraph)
library(ggraph)
library(SingleCellExperiment)
library(monocle)

if (!requireNamespace("scCancer", quietly = TRUE)) {
  stop(
    "Package 'scCancer' is required for the original cell-type color palette. ",
    "Install it with remotes::install_github('wguo-research/scCancer')."
  )
}

# Additional packages are called with package::function where practical.
# Some modules require optional packages such as CytoTRACE2, GSVA, limma,
# loomR, ggrastr, ggrepel, circlize, and RColorBrewer.

if (!file.exists(data_file)) {
  stop(
    "Input Seurat object not found. Set data_file in config/local_paths.R ",
    "or set the PANKBASE_RDS environment variable. Current path: ", data_file
  )
}

sr <- readRDS(data_file)
sr <- sr[, sr$treatments == "no_treatment"]

# Dataset-specific metadata harmonization used in the original analysis
sr$description_of_diabetes_status[
  sr$samples %in% c("HP-23135-01__Untreated", "HP-22234-01__Untreated")
] <- "non-diabetic"

sr$diabetes_group <- factor(
  sr$description_of_diabetes_status,
  levels = c("non-diabetic", "type 1 diabetes", "type 2 diabetes"),
  labels = c("Control", "T1D", "T2D")
)

sr$hba1c_group <- cut(
  sr$`hba1c_(percentage)`,
  breaks = c(-Inf, 5.7, 6.5, Inf),
  right = FALSE,
  labels = c("Normal", "Prediabetes", "Diabetes")
)

sr$hba1c_group[sr$samples == "HP-22234-01__Untreated"] <- "Prediabetes"
sr$hba1c_group[sr$samples == "HP-23135-01__Untreated"] <- "Normal"
sr$donor_accession[sr$samples == "HP-22234-01__Untreated"] <- "PKBDO5181VXNQ"
sr$donor_accession[sr$samples == "HP-23135-01__Untreated"] <- "PKBDO3790TCBJ"
