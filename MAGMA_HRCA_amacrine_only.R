# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(tidyverse)
library(MAGMA.Celltyping)


# Amacrine only (HRCA) ----------------------------------------------------


# Load CTD data
load("./ctd/hrca_snrna_amacrineonly/ctd_hrca_snrna_amacrineonly.rda")

# Run cell type associations pipeline
MAGMA_results <- celltype_associations_pipeline(
  magma_dirs = "./munged/SCZ/MAGMA_Files/SCZ_Munged.tsv.gz.35UP.10DOWN",
  ctd = ctd,
  ctd_species = "human",
  ctd_name = paste0("hrca_snrna_amacrineonly.", "SCZ"),
  force_new = TRUE,
  run_linear = TRUE,
  save_dir = NULL
)

t1 <- merge_results(MAGMA_results = MAGMA_results, level = 1, species = "human", save_dir = NULL)
t1$GWAS <- "SCZ"

t2 <- merge_results(MAGMA_results = MAGMA_results, level = 2, species = "human", save_dir = NULL)
t2$GWAS <- "SCZ"

write.table(t1, "Merged_results_hrca_snrna_amacrineonly_l1.tsv", quote = F, row.names = F, sep = "\t")
write.table(t2, "Merged_results_hrca_snrna_amacrineonly_l2.tsv", quote = F, row.names = F, sep = "\t")
