## Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(tidyverse)
library(MAGMA.Celltyping)
library(stats)


# Cowan Fovea -------------------------------------------------------------


# Load CTD data
load("./ctd/cowan_2020/ctd_human_fovea_cowan_2020.rda")

# Run cell type associations pipeline
MAGMA_results <- calculate_conditional_celltype_associations(
  ctd = ctd,
  ctd_species = "human",
  magma_dir = "./munged/SCZ/MAGMA_Files/SCZ_Munged.tsv.gz.35UP.10DOWN",
  force_new = T,
  EnrichmentMode = "Top 10%",
  controlledCTs = c("amacrine")
)

t <- MAGMA_results$level1$results %>%
  unique() %>%
  filter(CONTROL == "amacrine")
t$GWAS <- "SCZ"
t$FDR <- stats::p.adjust(t$P, method = "fdr")

# Clean environment
rm(list = (ls()[ls() != "magma_dirs"]))


# Cowan Periphery ---------------------------------------------------------


# Load CTD data
load("./ctd/cowan_2020/ctd_human_periphery_cowan_2020.rda")

# Run cell type associations pipeline
MAGMA_results <- calculate_conditional_celltype_associations(
  ctd = ctd,
  ctd_species = "human",
  magma_dir = "./munged/SCZ/MAGMA_Files/SCZ_Munged.tsv.gz.35UP.10DOWN",
  force_new = T,
  EnrichmentMode = "Top 10%",
  controlledCTs = c("amacrine")
)

t <- MAGMA_results$level1$results %>%
  unique() %>%
  filter(CONTROL == "amacrine")
t$GWAS <- "SCZ"
t$FDR <- stats::p.adjust(t$P, method = "fdr")

write.table(t, "Merged_results_human_periphery_cowan_2020_conditional.tsv", quote = F, row.names = F, sep = "\t")
