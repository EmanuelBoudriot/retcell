# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(tidyverse)
library(MAGMA.Celltyping)

# MAGMA file reference
magma_dirs <- data.frame(
  dir = c(
    "./munged/SCZ/MAGMA_Files/SCZ_Munged.tsv.gz.35UP.10DOWN",
    "./munged/BD/MAGMA_Files/BD_Munged.tsv.gz.35UP.10DOWN",
    "./munged/MDD/MAGMA_Files/MDD_Munged.tsv.gz.35UP.10DOWN",
    "./munged/PD/MAGMA_Files/PD_Munged.tsv.gz.35UP.10DOWN",
    "./munged/AD/MAGMA_Files/AD_Munged.tsv.gz.35UP.10DOWN",
    "./munged/MS/MAGMA_Files/MS_Munged.tsv.gz.35UP.10DOWN",
    "./munged/Stroke/MAGMA_Files/Stroke_Munged.tsv.gz.35UP.10DOWN"
  ),
  name = c("SCZ", "BD", "MDD", "PD", "AD", "MS", "Stroke"),
  stringsAsFactors = FALSE
)


# Developing Retina (Chen 2023) -------------------------------------------


# Load CTD data
load("./ctd/DevelopingRetina/ctd_DevelopingRetina.rda")

# Run celltype associations pipeline

### Loop

merged_results <- c()

for (i in 1:7) {
  print("#########################")
  print(magma_dirs$name[i])
  print("#########################")
  print("#########################")

  MAGMA_results <- celltype_associations_pipeline(
    magma_dirs = magma_dirs$dir[i],
    ctd = ctd,
    ctd_species = "human",
    ctd_name = paste0("ctd_DevelopingRetina.", magma_dirs$name[i]),
    run_linear = TRUE,
    save_dir = NULL
  )

  t <- merge_results(MAGMA_results = MAGMA_results, level = 1, species = "human", save_dir = NULL)
  t$GWAS <- magma_dirs$name[i]
  merged_results <- rbind(merged_results, t)
}

write.table(merged_results, "Merged_results_DevelopingRetina.tsv", quote = F, row.names = F, sep = "\t")
