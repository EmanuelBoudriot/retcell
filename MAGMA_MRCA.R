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


# MRCA --------------------------------------------------------------------


# Load CTD data
load("./ctd/mrca/ctd_mrca.rda")

# Run cell type associations pipeline

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
    ctd_species = "mouse",
    ctd_name = paste0("ctd_mrca.", magma_dirs$name[i]),
    force_new = T,
    run_linear = TRUE,
    save_dir = NULL
  )

  t <- merge_results(MAGMA_results = MAGMA_results, level = 1, species = "mouse", save_dir = NULL)
  t$GWAS <- magma_dirs$name[i]
  merged_results <- rbind(merged_results, t)
}

write.table(merged_results, "Merged_results_MRCA.tsv", quote = F, row.names = F, sep = "\t")
