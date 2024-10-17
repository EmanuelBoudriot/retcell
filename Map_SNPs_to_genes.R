# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(tidyverse)
library(MAGMA.Celltyping)

# Map SNPs to genes
traits <- c("AD", "BD", "MDD", "MS", "PD", "SCZ", "Stroke")
genesOutPath <- c()
for (trait in traits) {
  # Define input path
  inpath <- paste0("./munged/", trait, "/", trait, "_Munged.tsv.gz")

  # Do the mapping
  outs <- map_snps_to_genes(
    path_formatted = inpath,
    genome_build = "GRCh37",
    storage_dir = "."
  )

  # Add output path to the list
  genesOutPath <- c(genesOutPath, outs)
}
