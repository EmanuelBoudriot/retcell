# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(scKirby)
library(EWCE)
library(tidyverse)

# Load data
sce <- ingest_data("./data/mrca/scrna/all/d0183df5-815d-48c2-bcfe-fbf9b716505c.rds")

# CTD
annotLevels <- list(l1 = sce$majorclass)

fNames_counts <- generate_celltype_data(
  exp = sce,
  annotLevels = annotLevels,
  input_species = "mouse",
  groupName = "mrca",
  convert_orths = F,
  savePath = "./ctd/mrca"
)
