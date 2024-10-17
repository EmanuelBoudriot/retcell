# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
require(scKirby)
require(EWCE)

# Load data
sce <- ingest_data("./data/hrca/snrna/amacrine/5ae56a42-4f42-4a6d-99ba-a65ea100ddde.rds")

# CTD
annotLevels <- list(
  l1 = sce$AC_group1,
  l2 = sce$author_cell_type
)

fNames_counts <- generate_celltype_data(
  exp = sce,
  annotLevels = annotLevels,
  input_species = "human",
  groupName = "hrca_snrna_amacrineonly",
  savePath = "./ctd/hrca_snrna_amacrineonly"
)
