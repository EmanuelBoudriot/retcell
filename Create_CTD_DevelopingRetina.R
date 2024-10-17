# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
require(scKirby)
require(EWCE)

# Load data
sce <- ingest_data("./data/Chen_Developing_Retina_2023/92f2a709-e108-40b0-9879-49d0081b9eef.rds")

annotLevels <- list(l1 = sce$majorclass)

fNames_counts <- generate_celltype_data(
  exp = sce,
  annotLevels = annotLevels,
  input_species = "human",
  groupName = "DevelopingRetina",
  savePath = "./ctd/DevelopingRetina"
)
