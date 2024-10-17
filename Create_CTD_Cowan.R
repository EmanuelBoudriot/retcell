# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(tidyverse)
library(EWCE)
library(scKirby)


# Fovea -------------------------------------------------------------------


sce <- ingest_data("./data/cowan_2020/fovea/74a8e0fe-08ee-4a93-9c2f-9cb70acd714b.rds")

# CTD
annot <- list(l1 = sce$cell_type_group)

fNames_ALLCELLS <- generate_celltype_data(
  exp = sce,
  annotLevels = annot,
  input_species = "human",
  groupName = "human_fovea_cowan_2020",
  savePath = "./ctd/cowan_2020/"
)


# Periphery ---------------------------------------------------------------


sce <- ingest_data("./data/cowan_2020/periphery/20500dc8-c26d-4bd7-af60-8a0dcbe14656.rds")

# CTD
annot <- list(l1 = sce$cell_type_group)

fNames_ALLCELLS <- generate_celltype_data(
  exp = sce,
  annotLevels = annot,
  input_species = "human",
  groupName = "human_periphery_cowan_2020",
  savePath = "./ctd/cowan_2020/"
)
