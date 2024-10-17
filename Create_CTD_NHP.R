# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(EWCE)
library(tidyverse)

# Load expression matrix
mat <- readRDS("./data/nhp/FullExpressionMatrix.rds")

# Load metadata
meta <- read.table("./data/nhp/Macaque_NN_RGC_AC_BC_HC_PR_metadata_4_expanded.txt.gz", sep = ",", header = T)
meta <- meta %>% mutate(Cluster_Level3 = recode(Cluster_Level3,
  "Amacrine" = "AC",
  "Bipolar" = "BC",
  "Cones" = "Cone",
  "Ganglion Cells" = "RGC",
  "Horizontal" = "HC",
  "Microglia" = "Microglia",
  "Mueller Glia" = "MG",
  "Rods" = "Rod",
  "Vascular" = "Vascular"
))


# Fovea -------------------------------------------------------------------


metafovea <- meta %>% filter(Subcluster == "Fovea")

matfovea <- mat[, metafovea$NAME]

# CTD
annot <- list(l1 = metafovea$Cluster_Level3)

fNames_ALLCELLS <- generate_celltype_data(
  exp = matfovea,
  annotLevels = annot,
  input_species = "Crab-eating macaque",
  method = "gprofiler",
  convert_orths = TRUE,
  groupName = "nhp_fovea",
  savePath = "./ctd/nhp/"
)


# Periphery ---------------------------------------------------------------


metaperiphery <- meta %>% filter(Subcluster == "Periphery")

matperiphery <- mat[, metaperiphery$NAME]

# CTD
annot <- list(l1 = metaperiphery$Cluster_Level3)

fNames_ALLCELLS <- generate_celltype_data(
  exp = matperiphery,
  annotLevels = annot,
  input_species = "Crab-eating macaque",
  method = "gprofiler",
  convert_orths = TRUE,
  groupName = "nhp_periphery",
  savePath = "./ctd/nhp/"
)
