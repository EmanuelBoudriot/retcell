# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(BPCells)
library(Azimuth)
library(Seurat)
library(EWCE)

# Load data
# parse.data <- open_matrix_anndata_hdf5(
#   "./data/hrca/snrna/all/89f6a640-0537-4fd1-bdf9-540db9dd0b7d.h5ad"
# )
# write_matrix_dir(mat = parse.data, dir = "./bpcells/hrca")
hrca.mat <- open_matrix_dir(dir = "./bpcells/hrca")
metadata <- LoadH5ADobs("./data/hrca/snrna/all/89f6a640-0537-4fd1-bdf9-540db9dd0b7d.h5ad")
hrca.seurat <- CreateSeuratObject(counts = hrca.mat, meta.data = metadata)

# Downsample per majorclass
Idents(hrca.seurat) <- "majorclass"
hrca.seurat.downsampled <- subset(x = hrca.seurat, downsample = 5000)

# Convert to SCE
sce <- as.SingleCellExperiment(hrca.seurat.downsampled)

# CTD
annotLevels <- list(l1 = sce$majorclass)

fNames_counts <- generate_celltype_data(
  exp = sce,
  annotLevels = annotLevels,
  input_species = "human",
  groupName = "hrca_snrna_subset_majorclass_5000",
  savePath = "./ctd/hrca_snrna_subset_majorclass_5000"
)
