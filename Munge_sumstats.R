# Clear Workspace
rm(list = ls())

# setwd("...")

# Load packages
library(tidyverse)
library(MungeSumstats)


# Schizophrenia -----------------------------------------------------------


# Import file
sumstats <- read.table("./data/pgc_scz_2022/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz", sep = "\t", header = T)
sumstats$Neff <- sumstats$NEFFDIV2 * 2

# Run reformatting pipeline
sumstats_reformatted <- format_sumstats(path = sumstats, ref_genome = "GRCh37", rmv_chr = c("X", "Y", "MT"), save_format = "LDSC", compute_z = TRUE, sort_coordinates = TRUE)


# Bipolar Disorder --------------------------------------------------------


# Import file
sumstats <- read.table("./data/pgc_bd/daner_bip_pgc3_nm.gz", sep = "\t", header = TRUE)
sumstats$Neff <- sumstats$Neff_half * 2

# Run reformatting pipeline
sumstats_reformatted <- format_sumstats(path = sumstats, ref_genome = "GRCh37", rmv_chr = c("X", "Y", "MT"), save_format = "LDSC", compute_z = TRUE, sort_coordinates = TRUE, impute_beta = TRUE)


# Alzheimer's -------------------------------------------------------------


# Import file
sumstats <- read.table("./data/pgc_alz/PGCALZ2sumstatsExcluding23andMe.txt.gz", sep = "\t", header = T)

# Run reformatting pipeline
sumstats_reformatted <- format_sumstats(path = sumstats, ref_genome = "GRCh37", rmv_chr = c("X", "Y", "MT"), save_format = "LDSC", compute_z = TRUE, sort_coordinates = TRUE)


# MDD ---------------------------------------------------------------------


# Import file
sumstats <- read.table("./data/pgc_mdd_2023/mdd2023diverse_AllAncestry_Neff.csv", sep = ",", header = T)

# Run reformatting pipeline
sumstats_reformatted <- format_sumstats(path = sumstats, ref_genome = "GRCh37", rmv_chr = c("X", "Y", "MT"), save_format = "LDSC", impute_beta = T, compute_z = TRUE, sort_coordinates = TRUE)


# MS ----------------------------------------------------------------------


# Import file
sumstats <- read.table("./data/ms_2019/discovery_metav3.0.meta.gz", sep = " ", header = T)
sumstats$N <- 41505

# Run reformatting pipeline
sumstats_reformatted <- format_sumstats(path = sumstats, ref_genome = NULL, rmv_chr = c("X", "Y", "MT"), save_format = "LDSC", impute_beta = T, compute_z = TRUE, sort_coordinates = TRUE)


# Parkinson's -------------------------------------------------------------


# Import file
sumstats <- read.table("./data/pd_2024/GCST90275127.tsv", sep = "\t", header = T)
# Summary statistics include 44,358 cases, 18,618 proxy cases, and 966,017 controls (take total N as proxy)
sumstats$N <- 1028993
sumstats <- sumstats %>% select(!c(9:12, 15, 16))

# Run reformatting pipeline
sumstats_reformatted <- format_sumstats(path = sumstats, ref_genome = "GRCh37", rmv_chr = c("X", "Y", "MT"), save_format = "LDSC", compute_z = TRUE, sort_coordinates = TRUE)


# Stroke ------------------------------------------------------------------


# Import file
sumstats <- read.table("./data/stroke_2022/GCST90104534_buildGRCh37.tsv.gz", sep = "\t", header = T)
sumstats$N <- 1614080

# Run reformatting pipeline
sumstats_reformatted <- format_sumstats(path = sumstats, ref_genome = "GRCh37", rmv_chr = c("X", "Y", "MT"), save_format = "LDSC", compute_z = TRUE, sort_coordinates = TRUE)
