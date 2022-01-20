# All literature reanalysis was performed using (unless otherwise specified):
    # R 3.6.1
    # Seurat 3.1.5
    # tidyverse 1.3.0
    # loomR 0.2.1.9000

# 1.0 Load Packages & Scripts ---------------------------------------------
library(tidyverse)
library(readxl)
library(Seurat)
library(viridis)
library(beepr)
library(scCustomize)
library(patchwork)
library(loomR)
library(future) # Not necessary but speeds but larger dataset analysis
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 5000 * 1024^2)

# 2.0 Load Gene Lists for Module Scoring ----------------------------------
    # Gene lists available in SI Tables 04
# Column 3 (Microglia Myeloid Shared Act Score)
shared_sig <- "Load Microglia Meta Cell Score"
shared_sig <- list(shared_sig)

# Column 1 (Microglia Identity Score)
homeostatic_mg <- "Homeostatic microglia gene list"
homeostatic_mg <- list(homeostatic_mg)

# 2.1 Create gene lists compatible with Mizrak et al dataset
    # Mizrak et al dataset are in form of ensembl ID #s and not external gene names 
    # gene IDs also have suffix included

# Create compatible gene list for Mizrak module scoring
gene_of_interest <- gene_labels_mizrak %>% 
    filter(str_detect(gene_name, 'gene_of_interest')) %>% 
    pull(ensembl_ID)
        # see script 05_Mizrak et al Reanalysis for creation of 'gene_labels_mizrak' object

shared_sig_ensembl <- "Converted gene names to ensembl"
homeostatic_mg_ensembl <- "Converted gene names to ensembl"