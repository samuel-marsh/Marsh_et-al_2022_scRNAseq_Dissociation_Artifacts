# Environment Setup ---------------------------------------------------------------------------
# R 3.6.1
# Seurat 3.2.3
# tidyverse 1.3.0
# loomR 0.2.1.9000
# Azimuth reference mapping was performed using online portal

# Load Gene Lists for Module Scoring ----------------------------------
# Gene lists available in SI Tables 04
# Column 3 (Microglia Myeloid Shared Act Score)
shared_sig <- "Load Microglia Meta Cell Score"
shared_sig <- list(shared_sig)

# Convert to human gene symbols (homologous gene names previously confirmed so just using simple gene name capitalization conversion)
human_shared_sig <- lapply(shared_sig, function(x){str_to_upper(x)})

# Setup project
Setup_scRNAseq_Project()

# Load Libraries & Data -----------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(scCustomize)
library(marsh.utils)
library(patchwork)
library(viridis)

# Load Data 
control <- Read10X(data.dir = ".../control/filtered_feature_bc_matrix/")

inhib <- Read10X(data.dir = ".../inhib/filtered_feature_bc_matrix/")

merged <- Merge_Sparse_Data_All(matrix_list = c(control, inhib), add_cell_ids = c("control", "inhib"))

pbmc_mock_raw <- CreateSeuratObject(counts = merged, min.cells = 5, min.features = 200)

# Add meta data
pbmc_mock_raw <- Add_Mito_Ribo_Seurat(seurat_object = pbmc_mock_raw, species = "human")

# QC Filter and Subset
pbmc_mock_qc_filtered <- subset(x = pbmc_mock_raw, subset = nFeature_RNA > 400 & nFeature_RNA < 4500 & nCount_RNA < 22500 & percent_mito < 10)