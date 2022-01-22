# R 3.6.1
# Seurat 3.2.3
# tidyverse 1.3.0
# loomR 0.2.1.9000

# Load Gene Lists for Module Scoring ----------------------------------
# Gene lists available in SI Tables 04
# Column 3 (Microglia Myeloid Shared Act Score)
shared_sig <- "Load Microglia Meta Cell Score"
shared_sig <- list(shared_sig)

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(scCustomize)
library(patchwork)
library(viridis)
library(beepr)

# Load Full Data Set & Add Meta Data ------------------------------------------------------
full_list <- Read10X()

full_seurat <- CreateSeuratObject(counts = full_list, names.field = 2, names.delim = "-", min.cells = 5, min.features = 200)

# Update meta data
meta_data <- full_seurat@meta.data %>%
  rownames_to_column("barcodes")

# All meta data is provided on GEO
sample_meta <- read.csv(...)

joined_meta <- left_join(meta_data, sample_meta, by = c("orig.ident" = "cell_ranger_suffix"))

joined_meta <- joined_meta %>%
  mutate(Gene_UMI_Ratio = nFeature_RNA / nCount_RNA)


joined_meta <- joined_meta %>%
  select(-orig.ident, -nFeature_RNA, -nCount_RNA) %>%
  column_to_rownames("barcodes")

# Add meta data to Seurat Object
full_seurat <- AddMetaData(object = full_seurat, metadata = joined_meta)

# Add mito ribo percents
full_seurat <- Add_Mito_Ribo_Seurat(seurat_object = full_seurat, species = "mouse")

# QC Filter and subset
full_seurat_subset <- subset(x = full_seurat, subset = nCount_RNA < 15000 & nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent_mito < 20 & percent_ribo < 20 & percent_mito_ribo < 33)


write_rds(full_seurat_subset, "qc_filtered_seurat.RDS")