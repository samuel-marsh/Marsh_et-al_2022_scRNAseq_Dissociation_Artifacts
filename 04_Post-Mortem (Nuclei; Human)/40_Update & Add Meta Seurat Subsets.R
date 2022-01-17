# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) #v3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# Convert Objects to V3 and Add Sample ID ------------------------------------------------------------
opc_seurat <-  read_rds("RDS_Subcluster_Seurat/opc_rd2_seurat.RDS")
opc_seurat_v3 <- UpdateSeuratObject(opc_seurat)
opc_seurat_v3 <- PercentageFeatureSet(object = opc_seurat_v3, pattern = "^MT", col.name = "percent_mito")
opc_sample_names <- read_rds("RDS_Subcluster_Seurat/opc_sample_names.RDS")
opc_seurat_v3[["sample_id"]] <- opc_sample_names
write_rds(opc_seurat_v3, "RDS_Subcluster_Seurat/V3_Objects/opc_rd2_seurat_V3.RDS")

View(opc_seurat_v3@meta.data)

# micro
micro_seurat <-  read_rds("RDS_Subcluster_Seurat/micro_rd2_seurat.RDS")
micro_seurat_v3 <- UpdateSeuratObject(micro_seurat)
micro_seurat_v3 <- PercentageFeatureSet(object = micro_seurat_v3, pattern = "^MT", col.name = "percent_mito")
micro_sample_names <- read_rds("RDS_Subcluster_Seurat/micro_sample_names.RDS")
micro_seurat_v3[["sample_id"]] <- micro_sample_names
write_rds(micro_seurat_v3, "RDS_Subcluster_Seurat/V3_Objects/micro_rd2_seurat_V3.RDS")

View(micro_seurat_v3@meta.data)

# excit
excit_seurat <-  read_rds("RDS_Subcluster_Seurat/excit_rd3_seurat.RDS")
excit_seurat_v3 <- UpdateSeuratObject(excit_seurat)
excit_seurat_v3 <- PercentageFeatureSet(object = excit_seurat_v3, pattern = "^MT", col.name = "percent_mito")
excit_sample_names <- read_rds("RDS_Subcluster_Seurat/excit_sample_names.RDS")
excit_seurat_v3[["sample_id"]] <- excit_sample_names
write_rds(excit_seurat_v3, "RDS_Subcluster_Seurat/V3_Objects/excit_rd3_seurat_V3.RDS")

View(excit_seurat_v3@meta.data)

# inhib
inhib_seurat <-  read_rds("RDS_Subcluster_Seurat/inhib_rd4_seurat.RDS")
inhib_seurat_v3 <- UpdateSeuratObject(inhib_seurat)
inhib_seurat_v3 <- PercentageFeatureSet(object = inhib_seurat_v3, pattern = "^MT", col.name = "percent_mito")
inhib_sample_names <- read_rds("RDS_Subcluster_Seurat/inhib_sample_names.RDS")
inhib_seurat_v3[["sample_id"]] <- inhib_sample_names
write_rds(inhib_seurat_v3, "RDS_Subcluster_Seurat/V3_Objects/inhib_rd4_seurat_V3.RDS")

View(inhib_seurat_v3@meta.data)

# astro
astro_seurat <-  read_rds("RDS_Subcluster_Seurat/astro_rd3_seurat_B.RDS")
astro_seurat_v3 <- UpdateSeuratObject(astro_seurat)
astro_seurat_v3 <- PercentageFeatureSet(object = astro_seurat_v3, pattern = "^MT", col.name = "percent_mito")
astro_sample_names <- read_rds("RDS_Subcluster_Seurat/astro_sample_names_B.RDS")
astro_seurat_v3[["sample_id"]] <- astro_sample_names
write_rds(astro_seurat_v3, "RDS_Subcluster_Seurat/V3_Objects/astro_rd2_seurat_V3_B.RDS")

View(astro_seurat_v3@meta.data)

# oligo
oligo_seurat <-  read_rds("RDS_Subcluster_Seurat/oligo_rd2_seurat_b.RDS")
oligo_seurat_v3 <- UpdateSeuratObject(oligo_seurat)
oligo_seurat_v3 <- PercentageFeatureSet(object = oligo_seurat_v3, pattern = "^MT", col.name = "percent_mito")
oligo_sample_names <- read_rds("RDS_Subcluster_Seurat/oligo_sample_names_b.RDS")
oligo_seurat_v3[["sample_id"]] <- oligo_sample_names
write_rds(oligo_seurat_v3, "RDS_Subcluster_Seurat/V3_Objects/oligo_rd2_seurat_V3_B.RDS")

View(oligo_seurat_v3@meta.data)