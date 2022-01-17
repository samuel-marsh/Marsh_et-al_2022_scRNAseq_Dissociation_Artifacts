# All analysis performed with R 3.6.1
# Seurat version specified in code

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Starting with Seurat V3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(scCustomize)
library(beepr)

install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_3.1.5.tar.gz", repos = NULL, type = "source")

# Load Data ---------------------------------------------------------------
# Specify file paths
# Sample donorA PM
pm_donorA_01_path <- "~/Desktop/Transfer_to_O2/exp39_human_nucseq/pEXP27sHSrCTXPMidonorAd20191001dapi/filtered_feature_bc_matrix/"

pm_donorA_02_path <- "~/Desktop/Transfer_to_O2/exp39_human_nucseq/pEXP35sHSrCTXPMidonorAd20191001dapi2/filtered_feature_bc_matrix/"
# Sample donorB PM
pm_donorB_01_path <- "~/Desktop/Transfer_to_O2/exp39_human_nucseq/pEXP27sHSrCTXPMidonorBd20191001dapi/filtered_feature_bc_matrix/"

pm_donorB_02_path <- "~/Desktop/Transfer_to_O2/exp39_human_nucseq/pEXP35sHSrCTXPMidonorBd20191001dapi2/filtered_feature_bc_matrix/"

# Sample donorC PM
pm_donorC_01_path <- "~/Desktop/Transfer_to_O2/exp39_human_nucseq/pEXP27sHSrCTXPMidonorCd20191001dapi/filtered_feature_bc_matrix/"

pm_donorC_02_path <- "~/Desktop/Transfer_to_O2/exp39_human_nucseq/pEXP35sHSrCTXPMidonorCd20191001dapi2/filtered_feature_bc_matrix/"

# Read libraries into R
pmdonorA_01 <- Read10X(data.dir = pm_donorA_01_path, strip.suffix = FALSE)
pmdonorA_02 <- Read10X(data.dir = pm_donorA_02_path, strip.suffix = FALSE)

pmdonorB_01 <- Read10X(data.dir = pm_donorB_01_path, strip.suffix = FALSE)
pmdonorB_02 <- Read10X(data.dir = pm_donorB_02_path, strip.suffix = FALSE)

pmdonorC_01 <- Read10X(data.dir = pm_donorC_01_path, strip.suffix = FALSE)
pmdonorC_02 <- Read10X(data.dir = pm_donorC_02_path, strip.suffix = FALSE)

# Create individual Seurat Objects
pmdonorA_01 <- CreateSeuratObject(counts = pmdonorA_01, min.cells = 5, min.features = 200, project = "pmdonorA_01")

pmdonorA_02 <- CreateSeuratObject(counts = pmdonorA_02, min.cells = 5, min.features = 200, project = "pmdonorA_02")

pmdonorB_01 <- CreateSeuratObject(counts = pmdonorB_01, min.cells = 5, min.features = 200, project = "pmdonorB_01")

pmdonorB_02 <- CreateSeuratObject(counts = pmdonorB_02, min.cells = 5, min.features = 200, project = "pmdonorB_02")

pmdonorC_01 <- CreateSeuratObject(counts = pmdonorC_01, min.cells = 5, min.features = 200, project = "pmdonorC_01")

pmdonorC_02 <- CreateSeuratObject(counts = pmdonorC_02, min.cells = 5, min.features = 200, project = "pmdonorC_02")

# Merge Seurat objects into single object
marsh_post_seurat <- merge(x = pmdonorA_01, y = c(pmdonorA_02, pmdonorB_01, pmdonorB_02, pmdonorC_01, pmdonorC_02))

# Create combination column in metadata to combine libraries from the same donor
marsh_post_seurat@meta.data$sample_id[marsh_post_seurat@meta.data$orig.ident == "pmdonorA_01" | marsh_post_seurat@meta.data$orig.ident == "pmdonorA_02"] <- "pmdonorA"

marsh_post_seurat@meta.data$sample_id[marsh_post_seurat@meta.data$orig.ident == "pmdonorB_01" | marsh_post_seurat@meta.data$orig.ident == "pmdonorB_02"] <- "pmdonorB"

marsh_post_seurat@meta.data$sample_id[marsh_post_seurat@meta.data$orig.ident == "pmdonorC_01" | marsh_post_seurat@meta.data$orig.ident == "pmdonorC_02"] <- "pmdonorC"

# Add Mito %
marsh_post_seurat <- PercentageFeatureSet(object = marsh_post_seurat, pattern = "^MT", col.name = "percent_mito")

# QC ----------------------------------------------------------------------
# QC Data
QC_Plots_Genes(marsh_post_seurat, pt.size = 0.1, low_cutoff = 600)
QC_Plots_UMI(marsh_post_seurat, pt.size = 0.1, high_cutoff = 75000)
QC_Plots_Mito(marsh_post_seurat, high_cutoff = 20)

# Filter
marsh_post_seurat <- subset(x = marsh_post_seurat, subset = nFeature_RNA > 600 & nCount_RNA < 75000 & percent_mito < 20)

# Convert to Liger & Save -------------------------------------------------
marsh_post_seurat <- seuratToLiger(objects = marsh_post_seurat, combined.seurat = TRUE, meta.var = "sample_id", remove.missing = FALSE)

# Save to RDS
write_rds(marsh_post_seurat, "RDS_Objects/marsh_post_mortem_seuratV3_filtered_liger.RDS")