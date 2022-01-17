# Start with Seurat v3.1.5
# Analysis performed with R3.5.3

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_3.1.5.tar.gz", repos = NULL, type = "source")
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")

# Load Data
leng_sfg_all <- Read10X("~/Desktop/Literature Datasets/Human DataSets/Post-Mortem/Leng 2020 Biorxiv (Ctrl & AD)/SFG_ALL/")

# Create object
leng_sfg_all <- CreateSeuratObject(counts = leng_sfg_all, min.cells = 5, min.features = 200)
View(leng_sfg_all@meta.data)

# Add mito%
leng_sfg_all <- PercentageFeatureSet(object = leng_sfg_all, pattern = "^MT", col.name = "percent_mito")
write_rds(leng_sfg_all, "RDS_Objects/leng_sfg_full_seurat.RDS")

beep(sound = 2)

# QC Filtering ------------------------------------------------------------
# SFG Dataset
stats_sfg <- Cluster_Stats_All_Samples(leng_sfg_all)

QC_Plots_Genes(leng_sfg_all)
QC_Plots_UMI(leng_sfg_all, high_cutoff = 40000)
QC_Plots_Mito(leng_sfg_all, high_cutoff = 30)

leng_sfg_all <- subset(x = leng_sfg_all, subset = nCount_RNA < 40000)
leng_sfg_all <- subset(x = leng_sfg_all, subset = percent_mito < 30 & nFeature_RNA < 9000)
write_rds(leng_sfg_all, "RDS_Objects/leng_sfg_filtered_seurat.RDS")

# Convert to Liger
leng_sfg_all_liger <- seuratToLiger(objects = leng_sfg_all, combined.seurat = TRUE, meta.var = "orig.ident", remove.missing = FALSE)
write_rds(leng_sfg_all_liger, "RDS_Objects/leng_sfg_liger_RAW.RDS")

beep(sound = 2)

# Restart R session and Install Seurat v2.3.4 -----------------------------
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

# Round01 Clustering ------------------------------------------------------
# Load liger objects
leng_sfg_liger <- read_rds("RDS_Objects/leng_sfg_liger_RAW.RDS")

leng_sfg_liger <- liger::normalize(leng_sfg_liger)
leng_sfg_liger <- selectGenes(leng_sfg_liger, do.plot = T, num.genes = 400)
leng_sfg_liger <- scaleNotCenter(leng_sfg_liger)
leng_sfg_liger <- online_iNMF(leng_sfg_liger, k = 40, lambda = 5, miniBatch_size = 2000)
leng_sfg_liger <- quantile_norm(leng_sfg_liger, knn_k = 15)
leng_sfg_liger <- clusterLouvainJaccard(leng_sfg_liger,resolution = 0.7)
leng_sfg_liger <- runUMAP(leng_sfg_liger)

umap_dim <- plotByDatasetAndCluster(leng_sfg_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Save round01 object
write_rds(leng_sfg_liger, "RDS_Objects/leng_sfg_liger_Round01_clustered.RDS")

# Annotation --------------------------------------------------------------
leng_sfg_annotation <- tibble::tribble(
  ~cluster, ~cell_type,        ~color,
  0L,    "oligo",      "orange",
  1L,    "astro", "forestgreen",
  2L,    "oligo",      "orange",
  3L,    "excit",  "dodgerblue",
  4L,    "inhib",        "navy",
  5L,    "oligo",      "orange",
  6L,    "micro",        "gold",
  7L,      "OPC",  "darkorange",
  8L,    "excit",  "dodgerblue",
  9L,    "inhib",        "navy",
  10L,    "astro", "forestgreen",
  11L,    "oligo",      "orange",
  12L,    "excit",  "dodgerblue",
  13L,    "excit",  "dodgerblue",
  14L,    "excit",  "dodgerblue",
  15L,    "excit",  "dodgerblue",
  16L,    "excit",  "dodgerblue",
  17L,    "excit",  "dodgerblue",
  18L,    "excit",  "dodgerblue",
  19L,     "fibro",       "pink",
  20L,     "endo", "darkorchid3",
  21L,    "excit",  "dodgerblue",
  22L,    "excit",  "dodgerblue",
  23L,    "oligo",      "orange",
  24L,    "astro", "forestgreen",
  25L,    "excit",  "dodgerblue"
)

# Colors round01 sfg
round01_cluster_colors <- leng_sfg_annotation %>% 
  pull(color)

# Pull Clusters for Cleaning ----------------------------------------------
# Pull cell type cluster IDs
oligo_clu <- leng_sfg_annotation %>% 
  filter(cell_type == "oligo") %>% 
  pull(cluster)

excit_clu <- leng_sfg_annotation %>% 
  filter(cell_type == "excit") %>% 
  pull(cluster)

inhib_clu <- leng_sfg_annotation %>% 
  filter(cell_type == "inhib") %>% 
  pull(cluster)

micro_clu <- leng_sfg_annotation %>% 
  filter(cell_type == "micro") %>% 
  pull(cluster)

astro_clu <- leng_sfg_annotation %>% 
  filter(cell_type == "astro") %>% 
  pull(cluster)

opc_clu <- leng_sfg_annotation %>% 
  filter(cell_type == "OPC") %>% 
  pull(cluster)

endo_clu <- leng_sfg_annotation %>% 
  filter(cell_type == "endo") %>% 
  pull(cluster)

fibro_clu <- leng_sfg_annotation %>% 
  filter(cell_type == "fibro") %>% 
  pull(cluster)

# Microglia ---------------------------------------------------------------
micro_sub <- subsetLiger(object = leng_sfg_liger, clusters.use = micro_clu, remove.missing = FALSE)

micro_sub <- normalize(micro_sub)
micro_sub <- selectGenes(micro_sub, do.plot = T, num.genes = 200)
micro_sub <- scaleNotCenter(micro_sub)
micro_sub <- optimizeALS(micro_sub, k = 25, lambda = 5, nrep = 3)
micro_sub <- quantile_norm(micro_sub, knn_k = 15)
micro_sub <- clusterLouvainJaccard(micro_sub,resolution = 3)
micro_sub <- runUMAP(micro_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(micro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
# Really had to bump resolution to get some clusters to subdivde for removal
micro_doubl <- c("1", "20", "30")
doublet_micro <- data.frame(micro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_micro)[1] <- "cluster"
doublet_micro <- doublet_micro %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% micro_doubl) %>%
  pull(barcodes)

# Oligo ---------------------------------------------------------------
oligo_sub <- subsetLiger(object = leng_sfg_liger, clusters.use = oligo_clu, remove.missing = FALSE)

oligo_sub <- normalize(oligo_sub)
oligo_sub <- selectGenes(oligo_sub, do.plot = T, num.genes = 300)
oligo_sub <- scaleNotCenter(oligo_sub)
oligo_sub <- online_iNMF(oligo_sub, k = 25, lambda = 5, miniBatch_size = 200)
oligo_sub <- quantile_norm(oligo_sub, knn_k = 15)
oligo_sub <- clusterLouvainJaccard(oligo_sub,resolution = 1)
oligo_sub <- runUMAP(oligo_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(oligo_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
oligo_doubl <- c("6")
doublet_oligo <- data.frame(oligo_sub@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)

# Excit---------------------------------------------------------------
excit_sub <- subsetLiger(object = leng_sfg_liger, clusters.use = excit_clu, remove.missing = FALSE)
excit_sub <- normalize(excit_sub)
excit_sub <- selectGenes(excit_sub, do.plot = T, num.genes = 200)
excit_sub <- scaleNotCenter(excit_sub)
excit_sub <- online_iNMF(excit_sub, k = 25, lambda = 5, miniBatch_size = 700)
excit_sub <- quantile_norm(excit_sub, knn_k = 15)
excit_sub <- clusterLouvainJaccard(excit_sub,resolution = 0.6)
excit_sub <- runUMAP(excit_sub, min_dist = 0.3, n_neighbors = 40)

umap_dim <- plotByDatasetAndCluster(excit_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# # Elim Doublets
excit_doubl <- c("19", "7", "20")
doublet_excit <- data.frame(excit_sub@clusters, stringsAsFactors = FALSE)
names(doublet_excit)[1] <- "cluster"
doublet_excit <- doublet_excit %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% excit_doubl) %>%
  pull(barcodes)

# inhib ---------------------------------------------------------------
inhib_sub <- subsetLiger(object = leng_sfg_liger, clusters.use = inhib_clu, remove.missing = FALSE)
inhib_sub <- normalize(inhib_sub)
inhib_sub <- selectGenes(inhib_sub, do.plot = T, num.genes = 200)
inhib_sub <- scaleNotCenter(inhib_sub)
inhib_sub <- online_iNMF(inhib_sub, k = 25, lambda = 5, miniBatch_size = 250)
inhib_sub <- quantile_norm(inhib_sub, knn_k = 15)
inhib_sub <- clusterLouvainJaccard(inhib_sub,resolution = 1.2)
inhib_sub <- runUMAP(inhib_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(inhib_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
inhib_doubl <- c("13", "12", "15")
doublet_inhib <- data.frame(inhib_sub@clusters, stringsAsFactors = FALSE)
names(doublet_inhib)[1] <- "cluster"
doublet_inhib <- doublet_inhib %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% inhib_doubl) %>%
  pull(barcodes)

# Astro ---------------------------------------------------------------
astro_sub <- subsetLiger(object = leng_sfg_liger, clusters.use = astro_clu, remove.missing = FALSE)

astro_sub <- normalize(astro_sub)
astro_sub <- selectGenes(astro_sub, do.plot = T, num.genes = 200)
astro_sub <- scaleNotCenter(astro_sub)
astro_sub <- online_iNMF(astro_sub, k = 25, lambda = 5, miniBatch_size = 250)
astro_sub <- quantile_norm(astro_sub)
astro_sub <- clusterLouvainJaccard(astro_sub,resolution = 1.2)
astro_sub <- runUMAP(astro_sub, min_dist = 0.3, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(astro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
astro_doubl <- c("8", "7")
doublet_astro <- data.frame(astro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_astro)[1] <- "cluster"
doublet_astro <- doublet_astro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_doubl) %>% 
  pull(barcodes)

# OPC ---------------------------------------------------------------
opc_sub <- subsetLiger(object = leng_sfg_liger, clusters.use = opc_clu, remove.missing = FALSE)
opc_sub <- normalize(opc_sub)
opc_sub <- selectGenes(opc_sub, do.plot = T, num.genes = 200)
opc_sub <- scaleNotCenter(opc_sub)
opc_sub <- online_iNMF(opc_sub, k = 25, lambda = 5, miniBatch_size = 100)
opc_sub <- quantile_norm(opc_sub, knn_k = 15)
opc_sub <- clusterLouvainJaccard(opc_sub,resolution = 0.7)
opc_sub <- runUMAP(opc_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(opc_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
opc_doubl <- c("5")
doublet_opc <- data.frame(opc_sub@clusters, stringsAsFactors = FALSE)
names(doublet_opc)[1] <- "cluster"
doublet_opc <- doublet_opc %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% opc_doubl) %>% 
  pull(barcodes)

# Endo ---------------------------------------------------------------
endo_sub <- subsetLiger(object = leng_sfg_liger, clusters.use = endo_clu, remove.missing = FALSE)
endo_sub <- normalize(endo_sub)
endo_sub <- selectGenes(endo_sub, do.plot = T, num.genes = 100)
endo_sub <- scaleNotCenter(endo_sub)
endo_sub <- optimizeALS(endo_sub, k = 2, lambda = 5, nrep = 3)
endo_sub <- quantile_norm(endo_sub)
endo_sub <- clusterLouvainJaccard(endo_sub,resolution = 1.5)
endo_sub <- runUMAP(endo_sub)

umap_dim <- plotByDatasetAndCluster(endo_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
endo_doubl <- c("4")
doublet_endo <- data.frame(endo_sub@clusters, stringsAsFactors = FALSE)
names(doublet_endo)[1] <- "cluster"
doublet_endo <- doublet_endo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% endo_doubl) %>% 
  pull(barcodes)

# fibro ---------------------------------------------------------------
fibro_sub <- subsetLiger(object = leng_sfg_liger, clusters.use = fibro_clu, remove.missing = FALSE)

fibro_sub <- normalize(fibro_sub)
fibro_sub <- selectGenes(fibro_sub, do.plot = T, num.genes = 100)
fibro_sub <- scaleNotCenter(fibro_sub)
fibro_sub <- optimizeALS(fibro_sub, k = 5, lambda = 5, nrep = 3)
fibro_sub <- quantile_norm(fibro_sub, k=8)
fibro_sub <- clusterLouvainJaccard(fibro_sub,resolution = 1.5)
fibro_sub <- runUMAP(fibro_sub)

umap_dim <- plotByDatasetAndCluster(fibro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
fibro_doubl <- c("0", "1", "5", "6", "7", "2")
doublet_fibro <- data.frame(fibro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_fibro)[1] <- "cluster"
doublet_fibro <- doublet_fibro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% fibro_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete <- c(doublet_astro, doublet_endo, doublet_excit, doublet_fibro, doublet_inhib, doublet_micro, doublet_oligo, doublet_opc)

# Remove doublets
# Add barcodes to cell data
leng_sfg_liger@cell.data$barcodes <- rownames(leng_sfg_liger@cell.data)

# Remove doublet cells
leng_sfg_liger_cleaned <- subsetLiger(object = leng_sfg_liger, cells.use = setdiff(unique(leng_sfg_liger@cell.data$barcodes), doublets_complete), remove.missing = FALSE)

# Save object
write_rds(leng_sfg_liger_cleaned, "RDS_Objects/Leng_SFG_Cleaned_Round01.RDS")


# Round02 Clustering ------------------------------------------------------
# Load liger objects
leng_sfg_liger_cleaned <- liger::normalize(leng_sfg_liger_cleaned)
leng_sfg_liger_cleaned <- selectGenes(leng_sfg_liger_cleaned, do.plot = T, num.genes = 350)
leng_sfg_liger_cleaned <- scaleNotCenter(leng_sfg_liger_cleaned)
leng_sfg_liger_cleaned <- online_iNMF(leng_sfg_liger_cleaned, k = 33, lambda = 5, miniBatch_size = 2000)
leng_sfg_liger_cleaned <- quantile_norm(leng_sfg_liger_cleaned, knn_k = 15)
leng_sfg_liger_cleaned <- clusterLouvainJaccard(leng_sfg_liger_cleaned,resolution = 0.4)
leng_sfg_liger_cleaned <- runUMAP(leng_sfg_liger_cleaned, n_neighbors = 40)

umap_dim <- plotByDatasetAndCluster(leng_sfg_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Annotation --------------------------------------------------------------
leng_sfg_annotation_rd2 <- tibble::tribble(
  ~cluster, ~cell_type,        ~color,
  0L,    "excit",  "dodgerblue",
  1L,    "oligo",      "orange",
  2L,    "excit",  "dodgerblue",
  3L,    "oligo",      "orange",
  4L,    "astro", "forestgreen",
  5L,    "inhib",        "navy",
  6L,    "micro",        "gold",
  7L,    "oligo",      "orange",
  8L,    "astro", "forestgreen",
  9L,      "OPC",  "darkorange",
  10L,    "excit",  "dodgerblue",
  11L,    "inhib",        "navy",
  12L,    "excit",  "dodgerblue",
  13L,    "oligo",      "orange",
  14L,      "OPC",  "darkorange",
  15L,    "oligo",      "orange",
  16L,     "endo",        "pink",
  17L,    "oligo",      "orange",
  18L,    "excit",  "dodgerblue",
  19L,    "excit",  "dodgerblue"
)

# Plot recolored
# pull colors
sfg_cluster_colors <- leng_sfg_annotation_rd2 %>% 
  pull(color)


# 2nd Round Clean -------------------------------------------------
# Final mop up.  Looks to be mostly excitatory neurons that have other cell doublets
excit_clu <- leng_sfg_annotation_rd2 %>% 
  filter(cell_type == "excit") %>% 
  pull(cluster)

excit_sub_rd2 <- subsetLiger(object = leng_sfg_liger_rd02, clusters.use = excit_clu, remove.missing = FALSE)

excit_sub_rd2 <- normalize(excit_sub_rd2)
excit_sub_rd2 <- selectGenes(excit_sub_rd2, do.plot = T, num.genes = 150)
excit_sub_rd2 <- scaleNotCenter(excit_sub_rd2)
excit_sub_rd2 <- online_iNMF(excit_sub_rd2, k = 25, lambda = 5, miniBatch_size = 250)
#excit_sub_rd2 <- optimizeALS(excit_sub_rd2, k = 15, lambda = 1, nrep = 1)
excit_sub_rd2 <- quantile_norm(excit_sub_rd2, knn_k = 15)
excit_sub_rd2 <- clusterLouvainJaccard(excit_sub_rd2,resolution = 1)
excit_sub_rd2 <- runUMAP(excit_sub_rd2, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(excit_sub_rd2, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Extract doublet barcodes
excit_doubl_rd2 <- c("17", "18", "19")
doublet_excit_rd2 <- data.frame(excit_sub_rd2@clusters, stringsAsFactors = FALSE)
names(doublet_excit_rd2)[1] <- "cluster"
doublet_excit_rd2 <- doublet_excit_rd2 %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% excit_doubl_rd2) %>%
  pull(barcodes)

# Remove doublets ---------------------------------------------------------
# Remove doublet cells
leng_sfg_liger_rd03 <- subsetLiger(object = leng_sfg_liger_rd02, cells.use = setdiff(unique(leng_sfg_liger_rd02@cell.data$barcodes), doublet_excit_rd2), remove.missing = FALSE)

# Save object
write_rds(leng_sfg_liger_rd03, "RDS_Objects/leng_sfg_liger_Cleaned_Round02.RDS")

# Round 03 Clustering -----------------------------------------------------
leng_sfg_liger_rd03 <- read_rds("RDS_Objects/leng_sfg_liger_Cleaned_Round02.RDS")

leng_sfg_liger_rd03 <- liger::normalize(leng_sfg_liger_rd03)
leng_sfg_liger_rd03 <- selectGenes(leng_sfg_liger_rd03, do.plot = T, num.genes = 350)
leng_sfg_liger_rd03 <- scaleNotCenter(leng_sfg_liger_rd03)
leng_sfg_liger_rd03 <- online_iNMF(leng_sfg_liger_rd03, k = 20, lambda = 5, miniBatch_size = 2000)
leng_sfg_liger_rd03 <- quantile_norm(leng_sfg_liger_rd03)
leng_sfg_liger_rd03 <- clusterLouvainJaccard(leng_sfg_liger_rd03,resolution = 0.2)
leng_sfg_liger_rd03 <- runUMAP(leng_sfg_liger_rd03, n_neighbors = 40)

umap_dim <- plotByDatasetAndCluster(leng_sfg_liger_rd03, return.plots = TRUE, do.legend = TRUE, text.size = 0)

beep(sound = 2)

leng_sfg_annotation_final <- tibble::tribble(
  ~cluster, ~cell_type,        ~color,
  0L,    "excit",  "dodgerblue",
  1L,    "oligo",      "orange",
  2L,    "inhib",        "navy",
  3L,    "astro", "forestgreen",
  4L,    "oligo",      "orange",
  5L,    "excit",  "dodgerblue",
  6L,    "micro",        "gold",
  7L,      "OPC",  "darkorange",
  8L,     "endo",        "pink",
  9L,    "astro", "forestgreen"
)

write_rds(leng_sfg_liger_rd03, "RDS_Objects/Leng_SFG_Round03_Clustered_FINAL.RDS")

# Convert to Seurat and Update --------------------------------------------
# convert to seurat V2
leng_sfg_seurat <- ligerToSeurat(object = leng_sfg_liger_rd03)
write_rds(leng_sfg_seurat, "RDS_Objects/leng_sfg_seuratV2_converted.RDS")
beep(sound = 2)

# Restart R session and Install Seurat v3.1.5 -----------------------------
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_3.1.5.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(Seurat)

# Switch to Seurat v3 and update object and save
leng_sfg_seurat <- read_rds("RDS_Objects/leng_sfg_seuratV2_converted.RDS")
leng_sfg_seuratV3 <- UpdateSeuratObject(leng_sfg_seurat)
View(leng_sfg_seuratV3@meta.data)
leng_sfg_seuratV3 <- PercentageFeatureSet(object = leng_sfg_seuratV3, pattern = "^MT", col.name = "percent_mito")
write_rds(leng_sfg_seuratV3, "RDS_Objects/leng_sfg_seurat_FINAL.RDS")