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
leng_ec_all <- Read10X("~/Desktop/Literature Datasets/Human DataSets/Post-Mortem/Leng 2020 Biorxiv (Ctrl & AD)/EC_ALL/")

# Create object
leng_ec_all <- CreateSeuratObject(counts = leng_ec_all, min.cells = 5, min.features = 200)
View(leng_ec_all@meta.data)

# Add mito%
leng_ec_all <- PercentageFeatureSet(object = leng_ec_all, pattern = "^MT", col.name = "percent_mito")
write_rds(leng_ec_all, "RDS_Objects/leng_ec_full_seurat.RDS")

beep(sound = 2)

# QC Filtering ------------------------------------------------------------
# EC Dataset
stats_ec <- Cluster_Stats_All_Samples(leng_ec_all)

QC_Plots_Genes(leng_ec_all, high_cutoff = 7000, low_cutoff = 400)
QC_Plots_UMI(leng_ec_all, low_cutoff = 500) 
QC_Plots_Mito(leng_ec_all, high_cutoff = 30)

leng_ec_all <- subset(x = leng_ec_all, subset = nCount_RNA < 20000 & percent_mito < 30 & nFeature_RNA < 7000)
write_rds(leng_ec_all, "RDS_Objects/leng_ec_filtered_seurat.RDS")

# Convert to Liger
leng_ec_all_liger <- seuratToLiger(objects = leng_ec_all, combined.seurat = TRUE, meta.var = "orig.ident", remove.missing = FALSE)
write_rds(leng_ec_all_liger, "RDS_Objects/leng_ec_liger_RAW.RDS")

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
leng_ec_liger <- read_rds("RDS_Objects/leng_ec_liger_RAW.RDS")

leng_ec_liger <- liger::normalize(leng_ec_liger)
leng_ec_liger <- selectGenes(leng_ec_liger, do.plot = T, num.genes = 500)
leng_ec_liger <- scaleNotCenter(leng_ec_liger)
leng_ec_liger <- online_iNMF(leng_ec_liger, k = 35, lambda = 5, miniBatch_size = 1000)
leng_ec_liger <- quantile_norm(leng_ec_liger, knn_k = 15)
leng_ec_liger <- clusterLouvainJaccard(leng_ec_liger,resolution = 1.2)
leng_ec_liger <- runUMAP(leng_ec_liger)

umap_dim <- plotByDatasetAndCluster(leng_ec_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")

beep(sound = 2)

write_rds(leng_ec_liger, "RDS_Objects/leng_ec_liger_Round01_clustered.RDS")

# Annotation --------------------------------------------------------------
leng_ec_annotation <- tibble::tribble(
  ~cluster,   ~cell_type,                                    ~notes,
  0L,      "oligo",                                        NA,
  1L,        "OPC",                                        NA,
  2L,      "micro",                                        NA,
  3L,      "oligo",                                        NA,
  4L,      "astro",                                        NA,
  5L,      "astro",                                        NA,
  6L,      "excit",                                        NA,
  7L,      "inhib",                                        NA,
  8L,      "oligo",                                        NA,
  9L,      "inhib",                                        NA,
  10L,      "excit",                                        NA,
  11L,      "excit",                                        NA,
  12L,      "micro",                                        NA,
  13L,      "excit",                                        NA,
  14L,      "excit",                                        NA,
  15L,      "inhib",                                        NA,
  16L,      "micro",                                        NA,
  17L,      "excit",                                        NA,
  18L, "endo_fibro",                                        NA,
  19L,      "excit",                                        NA,
  20L,      "excit",                                        NA,
  21L,    "doublet",                   "micro neuron doublets",
  22L,      "astro",                                        NA,
  23L,    "doublet",                   "oligo neuron doublets",
  24L,      "astro",                                        NA,
  25L,        "OPC",                                        NA,
  26L,        "OPC",                                        NA,
  27L,      "excit",                                        NA,
  28L,      "micro", "could be PBMCs.  Check after subcluster",
  29L,    "doublet",                   "OPC with neuron genes",
  30L,      "astro",                                        NA
)


# Pull Clusters for Cleaning ----------------------------------------------

# Pull cell type cluster IDs
oligo_clu <- leng_ec_annotation %>% 
  filter(cell_type == "oligo") %>% 
  pull(cluster)

excit_clu <- leng_ec_annotation %>% 
  filter(cell_type == "excit") %>% 
  pull(cluster)

inhib_clu <- leng_ec_annotation %>% 
  filter(cell_type == "inhib") %>% 
  pull(cluster)

micro_clu <- leng_ec_annotation %>% 
  filter(cell_type == "micro") %>% 
  pull(cluster)

astro_clu <- leng_ec_annotation %>% 
  filter(cell_type == "astro") %>% 
  pull(cluster)

opc_clu <- leng_ec_annotation %>% 
  filter(cell_type == "OPC") %>% 
  pull(cluster)

endo_clu <- leng_ec_annotation %>% 
  filter(cell_type == "endo_fibro") %>% 
  pull(cluster)

doublet_clu <- leng_ec_annotation %>% 
  filter(cell_type == "doublet") %>% 
  pull(cluster)


# Doublets ----------------------------------------------------------------
doublets <- data.frame(leng_ec_liger@clusters, stringsAsFactors = FALSE)
names(doublets)[1] <- "cluster"
doublets <- doublets %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% doublet_clu) %>% 
  pull(barcodes)

# Microglia ---------------------------------------------------------------
micro_sub <- subsetLiger(object = leng_ec_liger, clusters.use = micro_clu, remove.missing = FALSE)

# Remove EC3 because only 26 cells and confounds downstream subclustering
micro_minus_EC3 <- data.frame(micro_sub@cell.data, stringsAsFactors = FALSE)
micro_minus_EC3 <- micro_minus_EC3 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(dataset != "EC3") %>% 
  pull(barcodes)

micro_minus_EC3 <- subsetLiger(micro_sub, cells.use = micro_minus_EC3, remove.missing = FALSE)

micro_minus_EC3 <- normalize(micro_minus_EC3)
micro_minus_EC3 <- selectGenes(micro_minus_EC3, do.plot = T, num.genes = 200)
micro_minus_EC3 <- scaleNotCenter(micro_minus_EC3)
micro_minus_EC3 <- online_iNMF(micro_minus_EC3, k = 25, lambda = 5, miniBatch_size = 100)
micro_minus_EC3 <- quantile_norm(micro_minus_EC3, knn_k = 15)
micro_minus_EC3 <- clusterLouvainJaccard(micro_minus_EC3,resolution = 1)
micro_minus_EC3 <- runUMAP(micro_minus_EC3, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(micro_minus_EC3, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Oligo ---------------------------------------------------------------
oligo_sub <- subsetLiger(object = leng_ec_liger, clusters.use = oligo_clu, remove.missing = FALSE)
# Remove EC3 because only 26 cells and confounds downstream subclustering
oligo_minus_EC3 <- data.frame(oligo_sub@cell.data, stringsAsFactors = FALSE)
oligo_minus_EC3 <- oligo_minus_EC3 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(dataset != "EC3") %>% 
  pull(barcodes)

oligo_minus_EC3 <- subsetLiger(oligo_sub, cells.use = oligo_minus_EC3, remove.missing = FALSE)

oligo_sub <- normalize(oligo_sub)
oligo_sub <- selectGenes(oligo_sub, do.plot = T, num.genes = 150)
oligo_sub <- scaleNotCenter(oligo_sub)
oligo_sub <- optimizeALS(object = oligo_sub, k = 25, lambda = 5, nrep = 1)
oligo_sub <- quantile_norm(oligo_sub, knn_k = 15)
oligo_sub <- clusterLouvainJaccard(oligo_sub,resolution = 2.5)
oligo_sub <- runUMAP(oligo_sub, min_dist = 0.3, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(oligo_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
oligo_doubl <- c("33", "16")
doublet_oligo <- data.frame(oligo_sub@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)

# Excit---------------------------------------------------------------
excit_sub <- subsetLiger(object = leng_ec_liger, clusters.use = excit_clu, remove.missing = FALSE)
excit_sub <- normalize(excit_sub)
excit_sub <- selectGenes(excit_sub, do.plot = T, num.genes = 300)
excit_sub <- scaleNotCenter(excit_sub)
excit_sub <- online_iNMF(excit_sub, k = 25, lambda = 5, miniBatch_size = 100)
excit_sub <- quantile_norm(excit_sub)
excit_sub <- clusterLouvainJaccard(excit_sub,resolution = 1.2)
excit_sub <- runUMAP(excit_sub)

umap_dim <- plotByDatasetAndCluster(excit_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# inhib ---------------------------------------------------------------
inhib_sub <- subsetLiger(object = leng_ec_liger, clusters.use = inhib_clu, remove.missing = FALSE)
inhib_sub <- normalize(inhib_sub)
inhib_sub <- selectGenes(inhib_sub, do.plot = T, num.genes = 200)
inhib_sub <- scaleNotCenter(inhib_sub)
inhib_sub <- online_iNMF(inhib_sub, k = 25, lambda = 5, miniBatch_size = 100)
inhib_sub <- quantile_norm(inhib_sub, knn_k = 15)
inhib_sub <- clusterLouvainJaccard(inhib_sub,resolution = 1.2)
inhib_sub <- runUMAP(inhib_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(inhib_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
inhib_doubl <- c("13")
doublet_inhib <- data.frame(inhib_sub@clusters, stringsAsFactors = FALSE)
names(doublet_inhib)[1] <- "cluster"
doublet_inhib <- doublet_inhib %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% inhib_doubl) %>%
  pull(barcodes)

# Astro ---------------------------------------------------------------
astro_sub <- subsetLiger(object = leng_ec_liger, clusters.use = astro_clu, remove.missing = FALSE)

# Remove EC3 because only 26 cells and confounds downstream subclustering
astro_minus_EC3 <- data.frame(astro_sub@cell.data, stringsAsFactors = FALSE)
astro_minus_EC3 <- astro_minus_EC3 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(dataset != "EC3") %>% 
  pull(barcodes)

astro_minus_EC3 <- subsetLiger(astro_sub, cells.use = astro_minus_EC3, remove.missing = FALSE)

astro_minus_EC3 <- normalize(astro_minus_EC3)
astro_minus_EC3 <- selectGenes(astro_minus_EC3, do.plot = T, num.genes = 100)
astro_minus_EC3 <- scaleNotCenter(astro_minus_EC3)
astro_minus_EC3 <- online_iNMF(astro_minus_EC3, k = 25, lambda = 5, miniBatch_size = 100)
#astro_minus_EC3 <- optimizeALS(astro_minus_EC3, k = 30, lambda = 5, nrep = 3)
astro_minus_EC3 <- quantile_norm(astro_minus_EC3,knn_k = 15)
astro_minus_EC3 <- clusterLouvainJaccard(astro_minus_EC3,resolution = 1.2)
astro_minus_EC3 <- runUMAP(astro_minus_EC3, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(astro_minus_EC3, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

# Elim Doublets
astro_doubl <- c("14")
doublet_astro <- data.frame(astro_minus_EC3@clusters, stringsAsFactors = FALSE)
names(doublet_astro)[1] <- "cluster"
doublet_astro <- doublet_astro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_doubl) %>% 
  pull(barcodes)

# OPC ---------------------------------------------------------------
opc_sub <- subsetLiger(object = leng_ec_liger, clusters.use = opc_clu, remove.missing = FALSE)
opc_sub <- normalize(opc_sub)
opc_sub <- selectGenes(opc_sub, do.plot = T, num.genes = 100)
opc_sub <- scaleNotCenter(opc_sub)
opc_sub <- online_iNMF(opc_sub, k = 25, lambda = 5, miniBatch_size = 100)
opc_sub <- quantile_norm(opc_sub, knn_k = 15)
opc_sub <- clusterLouvainJaccard(opc_sub,resolution = 0.5)
opc_sub <- runUMAP(opc_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(opc_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

# Elim Doublets
opc_doubl <- c("9", "3")
doublet_opc <- data.frame(opc_sub@clusters, stringsAsFactors = FALSE)
names(doublet_opc)[1] <- "cluster"
doublet_opc <- doublet_opc %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% opc_doubl) %>% 
  pull(barcodes)

# Endo ---------------------------------------------------------------
endo_sub <- subsetLiger(object = leng_ec_liger, clusters.use = endo_clu, remove.missing = FALSE)
endo_sub <- normalize(endo_sub)
endo_sub <- selectGenes(endo_sub, do.plot = T, num.genes = 200)
endo_sub <- scaleNotCenter(endo_sub)
endo_sub <- online_iNMF(endo_sub, k = 20, lambda = 5, miniBatch_size = 200)
endo_sub <- quantile_norm(endo_sub)
endo_sub <- clusterLouvainJaccard(endo_sub,resolution = 1.5)
endo_sub <- runUMAP(endo_sub)

umap_dim <- plotByDatasetAndCluster(endo_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)

# Elim Doublets
endo_doubl <- c("3", "8")
doublet_endo <- data.frame(endo_sub@clusters, stringsAsFactors = FALSE)
names(doublet_endo)[1] <- "cluster"
doublet_endo <- doublet_endo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% endo_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete <- c(doublet_astro, doublet_inhib, doublet_oligo, doublet_opc, doublets)

# Remove doublets
# Add barcodes to cell data
leng_ec_liger@cell.data$barcodes <- rownames(leng_ec_liger@cell.data)

# Remove doublet cells
leng_ec_liger_cleaned <- subsetLiger(object = leng_ec_liger, cells.use = setdiff(unique(leng_ec_liger@cell.data$barcodes), doublets_complete), remove.missing = FALSE)

# Round02 Clustering ------------------------------------------------------
# Load liger objects
leng_ec_liger_cleaned <- liger::normalize(leng_ec_liger_cleaned)
leng_ec_liger_cleaned <- selectGenes(leng_ec_liger_cleaned, do.plot = T, num.genes = 300)
leng_ec_liger_cleaned <- scaleNotCenter(leng_ec_liger_cleaned)
leng_ec_liger_cleaned <- online_iNMF(leng_ec_liger_cleaned, k = 30, lambda = 5, miniBatch_size = 1000)
leng_ec_liger_cleaned <- quantile_norm(leng_ec_liger_cleaned)
leng_ec_liger_cleaned <- clusterLouvainJaccard(leng_ec_liger_cleaned,resolution = 0.4)
leng_ec_liger_cleaned <- runUMAP(leng_ec_liger_cleaned, n_neighbors = 40)

umap_dim <- plotByDatasetAndCluster(leng_ec_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)

beep(sound = 2)

write_rds(leng_ec_liger_cleaned, "RDS_Objects/leng_ec_liger_cleaned_Round02_clustered.RDS")

beep(sound = 2)

# Annotation --------------------------------------------------------------
final_annotation <- tibble::tribble(
  ~cluster,        ~cell_type,        ~color,
  0L,           "oligo",      "orange",
  1L,           "micro",        "gold",
  2L,           "excit",  "dodgerblue",
  3L,             "OPC",  "darkorange",
  4L,           "astro", "forestgreen",
  5L,           "inhib",        "navy",
  6L,           "oligo",      "orange",
  7L,           "inhib",        "navy",
  8L,           "astro", "forestgreen",
  9L,           "excit",  "dodgerblue",
  10L,           "excit",  "dodgerblue",
  11L,           "excit",  "dodgerblue",
  12L, "endo_fibro_peri",        "pink"
)

# Plot recolored
# pull colors
ec_cluster_colors <- final_annotation %>% 
  pull(color)

# plot
umap_dim <- plotByDatasetAndCluster(leng_ec_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)

ggsave("plots/Round02/Leng_ec_round02_cluster_colored_FINAL.pdf")

# Convert to Seurat and Update --------------------------------------------
# convert to seurat V2
leng_ec_seurat <- ligerToSeurat(object = leng_ec_liger_cleaned)
write_rds(leng_ec_seurat, "RDS_Objects/leng_ec_seuratV2_converted.RDS")
beep(sound = 2)

# Restart R session and Install Seurat v3.1.5 -----------------------------
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_3.1.5.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(Seurat)

# Load and update object
leng_ec_seurat <- read_rds("RDS_Objects/leng_ec_seuratV2_converted.RDS")
leng_ec_seuratV3 <- UpdateSeuratObject(leng_ec_seurat)
View(leng_ec_seuratV3@meta.data)

# Add mito percent back
leng_ec_seuratV3 <- PercentageFeatureSet(object = leng_ec_seuratV3, pattern = "^MT", col.name = "percent_mito")

# Save seurat V3 object
write_rds(leng_ec_seuratV3, "RDS_Objects/Leng_EC_Seurat_FINAL.RDS")