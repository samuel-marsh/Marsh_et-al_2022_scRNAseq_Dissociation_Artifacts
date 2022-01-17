# Load Packages
library(tidyverse)
library(Seurat)
  # Starting with Seurat V3
library(liger)
  # Using "online" branch of liger
library(marsh.utils)
library(Matrix)
library(scCustomize)


# Morabito_RAW -----------------------------------------------------------
# Read raw data Seurat
morabito <- Read10X("lib_file_path")
    # Raw data downloaded from syn18915937

# Create Seurat
morabito <- CreateSeuratObject(counts = morabito, min.cells = 5, min.features = 200, names.field = 2, names.delim = "-")

# Mito %
morabito <- PercentageFeatureSet(object = morabito, pattern = "^MT", col.name = "percent_mito")

# Annotated libraries with metadata
sample23 <- "1"
sample187 <- c("2", "5")
sample280 <- c("3", "6")
sample121 <- c("4", "7")
sample165 <- c("8", "9")


# Add metadata to match library to sample ID
morabito@meta.data$sample_id[morabito@meta.data$orig.ident %in% sample23] <- "23"
morabito@meta.data$sample_id[morabito@meta.data$orig.ident %in% sample187] <- "187"
morabito@meta.data$sample_id[morabito@meta.data$orig.ident %in% sample280] <- "280"
morabito@meta.data$sample_id[morabito@meta.data$orig.ident %in% sample121] <- "121"
morabito@meta.data$sample_id[morabito@meta.data$orig.ident %in% sample165] <- "165"


# QC Filter 
# major outliers filter above 50K.  Filter and replot
morabito <- subset(x = morabito, subset = nCount_RNA < 50000)

# Subset the data agaain
morabito <- subset(x = morabito, subset = percent_mito < 30 & nFeature_RNA > 400)

# Convert to Liger & Save
morabito_liger <- seuratToLiger(objects = morabito, combined.seurat = TRUE, meta.var = "sample_id")
write_rds(morabito_liger, "RDS_Objects/RAW_liger/morabito_liger_RAW.RDS")



# Restart R Session and install Seurat V2 ---------------------------------
install.packages("~/seurat/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")
library(tidyverse)
library(Seurat)
library(liger)
library(marsh.utils)
library(patchwork)

# Read in liger object
morabito <- read_rds("RDS_Objects/RAW_liger/morabito_liger_RAW.RDS")

# Run first round analysis
morabito <-  normalize(morabito)
morabito <-  selectGenes(morabito, do.plot = T, num.genes = 900)
morabito <-  scaleNotCenter(morabito)
morabito <-  online_iNMF(morabito, k = 40, lambda = 10, miniBatch_size = 4000)
morabito <-  quantile_norm(morabito)
morabito <-  clusterLouvainJaccard(morabito,resolution = 1)
morabito <-  runUMAP(morabito)

umap_dim <- plotByDatasetAndCluster(morabito, return.plots = TRUE, do.legend = TRUE, text.size = 4)

# Annotate clusters using marker genes and factor loadings from liger for top level cell types for subclustering and doublet cleaning
morabito_annotation <- tibble::tribble(
  ~cluster,  ~cell_type,        ~color,
  0L,     "oligo",      "orange",
  1L,     "oligo",      "orange",
  2L,     "oligo",      "orange",
  3L,     "excit",  "dodgerblue",
  4L,     "excit",  "dodgerblue",
  5L,     "excit",  "dodgerblue",
  6L,       "OPC",  "darkorange",
  7L,     "astro", "forestgreen",
  8L,     "excit",  "dodgerblue",
  9L,     "astro", "forestgreen",
  10L,     "inhib",        "navy",
  11L,     "micro",        "gold",
  12L,     "inhib",        "navy",
  13L,     "inhib",        "navy",
  14L,     "excit",  "dodgerblue",
  15L,     "excit",  "dodgerblue",
  16L,     "inhib",        "navy",
  17L,     "oligo",      "orange",
  18L,     "astro", "forestgreen",
  19L,     "inhib",        "navy",
  20L,     "excit",  "dodgerblue",
  21L,     "oligo",      "orange",
  22L,     "astro", "forestgreen",
  23L,     "excit",  "dodgerblue",
  24L,     "excit",  "dodgerblue",
  25L,     "excit",  "dodgerblue",
  26L,     "inhib",        "navy",
  27L,     "oligo",      "orange",
  28L, "pericytes", "mediumorchid1",
  29L,     "Fibro",        "pink",
  30L,      "endo", "darkorchid3",
  31L,  "doublets",   "lightgray",
  32L,     "oligo",      "orange"
)


# Pull cell type cluster IDs
oligo_clu <- morabito_annotation %>% 
  filter(cell_type == "oligo") %>% 
  pull(cluster)

excit_clu <- morabito_annotation %>% 
  filter(cell_type == "excit") %>% 
  pull(cluster)

inhib_clu <- morabito_annotation %>% 
  filter(cell_type == "inhib") %>% 
  pull(cluster)

micro_clu <- morabito_annotation %>% 
  filter(cell_type == "micro") %>% 
  pull(cluster)

astro_clu <- morabito_annotation %>% 
  filter(cell_type == "astro") %>% 
  pull(cluster)

opc_clu <- morabito_annotation %>% 
  filter(cell_type == "OPC") %>% 
  pull(cluster)


# oligo cleaning
oligo_sub <- subsetLiger(object = morabito, clusters.use = oligo_clu)

# Reanalyze
oligo_sub <-  normalize(oligo_sub)
oligo_sub <-  selectGenes(oligo_sub, do.plot = T, num.genes = 400)
oligo_sub <-  scaleNotCenter(oligo_sub)
oligo_sub <-  online_iNMF(object = oligo_sub, k = 25, lambda = 5, miniBatch_size = 800)
oligo_sub <-  quantile_norm(oligo_sub)
oligo_sub <-  clusterLouvainJaccard(oligo_sub,resolution = 0.5)
oligo_sub <-  runUMAP(oligo_sub, min_dist = 0.3)

oligo_dim <- plotByDatasetAndCluster(oligo_sub, return.plots = TRUE, do.legend = TRUE, text.size = 4)

# Pull doublet barcodes
  # clu7 are likely neuronal and micro doublets
  # clu8 are endothelial doublets
oligo_doubl <- c("7", "8")
doublet_oligo <- data.frame(oligo_sub@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)


# Excitatory Neurons Cleaning
excit_sub <- subsetLiger(object = morabito, clusters.use = excit_clu)

# Reanalyze
excit_sub <-  normalize(excit_sub)
excit_sub <-  selectGenes(excit_sub, do.plot = T, num.genes = 500)
excit_sub <-  scaleNotCenter(excit_sub)
excit_sub <-  online_iNMF(object = excit_sub, k = 25, lambda = 5, miniBatch_size = 1000)
excit_sub <-  quantile_norm(excit_sub)
excit_sub <-  clusterLouvainJaccard(excit_sub,resolution = 1.2)
excit_sub <-  runUMAP(excit_sub, min_dist = 0.3, n_neighbors = 15)

excit_dim <- plotByDatasetAndCluster(excit_sub, return.plots = TRUE, do.legend = TRUE, text.size = 4)

# Pull doublet barcodes
  # clu17 are likely inhibitory neuron doublets
  # clu18 are oligo doublets
excit_doubl <- c("17", "18")
doublet_excit <- data.frame(excit_sub@clusters, stringsAsFactors = FALSE)
names(doublet_excit)[1] <- "cluster"
doublet_excit <- doublet_excit %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% excit_doubl) %>% 
  pull(barcodes)


# Microglia cleaning
micro_sub <- subsetLiger(object = morabito, clusters.use = micro_clu)

# Reanalyze
micro_sub <-  normalize(micro_sub)
micro_sub <-  selectGenes(micro_sub, do.plot = T, num.genes = 300)
micro_sub <-  scaleNotCenter(micro_sub)
micro_sub <-  online_iNMF(object = micro_sub, k = 25, lambda = 5, miniBatch_size = 50)
micro_sub <-  quantile_norm(micro_sub)
micro_sub <-  clusterLouvainJaccard(micro_sub,resolution = 1)
micro_sub <-  runUMAP(micro_sub, min_dist = 0.3, n_neighbors = 15)

# Pull doublet barcodes
  # clu6 are Oligo doublets
micro_doubl <- "6"
doublet_micro <- data.frame(micro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_micro)[1] <- "cluster"
doublet_micro <- doublet_micro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% micro_doubl) %>% 
  pull(barcodes)


# Inhibitory Neurons
inhib_sub <- subsetLiger(object = morabito, clusters.use = inhib_clu)

# Reanalyze
inhib_sub <-  normalize(inhib_sub)
inhib_sub <-  selectGenes(inhib_sub, do.plot = T, num.genes = 500)
inhib_sub <-  scaleNotCenter(inhib_sub)
inhib_sub <-  online_iNMF(object = inhib_sub, k = 25, lambda = 5, miniBatch_size = 500)
inhib_sub <-  quantile_norm(inhib_sub)
inhib_sub <-  clusterLouvainJaccard(inhib_sub,resolution = 0.8)
inhib_sub <-  runUMAP(inhib_sub, min_dist = 0.3, n_neighbors = 15)

# Pull doublet barcodes
  # clu15 are Oligo doublets
  # clu8 are excit neuron doublets
inhib_doubl <- c("8", "15")
doublet_inhib <- data.frame(inhib_sub@clusters, stringsAsFactors = FALSE)
names(doublet_inhib)[1] <- "cluster"
doublet_inhib <- doublet_inhib %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% inhib_doubl) %>% 
  pull(barcodes)

# Astrocytes
astro_sub <- subsetLiger(object = morabito, clusters.use = astro_clu)

# Reanalyze
astro_sub <-  normalize(astro_sub)
astro_sub <-  selectGenes(astro_sub, do.plot = T, num.genes = 300)
astro_sub <-  scaleNotCenter(astro_sub)
astro_sub <-  online_iNMF(object = astro_sub, k = 25, lambda = 5, miniBatch_size = 500)
astro_sub <-  quantile_norm(astro_sub)
astro_sub <-  clusterLouvainJaccard(astro_sub,resolution = 1.2)
astro_sub <-  runUMAP(astro_sub, min_dist = 0.3, n_neighbors = 15)

# Pull doublet barcodes
  # clu11, clu17, clu5 are oligo and OPC doublets
  # Clu10 are are oligo and or neuron doublets
astro_doubl <- c("11", "17", "5", "10")
doublet_astro <- data.frame(astro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_astro)[1] <- "cluster"
doublet_astro <- doublet_astro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_doubl) %>% 
  pull(barcodes)

# OPC Cleaning
opc_sub <- subsetLiger(object = morabito, clusters.use = opc_clu)

# Reanalyze
opc_sub <-  normalize(opc_sub)
opc_sub <-  selectGenes(opc_sub, do.plot = T, num.genes = 300)
opc_sub <-  scaleNotCenter(opc_sub)
opc_sub <-  online_iNMF(object = opc_sub, k = 25, lambda = 5, miniBatch_size = 200)
opc_sub <-  quantile_norm(opc_sub)
opc_sub <-  clusterLouvainJaccard(opc_sub,resolution = 0.8)
opc_sub <-  runUMAP(opc_sub, min_dist = 0.3, n_neighbors = 15)

# Pull doublet barcodes
  # clu5 are mature oligo doublets
  # Clu7 are neuron doublets
opc_doubl <- c("5", "7")
doublet_opc <- data.frame(opc_sub@clusters, stringsAsFactors = FALSE)
names(doublet_opc)[1] <- "cluster"
doublet_opc <- doublet_opc %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% opc_doubl) %>%
  pull(barcodes)


# Compile and remove Doublets ---------------------------------------------
doublets_complete <- c(doublet_astro, doublet_excit, doublet_inhib, doublet_micro, doublet_oligo, doublet_opc)

# Add barcodes to cell data
morabito@cell.data$barcodes <- rownames(morabito@cell.data)

# Subset to remove doublet cells
morabito_cleaned <- subsetLiger(object = morabito, cells.use = setdiff(unique(morabito@cell.data$barcodes), doublets_complete))




# Recluster ---------------------------------------------------------------
morabito_cleaned <-  normalize(morabito_cleaned)
morabito_cleaned <-  selectGenes(morabito_cleaned, do.plot = T, num.genes = 1000)
morabito_cleaned <-  scaleNotCenter(morabito_cleaned)
morabito_cleaned <-  online_iNMF(morabito_cleaned, k = 40, lambda = 10, miniBatch_size = 4000)
morabito_cleaned <-  quantile_norm(morabito_cleaned, knn_k = 15)
morabito_cleaned <-  clusterLouvainJaccard(morabito_cleaned,resolution = 1.2)
morabito_cleaned <-  runUMAP(morabito_cleaned, n_neighbors = 15)


morabito_cleaned_annotation <- tibble::tribble(
  ~cluster,  ~cell_type,          ~color,
  0L,     "excit",    "dodgerblue",
  1L,     "oligo",        "orange",
  2L,     "oligo",        "orange",
  3L,     "excit",    "dodgerblue",
  4L,     "astro",   "forestgreen",
  5L,     "oligo",        "orange",
  6L,       "OPC",    "darkorange",
  7L,     "oligo",        "orange",
  8L,     "oligo",        "orange",
  9L,     "excit",    "dodgerblue",
  10L,     "inhib",          "navy",
  11L,     "inhib",          "navy",
  12L,     "oligo",        "orange",
  13L,     "micro",          "gold",
  14L,     "astro",   "forestgreen",
  15L,     "astro",   "forestgreen",
  16L,     "excit",    "dodgerblue",
  17L,     "inhib",          "navy",
  18L,     "inhib",          "navy",
  19L,     "excit",    "dodgerblue",
  20L,     "inhib",          "navy",
  21L,     "excit",    "dodgerblue",
  22L,     "excit",    "dodgerblue",
  23L,     "excit",    "dodgerblue",
  24L,     "excit",    "dodgerblue",
  25L,     "excit",    "dodgerblue",
  26L,      "endo",          "pink",
  27L,     "excit",    "dodgerblue",
  28L,     "fibro", "mediumorchid1",
  29L, "pericytes",   "darkorchid3",
  30L,     "excit",    "dodgerblue",
  31L,     "astro",   "forestgreen"
)

