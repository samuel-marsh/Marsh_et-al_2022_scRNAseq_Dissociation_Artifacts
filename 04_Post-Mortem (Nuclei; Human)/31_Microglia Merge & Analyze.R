# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Seurat V3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# Load microglia subclustered Seurat Objects -------------------------------------------
marsh_seurat3_micro <- read_rds("RDS_SeuratV3/cell_type_subset/marsh_micro_seurat.RDS")
zhou_seurat3_micro <- read_rds("RDS_SeuratV3/cell_type_subset/zhou_micro_seurat.RDS")
morabito_seurat3_micro <- read_rds("RDS_SeuratV3/cell_type_subset/morabito_micro_seurat.RDS")
leng_ec_seurat3_micro <- read_rds("RDS_SeuratV3/cell_type_subset/leng_ec_micro_seurat.RDS")
leng_sfg_seurat3_micro <- read_rds("RDS_SeuratV3/cell_type_subset/leng_sfg_micro_seurat.RDS")

# merge the cells & save
micro_merged <- merge(x = marsh_seurat3_micro, y = c(morabito_seurat3_micro, zhou_seurat3_micro, leng_ec_seurat3_micro, leng_sfg_seurat3_micro))

# Write meta data
marsh_micro_meta <- micro_merged@meta.data
write_rds(marsh_micro_meta, "03_Meta Data/marsh_micro_seurat_meta.RDS")

# convert to liger
micro_merged_liger <- seuratToLiger(objects = micro_merged, combined.seurat = T, meta.var = "orig.ident", remove.missing = FALSE)
write_rds(micro_merged_liger, "RDS_subset_merge_liger/micro/micro_merged_liger_RAW.RDS")

# Restart R and install Seurat V2.3.4 -------------------------------------
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Seurat V2.3.4
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

# Liger Analyze -----------------------------------------------------------
micro_merged_liger <- read_rds("RDS_subset_merge_liger/micro/micro_merged_liger_RAW.RDS")
micro_meta_seurat <- read_rds("03_Meta Data/marsh_micro_seurat_meta.RDS")

micro_merged_liger@cell.data$source <- as.factor(micro_meta_seurat$Dataset)

micro_merged_liger <- reorganizeLiger(object = micro_merged_liger, by.feature = "source", remove.missing = FALSE)

micro_merged_liger <- normalize(micro_merged_liger)
micro_merged_liger <- selectGenes(micro_merged_liger, do.plot = T, num.genes = 300)
micro_merged_liger <- scaleNotCenter(micro_merged_liger)
micro_merged_liger <- online_iNMF(micro_merged_liger, k = 25, lambda = 5, miniBatch_size = 1000)
micro_merged_liger <- quantile_norm(micro_merged_liger)
micro_merged_liger <- clusterLouvainJaccard(micro_merged_liger,resolution = 2)
micro_merged_liger <- runUMAP(micro_merged_liger, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(micro_merged_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/micro_merged.pdf")
plotFactors(micro_merged_liger, num.genes = 8, plot.tsne = T)
dev.off()

# Doublets
# Pulling barcodes for cleaning
NK_doubl <- c("23")
doublet_NK <- data.frame(micro_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_NK)[1] <- "cluster"
doublet_NK <- doublet_NK %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% NK_doubl) %>% 
  pull(barcodes)

# Neurons
neuron_doubl <- c("4")
doublet_neuron <- data.frame(micro_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_neuron)[1] <- "cluster"
doublet_neuron <- doublet_neuron %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% neuron_doubl) %>% 
  pull(barcodes)

# Monocytes
mono_doubl <- c("13")
doublet_mono <- data.frame(micro_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_mono)[1] <- "cluster"
doublet_mono <- doublet_mono %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% mono_doubl) %>% 
  pull(barcodes)

# Neuron and oligos 
neuron_oligo_doubl <- c("19")
doublet_neuron_oligo <- data.frame(micro_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_neuron_oligo)[1] <- "cluster"
doublet_neuron_oligo <- doublet_neuron_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% neuron_oligo_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete <- c(doublet_mono, doublet_NK, doublet_neuron, doublet_neuron_oligo)

# Remove doublets
# Add barcodes to cell data
micro_merged_liger@cell.data$barcodes <- rownames(micro_merged_liger@cell.data)

# Remove doublet cells
micro_merged_liger_cleaned <- subsetLiger(object = micro_merged_liger, cells.use = setdiff(unique(micro_merged_liger@cell.data$barcodes), doublets_complete), remove.missing = FALSE)

# round02 -----------------------------------------------------------------
micro_merged_liger_cleaned <- normalize(micro_merged_liger_cleaned)
micro_merged_liger_cleaned <- selectGenes(micro_merged_liger_cleaned, do.plot = T, num.genes = 300)
micro_merged_liger_cleaned <- scaleNotCenter(micro_merged_liger_cleaned)
micro_merged_liger_cleaned <- online_iNMF(micro_merged_liger_cleaned, k = 20, lambda = 5, miniBatch_size = 1000)
micro_merged_liger_cleaned <- quantile_norm(micro_merged_liger_cleaned)
micro_merged_liger_cleaned <- clusterLouvainJaccard(micro_merged_liger_cleaned,resolution = 0.6)
micro_merged_liger_cleaned <- runUMAP(micro_merged_liger_cleaned, min_dist = 0.3, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(micro_merged_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/micro_merged_rd2_factors.pdf")
plotFactors(micro_merged_liger_cleaned, num.genes = 8, plot.tsne = T)
dev.off()

# Save Liger Objects
write_rds(micro_merged_liger_cleaned, "RDS_subset_merge_liger/micro_subset_cleaned.RDS")

micro_samples <- micro_merged_liger_cleaned@cell.data$orig.dataset
write_rds(micro_samples, "RDS_Subcluster_Seurat/micro_sample_names.RDS")

# Convert to Seurat
micro_seurat <- ligerToSeurat(micro_merged_liger_cleaned)
write_rds(micro_seurat, "RDS_Subcluster_Seurat/micro_rd2_seurat.RDS")
