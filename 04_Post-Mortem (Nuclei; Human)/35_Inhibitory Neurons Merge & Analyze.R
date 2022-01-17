# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Seurat V3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

# Load inhibglia subclustered -------------------------------------------
marsh_seurat3_inhib <- read_rds("RDS_SeuratV3/cell_type_subset/marsh_inhib_seurat.RDS")
zhou_seurat3_inhib <- read_rds("RDS_SeuratV3/cell_type_subset/zhou_inhib_seurat.RDS")
morabito_seurat3_inhib <- read_rds("RDS_SeuratV3/cell_type_subset/morabito_inhib_seurat.RDS")
leng_ec_seurat3_inhib <- read_rds("RDS_SeuratV3/cell_type_subset/leng_ec_inhib_seurat.RDS")
leng_sfg_seurat3_inhib <- read_rds("RDS_SeuratV3/cell_type_subset/leng_sfg_inhib_seurat.RDS")

# merge the cells & save
inhib_merged <- merge(x = marsh_seurat3_inhib, y = c(morabito_seurat3_inhib, zhou_seurat3_inhib, leng_ec_seurat3_inhib, leng_sfg_seurat3_inhib))

# meta data
marsh_inhib_meta <- inhib_merged@meta.data
write_rds(marsh_inhib_meta, "03_Meta Data/marsh_inhib_seurat_meta.RDS")

# convert to liger and save
inhib_merged_liger <- seuratToLiger(objects = inhib_merged, combined.seurat = T, meta.var = "orig.ident", remove.missing = FALSE)
write_rds(inhib_merged_liger, "RDS_subset_merge_liger/inhib/inhib_merged_liger_RAW.RDS")

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

# round01 Liger Analyze -----------------------------------------------------------
inhib_merged_liger <- read_rds("RDS_subset_merge_liger/inhib/inhib_merged_liger_RAW.RDS")
inhib_meta_seurat <- read_rds("03_Meta Data/marsh_inhib_seurat_meta.RDS")

inhib_merged_liger@cell.data$source <- as.factor(inhib_meta_seurat$Dataset)

inhib_merged_liger <- reorganizeLiger(object = inhib_merged_liger, by.feature = "source", remove.missing = FALSE)
View(inhib_merged_liger@cell.data)

inhib_merged_liger <- normalize(inhib_merged_liger)
inhib_merged_liger <- selectGenes(inhib_merged_liger, do.plot = T, num.genes = 600)
inhib_merged_liger <- scaleNotCenter(inhib_merged_liger)
inhib_merged_liger <- online_iNMF(inhib_merged_liger, k = 25, lambda = 5, miniBatch_size = 1000)
inhib_merged_liger <- quantile_norm(inhib_merged_liger)
inhib_merged_liger <- clusterLouvainJaccard(inhib_merged_liger,resolution = 2)
inhib_merged_liger <- runUMAP(inhib_merged_liger, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(inhib_merged_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/inhib_merged.pdf")
plotFactors(inhib_merged_liger, num.genes = 8, plot.tsne = T)
dev.off()

# Doublet Barcodes --------------------------------------------------------
# Pulling barcodes for cleaning
oligo_doubl <- c("23")
doublet_oligo <- data.frame(inhib_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)

endo_fibro_doubl <- c("26")
doublet_endo_fibro <- data.frame(inhib_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_endo_fibro)[1] <- "cluster"
doublet_endo_fibro <- doublet_endo_fibro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% endo_fibro_doubl) %>% 
  pull(barcodes)

micro_doubl <- c("25")
doublet_micro <- data.frame(inhib_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_micro)[1] <- "cluster"
doublet_micro <- doublet_micro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% micro_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_inhib <- c(doublet_oligo, doublet_endo_fibro, doublet_micro)

# Remove doublets
# Add barcodes to cell data
inhib_merged_liger@cell.data$barcodes <- rownames(inhib_merged_liger@cell.data)

# Remove doublet cells
inhib_merged_liger_cleaned <- subsetLiger(object = inhib_merged_liger, cells.use = setdiff(unique(inhib_merged_liger@cell.data$barcodes), doublets_complete_inhib), remove.missing = FALSE)

# round02 -----------------------------------------------------------------
inhib_merged_liger_cleaned <- normalize(inhib_merged_liger_cleaned)
inhib_merged_liger_cleaned <- selectGenes(inhib_merged_liger_cleaned, do.plot = T, num.genes = 625)
inhib_merged_liger_cleaned <- scaleNotCenter(inhib_merged_liger_cleaned)
inhib_merged_liger_cleaned <- online_iNMF(inhib_merged_liger_cleaned, k = 25, lambda = 5, miniBatch_size = 1000)
inhib_merged_liger_cleaned <- quantile_norm(inhib_merged_liger_cleaned)
inhib_merged_liger_cleaned <- clusterLouvainJaccard(inhib_merged_liger_cleaned,resolution = 2)
inhib_merged_liger_cleaned <- runUMAP(inhib_merged_liger_cleaned, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(inhib_merged_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/inhib_merged_cleaned.pdf")
plotFactors(inhib_merged_liger_cleaned, num.genes = 8, plot.tsne = T)
dev.off()

# Pulling barcodes for cleaning
astro_oligo_doubl <- c("27")
doublet_astro_oligo <- data.frame(inhib_merged_liger_cleaned@clusters, stringsAsFactors = FALSE)
names(doublet_astro_oligo)[1] <- "cluster"
doublet_astro_oligo <- doublet_astro_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_oligo_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_inhib2 <- c(doublet_astro_oligo)

# Remove doublets
# Add barcodes to cell data
inhib_merged_liger_cleaned@cell.data$barcodes <- rownames(inhib_merged_liger_cleaned@cell.data)

# Remove doublet cells
inhib_merged_liger_cleaned2 <- subsetLiger(object = inhib_merged_liger_cleaned, cells.use = setdiff(unique(inhib_merged_liger_cleaned@cell.data$barcodes), doublets_complete_inhib2), remove.missing = FALSE)

# round03 -----------------------------------------------------------------
inhib_merged_liger_cleaned2 <- normalize(inhib_merged_liger_cleaned2)
inhib_merged_liger_cleaned2 <- selectGenes(inhib_merged_liger_cleaned2, do.plot = T, num.genes = 700)
inhib_merged_liger_cleaned2 <- scaleNotCenter(inhib_merged_liger_cleaned2)
inhib_merged_liger_cleaned2 <- online_iNMF(inhib_merged_liger_cleaned2, k = 25, lambda = 5, miniBatch_size = 1000)
inhib_merged_liger_cleaned2 <- quantile_norm(inhib_merged_liger_cleaned2)
inhib_merged_liger_cleaned2 <- clusterLouvainJaccard(inhib_merged_liger_cleaned2,resolution = 2)
inhib_merged_liger_cleaned2 <- runUMAP(inhib_merged_liger_cleaned2)

umap_dim <- plotByDatasetAndCluster(inhib_merged_liger_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/inhib_merged_cleaned2.pdf")
plotFactors(inhib_merged_liger_cleaned2, num.genes = 8, plot.tsne = T)
dev.off()

# Pulling barcodes for cleaning
oligo_doubl2 <- c("26")
doublet_oligo2 <- data.frame(inhib_merged_liger_cleaned2@clusters, stringsAsFactors = FALSE)
names(doublet_oligo2)[1] <- "cluster"
doublet_oligo2 <- doublet_oligo2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl2) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_inhib3 <- c(doublet_oligo2)

# Remove doublets
# Add barcodes to cell data
inhib_merged_liger_cleaned@cell.data$barcodes <- rownames(inhib_merged_liger_cleaned@cell.data)

# Remove doublet cells
inhib_merged_liger_cleaned3 <- subsetLiger(object = inhib_merged_liger_cleaned2, cells.use = setdiff(unique(inhib_merged_liger_cleaned2@cell.data$barcodes), doublets_complete_inhib3), remove.missing = FALSE)

# round04 -----------------------------------------------------------------
inhib_merged_liger_cleaned3 <- normalize(inhib_merged_liger_cleaned3)
inhib_merged_liger_cleaned3 <- selectGenes(inhib_merged_liger_cleaned3, do.plot = T, num.genes = 600)
inhib_merged_liger_cleaned3 <- scaleNotCenter(inhib_merged_liger_cleaned3)
inhib_merged_liger_cleaned3 <- online_iNMF(inhib_merged_liger_cleaned3, k = 25, lambda = 5, miniBatch_size = 1000)
inhib_merged_liger_cleaned3 <- quantile_norm(inhib_merged_liger_cleaned3)
inhib_merged_liger_cleaned3 <- clusterLouvainJaccard(inhib_merged_liger_cleaned3,resolution = 2)
inhib_merged_liger_cleaned3 <- runUMAP(inhib_merged_liger_cleaned3, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(inhib_merged_liger_cleaned3, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/inhib_merged_cleaned3.pdf")
plotFactors(inhib_merged_liger_cleaned3, num.genes = 8, plot.tsne = T)
dev.off()

# save liger objects
write_rds(inhib_merged_liger_cleaned3, "RDS_subset_merge_liger/inhib_subset_cleaned3.RDS")

inhib_samples <- inhib_merged_liger_cleaned3@cell.data$orig.dataset
write_rds(inhib_samples, "RDS_Subcluster_Seurat/inhib_sample_names.RDS")

# Convert to Seurat and Save
inhib_seurat <- ligerToSeurat(inhib_merged_liger_cleaned3)
write_rds(inhib_seurat, "RDS_Subcluster_Seurat/inhib_rd4_seurat.RDS")