# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Seurat V3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# Load astrocytes subclustered -------------------------------------------
marsh_seurat3_astro <- read_rds("RDS_SeuratV3/cell_type_subset/marsh_astro_seurat.RDS")
zhou_seurat3_astro <- read_rds("RDS_SeuratV3/cell_type_subset/zhou_astro_seurat.RDS")
morabito_seurat3_astro <- read_rds("RDS_SeuratV3/cell_type_subset/morabito_astro_seurat.RDS")
leng_ec_seurat3_astro <- read_rds("RDS_SeuratV3/cell_type_subset/leng_ec_astro_seurat.RDS")
leng_sfg_b_seurat3_astro <- read_rds("RDS_SeuratV3/cell_type_subset/leng_sfg_astro_seurat_b.RDS")

# merge the cells & save
astro_merged_b <- merge(x = marsh_seurat3_astro, y = c(morabito_seurat3_astro, zhou_seurat3_astro, leng_ec_seurat3_astro, leng_sfg_b_seurat3_astro))

# Meta data
marsh_astro_meta <- astro_merged_b@meta.data
write_rds(marsh_astro_meta, "03_Meta Data/marsh_astro_seurat_meta_b.RDS")

# Convert to Liger
astro_merged_b_liger <- seuratToLiger(objects = astro_merged_b, combined.seurat = T, meta.var = "orig.ident", remove.missing = FALSE)
write_rds(astro_merged_b_liger, "RDS_subset_merge_liger/astro/astro_merged_b_liger_RAW_b.RDS")

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
astro_merged_b_liger <- read_rds("RDS_subset_merge_liger/astro/astro_merged_b_liger_RAW_b.RDS")
astro_meta_seurat <- read_rds("03_Meta Data/marsh_astro_seurat_meta_b.RDS")

astro_merged_b_liger@cell.data$source <- as.factor(astro_meta_seurat$Dataset)

astro_merged_b_liger <- reorganizeLiger(object = astro_merged_b_liger, by.feature = "source", remove.missing = FALSE)
View(astro_merged_b_liger@cell.data)

astro_merged_b_liger <- normalize(astro_merged_b_liger)
astro_merged_b_liger <- selectGenes(astro_merged_b_liger, do.plot = T, num.genes = 600)
astro_merged_b_liger <- scaleNotCenter(astro_merged_b_liger)
astro_merged_b_liger <- online_iNMF(astro_merged_b_liger, k = 25, lambda = 10, miniBatch_size = 1000)
astro_merged_b_liger <- quantile_norm(astro_merged_b_liger)
astro_merged_b_liger <- clusterLouvainJaccard(astro_merged_b_liger,resolution = 0.8)
astro_merged_b_liger <- runUMAP(astro_merged_b_liger, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(astro_merged_b_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

write_rds(astro_merged_b_liger, "RDS_subset_merge_liger/astro_subset_rd1.RDS")
astro_merged_b_liger <- read_rds("RDS_subset_merge_liger/astro_subset_rd1.RDS")

pdf("05_plots/subclustering/astro_merged_b_TESTING.pdf")
plotFactors(astro_merged_b_liger, num.genes = 8, plot.tsne = T)
dev.off()

# Pulling barcodes for cleaning
oligo_doubl <- c("9")
doublet_oligo <- data.frame(astro_merged_b_liger@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_astro <- c(doublet_oligo)

# Remove doublets
# Add barcodes to cell data
astro_merged_b_liger@cell.data$barcodes <- rownames(astro_merged_b_liger@cell.data)

# Remove doublet cells
astro_merged_b_liger_cleaned <- subsetLiger(object = astro_merged_b_liger, cells.use = setdiff(unique(astro_merged_b_liger@cell.data$barcodes), doublets_complete_astro), remove.missing = FALSE)

# round02 -----------------------------------------------------------------
astro_merged_b_liger_cleaned <- normalize(astro_merged_b_liger_cleaned)
astro_merged_b_liger_cleaned <- selectGenes(astro_merged_b_liger_cleaned, do.plot = T, num.genes = 400)
astro_merged_b_liger_cleaned <- scaleNotCenter(astro_merged_b_liger_cleaned)
astro_merged_b_liger_cleaned <- online_iNMF(astro_merged_b_liger_cleaned, k = 25, lambda = 10, miniBatch_size = 1000)
astro_merged_b_liger_cleaned <- quantile_norm(astro_merged_b_liger_cleaned)
astro_merged_b_liger_cleaned <- clusterLouvainJaccard(astro_merged_b_liger_cleaned,resolution = 1)
astro_merged_b_liger_cleaned <- runUMAP(astro_merged_b_liger_cleaned, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(astro_merged_b_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/astro_merged_b_cleaned_TESTING.pdf")
plotFactors(astro_merged_b_liger_cleaned, num.genes = 8, plot.tsne = T)
dev.off()

neuron_doubl_astro2 <- c("0")
doublet_neuron_astro2 <- data.frame(astro_merged_b_liger_cleaned@clusters, stringsAsFactors = FALSE)
names(doublet_neuron_astro2)[1] <- "cluster"
doublet_neuron_astro2 <- doublet_neuron_astro2 %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% neuron_doubl_astro2) %>%
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_astro2 <- c(doublet_neuron_astro2)

# Remove doublet cells
astro_merged_b_liger_cleaned2 <- subsetLiger(object = astro_merged_b_liger_cleaned, cells.use = setdiff(unique(astro_merged_b_liger_cleaned@cell.data$barcodes), doublets_complete_astro2), remove.missing = FALSE)

# round03 -----------------------------------------------------------------
astro_merged_b_liger_cleaned2 <- normalize(astro_merged_b_liger_cleaned2)
astro_merged_b_liger_cleaned2 <- selectGenes(astro_merged_b_liger_cleaned2, do.plot = T, num.genes = 350)
astro_merged_b_liger_cleaned2 <- scaleNotCenter(astro_merged_b_liger_cleaned2)
astro_merged_b_liger_cleaned2 <- online_iNMF(astro_merged_b_liger_cleaned2, k = 25, lambda = 10, miniBatch_size = 1000)
astro_merged_b_liger_cleaned2 <- quantile_norm(astro_merged_b_liger_cleaned2)
astro_merged_b_liger_cleaned2 <- clusterLouvainJaccard(astro_merged_b_liger_cleaned2,resolution = 1)
astro_merged_b_liger_cleaned2 <- runUMAP(astro_merged_b_liger_cleaned2, n_neighbors = 20)

umap_dim <- plotByDatasetAndCluster(astro_merged_b_liger_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

#   pull(barcodes)
opc_doubl <- c("14")
doublet_opc <- data.frame(astro_merged_b_liger_cleaned2@clusters, stringsAsFactors = FALSE)
names(doublet_opc)[1] <- "cluster"
doublet_opc <- doublet_opc %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% opc_doubl) %>%
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_astro3 <- c(doublet_opc)

# Remove doublet cells
astro_merged_b_liger_cleaned3 <- subsetLiger(object = astro_merged_b_liger_cleaned2, cells.use = setdiff(unique(astro_merged_b_liger_cleaned2@cell.data$barcodes), doublets_complete_astro3), remove.missing = FALSE)

# round04 -----------------------------------------------------------------
astro_merged_b_liger_cleaned3 <- normalize(astro_merged_b_liger_cleaned3)
astro_merged_b_liger_cleaned3 <- selectGenes(astro_merged_b_liger_cleaned3, do.plot = T, num.genes = 350)
astro_merged_b_liger_cleaned3 <- scaleNotCenter(astro_merged_b_liger_cleaned3)
astro_merged_b_liger_cleaned3 <- online_iNMF(astro_merged_b_liger_cleaned3, k = 25, lambda = 7, miniBatch_size = 1000)
astro_merged_b_liger_cleaned3 <- quantile_norm(astro_merged_b_liger_cleaned3)
astro_merged_b_liger_cleaned3 <- clusterLouvainJaccard(astro_merged_b_liger_cleaned3,resolution = 0.2)
astro_merged_b_liger_cleaned3 <- runUMAP(astro_merged_b_liger_cleaned3, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(astro_merged_b_liger_cleaned3_C, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

# factors
pdf("05_plots/subclustering/astro_merged_b_cleaned3_TESTING.pdf")
Customize_plotFactors(astro_merged_b_liger_cleaned3, num.genes = 8, plot.tsne = T)
dev.off()

# save liger object
write_rds(astro_merged_b_liger_cleaned3, "RDS_subset_merge_liger/astro_subset_cleaned3_C.RDS")

astro_samples <- astro_merged_b_liger_cleaned3@cell.data$orig.dataset
write_rds(astro_samples, "RDS_Subcluster_Seurat/astro_sample_names_C.RDS")

# Convert to Seurat
astro_seurat <- ligerToSeurat(astro_merged_b_liger_cleaned3)
write_rds(astro_seurat, "RDS_Subcluster_Seurat/astro_rd3_seurat_C.RDS")