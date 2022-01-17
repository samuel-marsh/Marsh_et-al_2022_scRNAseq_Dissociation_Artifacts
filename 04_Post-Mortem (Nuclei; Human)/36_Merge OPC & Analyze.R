# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Seurat V3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

# Load opcglia subclustered -------------------------------------------
marsh_seurat3_opc <- read_rds("RDS_SeuratV3/cell_type_subset/marsh_opc_seurat.RDS")
zhou_seurat3_opc <- read_rds("RDS_SeuratV3/cell_type_subset/zhou_opc_seurat.RDS")
morabito_seurat3_opc <- read_rds("RDS_SeuratV3/cell_type_subset/morabito_opc_seurat.RDS")
leng_ec_seurat3_opc <- read_rds("RDS_SeuratV3/cell_type_subset/leng_ec_opc_seurat.RDS")
leng_sfg_seurat3_opc <- read_rds("RDS_SeuratV3/cell_type_subset/leng_sfg_opc_seurat.RDS")

# merge the cells & save
opc_merged <- merge(x = marsh_seurat3_opc, y = c(morabito_seurat3_opc, zhou_seurat3_opc, leng_ec_seurat3_opc, leng_sfg_seurat3_opc))

# meta data
marsh_opc_meta <- opc_merged@meta.data
write_rds(marsh_opc_meta, "03_Meta Data/marsh_opc_seurat_meta.RDS")

# Convert to liger and save
opc_merged_liger <- seuratToLiger(objects = opc_merged, combined.seurat = T, meta.var = "orig.ident", remove.missing = FALSE)
write_rds(opc_merged_liger, "RDS_subset_merge_liger/opc/opc_merged_liger_RAW.RDS")

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
opc_merged_liger <- read_rds("RDS_subset_merge_liger/opc/opc_merged_liger_RAW.RDS")
opc_meta_seurat <- read_rds("03_Meta Data/marsh_opc_seurat_meta.RDS")

opc_merged_liger@cell.data$source <- as.factor(opc_meta_seurat$Dataset)

opc_merged_liger <- reorganizeLiger(object = opc_merged_liger, by.feature = "source", remove.missing = FALSE)
View(opc_merged_liger@cell.data)

opc_merged_liger <- normalize(opc_merged_liger)
opc_merged_liger <- selectGenes(opc_merged_liger, do.plot = T, num.genes = 800)
opc_merged_liger <- scaleNotCenter(opc_merged_liger)
opc_merged_liger <- online_iNMF(opc_merged_liger, k = 20, lambda = 5, miniBatch_size = 500)
opc_merged_liger <- quantile_norm(opc_merged_liger)
opc_merged_liger <- clusterLouvainJaccard(opc_merged_liger,resolution = 1.2)
opc_merged_liger <- runUMAP(opc_merged_liger, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(opc_merged_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

opc_merged_liger_mito <- AddMito(opc_merged_liger, species = "human")
View(opc_merged_liger_mito@cell.data)

pdf("05_plots/subclustering/opc_merged.pdf")
plotFactors(opc_merged_liger, num.genes = 8, plot.tsne = T)
dev.off()

# Pulling barcodes for cleaning
oligo_doubl <- c("11")
doublet_oligo <- data.frame(opc_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)

neuron_doubl <- c("1", "3")
doublet_neuron <- data.frame(opc_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_neuron)[1] <- "cluster"
doublet_neuron <- doublet_neuron %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% neuron_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_opc <- c(doublet_oligo, doublet_neuron)

# Remove doublets
# Add barcodes to cell data
opc_merged_liger@cell.data$barcodes <- rownames(opc_merged_liger@cell.data)

# Remove doublet cells
opc_merged_liger_cleaned <- subsetLiger(object = opc_merged_liger, cells.use = setdiff(unique(opc_merged_liger@cell.data$barcodes), doublets_complete_opc), remove.missing = FALSE)

# round02 -----------------------------------------------------------------
opc_merged_liger_cleaned <- normalize(opc_merged_liger_cleaned)
opc_merged_liger_cleaned <- selectGenes(opc_merged_liger_cleaned, do.plot = T, num.genes = 400)
opc_merged_liger_cleaned <- scaleNotCenter(opc_merged_liger_cleaned)
opc_merged_liger_cleaned <- online_iNMF(opc_merged_liger_cleaned, k = 20, lambda = 7, miniBatch_size = 500)
opc_merged_liger_cleaned <- quantile_norm(opc_merged_liger_cleaned)
opc_merged_liger_cleaned <- clusterLouvainJaccard(opc_merged_liger_cleaned,resolution = 0.8)
opc_merged_liger_cleaned <- runUMAP(opc_merged_liger_cleaned, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(opc_merged_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/opc_merged_cleaned.pdf")
plotFactors(opc_merged_liger_cleaned, num.genes = 8, plot.tsne = T)
dev.off()

# Pulling barcodes for cleaning
neuron_doubl2 <- c("9")
doublet_neuron2 <- data.frame(opc_merged_liger_cleaned@clusters, stringsAsFactors = FALSE)
names(doublet_neuron2)[1] <- "cluster"
doublet_neuron2 <- doublet_neuron2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% neuron_doubl2) %>% 
  pull(barcodes)

astro_doubl2 <- c("8")
doublet_astro2 <- data.frame(opc_merged_liger_cleaned@clusters, stringsAsFactors = FALSE)
names(doublet_astro2)[1] <- "cluster"
doublet_astro2 <- doublet_astro2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_doubl2) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_opc2 <- c(doublet_neuron2, doublet_astro2)

# Remove doublet cells
opc_merged_liger_cleaned2 <- subsetLiger(object = opc_merged_liger_cleaned, cells.use = setdiff(unique(opc_merged_liger_cleaned@cell.data$barcodes), doublets_complete_opc2), remove.missing = FALSE)

# round03 -----------------------------------------------------------------
opc_merged_liger_cleaned2 <- normalize(opc_merged_liger_cleaned2)
opc_merged_liger_cleaned2 <- selectGenes(opc_merged_liger_cleaned2, do.plot = T, num.genes = 400)
opc_merged_liger_cleaned2 <- scaleNotCenter(opc_merged_liger_cleaned2)
opc_merged_liger_cleaned2 <- online_iNMF(opc_merged_liger_cleaned2, k = 20, lambda = 7, miniBatch_size = 500)
opc_merged_liger_cleaned2 <- quantile_norm(opc_merged_liger_cleaned2)
opc_merged_liger_cleaned2 <- clusterLouvainJaccard(opc_merged_liger_cleaned2,resolution = 0.8)
opc_merged_liger_cleaned2 <- runUMAP(opc_merged_liger_cleaned2, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(opc_merged_liger_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/opc_merged_cleaned2.pdf")
plotFactors(opc_merged_liger_cleaned2, num.genes = 8, plot.tsne = T)
dev.off()

# save liger objects
write_rds(opc_merged_liger_cleaned2, "RDS_subset_merge_liger/opc_subset_cleaned2.RDS")

opc_samples <- opc_merged_liger_cleaned2@cell.data$orig.dataset
write_rds(opc_samples, "RDS_Subcluster_Seurat/opc_sample_names.RDS")
View(opc_merged_liger_cleaned2@cell.data)

# Convert to Seurat and save
opc_seurat <- ligerToSeurat(opc_merged_liger_cleaned2)
write_rds(opc_seurat, "RDS_Subcluster_Seurat/opc_rd2_seurat.RDS")