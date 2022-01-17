# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Seurat V3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# Load oligoglia subclustered -------------------------------------------
marsh_seurat3_oligo <- read_rds("RDS_SeuratV3/cell_type_subset/marsh_oligo_seurat.RDS")
zhou_seurat3_oligo <- read_rds("RDS_SeuratV3/cell_type_subset/zhou_oligo_seurat.RDS")
morabito_seurat3_oligo <- read_rds("RDS_SeuratV3/cell_type_subset/morabito_oligo_seurat.RDS")
leng_ec_seurat3_oligo <- read_rds("RDS_SeuratV3/cell_type_subset/leng_ec_oligo_seurat.RDS")
leng_sfg_seurat3_oligo_b <- read_rds("RDS_SeuratV3/cell_type_subset/leng_sfg_oligo_seurat_b.RDS")

# merge the cells & save
oligo_merged_b <- merge(x = marsh_seurat3_oligo, y = c(morabito_seurat3_oligo, zhou_seurat3_oligo, leng_ec_seurat3_oligo, leng_sfg_seurat3_oligo_b))

# meta data
marsh_oligo_meta <- oligo_merged_b@meta.data
write_rds(marsh_oligo_meta, "03_Meta Data/marsh_oligo_seurat_meta_b.RDS")

# convert to liger
oligo_merged_b_liger <- seuratToLiger(objects = oligo_merged_b, combined.seurat = T, meta.var = "orig.ident", remove.missing = FALSE)
write_rds(oligo_merged_b_liger, "RDS_subset_merge_liger/oligo/oligo_merged_b_liger_RAW_b.RDS")

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

# round01 -----------------------------------------------------------
oligo_merged_b_liger <- read_rds("RDS_subset_merge_liger/oligo/oligo_merged_b_liger_RAW_b.RDS")
oligo_meta_seurat <- read_rds("03_Meta Data/marsh_oligo_seurat_meta_b.RDS")

oligo_merged_b_liger@cell.data$source <- as.factor(oligo_meta_seurat$Dataset)

oligo_merged_b_liger <- reorganizeLiger(object = oligo_merged_b_liger, by.feature = "source", remove.missing = FALSE)
View(oligo_merged_b_liger@cell.data)

oligo_merged_b_liger <- normalize(oligo_merged_b_liger)
oligo_merged_b_liger <- selectGenes(oligo_merged_b_liger, do.plot = T, num.genes = 1000)
oligo_merged_b_liger <- scaleNotCenter(oligo_merged_b_liger)
oligo_merged_b_liger <- online_iNMF(oligo_merged_b_liger, k = 25, lambda = 5, miniBatch_size = 5000)
oligo_merged_b_liger <- quantile_norm(oligo_merged_b_liger)
oligo_merged_b_liger <- clusterLouvainJaccard(oligo_merged_b_liger,resolution = 1.5)
oligo_merged_b_liger <- runUMAP(oligo_merged_b_liger, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(oligo_merged_b_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6, pt.size = 0.05)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/oligo_merged_b_rd1_factors_b_1000_knn15.pdf")
plotFactors(oligo_merged_b_liger, num.genes = 8, plot.tsne = T)
dev.off()

# Pulling barcodes for cleaning
micro_doubl_oligo <- c("23")
doublet_micro_oligo <- data.frame(oligo_merged_b_liger@clusters, stringsAsFactors = FALSE)
names(doublet_micro_oligo)[1] <- "cluster"
doublet_micro_oligo <- doublet_micro_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% micro_doubl_oligo) %>% 
  pull(barcodes)

# Neurons
neuron_doubl_oligo <- c("7")
doublet_neuron_oligo <- data.frame(oligo_merged_b_liger@clusters, stringsAsFactors = FALSE)
names(doublet_neuron_oligo)[1] <- "cluster"
doublet_neuron_oligo <- doublet_neuron_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% neuron_doubl_oligo) %>% 
  pull(barcodes)

# Astros
astro_doubl_oligo <- c("22")
doublet_astro_oligo <- data.frame(oligo_merged_b_liger@clusters, stringsAsFactors = FALSE)
names(doublet_astro_oligo)[1] <- "cluster"
doublet_astro_oligo <- doublet_astro_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_doubl_oligo) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_oligo <- c(doublet_micro_oligo, doublet_neuron_oligo, doublet_astro_oligo)

# Remove doublets
# Add barcodes to cell data
oligo_merged_b_liger@cell.data$barcodes <- rownames(oligo_merged_b_liger@cell.data)

# Remove doublet cells
oligo_merged_b_liger_cleaned <- subsetLiger(object = oligo_merged_b_liger, cells.use = setdiff(unique(oligo_merged_b_liger@cell.data$barcodes), doublets_complete_oligo), remove.missing = FALSE)

# round02 -----------------------------------------------------------------

oligo_merged_b_liger_cleaned <- normalize(oligo_merged_b_liger_cleaned)
oligo_merged_b_liger_cleaned <- selectGenes(oligo_merged_b_liger_cleaned, do.plot = T, num.genes = 1000)
oligo_merged_b_liger_cleaned <- scaleNotCenter(oligo_merged_b_liger_cleaned)
oligo_merged_b_liger_cleaned <- online_iNMF(oligo_merged_b_liger_cleaned, k = 25, lambda = 5, miniBatch_size = 5000)
oligo_merged_b_liger_cleaned <- quantile_norm(oligo_merged_b_liger_cleaned)
oligo_merged_b_liger_cleaned <- clusterLouvainJaccard(oligo_merged_b_liger_cleaned,resolution = 0.3)
oligo_merged_b_liger_cleaned <- runUMAP(oligo_merged_b_liger_cleaned, n_neighbors = 20)

umap_dim <- plotByDatasetAndCluster(oligo_merged_b_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/oligo_merged_b_rd2_factors_b_1000_k25.pdf")
plotFactors(oligo_merged_b_liger_cleaned, num.genes = 8, plot.tsne = T)
dev.off()

# Save liger objects
write_rds(oligo_merged_b_liger_cleaned, "RDS_subset_merge_liger/oligo_merge_cleaned_b.RDS")
oligo_merged_b_liger_cleaned <- read_rds("RDS_subset_merge_liger/oligo_merge_cleaned.RDS")

oligo_samples <- oligo_merged_b_liger_cleaned@cell.data$orig.dataset
write_rds(oligo_samples, "RDS_Subcluster_Seurat/oligo_sample_names_b.RDS")

# Convert to Seurat
oligo_seurat <- ligerToSeurat(oligo_merged_b_liger_cleaned)
write_rds(oligo_seurat, "RDS_Subcluster_Seurat/oligo_rd2_seurat_b.RDS")