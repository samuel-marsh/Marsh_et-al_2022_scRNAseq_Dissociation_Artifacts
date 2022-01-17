# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat) # Seurat V3.1.5
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

# Load excitglia subclustered -------------------------------------------
marsh_seurat3_excit <- read_rds("RDS_SeuratV3/cell_type_subset/marsh_excit_seurat.RDS")
zhou_seurat3_excit <- read_rds("RDS_SeuratV3/cell_type_subset/zhou_excit_seurat.RDS")
morabito_seurat3_excit <- read_rds("RDS_SeuratV3/cell_type_subset/morabito_excit_seurat.RDS")
leng_ec_seurat3_excit <- read_rds("RDS_SeuratV3/cell_type_subset/leng_ec_excit_seurat.RDS")
leng_sfg_seurat3_excit <- read_rds("RDS_SeuratV3/cell_type_subset/leng_sfg_excit_seurat.RDS")

# merge the cells & save
excit_merged <- merge(x = marsh_seurat3_excit, y = c(morabito_seurat3_excit, zhou_seurat3_excit, leng_ec_seurat3_excit, leng_sfg_seurat3_excit))

# meta data
marsh_excit_meta <- excit_merged@meta.data
write_rds(marsh_excit_meta, "03_Meta Data/marsh_excit_seurat_meta.RDS")

# convert to liger and save
excit_merged_liger <- seuratToLiger(objects = excit_merged, combined.seurat = T, meta.var = "orig.ident", remove.missing = FALSE)
write_rds(excit_merged_liger, "RDS_subset_merge_liger/excit/excit_merged_liger_RAW.RDS")

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
excit_merged_liger <- read_rds("RDS_subset_merge_liger/excit/excit_merged_liger_RAW.RDS")
excit_meta_seurat <- read_rds("03_Meta Data/marsh_excit_seurat_meta.RDS")

excit_merged_liger@cell.data$source <- as.factor(excit_meta_seurat$Dataset)

excit_merged_liger <- reorganizeLiger(object = excit_merged_liger, by.feature = "source", remove.missing = FALSE)

excit_merged_liger <- normalize(excit_merged_liger)
excit_merged_liger <- selectGenes(excit_merged_liger, do.plot = T, num.genes = 600)
excit_merged_liger <- scaleNotCenter(excit_merged_liger)
excit_merged_liger <- online_iNMF(excit_merged_liger, k = 25, lambda = 7, miniBatch_size = 2500)
excit_merged_liger <- quantile_norm(excit_merged_liger)
excit_merged_liger <- clusterLouvainJaccard(excit_merged_liger,resolution = 1.2)
excit_merged_liger <- runUMAP(excit_merged_liger, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(excit_merged_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/excit_merged.pdf")
plotFactors(excit_merged_liger, num.genes = 8, plot.tsne = T)
dev.off()

# Pulling barcodes for cleaning
oligo_doubl <- c("29")
doublet_oligo <- data.frame(excit_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)

astro_endo_doubl <- c("25")
doublet_astro_endo <- data.frame(excit_merged_liger@clusters, stringsAsFactors = FALSE)
names(doublet_astro_endo)[1] <- "cluster"
doublet_astro_endo <- doublet_astro_endo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_endo_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublets_complete_excit <- c(doublet_oligo, doublet_astro_endo)

# Remove doublets
# Add barcodes to cell data
excit_merged_liger@cell.data$barcodes <- rownames(excit_merged_liger@cell.data)

# Remove doublet cells
excit_merged_liger_cleaned <- subsetLiger(object = excit_merged_liger, cells.use = setdiff(unique(excit_merged_liger@cell.data$barcodes), doublets_complete_excit), remove.missing = FALSE)

# round02 -----------------------------------------------------------------
excit_merged_liger_cleaned <- normalize(excit_merged_liger_cleaned)
excit_merged_liger_cleaned <- selectGenes(excit_merged_liger_cleaned, do.plot = T, num.genes = 600)
excit_merged_liger_cleaned <- scaleNotCenter(excit_merged_liger_cleaned)
excit_merged_liger_cleaned <- online_iNMF(excit_merged_liger_cleaned, k = 25, lambda = 7, miniBatch_size = 2500)
excit_merged_liger_cleaned <- quantile_norm(excit_merged_liger_cleaned)
excit_merged_liger_cleaned <- clusterLouvainJaccard(excit_merged_liger_cleaned,resolution = 1.2)
excit_merged_liger_cleaned <- runUMAP(excit_merged_liger_cleaned, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(excit_merged_liger_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

# astrocytes and endo
astro_doubl2 <- c("28")
doublet_astro2 <- data.frame(excit_merged_liger_cleaned@clusters, stringsAsFactors = FALSE)
names(doublet_astro2)[1] <- "cluster"
doublet_astro2 <- doublet_astro2 %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% astro_doubl2) %>%
  pull(barcodes)

# astrocytes and endo
micro_doubl2 <- c("29")
doublet_micro2 <- data.frame(excit_merged_liger_cleaned@clusters, stringsAsFactors = FALSE)
names(doublet_micro2)[1] <- "cluster"
doublet_micro2 <- doublet_micro2 %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% micro_doubl2) %>%
  pull(barcodes)

# astrocytes and endo
inhib_doubl2 <- c("27")
doublet_inhib2 <- data.frame(excit_merged_liger_cleaned@clusters, stringsAsFactors = FALSE)
names(doublet_inhib2)[1] <- "cluster"
doublet_inhib2 <- doublet_inhib2 %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster %in% inhib_doubl2) %>%
  pull(barcodes)

#  Compile Doublets and Clean ----------------------------------------------
doublets_complete_excit2 <- c(doublet_astro2, doublet_micro2, doublet_inhib2)

# Remove doublet cells
excit_merged_liger_cleaned2 <- subsetLiger(object = excit_merged_liger_cleaned, cells.use = setdiff(unique(excit_merged_liger_cleaned@cell.data$barcodes), doublets_complete_excit2), remove.missing = FALSE)

# round03 -----------------------------------------------------------------
excit_merged_liger_cleaned2 <- normalize(excit_merged_liger_cleaned2)
excit_merged_liger_cleaned2 <- selectGenes(excit_merged_liger_cleaned2, do.plot = T, num.genes = 500)
excit_merged_liger_cleaned2 <- scaleNotCenter(excit_merged_liger_cleaned2)
excit_merged_liger_cleaned2 <- online_iNMF(excit_merged_liger_cleaned2, k = 25, lambda = 7, miniBatch_size = 2500)
excit_merged_liger_cleaned2 <- quantile_norm(excit_merged_liger_cleaned2)
excit_merged_liger_cleaned2 <- clusterLouvainJaccard(excit_merged_liger_cleaned2,resolution = 1.2)
excit_merged_liger_cleaned2 <- runUMAP(excit_merged_liger_cleaned2, n_neighbors = 30)

umap_dim <- plotByDatasetAndCluster(excit_merged_liger_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("05_plots/subclustering/excit_merged_cleaned2.pdf")
plotFactors(excit_merged_liger_cleaned2, num.genes = 8, plot.tsne = T)
dev.off()

# save gene list
colnames(excit_rd2_factor9) <- "excit_rd2_factor9"
write_rds(excit_rd2_factor9, "subcluster_factor_gene_list/liger_excit_rd2_gene_list.RDS")

# save liger object
write_rds(excit_merged_liger_cleaned2, "RDS_subset_merge_liger/excit_subset_cleaned2.RDS")

excit_samples <- excit_merged_liger_cleaned2@cell.data$orig.dataset
write_rds(excit_samples, "RDS_Subcluster_Seurat/excit_sample_names.RDS")

# Convert to Seurat and Save
excit_seurat <- ligerToSeurat(excit_merged_liger_cleaned2)
write_rds(excit_seurat, "RDS_Subcluster_Seurat/excit_rd3_seurat.RDS")