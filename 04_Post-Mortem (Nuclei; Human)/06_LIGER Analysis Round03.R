# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")

# Round03 Analysis --------------------------------------------------------
marsh_post_liger_cleaned2 <- liger::normalize(marsh_post_liger_cleaned2)
marsh_post_liger_cleaned2 <- selectGenes(marsh_post_liger_cleaned2, do.plot = T, num.genes = 1000)
marsh_post_liger_cleaned2 <- scaleNotCenter(marsh_post_liger_cleaned2)
marsh_post_liger_cleaned2 <- online_iNMF(marsh_post_liger_cleaned2, k = 20, lambda = 5, miniBatch_size = 5000)
marsh_post_liger_cleaned2 <- quantile_norm(marsh_post_liger_cleaned2)
marsh_post_liger_cleaned2 <- clusterLouvainJaccard(marsh_post_liger_cleaned2,resolution = 0.2)
marsh_post_liger_cleaned2 <- runUMAP(marsh_post_liger_cleaned2, n_neighbors = 40, min_dist = 0.1)

umap_dim <- plotByDatasetAndCluster(marsh_post_liger_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

# remove clustering factor
marsh_post_liger_cleaned2 <- quantile_norm(marsh_post_liger_cleaned2, dims.use = setdiff(1:20, 19))
marsh_post_liger_cleaned2 <- clusterLouvainJaccard(marsh_post_liger_cleaned2,resolution = 0.2)
marsh_post_liger_cleaned2 <- runUMAP(marsh_post_liger_cleaned2, n_neighbors = 45, min_dist = 0.2)

umap_dim <- plotByDatasetAndCluster(marsh_post_liger_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

# Annotate
post_final_annotation <- tibble::tribble(
  ~cluster, ~cell_type,        ~color,
  0L,    "oligo",      "orange",
  1L,    "oligo",      "orange",
  2L,    "oligo",      "orange",
  3L,    "excit",  "dodgerblue",
  4L,    "excit",  "dodgerblue",
  5L,    "astro", "forestgreen",
  6L,    "micro",        "gold",
  7L,    "astro", "forestgreen",
  8L,    "inhib",        "navy",
  9L,    "inhib",        "navy",
  10L,      "OPC", "darkorange2",
  11L,     "endo",        "orchid",
  12L,     "PBMC",        "gray",
  13L,    "fibro", "darkorchid3"
)

cluster_colors_final <- post_final_annotation %>% 
  pull(color)

umap_dim <- plotByDatasetAndCluster(marsh_post_liger_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = cluster_colors_final) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

# Plot Marker genes
marker.genes <- c('OLIG1','CX3CR1','AQP4','TNR','THEMIS','SLC17A7','GAD2','CLDN5', "FLT1", "KCNJ8", "SLC6A12", "COL1A1", "COL1A2", "DCN", "KCNK17", "RBFOX3", "THY1", "SYT1", "RORB", "GAD1", "GAD2", "ADARB2", "LAMP5", "SST", "VIP", "AQP4", "FGFR3", "OPALIN", "COL18A1", "C1QA", "PTPRC")

pdf("plots/Round02/marsh_post-mortem_liger_round02_marker-genes_NEW.pdf")
lapply(marker.genes,function(x){print(plotGene_keep_scale(marsh_post_liger_cleaned2, gene = x, plot.by = "none"))})
dev.off()

# Convert to Seurat  ------------------------------------------------------
# Object will be in V2 structure
marsh_post_seurat <- ligerToSeurat(marsh_post_liger_cleaned2)
write_rds(marsh_post_seurat, "RDS_Objects/marsh_post_seurat_round03_FINAL_tempv2.RDS")

# Restart R Session and install Seurat V3 ---------------------------------
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_3.1.5.tar.gz", repos = NULL, type = "source")

# Load Packages 
library(tidyverse)
library(Seurat) # Starting with Seurat V3.1.5
library(patchwork)
library(marsh_post.utils)
library(viridis)
library(liger)
library(scCustomize)
library(beepr)

# Read in seurat v2 object
marsh_post_seurat <- read_rds("RDS_Objects/marsh_post_seurat_round03_FINAL_tempv2.RDS")
marsh_post_seurat_v3 <- UpdateSeuratObject(marsh_post_seurat)

write_rds(marsh_post_seurat_v3, "RDS_SeuratV3/marsh_post_seuratv3_final.RDS")