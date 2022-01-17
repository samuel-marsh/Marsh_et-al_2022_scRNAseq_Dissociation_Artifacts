# Load Libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)

install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")

# Load Liger Object -------------------------------------------------------
marsh_post_liger_cleaned <- read_rds("RDS_Objects/marsh_post_liger_round01_cleaned.RDS")

marsh_post_liger_cleaned <- liger::normalize(marsh_post_liger_cleaned)
marsh_post_liger_cleaned <- selectGenes(marsh_post_liger_cleaned, do.plot = T, num.genes = 1000)
marsh_post_liger_cleaned <- scaleNotCenter(marsh_post_liger_cleaned)
marsh_post_liger_cleaned <- online_iNMF(marsh_post_liger_cleaned, k = 35, lambda = 10, miniBatch_size = 5000)
marsh_post_liger_cleaned <- quantile_norm(marsh_post_liger_cleaned, knn_k = 15)
marsh_post_liger_cleaned <- clusterLouvainJaccard(marsh_post_liger_cleaned,resolution = 0.5)
marsh_post_liger_cleaned <- runUMAP(marsh_post_liger_cleaned, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(marsh_post_liger_cleaned, return.plots = TRUE, do.legend = TRUE)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

write_rds(marsh_post_liger_cleaned, "RDS_Objects/marsh_post_liger_cleaned_Round02_clustered.RDS")

ggsave("plots/Round02/marsh_post-mortem_full_liger_round02_UMAP.pdf")

# Check marker genes
marker.genes <- c('OLIG1','CX3CR1','AQP4','TNR','THEMIS','SLC17A7','GAD2','CLDN5', "FLT1", "KCNJ8", "SLC6A12", "COL1A1", "COL1A2", "DCN", "KCNK17", "RBFOX3", "THY1", "SYT1", "RORB", "GAD1", "GAD2", "ADARB2", "LAMP5", "SST", "VIP", "AQP4", "FGFR3", "OPALIN", "COL18A1", "C1QA", "PTPRC")

pdf("plots/Round02/marsh_post-mortem_liger_round02_marker-genes_NEW.pdf")
lapply(marker.genes,function(x){print(plotGene_keep_scale(marsh_post_liger_cleaned, gene = x, plot.by = "none"))})
dev.off()

plotGene_keep_scale(object = marsh_post_liger_cleaned, "PTPRC", plot.by = "none")

# Annotate clusters
marsh_post_mortem_round02_annotation <- tibble::tribble(
  ~cluster, ~cell_type,        ~color,
  0L,    "oligo",      "orange",
  1L,    "oligo",      "orange",
  2L,    "oligo",      "orange",
  3L,    "excit",  "dodgerblue",
  4L,    "oligo",      "orange",
  5L,    "excit",  "dodgerblue",
  6L,    "astro", "forestgreen",
  7L,    "astro", "forestgreen",
  8L,    "micro",        "gold",
  9L,      "OPC",  "darkorange",
  10L,    "inhib",        "navy",
  11L,    "excit",  "dodgerblue",
  12L,    "excit",  "dodgerblue",
  13L,    "inhib",        "navy",
  14L,    "inhib",        "navy",
  15L,    "oligo",      "orange",
  16L,    "inhib",        "navy",
  17L,     "endo",        "pink",
  18L,  "doublet",        "gray",
  19L,  "doublet",        "gray",
  20L,     "PBMC",       "black",
  21L,    "oligo",      "orange",
  22L,    "fibro", "darkorchid3"
)

# Recolor UMAP plot
cluster_color_cleaned <- marsh_post_mortem_round02_annotation %>% 
  pull(color)

umap_dim <- plotByDatasetAndCluster(marsh_post_liger_cleaned, return.plots = TRUE, do.legend = TRUE)
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = cluster_color_cleaned) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

ggsave("plots/Round02/marsh_post-mortem_full_liger_round02_UMAP_colored.pdf")

# Pull cell type cluster IDs
oligo_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "oligo") %>% 
  pull(cluster)

excit_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "excit") %>% 
  pull(cluster)

inhib_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "inhib") %>% 
  pull(cluster)

micro_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "micro") %>% 
  pull(cluster)

astro_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "astro") %>% 
  pull(cluster)

opc_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "OPC") %>% 
  pull(cluster)

endo_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "endo") %>% 
  pull(cluster)

fibro_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "fibro") %>% 
  pull(cluster)

immune_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "PBMC") %>% 
  pull(cluster)

doublet_clu2 <- marsh_post_mortem_round02_annotation %>% 
  filter(cell_type == "doublet") %>% 
  pull(cluster)