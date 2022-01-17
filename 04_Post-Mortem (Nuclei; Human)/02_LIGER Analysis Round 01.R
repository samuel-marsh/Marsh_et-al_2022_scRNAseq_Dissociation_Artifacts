# Restart R session following QC and install Seurat v2.3.4
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")

# Load Packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

# Liger Round01 -----------------------------------------------------------
marsh_post_liger <- read_rds("RDS_Objects/marsh_post_mortem_seuratV3_filtered_liger.RDS")

marsh_post_liger <- liger::normalize(marsh_post_liger)
marsh_post_liger <- selectGenes(marsh_post_liger, do.plot = T, num.genes = 1200)
marsh_post_liger <- scaleNotCenter(marsh_post_liger)
marsh_post_liger <- online_iNMF(marsh_post_liger, k = 40, lambda = 12, miniBatch_size = 5000)
marsh_post_liger <- quantile_norm(marsh_post_liger)
marsh_post_liger <- clusterLouvainJaccard(marsh_post_liger,resolution = 0.8)
marsh_post_liger <- runUMAP(marsh_post_liger)

umap_dim <- plotByDatasetAndCluster(marsh_post_liger, return.plots = TRUE, do.legend = TRUE)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

marker.genes <- c('OLIG1','CX3CR1','AQP4','TNR','THEMIS','SLC17A7','GAD2','CLDN5', "FLT1", "KCNJ8", "SLC6A12", "COL1A1", "COL1A2", "DCN", "KCNK17", "RBFOX3", "THY1", "SYT1", "RORB", "GAD1", "GAD2", "ADARB2", "LAMP5", "SST", "VIP", "AQP4", "FGFR3", "OPALIN", "COL18A1", "C1QA", "PTPRC")

pdf("plots/Round01/marsh_post-mortem_liger_round01_marker-genes.pdf")
lapply(marker.genes,function(x){print(plotGene_keep_scale(marsh_post_liger, gene = x, plot.by = "none"))})
dev.off()

# Annotate Clusters
marsh_post_mortem_annotation <- tibble::tribble(
  ~cluster,    ~cell_type,        ~color,
  0L,       "oligo",      "orange",
  1L,       "oligo",      "orange",
  2L,       "oligo",      "orange",
  3L,       "astro", "forestgreen",
  4L,       "micro",        "gold",
  5L,       "excit",  "dodgerblue",
  6L,       "oligo",      "orange",
  7L,       "oligo",      "orange",
  8L,       "astro", "forestgreen",
  9L,         "OPC",  "darkorange",
  10L,       "excit",  "dodgerblue",
  11L,       "excit",  "dodgerblue",
  12L,       "inhib",        "navy",
  13L,       "excit",  "dodgerblue",
  14L,       "excit",  "dodgerblue",
  15L,       "excit",  "dodgerblue",
  16L,       "inhib",        "navy",
  17L,       "excit",  "dodgerblue",
  18L,       "inhib",        "navy",
  19L,       "inhib",        "navy",
  20L,  "fibro_peri",        "pink",
  21L,        "endo",        "endo",
  22L,       "oligo",      "orange",
  23L,       "oligo",      "orange",
  24L,       "excit",  "dodgerblue",
  25L,        "PBMC",       "black",
  26L,       "oligo",      "orange",
  27L, "granulocyte",       "black",
  28L,     "doublet",        "gray"
)

# Pull cell type cluster IDs
oligo_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "oligo") %>% 
  pull(cluster)

excit_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "excit") %>% 
  pull(cluster)

inhib_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "inhib") %>% 
  pull(cluster)

micro_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "micro") %>% 
  pull(cluster)

astro_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "astro") %>% 
  pull(cluster)

opc_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "OPC") %>% 
  pull(cluster)

endo_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "endo") %>% 
  pull(cluster)

fibro_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "fibro_peri") %>% 
  pull(cluster)

immune_cells <- c("PBMC", "granulocyte")
immune_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type %in% immune_cells) %>% 
  pull(cluster)

unknown_clu <- marsh_post_mortem_annotation %>% 
  filter(cell_type == "doublet") %>% 
  pull(cluster)