# Start with Seurat v3.1.5
# Analysis performed with R3.5.3

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_3.1.5.tar.gz", repos = NULL, type = "source")
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")

# Load Data ---------------------------------------------------------------
# Read in full meta data table
snTrem2_meta <- read.csv("Zhou_snTrem2/snTrem2_Meta_Data.csv", stringsAsFactors = FALSE, header = 1)

# Filter for AD & Trem2 samples
snTrem2_ctrl_meta <- snTrem2_meta %>% 
  filter(Sample_Source == "Control - Rush")

snTrem2_AD_meta <- snTrem2_meta %>% 
  filter(Sample_Source == "AD - Rush")

snTrem2_TREM2_meta <- snTrem2_meta %>% 
  filter(Sample_Source == "TREM2 - R62H - Rush")

# Create Full meta --------------------------------------------------------
rosmap <- read.csv("ROSMAP_Clinical_2019-05_v3.csv", stringsAsFactors = FALSE)

sample_id_ctrl <- snTrem2_ctrl_meta %>% 
  pull(Sample_Identifier)

sample_id_ad <- snTrem2_AD_meta %>% 
  pull(Sample_Identifier)

rosmap_ctrl <- rosmap %>% 
  filter(projid %in% sample_id_ctrl)

rosmap_AD <- rosmap %>% 
  filter(projid %in% sample_id_ad)

rosmap_full_all_samples <- full_join(rosmap_ctrl, rosmap_AD)

write.csv(rosmap_full_all_samples, "ROSMAP_all_meta_filtered.csv")

# Create Ctrl & AD Dataset ------------------------------------------------
# Read in full meta data table
snTrem2_meta <- read.csv("Zhou_snTrem2/snTrem2_Meta_Data.csv", stringsAsFactors = FALSE, header = 1)

samples <- c("Control - Rush", "AD - Rush")

# Filter for AD & Trem2 samples
snTrem2_ctrl_meta <- snTrem2_meta %>% 
  filter(Sample_Source %in% samples)

work.dir <- "~/raw_data/sn_Trem2_10X/"
lib.list <- snTrem2_ctrl_meta %>% 
  pull(Sample_ID)
lib.names <- snTrem2_ctrl_meta %>% 
  pull(Sample_ID)
dge.full <- lapply(c(1:length(lib.list)), function(x){
  dge.1 <- readMM(paste0(work.dir,lib.list[x], '_matrix.mtx.gz') )
  bx2 <- read.table(paste0(work.dir,lib.list[x],'_barcodes.tsv.gz'))
  genes2 <- read.table(paste0(work.dir,lib.list[x],'_features.tsv.gz'))
  rownames(dge.1) <- genes2$V2
  bx.use <- as.vector(unlist(lapply(bx2$V1,function(y){paste0(lib.names[x],"_",y)})))
  colnames(dge.1) <- bx.use
  gc()
  print(paste0('Returning ',as.character(x)))
  return(dge.1)
})
names(dge.full) <- lib.names

A1 <- CreateSeuratObject(counts = dge.full[[1]], min.cells = 1, min.features = 200, lib.names[[1]])
A10 <- CreateSeuratObject(counts = dge.full[[2]], min.cells = 1, min.features = 200, lib.names[[2]])
A11 <- CreateSeuratObject(counts = dge.full[[3]], min.cells = 1, min.features = 200, lib.names[[3]])
A12 <- CreateSeuratObject(counts = dge.full[[4]], min.cells = 1, min.features = 200, lib.names[[4]])
A13 <- CreateSeuratObject(counts = dge.full[[5]], min.cells = 1, min.features = 200, lib.names[[5]])
A2 <- CreateSeuratObject(counts = dge.full[[6]], min.cells = 1, min.features = 200, lib.names[[6]])
A3 <- CreateSeuratObject(counts = dge.full[[7]], min.cells = 1, min.features = 200, lib.names[[7]])
A5 <- CreateSeuratObject(counts = dge.full[[8]], min.cells = 1, min.features = 200, lib.names[[8]])
A7 <- CreateSeuratObject(counts = dge.full[[9]], min.cells = 1, min.features = 200, lib.names[[9]])
A8 <- CreateSeuratObject(counts = dge.full[[10]], min.cells = 1, min.features = 200, lib.names[[10]])
A9 <- CreateSeuratObject(counts = dge.full[[11]], min.cells = 1, min.features = 200, lib.names[[11]])
C1 <- CreateSeuratObject(counts = dge.full[[12]], min.cells = 1, min.features = 200, lib.names[[12]])
C11 <- CreateSeuratObject(counts = dge.full[[13]], min.cells = 1, min.features = 200, lib.names[[13]])
C12 <- CreateSeuratObject(counts = dge.full[[14]], min.cells = 1, min.features = 200, lib.names[[14]])
C2 <- CreateSeuratObject(counts = dge.full[[15]], min.cells = 1, min.features = 200, lib.names[[15]])
C3 <- CreateSeuratObject(counts = dge.full[[16]], min.cells = 1, min.features = 200, lib.names[[16]])
C4 <- CreateSeuratObject(counts = dge.full[[17]], min.cells = 1, min.features = 200, lib.names[[17]])
C5 <- CreateSeuratObject(counts = dge.full[[18]], min.cells = 1, min.features = 200, lib.names[[18]])
C6 <- CreateSeuratObject(counts = dge.full[[19]], min.cells = 1, min.features = 200, lib.names[[19]])
C7 <- CreateSeuratObject(counts = dge.full[[20]], min.cells = 1, min.features = 200, lib.names[[20]])
C8 <- CreateSeuratObject(counts = dge.full[[21]], min.cells = 1, min.features = 200, lib.names[[21]])
C9 <- CreateSeuratObject(counts = dge.full[[22]], min.cells = 1, min.features = 200, lib.names[[22]])

# Create unified object
snTrem2_ctrlAD_seurat <- merge(x = C1, y = c(C2, C3, C4, C5, C6, C7, C8, C9, C11, C12, A1, A2, A3, A5, A7, A8, A9, A10, A11, A12, A13))

stats <- Cluster_Stats_All_Samples(snTrem2_ctrlAD_seurat)

# Add mito
snTrem2_ctrlAD_seurat <- PercentageFeatureSet(object = snTrem2_ctrlAD_seurat, pattern = "^MT", col.name = "percent_mito")

# QC Filter ---------------------------------------------------------------
QC_Plots_UMI(snTrem2_ctrlAD_seurat)

snTrem2_ctrlAD_seurat <- subset(x = snTrem2_ctrlAD_seurat, subset = nCount_RNA < 50000)

QC_Plots_UMI(snTrem2_ctrlAD_seurat)
QC_Plots_Genes(snTrem2_ctrlAD_seurat, low_cutoff = 300, high_cutoff = 9000) + ylim(0, 600)
QC_Plots_Mito(snTrem2_ctrlAD_seurat)

snTrem2_ctrlAD_seurat50 <- subset(x = snTrem2_ctrlAD_seurat, subset = percent_mito < 50)

# Save objects
write_rds(snTrem2_ctrlAD_seurat50, "RDS_Objects/RAW_SeuratV3_Filtered/snTrem2_ctrlAD_Seurat50_RAW.RDS")

# REMOVE DONOR C4
# ADD RATIONALE NOTES
snTrem2_ctrlAD_seurat50 <- subset(x = snTrem2_ctrlAD_seurat50, idents = "C4", invert = TRUE)

snTrem2_ctrlAD_liger50 <- seuratToLiger(objects = snTrem2_ctrlAD_seurat50, combined.seurat = TRUE, meta.var = "orig.ident")
write_rds(snTrem2_ctrlAD_liger50, "RDS_Objects/RAW_liger/snTrem2_ctrl_liger50_RAW.RDS")

# Restart R session and Install Seurat v2.3.4 -----------------------------
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_2.3.4.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)

# LIGER Round01 -----------------------------------------------------------
ctrl_AD_comb <- read_rds("RDS_Objects/RAW_liger/snTrem2_ctrl_liger50_RAW.RDS")
ctrl_AD_comb <-  normalize(ctrl_AD_comb)
ctrl_AD_comb <-  selectGenes(ctrl_AD_comb, do.plot = T, num.genes = 800)
ctrl_AD_comb <-  scaleNotCenter(ctrl_AD_comb)
ctrl_AD_comb <-  online_iNMF(ctrl_AD_comb, k = 40, lambda = 10, miniBatch_size = 900)
ctrl_AD_comb <-  quantile_norm(ctrl_AD_comb)
ctrl_AD_comb <-  clusterLouvainJaccard(ctrl_AD_comb,resolution = 0.7)
ctrl_AD_comb <-  runUMAP(ctrl_AD_comb)

umap_dim <- plotByDatasetAndCluster(ctrl_AD_comb, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

write_rds(ctrl_AD_comb, "RDS_Objects/Round01_Liger/snTrem2_ctrlAD_liger_ROUND01.RDS")

ggsave("plots/Round01_Clustering/snTrem2_ctrlAD_liger_round01_UMAP.pdf")

marker.genes <- c('OLIG1','CX3CR1','AQP4','TNR','THEMIS','SLC17A7','GAD2','CLDN5', "FLT1", "KCNJ8", "COL1A1", "RBFOX3", "THY1", "SYT1", "GAD1", "GAD2", "ADARB2", "SST", "FGFR3", "OPALIN", "COL18A1", "C1QA")

pdf("plots/Round01_Markers/snTrem2_ctrlAD_liger_round01_marker-genes.pdf")
lapply(marker.genes,function(x){print(plotGene_keep_scale(ctrl_AD_comb, gene = x, plot.by = "none"))})
dev.off()

pdf("plots/Round01_Factors/snTrem2_ctrlAD_liger_round01_factors.pdf")
plotFactors(liger_object = ctrl_AD_comb, num.genes = 8, plot.tsne = T)
dev.off()

# Cleaning Round01 --------------------------------------------------------
# Annotate
ctrl_AD_comb_annotation <- tibble::tribble(
  ~cluster, ~cell_type,        ~color,
  0L,    "oligo",      "orange",
  1L,    "oligo",      "orange",
  2L,    "excit",  "dodgerblue",
  3L,    "excit",  "dodgerblue",
  4L,    "oligo",      "orange",
  5L,    "excit",  "dodgerblue",
  6L,    "excit",  "dodgerblue",
  7L,    "inhib",        "navy",
  8L,    "excit",  "dodgerblue",
  9L,    "micro",        "gold",
  10L,    "astro", "forestgreen",
  11L,    "excit",  "dodgerblue",
  12L,    "astro", "forestgreen",
  13L,    "astro", "forestgreen",
  14L,    "inhib",        "navy",
  15L,      "OPC",  "darkorange",
  16L,    "excit",  "dodgerblue",
  17L,      "OPC",  "darkorange",
  18L,     "endo",        "pink",
  19L,    "inhib",        "navy",
  20L,    "excit",  "dodgerblue"
)

ctrl_AD_comb_cluster_colors <- ctrl_AD_comb_annotation %>% 
  pull(color)

ctrl_AD_comb_dim <- plotByDatasetAndCluster(ctrl_AD_comb, return.plots = TRUE, do.legend = TRUE, text.size = 4)
ctrl_AD_comb_dim <- ctrl_AD_comb_dim[[2]] + 
  scale_color_manual(values = ctrl_AD_comb_cluster_colors) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
ctrl_AD_comb_dim

ggsave("plots/Round01_Clustering/ctrl_AD_comb_liger_round01_colored.pdf")

# oligo -------------------------------------------------------------------
oligo_sub <- subsetLiger(object = ctrl_AD_comb, clusters.use = c("0", "1", "4"))

# Reanalyze
oligo_sub <-  normalize(oligo_sub)
oligo_sub <-  selectGenes(oligo_sub, do.plot = T, num.genes = 100)
oligo_sub <-  scaleNotCenter(oligo_sub)
#oligo_sub <- optimizeALS(object = oligo_sub, k = 25, lambda = 5, nrep = 3)
oligo_sub <-  online_iNMF(object = oligo_sub, k = 25, lambda = 5, miniBatch_size = 200)
oligo_sub <-  quantile_norm(oligo_sub)
oligo_sub <-  clusterLouvainJaccard(oligo_sub,resolution = 0.7)
oligo_sub <-  runUMAP(oligo_sub)

oligo_dim <- plotByDatasetAndCluster(oligo_sub, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "polychrome")
oligo_dim <- oligo_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
oligo_dim

pdf("plots/Cleaning/snTrem2_ctrlAD_cleaning/oligo_sub_factors2.pdf")
plotFactors(liger_object = oligo_sub, num.genes = 8, plot.tsne = TRUE)
dev.off()

plotGene_keep_scale(oligo_sub, "OPALIN", plot.by = "none")
plotGene_keep_scale(oligo_sub, "CLU", plot.by = "none")
plotGene_keep_scale(oligo_sub, "ZCCHC11", plot.by = "none")

seurat_sub <- ligerToSeurat(oligo_sub)
FeaturePlot(seurat_sub, "NR")

clu16 <- prestowrapper(seurat_sub, ident.1 = 16)

# Remove Doublets
oligo_doubl <- c("13", "15")
doublet_oligo <- data.frame(oligo_sub@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)

# excit -------------------------------------------------------------------
# Subcluster and run pipeline again
excit_sub <- subsetLiger(object = ctrl_AD_comb, clusters.use = c("2", "3", "5", "6", "8", "11", "16", "20"))

# Reanalyze
excit_sub <-  normalize(excit_sub)
excit_sub <-  selectGenes(excit_sub, do.plot = T, num.genes = 250)
excit_sub <-  scaleNotCenter(excit_sub)
excit_sub <- optimizeALS(object = excit_sub, k = 25, lambda = 5, nrep = 3)
#excit_sub <-  online_iNMF(object = excit_sub, k = 20, lambda = 5, miniBatch_size = 300)
excit_sub <-  quantile_norm(excit_sub)
excit_sub <-  clusterLouvainJaccard(excit_sub,resolution = 0.5)
excit_sub <-  runUMAP(excit_sub)

excit_dim <- plotByDatasetAndCluster(excit_sub, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
excit_dim <- excit_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
excit_dim

pdf("plots/Cleaning/snTrem2_ctrlAD_cleaning/excit_sub_factors.pdf")
plotFactors(liger_object = excit_sub, num.genes = 8, plot.tsne = TRUE)
dev.off()

seurat_sub <- ligerToSeurat(excit_sub)
clu14 <- prestowrapper(seurat_sub, ident.1 = "14")

plotGene_keep_scale(excit_sub, "GFAP", plot.by = "none")
plotGene_keep_scale(excit_sub, "DCN", plot.by = "none")
plotGene_keep_scale(excit_sub, "ST18", plot.by = "none")
plotGene_keep_scale(excit_sub, "CST3", plot.by = "none")

# Clu12, 4, 9, 8, 10, 11, 15, 18, 17 are doublets
excit_clu_doublets <- c("12", "4", "9", "8", "10", "11", "15", "18", "17")

doublet_excit <- data.frame(excit_sub@clusters, stringsAsFactors = FALSE)
names(doublet_excit)[1] <- "cluster"
doublet_excit <- doublet_excit %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% excit_clu_doublets) %>% 
  pull(barcodes)

# micro -------------------------------------------------------------------
# Subcluster and run pipeline again
micro_sub <- subsetLiger(object = ctrl_AD_comb, clusters.use = "9")

# Reanalyze
micro_sub <-  normalize(micro_sub)
micro_sub <-  selectGenes(micro_sub, do.plot = T, num.genes = 300)
micro_sub <-  scaleNotCenter(micro_sub)
micro_sub <- optimizeALS(object = micro_sub, k = 14, lambda = 5, nrep = 3)
#micro_sub <-  online_iNMF(object = micro_sub, k = 14, lambda = 5, miniBatch_size = 10)
micro_sub <-  quantile_norm(micro_sub,knn_k = 14)
micro_sub <-  clusterLouvainJaccard(micro_sub,resolution = 1)
micro_sub <-  runUMAP(micro_sub, min_dist = 0.3, n_neighbors = 15)

micro_dim <- plotByDatasetAndCluster(micro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "polychrome")
micro_dim <- micro_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
micro_dim

pdf("plots/Cleaning/snTrem2_ctrlAD_cleaning/micro_sub_factors.pdf")
plotFactors(liger_object = micro_sub, num.genes = 8, plot.tsne = TRUE)
dev.off()

plotGene_keep_scale(micro_sub, "FAM107A", plot.by = "none")
plotGene_keep_scale(micro_sub, "CX3CR1", plot.by = "none")
plotGene_keep_scale(micro_sub, "TMEM119", plot.by = "none")
plotGene_keep_scale(micro_sub, "SPP1", plot.by = "none")

seurat_sub <- ligerToSeurat(micro_sub)
FeaturePlot(object = seurat_sub, features.plot = "S100A8")

clu7_markers <- prestowrapper(object = seurat_sub, ident.1 = 7)
# Likely PBMCs Neutrophils

# Clu7 are astro/neuron
micro_doubl <- c("5", "8")
doublet_micro <- data.frame(micro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_micro)[1] <- "cluster"
doublet_micro <- doublet_micro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% micro_doubl) %>% 
  pull(barcodes)

# Inhib -------------------------------------------------------------------
# Subcluster and run pipeline again
inhib_sub <- subsetLiger(object = ctrl_AD_comb, clusters.use = c("7", "14", "19"))

# AD1 and AD8 have low cell numbers and pipeline will not run with them included.
AD1_AD8 <- c("AD1", "AD8")
rmv_AD8_AD1 <- data.frame(inhib_sub@cell.data, stringsAsFactors = FALSE)
rmv_AD8_AD1 <- rmv_AD8_AD1 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(dataset %in% AD1_AD8) %>% 
  pull(barcodes)

# Add barcodes to cell data
inhib_sub@cell.data$barcodes <- rownames(inhib_sub@cell.data)

# Remove doublet cells
inhib_sub_minus_AD1_AD8 <- subsetLiger(object = inhib_sub, cells.use = setdiff(unique(inhib_sub@cell.data$barcodes), rmv_AD8_AD1))

# Reanalyze
inhib_sub_minus_AD1_AD8 <-  normalize(inhib_sub_minus_AD1_AD8)
inhib_sub_minus_AD1_AD8 <-  selectGenes(inhib_sub_minus_AD1_AD8, do.plot = T, num.genes = 250)
inhib_sub_minus_AD1_AD8 <-  scaleNotCenter(inhib_sub_minus_AD1_AD8)
inhib_sub_minus_AD1_AD8 <- optimizeALS(object = inhib_sub_minus_AD1_AD8, k = 25, lambda = 5, nrep = 3)
#inhib_sub_minus_AD8 <-  online_iNMF(object = inhib_sub_minus_AD8, k = 25, lambda = 5, miniBatch_size = 50)
inhib_sub_minus_AD1_AD8 <-  quantile_norm(inhib_sub_minus_AD1_AD8)
inhib_sub_minus_AD1_AD8 <-  clusterLouvainJaccard(inhib_sub_minus_AD1_AD8,resolution = 1.2)
inhib_sub_minus_AD1_AD8 <-  runUMAP(inhib_sub_minus_AD1_AD8, min_dist = 0.3, n_neighbors = 15)

inhib_dim <- plotByDatasetAndCluster(inhib_sub_minus_AD1_AD8, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
inhib_dim <- inhib_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
inhib_dim

pdf("plots/Cleaning/snTrem2_ctrlAD_cleaning/inhib_sub_factors.pdf")
plotFactors(liger_object = inhib_sub_minus_AD1_AD8, num.genes = 8, plot.tsne = TRUE)
dev.off()

plotGene_keep_scale(object = inhib_sub_minus_AD1_AD8, "COL18A1", plot.by = "none")
plotGene_keep_scale(object = inhib_sub_minus_AD1_AD8, "GAD2", plot.by = "none")
plotGene_keep_scale(object = inhib_sub_minus_AD1_AD8, "FGFR3", plot.by = "none")

seurat_sub <- ligerToSeurat(inhib_sub_minus_AD1_AD8)
clu21 <- prestowrapper(seurat_sub, ident.1 = "21")

# Clu10 are oligo
inhib_doubl <- c("19", "16", "17", "14", "18")
doublet_inhib <- data.frame(inhib_sub_minus_AD1_AD8@clusters, stringsAsFactors = FALSE)
names(doublet_inhib)[1] <- "cluster"
doublet_inhib <- doublet_inhib %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% inhib_doubl) %>% 
  pull(barcodes)

# Astro -------------------------------------------------------------------
# Subcluster and run pipeline again
astro_sub <- subsetLiger(object = ctrl_AD_comb, clusters.use = c("10", "12", "13"))

# Reanalyze
astro_sub <-  normalize(astro_sub)
astro_sub <-  selectGenes(astro_sub, do.plot = T, num.genes = 300)
astro_sub <-  scaleNotCenter(astro_sub)
astro_sub <-  online_iNMF(object = astro_sub, k = 25, lambda = 5, miniBatch_size = 100)
astro_sub <-  quantile_norm(astro_sub)
astro_sub <-  clusterLouvainJaccard(astro_sub,resolution = 1.2)
astro_sub <-  runUMAP(astro_sub, min_dist = 0.3, n_neighbors = 15)

astro_dim <- plotByDatasetAndCluster(astro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
astro_dim <- astro_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
astro_dim

pdf("plots/Cleaning/snTrem2_ctrlAD_cleaning/astro_sub_factors.pdf")
lotFactors(liger_object = astro_sub, num.genes = 8, plot.tsne = TRUE)
dev.off()

plotGene_keep_scale(object = astro_sub, gene = "TNR", plot.by = "none")
plotGene_keep_scale(object = astro_sub, gene = "PRKCB", plot.by = "none")
plotGene_keep_scale(object = astro_sub, gene = "SLC17A7", plot.by = "none")

seurat_sub <- ligerToSeurat(astro_sub)

clu8 <- prestowrapper(seurat_sub, ident.1 = 8)

# Clu10 are oligo
astro_doubl <- c("18", "8", "17")
doublet_astro <- data.frame(astro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_astro)[1] <- "cluster"
doublet_astro <- doublet_astro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_doubl) %>% 
  pull(barcodes)

# OPC -------------------------------------------------------------------
# Subcluster and run pipeline again
opc_sub <- subsetLiger(object = ctrl_AD_comb, clusters.use = c("15", "17"))

# Reanalyze
opc_sub <-  normalize(opc_sub)
opc_sub <-  selectGenes(opc_sub, do.plot = T, num.genes = 200)
opc_sub <-  scaleNotCenter(opc_sub)
opc_sub <- optimizeALS(object = opc_sub, k = 25, lambda = 5, nrep = 3)
#opc_sub <-  online_iNMF(object = opc_sub, k = 25, lambda = 5, miniBatch_size = 25)
opc_sub <-  quantile_norm(opc_sub, knn_k = 15)
opc_sub <-  clusterLouvainJaccard(opc_sub,resolution = 0.8)
opc_sub <-  runUMAP(opc_sub, min_dist = 0.3, n_neighbors = 15)

opc_sub_dim <- plotByDatasetAndCluster(opc_sub, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
opc_sub_dim <- opc_sub_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
opc_sub_dim

pdf("plots/Cleaning/snTrem2_ctrlAD_cleaning/opc_sub_factors.pdf")
plotFactors(liger_object = opc_sub, num.genes = 8, plot.tsne = TRUE)
dev.off()

plotGene_keep_scale(object = opc_sub, "AQP4", plot.by = "none")

seurat_sub <- ligerToSeurat(object = opc_sub)
FeaturePlot(seurat_sub, "GPNMB")

# # Clu5 are oligos
doublet_opc_clu5_oligo <- data.frame(opc_sub@clusters, stringsAsFactors = FALSE)
names(doublet_opc_clu5_oligo)[1] <- "cluster"
doublet_opc_clu5_oligo <- doublet_opc_clu5_oligo %>%
  rownames_to_column(var = "barcodes") %>%
  filter(cluster == "5") %>%
  pull(barcodes)

# Endo -------------------------------------------------------------------
# Attempt without cleaning for now.  Cell number very low

# Compile and Clean -------------------------------------------------------
# compile doublets
doublets_complete <- c(doublet_astro, doublet_excit, doublet_inhib, doublet_micro, doublet_oligo)

# Add barcodes to cell data
ctrl_AD_comb@cell.data$barcodes <- rownames(ctrl_AD_comb@cell.data)

# Remove doublet cells
ctrl_AD_comb_cleaned <- subsetLiger(object = ctrl_AD_comb, cells.use = setdiff(unique(ctrl_AD_comb@cell.data$barcodes), doublets_complete))

# LIGER Round02 -----------------------------------------------------------
# Reanalyze and Cluster
ctrl_AD_comb_cleaned <-  normalize(ctrl_AD_comb_cleaned)
ctrl_AD_comb_cleaned <-  selectGenes(ctrl_AD_comb_cleaned, do.plot = T, num.genes = 800)
ctrl_AD_comb_cleaned <-  scaleNotCenter(ctrl_AD_comb_cleaned)
ctrl_AD_comb_cleaned <-  online_iNMF(ctrl_AD_comb_cleaned, k = 40, lambda = 10, miniBatch_size = 800)
ctrl_AD_comb_cleaned <-  quantile_norm(ctrl_AD_comb_cleaned, knn_k = 15)
ctrl_AD_comb_cleaned <-  clusterLouvainJaccard(ctrl_AD_comb_cleaned,resolution = 0.7)
ctrl_AD_comb_cleaned <-  runUMAP(ctrl_AD_comb_cleaned)

umap_dim <- plotByDatasetAndCluster(ctrl_AD_comb_cleaned, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

marker.genes <- c('OLIG1','CX3CR1','AQP4','TNR','THEMIS','SLC17A7','GAD2','CLDN5', "FLT1", "KCNJ8", "COL1A1", "RBFOX3", "THY1", "SYT1", "GAD1", "GAD2", "ADARB2", "SST", "FGFR3", "OPALIN", "COL18A1", "C1QA")

pdf("plots/Round02_Markers/snTrem2_ctrlAD_liger_round02_marker-genes.pdf")
lapply(marker.genes,function(x){print(plotGene_keep_scale(ctrl_AD_comb_cleaned, gene = x, plot.by = "none"))})
dev.off()

# Clean Round02 -----------------------------------------------------------
# Remove cluster 14 astro doublets and recluster
doublets_full_clu14_neuron <- data.frame(ctrl_AD_comb_cleaned@clusters, stringsAsFactors = FALSE)
names(doublets_full_clu14_neuron)[1] <- "cluster"
doublets_full_clu14_neuron <- doublets_full_clu14_neuron %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster == "14") %>% 
  pull(barcodes)

ctrl_AD_comb_cleaned2 <- subsetLiger(object = ctrl_AD_comb_cleaned, cells.use = setdiff(unique(ctrl_AD_comb_cleaned@cell.data$barcodes), doublets_full_clu14_neuron))

# LIGER Round03 -----------------------------------------------------------
ctrl_AD_comb_cleaned2 <-  normalize(ctrl_AD_comb_cleaned2)
ctrl_AD_comb_cleaned2 <-  selectGenes(ctrl_AD_comb_cleaned2, do.plot = T, num.genes = 800)
ctrl_AD_comb_cleaned2 <-  scaleNotCenter(ctrl_AD_comb_cleaned2)
ctrl_AD_comb_cleaned2 <-  online_iNMF(ctrl_AD_comb_cleaned2, k = 40, lambda = 10, miniBatch_size = 800)
ctrl_AD_comb_cleaned2 <-  quantile_norm(ctrl_AD_comb_cleaned2)
ctrl_AD_comb_cleaned2 <-  clusterLouvainJaccard(ctrl_AD_comb_cleaned2,resolution = 0.5)
ctrl_AD_comb_cleaned2 <-  runUMAP(ctrl_AD_comb_cleaned2)

umap_dim <- plotByDatasetAndCluster(ctrl_AD_comb_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 4)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

# Annotate
final_annotation <- tibble::tribble(
  ~cluster,  ~cell_type,        ~color,
  0L,     "oligo",      "orange",
  1L,     "oligo",      "orange",
  2L,     "oligo",      "orange",
  3L,     "Excit",  "dodgerblue",
  4L,     "Excit",  "dodgerblue",
  5L,     "Inhib",        "navy",
  6L,     "Astro", "forestgreen",
  7L,     "Excit",  "dodgerblue",
  8L,     "Excit",  "dodgerblue",
  9L,     "Excit",  "dodgerblue",
  10L,     "Excit",  "dodgerblue",
  11L,       "OPC",  "darkorange",
  12L,     "micro",        "gold",
  13L,     "Astro", "forestgreen",
  14L,     "Excit",  "dodgerblue",
  15L,     "Excit",  "dodgerblue",
  16L,     "Inhib",        "navy",
  17L, "Endo/Peri",        "pink",
  18L,     "oligo",      "orange",
  19L,     "Excit",  "dodgerblue"
)

final_clu_colors <- final_annotation %>% 
  pull(color)

umap_dim <- plotByDatasetAndCluster(ctrl_AD_comb_cleaned2, return.plots = TRUE, do.legend = TRUE, text.size = 0)
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = final_clu_colors) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

ctrl_AD_comb_cleaned2_seurat <- ligerToSeurat(ctrl_AD_comb_cleaned2)
write_rds(ctrl_AD_comb_cleaned2_seurat, "zhou_ctrl_AD_seuratv2_temp.RDS")

# Restart R and Install Seurat v3.1.5 -------------------------------------
install.packages("~/Desktop/Bioinformatics Tools/R Source Packages/Seurat_3.1.5.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# Load and Update Object to V3
ctrl_AD_comb_cleaned2_seuratv2 <- read_rds("zhou_ctrl_AD_seuratv2_temp.RDS")
ctrl_AD_comb_cleaned2_seuratv3 <- UpdateSeuratObject(object = ctrl_AD_comb_cleaned2_seuratv2)

# Save
write_rds(ctrl_AD_comb_cleaned2_seuratv3, "zhou_ctrl_AD_seurat_v3_final.RDS")