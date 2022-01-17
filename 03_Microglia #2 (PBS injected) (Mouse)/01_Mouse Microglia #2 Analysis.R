# Code for scRNA-seq Mouse microglia (PBS-Injected) ---------------------------------
    # Analysis was run with
        # R 3.6.1
        # Seurat 3.1.5

# Load packages
library(tidyverse)
library(Seurat)
library(reticulate)
library(viridis)
library(scCustomize)
library(patchwork)
library(beepr)

# Read Raw data
wt_pbs_inj <- Read10X("file_path_here")

# Create Seurat Object with low thresholds
wt_pbs_inj <- CreateSeuratObject(counts = wt_pbs_inj, project = "wt_pbs_injPBS", min.cells = 10, min.features = 200, names.delim = "-", names.field = 2)

# Add mitochondrial gene percentage to meta data
wt_pbs_inj[["percent.mt"]] <- PercentageFeatureSet(wt_pbs_inj, pattern = "^mt-")

# Filter object by thresholds
wt_pbs_inj <- subset(x = wt_pbs_inj, subset = nFeature_RNA > 600 & nFeature_RNA < 3000 & nCount_RNA < 8000 & percent.mt < 10)

# Normalize data
wt_pbs_inj <- NormalizeData(wt_pbs_inj)

# Find Variable Genes
wt_pbs_inj <- FindVariableFeatures(object = wt_pbs_inj, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

# Scale data and regress UMIs and mito%
all.genes <- rownames(wt_pbs_inj)
wt_pbs_inj <- ScaleData(wt_pbs_inj, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mt")); beep(sound = 2)

# Run PCA
Var_Feat <- VariableFeatures(wt_pbs_inj)
wt_pbs_inj <- RunPCA(wt_pbs_inj, features = Var_Feat, npcs = 50)
ElbowPlot(wt_pbs_inj, ndims = 25)

wt_pbs_inj <- JackStraw(wt_pbs_inj, num.replicate = 100)
wt_pbs_inj <- ScoreJackStraw(wt_pbs_inj, dims = 1:20)
JackStrawPlot(wt_pbs_inj, dims = 1:20)
beep(sound = 2)

# Cluster & Plot
wt_pbs_inj <- FindNeighbors(wt_pbs_inj, dims = 1:14)
wt_pbs_inj <- FindClusters(wt_pbs_inj, resolution = 0.2)
wt_pbs_inj <- RunTSNE(wt_pbs_inj, dims = 1:14)

# Add Activation & Microglial Scores
# Features used for microglial and activation scores can be found in SI Tables XX and XX
wt_pbs_inj <- AddModuleScore(wt_pbs_inj, features = microglia_features_list, name = "mg_score")
wt_pbs_inj <- AddModuleScore(wt_pbs_inj, features = shared_sig_list, name = "act_score")

# Plot Act Score on tSNE
FeaturePlot(wt_pbs_inj, features = "act_score1", cols = c("navy", "gold"), pt.size = 4)

# Subset Object by sample to enable individual FeatureScatter Plots
wt_pbs_inj_01 <- subset(wt_pbs_inj, subset = orig.ident == "2")
wt_pbs_inj_02 <- subset(wt_pbs_inj, subset = orig.ident == "3")
wt_pbs_inj_03 <- subset(wt_pbs_inj, subset = orig.ident == "6")
wt_pbs_inj_04 <- subset(wt_pbs_inj, subset = orig.ident == "8")

# Extract min and max values to preserve axes across replicates when split
y_min <- min(wt_pbs_inj@meta.data$mg_score1) - 0.01
y_max <- max(wt_pbs_inj@meta.data$mg_score1) + 0.01
x_min <- min(wt_pbs_inj@meta.data$act_score1) - 0.01

# Create FeatureScatter Plots
p1 <- FeatureScatter(wt_pbs_inj_01, feature1 = "act_score1", feature2 = "mg_score1", cols = "dodgerblue", group.by = "orig.ident", pt.size = 4) + xlim(x_min, 1.5) + ylim(y_min, y_max)

p2 <- FeatureScatter(wt_pbs_inj_02, feature1 = "act_score1", feature2 = "mg_score1", cols = "palegreen4", group.by = "orig.ident", pt.size = 4) + xlim(x_min, 1.5) + ylim(y_min, y_max)

p3 <- FeatureScatter(wt_pbs_inj_03, feature1 = "act_score1", feature2 = "mg_score1", cols = "orchid3", group.by = "orig.ident", pt.size = 4) + xlim(x_min, 1.5) + ylim(y_min, y_max)

p4 <- FeatureScatter(wt_pbs_inj_04, feature1 = "act_score1", feature2 = "mg_score1", cols = "firebrick", group.by = "orig.ident", pt.size = 4) + xlim(x_min, 1.5) + ylim(y_min, y_max)

# View plots before saving
(p1 | p2)/(p3 | p4)
