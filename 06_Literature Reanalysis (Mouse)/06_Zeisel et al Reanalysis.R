# Load Data & Create Object -------------------------------------------
zeisel <- connect(filename = "~/Desktop/Literature Reanalysis/Zeisel et al (Papain 10X)/l6_r4_microglia.loom", mode = "r")

zeisel

zeisel_seurat <- as.Seurat(zeisel)

zeisel$close_all()

# save raw seurat object after conversion from loom
write_rds(zeisel_seurat, "Final_RDS_Objects/zeisel_raw_seurat.RDS")

zeisel_seurat <- read_rds("Final_RDS_Objects/zeisel_raw_seurat.RDS")

# Add mito percentage to data
zeisel_seurat <- PercentageFeatureSet(zeisel_seurat, pattern = "^mt-", col.name = "percent_mito")

# QC Filtering
QC_Plots_Genes(zeisel_seurat, low_cutoff = 400)
QC_Plots_UMI(zeisel_seurat, high_cutoff = 5000)
QC_Plots_Mito(zeisel_seurat)

# QC Filter
zeisel_seurat <- subset(x = zeisel_seurat, subset = nCount_RNA < 5000 & percent_mito < 10 & nFeature_RNA > 400)

# Normalize Data ------------------------------------------------------
zeisel_seurat <- NormalizeData(zeisel_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

zeisel_seurat <- FindVariableFeatures(zeisel_seurat, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))

# Scale the data
all_genes <- rownames(zeisel_seurat)
zeisel_seurat <- ScaleData(zeisel_seurat, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))

zeisel_seurat <- RunPCA(zeisel_seurat, features = VariableFeatures(object = zeisel_seurat), npcs = 100)

ElbowPlot(zeisel_seurat, ndims = 40)

beep(sound = 2)

zeisel_seurat <- JackStraw(zeisel_seurat)
zeisel_seurat <- ScoreJackStraw(zeisel_seurat, dims = 1:20)
JackStrawPlot(zeisel_seurat, dims = 1:20)

beep(sound = 2)

DimHeatmap(zeisel_seurat, dims = 12)

# Clustering ----------------------------------------------------------
zeisel_seurat <- FindNeighbors(zeisel_seurat, dims = 1:11)
zeisel_seurat <- FindClusters(zeisel_seurat, resolution = 0.2)

zeisel_seurat <- RunTSNE(zeisel_seurat, dims = 1:11)

DimPlot(zeisel_seurat, label = TRUE, reduction = "tsne", label.size = 5)

# Examine Clusters
markers <- FindAllMarkers(zeisel_seurat)
    # Cluster 5 are endothelial cells

# Save Object ---------------------------------------------------------
write_rds(zeisel_seurat, "Final_RDS_Objects/zeisel_reanalyzed_clustered.RDS")

# Subset the endothelial cells out of the object
zeisel_micro <- subset(zeisel_seurat, idents = 5, invert = TRUE)

write_rds(zeisel_micro, "Final_RDS_Objects/zeisel_micro_subset_raw.RDS")

# QC Filtering
QC_Plots_Genes(zeisel_micro, low_cutoff = 400)
QC_Plots_UMI(zeisel_micro, high_cutoff = 5000)
QC_Plots_Mito(zeisel_micro)

# Normalize Data ------------------------------------------------------
zeisel_micro <- NormalizeData(zeisel_micro, normalization.method = "LogNormalize", scale.factor = 10000)

zeisel_micro <- FindVariableFeatures(zeisel_micro, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))

# Scale the data
all_genes <- rownames(zeisel_micro)
zeisel_micro <- ScaleData(zeisel_micro, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))

zeisel_micro <- RunPCA(zeisel_micro, features = VariableFeatures(object = zeisel_micro), npcs = 100)

ElbowPlot(zeisel_micro, ndims = 40)

beep(sound = 2)

zeisel_micro <- JackStraw(zeisel_micro)
zeisel_micro <- ScoreJackStraw(zeisel_micro, dims = 1:20)
JackStrawPlot(zeisel_micro, dims = 1:20)

beep(sound = 2)

DimHeatmap(zeisel_micro, dims = 12)

# Clustering ----------------------------------------------------------
zeisel_micro <- FindNeighbors(zeisel_micro, dims = 1:12)
zeisel_micro <- FindClusters(zeisel_micro, resolution = 0.2)

zeisel_micro <- RunTSNE(zeisel_micro, dims = 1:12)

DimPlot(zeisel_micro, label = TRUE, reduction = "tsne", label.size = 5)

# Examine Clusters
markers <- FindAllMarkers(zeisel_micro)

# Save Object ---------------------------------------------------------
write_rds(zeisel_micro, "Final_RDS_Objects/zeisel_micro_reanalyzed_clustered_FINAL.RDS")

# Plot Activation Score -----------------------------------------------
zeisel_micro <- read_rds("Final_RDS_Objects/zeisel_micro_reanalyzed_clustered_FINAL.RDS")

zeisel_micro <- AddModuleScore(zeisel_micro, features = shared_sig, name = "sg")
# One gene not found Hist2h2aa1.  Old synonym not found.  Excluded from score.
zeisel_micro <- AddModuleScore(zeisel_micro, features = homeostatic_mg, name = "mg")

# Plot Scores
p <- FeaturePlot(object = zeisel_micro, features = "sg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/zeisel_micro_sg1.pdf", height = 8, width = 9.2)

p <- FeaturePlot(object = zeisel_micro, features = "mg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/zeisel_micro_mg1.pdf", height = 8, width = 9.2)

# Save Module Scored Object
write_rds(zeisel_micro, "Final_RDS_Objects/zeisel_micro_module_scored_FINAL.RDS")

# Check cell number
stats_obj <- read_rds("Final_RDS_Objects/zeisel_micro_module_scored_FINAL.RDS")

stats <- Cluster_Stats_All_Samples(stats_obj)