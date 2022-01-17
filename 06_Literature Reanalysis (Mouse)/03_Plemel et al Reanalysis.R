# Load Data & Create Object -------------------------------------------
raw_data <- read("file_path")

plemel <- CreateSeuratObject(counts = raw_data, project = "plemel_reanalysis", min.cells = 3, min.features = 200)

# Add mito percentage to data
plemel <- PercentageFeatureSet(plemel, pattern = "^mt-", col.name = "percent_mito")

QC_Plots_Genes(plemel)
QC_Plots_UMI(plemel)
QC_Plots_Mito(plemel)

# QC Filter on Counts
plemel <- subset(x = plemel, subset = nCount_RNA < 6000 & percent_mito < 10 & nFeature_RNA > 600)

# Normalize Data ------------------------------------------------------
plemel <- NormalizeData(plemel, normalization.method = "LogNormalize", scale.factor = 10000)

plemel <- FindVariableFeatures(plemel, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

# Scale the data
all_genes <- rownames(plemel)
plemel <- ScaleData(plemel, features = all_genes, vars.to.regress = c("percent_mito", "nCount_RNA"))

plemel <- RunPCA(plemel, features = VariableFeatures(object = plemel), npcs = 50)

ElbowPlot(plemel, ndims = 20)

plemel <- JackStraw(plemel)
plemel <- ScoreJackStraw(plemel, dims = 1:20)
JackStrawPlot(plemel, dims = 1:10)
beep(sound = 2)

# Clustering ----------------------------------------------------------
plemel <- FindNeighbors(plemel, dims = 1:3)
plemel <- FindClusters(plemel, resolution = 0.2)

plemel <- RunTSNE(plemel, dims = 1:3)

DimPlot(plemel, label = TRUE, reduction = "tsne")

# Save Object ---------------------------------------------------------
write_rds(plemel, "Final_RDS_Objects/plemel_reanalyzed_clustered_FINAL.RDS")

# Module Scoring
plemel <- AddModuleScore(plemel, features = shared_sig, name = "sg")
# One gene not found Hist2h2aa1.  Old synonym not found.  Excluded from score.
plemel <- AddModuleScore(plemel, features = homeostatic_mg, name = "mg")
# One gene not found Adgrg1.  Dataset has old gene name Gpr56.  Use altered list.

# Plot Scores
p <- FeaturePlot(object = plemel, features = "sg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/plemel_sg1.pdf", height = 8, width = 9.2)

p <- FeaturePlot(object = plemel, features = "mg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/plemel_mg1.pdf", height = 8, width = 9.2)

# Save Module Scored Object
write_rds(plemel, "Final_RDS_Objects/plemel_module_scored_FINAL.RDS")

# Check final cell number
stats_obj <- read_rds("Final_RDS_Objects/plemel_module_scored_FINAL.RDS")

stats <- Cluster_Stats_All_Samples(stats_obj)