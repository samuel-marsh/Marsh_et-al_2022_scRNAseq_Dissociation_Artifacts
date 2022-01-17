# Load Data & Create Object -------------------------------------------
raw_data <- read("file_path")

mathys <- CreateSeuratObject(counts = raw_data, project = "mathys_reanalysis", min.cells = 3, min.features = 200)

# Add mito percentage to data
mathys <- PercentageFeatureSet(mathys, pattern = "^mt-", col.name = "percent_mito")

QC_Plots_Genes(mathys)
QC_Plots_UMI(mathys)
QC_Plots_Mito(mathys)

# No mito genes present in gene list
gene_list <- data.frame(mathys@assays$RNA@counts@Dimnames[[1]])

# QC Filter on Counts
mathys <- subset(x = mathys, subset = nCount_RNA < 500000)

# Normalize Data ------------------------------------------------------
mathys <- NormalizeData(mathys, normalization.method = "LogNormalize", scale.factor = 10000)

mathys <- FindVariableFeatures(mathys, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))

# Scale the data
all_genes <- rownames(mathys)
mathys <- ScaleData(mathys, features = all_genes, vars.to.regress = "nCount_RNA")

mathys <- RunPCA(mathys, features = VariableFeatures(object = mathys), npcs = 50)

ElbowPlot(mathys, ndims = 30)

mathys <- JackStraw(mathys)
mathys <- ScoreJackStraw(mathys, dims = 1:20)
JackStrawPlot(mathys, dims = 1:20)
beep(sound = 2)

# Clustering ----------------------------------------------------------
mathys <- FindNeighbors(mathys, dims = 1:5)
mathys <- FindClusters(mathys, resolution = 0.2)

mathys <- RunTSNE(mathys, dims = 1:5)

DimPlot(mathys, label = TRUE, reduction = "tsne")

# Check clustering
markers <- FindAllMarkers(mathys)

# Save Object ---------------------------------------------------------
write_rds(mathys, "Final_RDS_Objects/mathys_reanalyzed_clustered_FINAL.RDS")

# Plot Activation Score -----------------------------------------------
mathys <- AddModuleScore(mathys, features = shared_sig, name = "sg")
# One gene not found Hist2h2aa1.  Old synonym not found.  Excluded from score.
mathys <- AddModuleScore(mathys, features = homeostatic_mg_mathys, name = "mg")
# One gene not found Adgrg1.  Dataset has old gene name Gpr56.  Use altered list.

# Plot Scores
p <- FeaturePlot(object = mathys, features = "sg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/mathys_sg1.pdf", height = 8, width = 9.2)

p <- FeaturePlot(object = mathys, features = "mg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/mathys_mg1.pdf", height = 8, width = 9.2)

# Save Module Scored Object
write_rds(mathys, "Final_RDS_Objects/mathys_module_scored_FINAL.RDS")

# Get final cell number
stats_obj <- read_rds("Final_RDS_Objects/mathys_module_scored_FINAL.RDS")

stats <- Cluster_Stats_All_Samples(stats_obj)