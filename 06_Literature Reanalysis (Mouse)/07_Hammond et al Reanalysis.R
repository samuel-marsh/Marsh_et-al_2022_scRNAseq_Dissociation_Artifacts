# Load Data & Create Object -------------------------------------------
raw_data1 <- read.table("~/Desktop/Literature Reanalysis/Hammond et al (Dounce 10X)/GSE121654_RAW/GSM3442026_P100_Male_1.dge.txt.gz", row.names = 1, header = TRUE)

hammond1 <- CreateSeuratObject(counts = raw_data1, project = "hammond_reanalysis1", min.cells = 3, min.features = 200)

hammond1 <- RenameCells(hammond1, add.cell.id = "A")

# Sample2
raw_data2 <- read.table("~/Desktop/Literature Reanalysis/Hammond et al (Dounce 10X)/GSE121654_RAW/GSM3442027_P100_Male_2.dge.txt.gz", row.names = 1, header = TRUE)

hammond2 <- CreateSeuratObject(counts = raw_data2, project = "hammond_reanalysis2", min.cells = 3, min.features = 200)

hammond2 <- RenameCells(hammond2, add.cell.id = "B")

# Sample3
raw_data3 <- read.table("~/Desktop/Literature Reanalysis/Hammond et al (Dounce 10X)/GSE121654_RAW/GSM3442028_P100_female_1.dge.txt.gz", row.names = 1, header = TRUE)

hammond3 <- CreateSeuratObject(counts = raw_data3, project = "hammond_reanalysis3", min.cells = 3, min.features = 200)

hammond3 <- RenameCells(hammond3, add.cell.id = "C")

# Sample4
raw_data4 <- read.table("~/Desktop/Literature Reanalysis/Hammond et al (Dounce 10X)/GSE121654_RAW/GSM3442029_P100_female_2.dge.txt.gz", row.names = 1, header = TRUE)

hammond4 <- CreateSeuratObject(counts = raw_data4, project = "hammond_reanalysis4", min.cells = 3, min.features = 200)

hammond4 <- RenameCells(hammond4, add.cell.id = "D")

# Merge Samples into Single Object
hammond_comb <- merge(x = hammond1, y = c(hammond2, hammond3, hammond4))

# Add mito percentage to data
hammond_comb <- PercentageFeatureSet(hammond_comb, pattern = "^mt.", col.name = "percent_mito")

# save raw object
write_rds(hammond_comb, "Final_RDS_Objects/hammond_comb_raw.RDS")

# QC Filter
QC_Plots_Genes(hammond_comb)
QC_Plots_UMI(hammond_comb)
QC_Plots_Mito(hammond_comb)

# QC Filter
hammond_comb <- subset(x = hammond_comb, subset = nCount_RNA < 6000 & percent.mt < 10 & nFeature_RNA > 600 & nFeature_RNA < 3000)

# Normalize Data ------------------------------------------------------
hammond_comb <- NormalizeData(hammond_comb, normalization.method = "LogNormalize", scale.factor = 10000)

hammond_comb <- FindVariableFeatures(hammond_comb, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))

# Scale the data
all_genes <- rownames(hammond_comb)
hammond_comb <- ScaleData(hammond_comb, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))

hammond_comb <- RunPCA(hammond_comb, features = VariableFeatures(object = hammond_comb), npcs = 100)

ElbowPlot(hammond_comb, ndims = 30)
beep(sound = 2)

hammond_comb <- JackStraw(hammond_comb)
hammond_comb <- ScoreJackStraw(hammond_comb, dims = 1:20)
JackStrawPlot(hammond_comb, dims = 1:20)

beep(sound = 2)

DimHeatmap(hammond_comb, dims = 7)

# Clustering ----------------------------------------------------------
hammond_comb <- FindNeighbors(hammond_comb, dims = 1:7)
hammond_comb <- FindClusters(hammond_comb, resolution = 0.2)

hammond_comb <- RunTSNE(hammond_comb, dims = 1:7)

DimPlot(hammond_comb, label = TRUE, reduction = "tsne")

# Check clustering
markers <- FindAllMarkers(hammond_comb)

# Save Object ---------------------------------------------------------
write_rds(hammond_comb, "Final_RDS_Objects/hammond_comb_reanalyzed_clustered.RDS")

# Plot Activation Score -----------------------------------------------
hammond_comb <- AddModuleScore(hammond_comb, features = shared_sig, name = "sg")
# One gene not found Hist2h2aa1.  Old synonym not found.  Excluded from score.
hammond_comb <- AddModuleScore(hammond_comb, features = homeostatic_mg, name = "mg")

# Plot Scores
p <- FeaturePlot(object = hammond_comb, features = "sg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/hammond_comb_sg1.pdf", height = 8, width = 9.2)

p <- FeaturePlot(object = hammond_comb, features = "mg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/hammond_comb_mg1.pdf", height = 8, width = 9.2)

# Save Module Scored Object
write_rds(hammond_comb, "Final_RDS_Objects/hammond_comb_module_scored_FINAL.RDS")

test <- readRDS("Final_RDS_Objects/hammond_comb_module_scored_FINAL.RDS")

# Check final cell number
stats_obj <- read_rds("Final_RDS_Objects/hammond_comb_module_scored_FINAL.RDS")

stats <- Cluster_Stats_All_Samples(stats_obj)