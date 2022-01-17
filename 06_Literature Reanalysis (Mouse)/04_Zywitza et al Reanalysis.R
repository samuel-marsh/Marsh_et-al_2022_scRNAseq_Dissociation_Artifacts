# Load Data & Create Object -------------------------------------------
raw_data1 <- read.table("~/Desktop/Literature Reanalysis/Zywitza et al (Papain Drop-seq)/GSE111527_RAW/GSM3032675_an002_dge.txt.gz", row.names = 1, header = TRUE)

zywitza1 <- CreateSeuratObject(counts = raw_data1, project = "an002", min.cells = 3, min.features = 200)

zywitza1 <- RenameCells(zywitza1, add.cell.id = "A")

# Sample2
raw_data2 <- read.table("~/Desktop/Literature Reanalysis/Zywitza et al (Papain Drop-seq)/GSE111527_RAW/GSM3032676_an003F_dge.txt.gz", row.names = 1, header = TRUE)

zywitza2 <- CreateSeuratObject(counts = raw_data2, project = "an003F", min.cells = 3, min.features = 200)

zywitza2 <- RenameCells(zywitza2, add.cell.id = "B")

# Sample3
raw_data3 <- read.table("~/Desktop/Literature Reanalysis/Zywitza et al (Papain Drop-seq)/GSE111527_RAW/GSM3032677_an003L_dge.txt.gz", row.names = 1, header = TRUE)

zywitza3 <- CreateSeuratObject(counts = raw_data3, project = "an003L", min.cells = 3, min.features = 200)

zywitza3 <- RenameCells(zywitza3, add.cell.id = "C")

# Sample4
raw_data4 <- read.table("~/Desktop/Literature Reanalysis/Zywitza et al (Papain Drop-seq)/GSE111527_RAW/GSM3032678_an008_dge.txt.gz", row.names = 1, header = TRUE)

zywitza4 <- CreateSeuratObject(counts = raw_data4, project = "an008", min.cells = 3, min.features = 200)

zywitza4 <- RenameCells(zywitza4, add.cell.id = "D")

# Sample5
raw_data5 <- read.table("~/Desktop/Literature Reanalysis/Zywitza et al (Papain Drop-seq)/GSE111527_RAW/GSM3032679_an009_dge.txt.gz", row.names = 1, header = TRUE)

zywitza5 <- CreateSeuratObject(counts = raw_data5, project = "an009", min.cells = 3, min.features = 200)

zywitza5 <- RenameCells(zywitza5, add.cell.id = "E")

# Merge Samples into Single Object
zywitza_comb <- merge(x = zywitza1, y = c(zywitza2, zywitza3, zywitza4, zywitza5))

# Add mito percentage to data
zywitza_comb <- PercentageFeatureSet(zywitza_comb, pattern = "^mt-", col.name = "percent_mito")

QC_Plots_Genes(zywitza_comb)
QC_Plots_UMI(zywitza_comb)
QC_Plots_Mito(zywitza_comb)

# QC Filter
zywitza_comb <- subset(x = zywitza_comb, subset = nCount_RNA < 10000 & percent_mito < 15 & nFeature_RNA > 200 & nFeature_RNA < 4000)

# Normalize Data ------------------------------------------------------
zywitza_comb <- NormalizeData(zywitza_comb, normalization.method = "LogNormalize", scale.factor = 10000)

zywitza_comb <- FindVariableFeatures(zywitza_comb, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))

# Scale the data
all_genes <- rownames(zywitza_comb)
zywitza_comb <- ScaleData(zywitza_comb, features = all_genes, vars.to.regress = c("percent_mito", "nCount_RNA"))

zywitza_comb <- RunPCA(zywitza_comb, features = VariableFeatures(object = zywitza_comb), npcs = 100)

ElbowPlot(zywitza_comb, ndims = 40)

zywitza_comb <- JackStraw(zywitza_comb)
zywitza_comb <- ScoreJackStraw(zywitza_comb, dims = 1:20)
JackStrawPlot(zywitza_comb, dims = 1:20)

beep(sound = 2)

# Clustering ----------------------------------------------------------
zywitza_comb <- FindNeighbors(zywitza_comb, dims = 1:25)
zywitza_comb <- FindClusters(zywitza_comb, resolution = 0.6)

zywitza_comb <- RunTSNE(zywitza_comb, dims = 1:25)

DimPlot(zywitza_comb, label = TRUE, reduction = "tsne", label.size = 5)

# Identify myeloid clusters
markers_zywitza_comb <- FindAllMarkers(zywitza_comb)
# Cluster 4 microglia
# Cluster 14 Macro

# Subset myeloid cells ------------------------------------------------
zywitza_micro <- subset(zywitza_comb, idents = c(4, 14))

QC_Plots_Genes(zywitza_micro)
QC_Plots_UMI(zywitza_micro)
QC_Plots_Mito(zywitza_micro)

# QC Filter on Counts
zywitza_micro <- subset(x = zywitza_micro, subset = nCount_RNA < 4000)

# Normalize Data ------------------------------------------------------
zywitza_micro <- NormalizeData(zywitza_micro, normalization.method = "LogNormalize", scale.factor = 10000)

zywitza_micro <- FindVariableFeatures(zywitza_micro, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))

# Scale the data
all_genes <- rownames(zywitza_micro)
zywitza_micro <- ScaleData(zywitza_micro, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))

zywitza_micro <- RunPCA(zywitza_micro, features = VariableFeatures(object = zywitza_micro), npcs = 100)

ElbowPlot(zywitza_micro, ndims = 20)

zywitza_micro <- JackStraw(zywitza_micro)
zywitza_micro <- ScoreJackStraw(zywitza_micro, dims = 1:20)
JackStrawPlot(zywitza_micro, dims = 1:20)

beep(sound = 2)

# Examining PCs shows that PC3 loads on astrocyte genes and is most likely doublets.  Check after clustering and remove.
DimHeatmap(zywitza_micro, dims = 3)

# Clustering ----------------------------------------------------------
zywitza_micro <- FindNeighbors(zywitza_micro, dims = 1:5)
zywitza_micro <- FindClusters(zywitza_micro, resolution = 0.2)

zywitza_micro <- RunTSNE(zywitza_micro, dims = 1:5)

DimPlot(zywitza_micro, label = TRUE, reduction = "tsne")

Customize_SeuratV3_FeaturePlot(zywitza_micro, features = "Fcrls", reduction = "tsne")
Customize_SeuratV3_FeaturePlot(zywitza_micro, features = "Ms4a7", reduction = "tsne")
Customize_SeuratV3_FeaturePlot(zywitza_micro, features = "Aldoc", reduction = "tsne")
# Cluster 3 are astrocytes or astro-microglia doublets

# Remove Astrocyte Doublets and Rerun ------------------------------------------------
zywitza_micro_rd2 <- subset(zywitza_micro, idents = 3, invert = TRUE)

# Normalize Data ------------------------------------------------------
zywitza_micro_rd2 <- NormalizeData(zywitza_micro_rd2, normalization.method = "LogNormalize", scale.factor = 10000)

zywitza_micro_rd2 <- FindVariableFeatures(zywitza_micro_rd2, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))

# Scale the data
all_genes <- rownames(zywitza_micro_rd2)
zywitza_micro_rd2 <- ScaleData(zywitza_micro_rd2, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))

zywitza_micro_rd2 <- RunPCA(zywitza_micro_rd2, features = VariableFeatures(object = zywitza_micro_rd2), npcs = 100)

ElbowPlot(zywitza_micro_rd2, ndims = 20)

zywitza_micro_rd2 <- JackStraw(zywitza_micro_rd2)
zywitza_micro_rd2 <- ScoreJackStraw(zywitza_micro_rd2, dims = 1:20)
JackStrawPlot(zywitza_micro_rd2, dims = 1:20)

beep(sound = 2)

# Examine PCs
DimHeatmap(zywitza_micro_rd2, dims = 4)

# Clustering ----------------------------------------------------------
zywitza_micro_rd2 <- FindNeighbors(zywitza_micro_rd2, dims = 1:3)
zywitza_micro_rd2 <- FindClusters(zywitza_micro_rd2, resolution = 0.2)

zywitza_micro_rd2 <- RunTSNE(zywitza_micro_rd2, dims = 1:3)

DimPlot(zywitza_micro_rd2, label = TRUE, reduction = "tsne")

# Save Object ---------------------------------------------------------
write_rds(zywitza_micro_rd2, "Final_RDS_Objects/zywitza_micro_rd3_reanalyzed_clustered_FINAL.RDS")

# Plot Activation Score -----------------------------------------------
zywitza_micro_rd2 <- read_rds("Final_RDS_Objects/zywitza_micro_rd3_reanalyzed_clustered_FINAL.RDS")

zywitza_micro_rd2 <- AddModuleScore(zywitza_micro_rd2, features = shared_sig, name = "sg")
# One gene not found Hist2h2aa1.  Old synonym not found.  Excluded from score.
zywitza_micro_rd2 <- AddModuleScore(zywitza_micro_rd2, features = homeostatic_mg, name = "mg")

# Plot Scores
p <- FeaturePlot(object = zywitza_micro_rd2, features = "sg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/zywitza_micro_rd2_sg1.pdf", height = 8, width = 9.2)

p <- FeaturePlot(object = zywitza_micro_rd2, features = "mg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/zywitza_micro_rd2_mg1.pdf", height = 8, width = 9.2)

# Save Module Scored Object
write_rds(zywitza_micro_rd2, "Final_RDS_Objects/zywitza_micro_rd3_module_scored_FINAL.RDS")

# Check final Cell Number
stats_obj <- read_rds("Final_RDS_Objects/zywitza_micro_rd3_module_scored_FINAL.RDS")

stats <- Cluster_Stats_All_Samples(stats_obj)