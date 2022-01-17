# Load Data & Create Object -------------------------------------------
raw_data1 <- read.table("~/Desktop/Literature Reanalysis/Mizrak et al (Papain Microwell)/GSE109447_13055_cells_matrix.txt.gz", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
beep(sound = 2)

# Extract gene name list for mizrak dataset 1
gene_list1 <- rownames(raw_data1)

# 2 more genes expressed in dataset 1 than dataset 2.  Pull gene list from dataset one to create mito and score gene list conversions to ensembl IDs
gene_labels_mizrak <- raw_data1 %>% 
  rownames_to_column(var = "ensembl_ID") %>% 
  as_tibble() %>% 
  rename(gene_name = CRERECOMB.1) %>% 
  select(ensembl_ID ,gene_name)

raw_data1 <- raw_data1 %>% 
  select(-CRERECOMB.1)

# Create Seurat Object
mizrak1 <- CreateSeuratObject(counts = raw_data1, project = "mizrak_reanalysis", min.cells = 3, min.features = 200)

# Add identifier to the cells from dataset 1
mizrak1 <- RenameCells(mizrak1, add.cell.id = "A")
beep(sound = 2)

# Load Data from second dataset
raw_data2 <- read.table("~/Desktop/Literature Reanalysis/Mizrak et al (Papain Microwell)/GSE109447_29319_cells.matrix.txt.gz", row.names = 1, stringsAsFactors = FALSE)

# Extract gene name list for mizrak dataset 2
gene_list2 <- rownames(raw_data2)

raw_data2 <- raw_data2 %>% 
  select(-V2)

mizrak2 <- CreateSeuratObject(counts = raw_data2, project = "mizrak_reanalysis1", min.cells = 3, min.features = 200)

mizrak2 <- RenameCells(mizrak2, add.cell.id = "B")

# Identify Ensembl IDs for mito genes for adding mito percent
mito <- gene_labels_mizrak %>% 
  filter(str_detect(gene_name, '^mt-'))

mito_list <- mito %>% 
  pull(ensembl_ID)

# Merge Samples into Single Object
mizrak_comb <- merge(x = mizrak1, y = mizrak2)

# Add mito percentage to data
mizrak_comb <- PercentageFeatureSet(mizrak_comb, features = mito_list, col.name = "percent_mito")

# Plot QC Metrics
QC_Plots_Genes(mizrak_comb, low_cutoff = 400, high_cutoff = 4000)
QC_Plots_UMI(mizrak_comb)
QC_Plots_Mito(mizrak_comb)

# QC Filter 
mizrak_comb <- subset(x = mizrak_comb, subset = nCount_RNA < 10000 & nFeature_RNA > 200 & nFeature_RNA < 4000 & percent_mito < 15)

# Normalize Data ------------------------------------------------------
mizrak_comb <- NormalizeData(mizrak_comb, normalization.method = "LogNormalize", scale.factor = 10000)

mizrak_comb <- FindVariableFeatures(mizrak_comb, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

# Scale the data
all_genes <- rownames(mizrak_comb)
mizrak_comb <- ScaleData(mizrak_comb, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))

mizrak_comb <- RunPCA(mizrak_comb, features = VariableFeatures(object = mizrak_comb), npcs = 100)

ElbowPlot(mizrak_comb, ndims = 60)

# Clustering ----------------------------------------------------------
mizrak_comb <- FindNeighbors(mizrak_comb, dims = 1:40)
mizrak_comb <- FindClusters(mizrak_comb, resolution = 0.6)

mizrak_comb <- RunTSNE(mizrak_comb, dims = 1:40)

DimPlot(mizrak_comb, label = TRUE, reduction = "tsne", label.size = 5)

# Identify myeloid clusters
markers <- FindAllMarkers(mizrak_comb)
# cluster 2 is microglia
# cluster 16 is macro
# cluster 21 is endothelial (maybe with monocyte doublets)

#  Save Object ---------------------------------------------------------
write_rds(mizrak_comb, "Final_RDS_Objects/mizrak_comb_reanalyzed_clustered.RDS")
mizrak_comb <- read_rds("Final_RDS_Objects/mizrak_comb_reanalyzed_clustered.RDS")

#  Subclustering Myeloid ---------------------------------------------
mizrak_micro <- subset(mizrak_comb, idents = c(2, 16))

write_rds(mizrak_micro, "Final_RDS_Objects/mizrak_micro_raw_subset.RDS")

mizrak_micro <- read_rds("Final_RDS_Objects/mizrak_micro_raw_subset.RDS")

# Plot QC Metrics
QC_Plots_Genes(mizrak_micro, low_cutoff = 400)
QC_Plots_UMI(mizrak_micro)
QC_Plots_Mito(mizrak_micro)

# QC Filter
mizrak_micro <- subset(x = mizrak_micro, subset = nCount_RNA < 5000 & nFeature_RNA < 2000)

# Normalize Data ------------------------------------------------------
mizrak_micro <- NormalizeData(mizrak_micro, normalization.method = "LogNormalize", scale.factor = 10000)

mizrak_micro <- FindVariableFeatures(mizrak_micro, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(.75, Inf))

# Scale the data
all_genes <- rownames(mizrak_micro)
mizrak_micro <- ScaleData(mizrak_micro, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))

mizrak_micro <- RunPCA(mizrak_micro, features = VariableFeatures(object = mizrak_micro), npcs = 100)

ElbowPlot(mizrak_micro, ndims = 30)
beep(sound = 2)

mizrak_micro <- JackStraw(mizrak_micro)
mizrak_micro <- ScoreJackStraw(mizrak_micro, dims = 1:20)
JackStrawPlot(mizrak_micro, dims = 1:20)

beep(sound = 2)

# Clustering ----------------------------------------------------------
mizrak_micro <- FindNeighbors(mizrak_micro, dims = 1:20)
mizrak_micro <- FindClusters(mizrak_micro, resolution = 0.3)

mizrak_micro <- RunTSNE(mizrak_micro, dims = 1:20)

DimPlot(mizrak_micro, label = TRUE, reduction = "tsne", label.size = 5)

# Check for doublets/contamination
markers <- FindAllMarkers(mizrak_micro)
# Cluster 7 are T cells
# Cluster 8 are B cells
# Cluster 3 are  astro doublets

write_rds(mizrak_micro, "Final_RDS_Objects/mizrak_micro_rd1_clustered.RDS")

# Subclustering Microglia ---------------------------------------------
mizrak_micro <- read_rds("Final_RDS_Objects/mizrak_micro_rd1_clustered.RDS")

mizrak_micro_rd2 <- subset(mizrak_micro, idents = c(3, 7, 8), invert = TRUE)

write_rds(mizrak_micro_rd2, "Final_RDS_Objects/mizrak_micro_rd2_raw_subset.RDS")

mizrak_micro_rd2 <- read_rds("Final_RDS_Objects/mizrak_micro_rd2_raw_subset.RDS")

# Plot QC Metrics
QC_Plots_Genes(mizrak_micro_rd2, low_cutoff = 400)
QC_Plots_UMI(mizrak_micro_rd2)
QC_Plots_Mito(mizrak_micro_rd2)

# Normalize Data ------------------------------------------------------
mizrak_micro_rd2 <- NormalizeData(mizrak_micro_rd2, normalization.method = "LogNormalize", scale.factor = 10000)

mizrak_micro_rd2 <- FindVariableFeatures(mizrak_micro_rd2, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))

# Scale the data
all_genes <- rownames(mizrak_micro_rd2)
mizrak_micro_rd2 <- ScaleData(mizrak_micro_rd2, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))

mizrak_micro_rd2 <- RunPCA(mizrak_micro_rd2, features = VariableFeatures(object = mizrak_micro_rd2), npcs = 50)

ElbowPlot(mizrak_micro_rd2, ndims = 30)
beep(sound = 2)

mizrak_micro_rd2 <- JackStraw(mizrak_micro_rd2)
mizrak_micro_rd2 <- ScoreJackStraw(mizrak_micro_rd2, dims = 1:20)
JackStrawPlot(mizrak_micro_rd2, dims = 1:20)

beep(sound = 2)

DimHeatmap(mizrak_micro_rd2, dims = 9)

# Clustering ----------------------------------------------------------
mizrak_micro_rd2 <- FindNeighbors(mizrak_micro_rd2, dims = 1:10)
mizrak_micro_rd2 <- FindClusters(mizrak_micro_rd2, resolution = 0.2)

mizrak_micro_rd2 <- RunTSNE(mizrak_micro_rd2, dims = 1:10)

DimPlot(mizrak_micro_rd2, label = TRUE, reduction = "tsne", label.size = 5)

# Save Clustered Object
write_rds(mizrak_micro_rd2, "Final_RDS_Objects/mizrak_micro_rd2_clustered.RDS")

mizrak_micro_rd2 <- read_rds("Final_RDS_Objects/mizrak_micro_rd2_clustered.RDS")

# Module Score on final object
mizrak_micro_rd2 <- AddModuleScore(mizrak_micro_rd2, features = shared_sig_ensembl, name = "sg")
# One gene not found Hist2h2aa1.  Old synonym not found.  Excluded from score.
mizrak_micro_rd2 <- AddModuleScore(mizrak_micro_rd2, features = homeostatic_mg_ensembl, name = "mg")

# Plot Scores
p <- FeaturePlot(object = mizrak_micro_rd2, features = "sg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/mizrak_micro_rd2_sg1.pdf", height = 8, width = 9.2)

p <- FeaturePlot(object = mizrak_micro_rd2, features = "mg1", cols = c("navy", "gold"), pt.size = 6, reduction = "tsne")
p
ggsave("final_module_plots/mizrak_micro_rd2_mg1.pdf", height = 8, width = 9.2)

# Save Module Scored Object
write_rds(mizrak_micro_rd2, "Final_RDS_Objects/mizrak_micro_rd2_module_scored_FINAL.RDS")

# Check cell number
stats_obj <- read_rds("Final_RDS_Objects/mizrak_micro_rd2_module_scored_FINAL.RDS")

stats <- Cluster_Stats_All_Samples(stats_obj)