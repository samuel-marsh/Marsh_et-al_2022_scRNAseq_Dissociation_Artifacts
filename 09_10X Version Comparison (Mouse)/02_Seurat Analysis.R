# Round01 Analysis -------------------------------------------------
# Normalize, HVG, and Scale
version_analysis_rd1 <- NormalizeData(object = version_analysis_rd1)

version_analysis_rd1 <- FindVariableFeatures(object = version_analysis_rd1, selection.method = "vst", nfeatures = 500)

version_analysis_rd1 <- ScaleData(object = version_analysis_rd1, features = rownames(version_analysis_rd1))

# PCA
version_analysis_rd1 <- RunPCA(object = version_analysis_rd1, npcs = 50)
version_analysis_rd1 <- JackStraw(object = version_analysis_rd1, num.replicate = 100)
version_analysis_rd1 <- ScoreJackStraw(object = version_analysis_rd1, dims = 1:20)
JackStrawPlot(object = version_analysis_rd1, 1:20)
ElbowPlot(version_analysis_rd1, ndims = 30)
beep(sound = 2)

# Clustering
version_analysis_rd1 <- FindNeighbors(object = version_analysis_rd1, dims = 1:15)
version_analysis_rd1 <- FindClusters(object = version_analysis_rd1, resolution = 0.2)
version_analysis_rd1 <- RunTSNE(object = version_analysis_rd1, dims = 1:15)

DimPlot(object = version_analysis_rd1, cols = polychrome, label = T, label.size = 6, reduction = "tsne", repel = T)
beep(sound = 5)

# Remove Contaminating Populations & Clean Object --------------------------------------------
# Remove granulocytes, mac/mono, and BAMs
version_analysis_rd2 <- subset(x = version_analysis_rd1, idents = c(3, 4, 5), invert = TRUE)
version_analysis_rd2 <- DietSeurat(version_analysis_rd2)

# Analysis Round 02 ---------------------------------------------------------------------------
# Normalise HVG and scale
version_analysis_rd2 <- NormalizeData(object = version_analysis_rd2)

version_analysis_rd2 <- FindVariableFeatures(object = version_analysis_rd2, selection.method = "vst", nfeatures = 300)

version_analysis_rd2 <- ScaleData(object = version_analysis_rd2, vars.to.regress = c("percent_mito", "version_10x"))

# PCA
version_analysis_rd2 <- RunPCA(object = version_analysis_rd2, npcs = 50)
version_analysis_rd2 <- JackStraw(object = version_analysis_rd2, num.replicate = 100)
version_analysis_rd2 <- ScoreJackStraw(object = version_analysis_rd2, dims = 1:20)
JackStrawPlot(object = version_analysis_rd2, 11:15)
ElbowPlot(version_analysis_rd2, ndims = 30)
beep(sound = 2)

# Clustering
version_analysis_rd2 <- FindNeighbors(object = version_analysis_rd2, dims = 1:13)
version_analysis_rd2 <- FindClusters(object = version_analysis_rd2, resolution = 0.2)
version_analysis_rd2 <- RunTSNE(object = version_analysis_rd2, dims = 1:13)
beep(sound = 5)

DimPlot(object = version_analysis_rd2, cols = polychrome, label = T, label.size = 6, reduction = "tsne", repel = T)

# Add module score
version_analysis_rd2 <- AddModuleScore(object = version_analysis_rd2, features = shared_sig, name = "exAM_score")

# Annotate Clusters ---------------------------------------------------------------------------
new_idents <- c("Homeostatic", "Homeostatic", "Ifn-Responsive", "Chemokine")

names(new_idents) <- levels(version_analysis_rd2)

version_analysis_rd2 <- RenameIdents(object = version_analysis_rd2, new_idents)

# Save Object ---------------------------------------------------------------------------------
write_rds(version_analysis_rd2, "10X_Version_Comparison_Microglia_Final.RDS")
