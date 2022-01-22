# Analysis Round 01 ---------------------------------------------------------------------------
# Normalize HVG and Scale
pbmc_mock_rd1 <- NormalizeData(object = pbmc_mock_qc_filtered)
pbmc_mock_rd1 <- FindVariableFeatures(object = pbmc_mock_rd1, selection.method = "vst", nfeatures = 1000)
pbmc_mock_rd1 <- ScaleData(object = pbmc_mock_rd1)

# PCA
pbmc_mock_rd1 <- RunPCA(object = pbmc_mock_rd1, npcs = 75)
pbmc_mock_rd1 <- JackStraw(object = pbmc_mock_rd1, num.replicate = 100, dims = 30)
pbmc_mock_rd1 <- ScoreJackStraw(object = pbmc_mock_rd1, dims = 1:30)
JackStrawPlot(object = pbmc_mock_rd1, dims = 1:20)
JackStrawPlot(object = pbmc_mock_rd1, dims = 20:30)
ElbowPlot(object = pbmc_mock_rd1, ndims = 30)
beep(sound = 5)

# Cluster
pbmc_mock_rd1 <- FindNeighbors(object = pbmc_mock_rd1, dims = 1:30)
pbmc_mock_rd1 <- FindClusters(object = pbmc_mock_rd1, resolution = 0.6)
pbmc_mock_rd1 <- RunUMAP(pbmc_mock_rd1, dims = 1:30)

DimPlot(object = pbmc_mock_rd1, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

# Azimuth Annotation Rd1 ------------------------------------------------------
azimuth_mito10 <- read.delim(".../azimuth_pred.tsv")

meta_azimuth_merge <- pbmc_mock_rd1@meta.data %>%
  rownames_to_column("barcodes")

joined_azimuth <- left_join(meta_azimuth_merge, azimuth_mito10, by = c("barcodes" = "cell"))

joined_azimuth <- joined_azimuth %>%
  select(barcodes, predicted.celltype.l2, predicted.celltype.l2.score, mapping.score) %>%
  column_to_rownames("barcodes")

pbmc_mock_rd1 <- AddMetaData(object = pbmc_mock_rd1, metadata = joined_azimuth)

DimPlot(pbmc_mock_rd1, cols = DiscretePalette(36, "polychrome"), group.by = "predicted.celltype.l2", label = T, label.size = 4)
DimPlot(pbmc_mock_rd1, cols = DiscretePalette(36, "polychrome"), group.by = "predicted.celltype.l1", label = T, label.size = 4)
DimPlot(pbmc_mock_rd1, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 4)

azimuth_mito10 <- read.delim("...//azimuth_pred_L1.tsv")

meta_azimuth_merge <- pbmc_mock_rd1@meta.data %>%
  rownames_to_column("barcodes")

joined_azimuth <- left_join(meta_azimuth_merge, azimuth_mito10, by = c("barcodes" = "cell"))

joined_azimuth <- joined_azimuth %>%
  select(barcodes, predicted.celltype.l1, predicted.celltype.l1.score) %>%
  column_to_rownames("barcodes")

pbmc_mock_rd1 <- AddMetaData(object = pbmc_mock_rd1, metadata = joined_azimuth)

# Annotate Round 01 --------------------------------------------------
pbmc_mock_rd1_annotation <- Pull_Cluster_Annotation(seurat_object = pbmc_mock_rd1, file_name = ".../pbmc_cluster_annotation_Round01.csv")

pbmc_mock_rd1_annotate <- pbmc_mock_rd1

new_ids_rd1 <- pbmc_mock_rd1_annotate[["new_cluster_idents"]]
names(new_ids_rd1) <- levels(pbmc_mock_rd1_annotate)

pbmc_mock_rd1_annotate <- RenameIdents(object = pbmc_mock_rd1_annotate, new_ids_rd1)

DimPlot(pbmc_mock_rd1_annotate, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 4)

# Subclustering: Filter & Clean Doublets ------------------------------------------------------
# T Cells -----------------------------------------------------------------
pbmc_rd1_tcell <- subset(pbmc_mock_rd1_annotate, idents = "T Cell")
pbmc_rd1_tcell <- NormalizeData(object = pbmc_rd1_tcell)
pbmc_rd1_tcell <- FindVariableFeatures(object = pbmc_rd1_tcell, selection.method = "vst", nfeatures = 750)
pbmc_rd1_tcell <- ScaleData(object = pbmc_rd1_tcell)
pbmc_rd1_tcell <- RunPCA(object = pbmc_rd1_tcell, npcs = 75)
pbmc_rd1_tcell <- JackStraw(object = pbmc_rd1_tcell, num.replicate = 100, dims = 30)
pbmc_rd1_tcell <- ScoreJackStraw(object = pbmc_rd1_tcell, dims = 1:30)
JackStrawPlot(object = pbmc_rd1_tcell, dims = 1:20)
JackStrawPlot(object = pbmc_rd1_tcell, dims = 20:30)
ElbowPlot(object = pbmc_rd1_tcell, ndims = 50)
beep(sound = 2)

pbmc_rd1_tcell <- FindNeighbors(object = pbmc_rd1_tcell, dims = 1:20)
pbmc_rd1_tcell <- FindClusters(object = pbmc_rd1_tcell, resolution = 0.6)
pbmc_rd1_tcell <- RunUMAP(pbmc_rd1_tcell, dims = 1:20)

DimPlot(object = pbmc_rd1_tcell, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

FeaturePlot_scCustom(pbmc_rd1_tcell, features = "CD3E")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "CD4")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "CD8A")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "MS4A1")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "CD14")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "FCGR3A")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "MS4A7")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "NKG7")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "GNLY")

FeaturePlot_scCustom(pbmc_rd1_tcell, features = "FCER1G")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = "TYROBP")
FeaturePlot_scCustom(pbmc_rd1_tcell, features = c("KLRB1", "NKG7", "GZMK", "IL7R", "SLC4A10", "GZMA", "CXCR6", "PRSS35", "RBM24", "NCR3")
)

# NK Cells -----------------------------------------------------------------
pbmc_rd1_nkcell <- subset(pbmc_mock_rd1_annotate, idents = "NK")
pbmc_rd1_nkcell <- NormalizeData(object = pbmc_rd1_nkcell)
pbmc_rd1_nkcell <- FindVariableFeatures(object = pbmc_rd1_nkcell, selection.method = "vst", nfeatures = 1000)
pbmc_rd1_nkcell <- ScaleData(object = pbmc_rd1_nkcell)
pbmc_rd1_nkcell <- RunPCA(object = pbmc_rd1_nkcell, npcs = 75)
pbmc_rd1_nkcell <- JackStraw(object = pbmc_rd1_nkcell, num.replicate = 100, dims = 30)
pbmc_rd1_nkcell <- ScoreJackStraw(object = pbmc_rd1_nkcell, dims = 1:30)
JackStrawPlot(object = pbmc_rd1_nkcell, dims = 1:20)
JackStrawPlot(object = pbmc_rd1_nkcell, dims = 20:30)
ElbowPlot(object = pbmc_rd1_nkcell, ndims = 30)
beep(sound = 2)

pbmc_rd1_nkcell <- FindNeighbors(object = pbmc_rd1_nkcell, dims = 1:14)
pbmc_rd1_nkcell <- FindClusters(object = pbmc_rd1_nkcell, resolution = 0.6)
pbmc_rd1_nkcell <- RunUMAP(pbmc_rd1_nkcell, dims = 1:14)

DimPlot(object = pbmc_rd1_nkcell, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "CD3E")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "CD4")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "CD8A")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "MS4A1")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "CD14")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "FCGR3A")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "MS4A7")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "NKG7")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "FCER1G")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = "TYROBP")
FeaturePlot_scCustom(pbmc_rd1_nkcell, features = c("CCL5", "GZMH", "CD8A", "TRAC", "KLRD1", "NKG7", "GZMK", "CST7", "CD8B", "TRGC2"))

# isolate doublets
doublet_nk <- WhichCells(object = pbmc_rd1_nkcell, idents = 5)
# Cluster 5 appears to me more NK/CD8 mix

# Mono -----------------------------------------------------------------
pbmc_rd1_mono <- subset(pbmc_mock_rd1_annotate, idents = "Mono")
pbmc_rd1_mono <- NormalizeData(object = pbmc_rd1_mono)
pbmc_rd1_mono <- FindVariableFeatures(object = pbmc_rd1_mono, selection.method = "vst", nfeatures = 1000)
pbmc_rd1_mono <- ScaleData(object = pbmc_rd1_mono)
pbmc_rd1_mono <- RunPCA(object = pbmc_rd1_mono, npcs = 75)
pbmc_rd1_mono <- JackStraw(object = pbmc_rd1_mono, num.replicate = 100, dims = 30)
pbmc_rd1_mono <- ScoreJackStraw(object = pbmc_rd1_mono, dims = 1:30)
JackStrawPlot(object = pbmc_rd1_mono, dims = 1:20)
JackStrawPlot(object = pbmc_rd1_mono, dims = 20:30)
ElbowPlot(object = pbmc_rd1_mono, ndims = 30)
beep(sound = 2)

pbmc_rd1_mono <- FindNeighbors(object = pbmc_rd1_mono, dims = 1:15)
pbmc_rd1_mono <- FindClusters(object = pbmc_rd1_mono, resolution = 0.6)
pbmc_rd1_mono <- RunUMAP(pbmc_rd1_mono, dims = 1:15)

DimPlot(object = pbmc_rd1_mono, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

FeaturePlot_scCustom(pbmc_rd1_mono, features = "CD3E")
FeaturePlot_scCustom(pbmc_rd1_mono, features = "CD4")
FeaturePlot_scCustom(pbmc_rd1_mono, features = "CD8A")
FeaturePlot_scCustom(pbmc_rd1_mono, features = "MS4A1")
FeaturePlot_scCustom(pbmc_rd1_mono, features = "CD14")
FeaturePlot_scCustom(pbmc_rd1_mono, features = "FCGR3A")
FeaturePlot_scCustom(pbmc_rd1_mono, features = "MS4A7")
FeaturePlot_scCustom(pbmc_rd1_mono, features = "NKG7")

# Isolate Doublets
doublet_mono <- WhichCells(object = pbmc_rd1_mono, idents = c(4, 7, 8))
# Cluster 4, 7 are T Cells/ NK Cells
# Cluster 8 are B Cells

# B Cells -----------------------------------------------------------------
pbmc_rd1_bcell <- subset(pbmc_mock_rd1_annotate, idents = "B Cell")
pbmc_rd1_bcell <- NormalizeData(object = pbmc_rd1_bcell)
pbmc_rd1_bcell <- FindVariableFeatures(object = pbmc_rd1_bcell, selection.method = "vst", nfeatures = 1000)
pbmc_rd1_bcell <- ScaleData(object = pbmc_rd1_bcell)
pbmc_rd1_bcell <- RunPCA(object = pbmc_rd1_bcell, npcs = 75)
pbmc_rd1_bcell <- JackStraw(object = pbmc_rd1_bcell, num.replicate = 100, dims = 30)
pbmc_rd1_bcell <- ScoreJackStraw(object = pbmc_rd1_bcell, dims = 1:30)
JackStrawPlot(object = pbmc_rd1_bcell, dims = 1:20)
JackStrawPlot(object = pbmc_rd1_bcell, dims = 20:30)
ElbowPlot(object = pbmc_rd1_bcell, ndims = 30)
beep(sound = 2)

pbmc_rd1_bcell <- FindNeighbors(object = pbmc_rd1_bcell, dims = 1:12)
pbmc_rd1_bcell <- FindClusters(object = pbmc_rd1_bcell, resolution = 0.6)
pbmc_rd1_bcell <- RunUMAP(pbmc_rd1_bcell, dims = 1:12)

DimPlot(object = pbmc_rd1_bcell, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

FeaturePlot_scCustom(pbmc_rd1_bcell, features = "CD3E")
FeaturePlot_scCustom(pbmc_rd1_bcell, features = "CD4")
FeaturePlot_scCustom(pbmc_rd1_bcell, features = "CD8A")
FeaturePlot_scCustom(pbmc_rd1_bcell, features = "MS4A1")
FeaturePlot_scCustom(pbmc_rd1_bcell, features = "CD14")
FeaturePlot_scCustom(pbmc_rd1_bcell, features = "FCGR3A")
FeaturePlot_scCustom(pbmc_rd1_bcell, features = "MS4A7")
FeaturePlot_scCustom(pbmc_rd1_bcell, features = "NKG7")

# Isolate Doublets
doublet_bcell <- WhichCells(object = pbmc_rd1_bcell, idents = 4)
# Cluster 4 is T/NK Cells

# DCs -----------------------------------------------------------------
pbmc_rd1_DC <- subset(pbmc_mock_rd1_annotate, idents = "DC")
pbmc_rd1_DC <- NormalizeData(object = pbmc_rd1_DC)
pbmc_rd1_DC <- FindVariableFeatures(object = pbmc_rd1_DC, selection.method = "vst", nfeatures = 1000)
pbmc_rd1_DC <- ScaleData(object = pbmc_rd1_DC)
pbmc_rd1_DC <- RunPCA(object = pbmc_rd1_DC, npcs = 75)
pbmc_rd1_DC <- JackStraw(object = pbmc_rd1_DC, num.replicate = 100, dims = 30)
pbmc_rd1_DC <- ScoreJackStraw(object = pbmc_rd1_DC, dims = 1:30)
JackStrawPlot(object = pbmc_rd1_DC, dims = 1:20)
JackStrawPlot(object = pbmc_rd1_DC, dims = 20:30)
ElbowPlot(object = pbmc_rd1_DC, ndims = 30)
beep(sound = 2)

pbmc_rd1_DC <- FindNeighbors(object = pbmc_rd1_DC, dims = 1:7)
pbmc_rd1_DC <- FindClusters(object = pbmc_rd1_DC, resolution = 0.6)
pbmc_rd1_DC <- RunUMAP(pbmc_rd1_DC, dims = 1:7)

DimPlot(object = pbmc_rd1_DC, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

FeaturePlot_scCustom(pbmc_rd1_DC, features = "CD3E")
FeaturePlot_scCustom(pbmc_rd1_DC, features = "CD4")
FeaturePlot_scCustom(pbmc_rd1_DC, features = "CD8A")
FeaturePlot_scCustom(pbmc_rd1_DC, features = "MS4A1")
FeaturePlot_scCustom(pbmc_rd1_DC, features = "CD14")
FeaturePlot_scCustom(pbmc_rd1_DC, features = "FCGR3A")
FeaturePlot_scCustom(pbmc_rd1_DC, features = "MS4A7")
FeaturePlot_scCustom(pbmc_rd1_DC, features = "NKG7")

# Platelet -----------------------------------------------------------------
pbmc_rd1_platelet <- subset(pbmc_mock_rd1_annotate, idents = "platelet")
pbmc_rd1_platelet <- NormalizeData(object = pbmc_rd1_platelet)
pbmc_rd1_platelet <- FindVariableFeatures(object = pbmc_rd1_platelet, selection.method = "vst", nfeatures = 500)
pbmc_rd1_platelet <- ScaleData(object = pbmc_rd1_platelet)
pbmc_rd1_platelet <- RunPCA(object = pbmc_rd1_platelet, npcs = 75)
pbmc_rd1_platelet <- JackStraw(object = pbmc_rd1_platelet, num.replicate = 100, dims = 30)
pbmc_rd1_platelet <- ScoreJackStraw(object = pbmc_rd1_platelet, dims = 1:30)
JackStrawPlot(object = pbmc_rd1_platelet, dims = 1:20)
JackStrawPlot(object = pbmc_rd1_platelet, dims = 20:30)
ElbowPlot(object = pbmc_rd1_platelet, ndims = 30)
beep(sound = 2)

pbmc_rd1_platelet <- FindNeighbors(object = pbmc_rd1_platelet, dims = 1:10)
pbmc_rd1_platelet <- FindClusters(object = pbmc_rd1_platelet, resolution = 0.6)
pbmc_rd1_platelet <- RunUMAP(pbmc_rd1_platelet, dims = 1:10)

DimPlot(object = pbmc_rd1_platelet, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

FeaturePlot_scCustom(pbmc_rd1_platelet, features = "CD3E")
FeaturePlot_scCustom(pbmc_rd1_platelet, features = "CD4")
FeaturePlot_scCustom(pbmc_rd1_platelet, features = "CD8A")
FeaturePlot_scCustom(pbmc_rd1_platelet, features = "MS4A1")
FeaturePlot_scCustom(pbmc_rd1_platelet, features = "CD14")
FeaturePlot_scCustom(pbmc_rd1_platelet, features = "FCGR3A")
FeaturePlot_scCustom(pbmc_rd1_platelet, features = "MS4A7")
FeaturePlot_scCustom(pbmc_rd1_platelet, features = "NKG7")
FeaturePlot_scCustom(pbmc_rd1_platelet, features = c("S100A9", "CTSS", "S100A8", "LYZ", "VCAN", "S100A12", "IL1B", "CD14", "G0S2", "FCN1"))

# Isolate Doublets
doublet_platelet <- WhichCells(object = pbmc_rd1_platelet, idents = 1)
# Cluster 1 are monocyte doublets

# HSPC -----------------------------------------------------------------
pbmc_rd1_HSPC <- subset(pbmc_mock_rd1_annotate, idents = "HSPC")
pbmc_rd1_HSPC <- NormalizeData(object = pbmc_rd1_HSPC)
pbmc_rd1_HSPC <- FindVariableFeatures(object = pbmc_rd1_HSPC, selection.method = "vst", nfeatures = 250)
pbmc_rd1_HSPC <- ScaleData(object = pbmc_rd1_HSPC)
pbmc_rd1_HSPC <- RunPCA(object = pbmc_rd1_HSPC, npcs = 25)
pbmc_rd1_HSPC <- JackStraw(object = pbmc_rd1_HSPC, num.replicate = 100, dims = 25)
pbmc_rd1_HSPC <- ScoreJackStraw(object = pbmc_rd1_HSPC, dims = 1:25)
JackStrawPlot(object = pbmc_rd1_HSPC, dims = 1:20)
JackStrawPlot(object = pbmc_rd1_HSPC, dims = 20:25)
ElbowPlot(object = pbmc_rd1_HSPC, ndims = 25)
beep(sound = 5)

pbmc_rd1_HSPC <- FindNeighbors(object = pbmc_rd1_HSPC, dims = 1:10)
pbmc_rd1_HSPC <- FindClusters(object = pbmc_rd1_HSPC, resolution = 0.6)
pbmc_rd1_HSPC <- RunUMAP(pbmc_rd1_HSPC, dims = 1:10)

DimPlot(object = pbmc_rd1_HSPC, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

FeaturePlot_scCustom(pbmc_rd1_HSPC, features = "CD3E")
FeaturePlot_scCustom(pbmc_rd1_HSPC, features = "CD4")
FeaturePlot_scCustom(pbmc_rd1_HSPC, features = "CD8A")
FeaturePlot_scCustom(pbmc_rd1_HSPC, features = "MS4A1")
FeaturePlot_scCustom(pbmc_rd1_HSPC, features = "CD14")
FeaturePlot_scCustom(pbmc_rd1_HSPC, features = "FCGR3A")
FeaturePlot_scCustom(pbmc_rd1_HSPC, features = "MS4A7")
FeaturePlot_scCustom(pbmc_rd1_HSPC, features = "NKG7")

# Doublet Clusters --------------------------------------------------------
doublet_doublet <- WhichCells(object = pbmc_mock_rd1_annotate, idents = "Doublet")

# Clean & Remove Doublets ---------------------------------------------------------
doublets_combined <- c(doublet_bcell, doublet_doublet, doublet_mono, doublet_nk, doublet_platelet)

# Remove doublets from original object
pbmc_mock_rd1_annotate_clean <- subset(x = pbmc_mock_rd1_annotate, cells = doublets_combined, invert = TRUE)
pbmc_mock_rd2 <- DietSeurat(pbmc_mock_rd1_annotate_clean)

# Analysis Rd2 --------------------------------------------------------------------------------
# Normalize HVG scale
pbmc_mock_rd2 <- NormalizeData(object = pbmc_mock_rd2)
pbmc_mock_rd2 <- FindVariableFeatures(object = pbmc_mock_rd2, selection.method = "vst", nfeatures = 1000)
pbmc_mock_rd2 <- ScaleData(object = pbmc_mock_rd2)

# PCA
pbmc_mock_rd2 <- RunPCA(object = pbmc_mock_rd2, npcs = 75)
pbmc_mock_rd2 <- JackStraw(object = pbmc_mock_rd2, num.replicate = 100, dims = 30)
pbmc_mock_rd2 <- ScoreJackStraw(object = pbmc_mock_rd2, dims = 1:30)
JackStrawPlot(object = pbmc_mock_rd2, dims = 1:20)
JackStrawPlot(object = pbmc_mock_rd2, dims = 20:30)
ElbowPlot(object = pbmc_mock_rd2, ndims = 40)
beep(sound = 5)

# Cluster
pbmc_mock_rd2 <- FindNeighbors(object = pbmc_mock_rd2, dims = 1:25)
pbmc_mock_rd2 <- FindClusters(object = pbmc_mock_rd2, resolution = 0.4)
pbmc_mock_rd2_annotated <- RunUMAP(object = pbmc_mock_rd2_annotated, dims = 1:25, n.neighbors = 20, min.dist = 0.5)

DimPlot(object = pbmc_mock_rd2, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

# Add Azimuth Mapping Info ------------------------------------------------
azimuth_rd2 <- read.delim(".../azimuth_pred_rd2.tsv")

meta_azimuth_merge <- pbmc_mock_rd2@meta.data %>%
  rownames_to_column("barcodes")

joined_azimuth <- left_join(meta_azimuth_merge, azimuth_mito10, by = c("barcodes" = "cell"))

joined_azimuth <- joined_azimuth %>%
  select(barcodes, predicted.celltype.l2_Rd2, predicted.celltype.l1_Rd2, predicted.celltype.l2.score_Rd2, predicted.celltype.l1.score_Rd2, mapping.score_Rd2) %>%
  column_to_rownames("barcodes")

pbmc_mock_rd2 <- AddMetaData(object = pbmc_mock_rd2, metadata = joined_azimuth)

p2 <- DimPlot(pbmc_mock_rd2, cols = DiscretePalette(36, "polychrome"), group.by = "predicted.celltype.l2_Rd2", label = T, label.size = 4)
DimPlot(pbmc_mock_rd2, cols = DiscretePalette(36, "polychrome"), group.by = "predicted.celltype.l1_Rd2", label = T, label.size = 4)
p1 <- DimPlot(pbmc_mock_rd2, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 4)
wrap_plots(p1, p2)

# Create Final Annotation
annotation_pbmc_rd2 <- Pull_Cluster_Annotation(seurat_object = pbmc_mock_rd2, file_name = ".../pbmc_cluster_annotation_Round02.csv")

pbmc_mock_rd2_annotated <- pbmc_mock_rd2

basic_idents <- annotation_pbmc_rd2[["basic_type"]]
names(basic_idents) <- levels(pbmc_mock_rd2_annotated)

DimPlot(pbmc_mock_rd2_annotated, label = T)

pbmc_mock_rd2_annotated <- RenameIdents(pbmc_mock_rd2_annotated, basic_idents)

pbmc_mock_rd2_annotated[["basic_idents"]] <- Idents(pbmc_mock_rd2_annotated)

store_basic_mapping <- pbmc_mock_rd2_annotated@meta.data %>%
  rownames_to_column("barcodes") %>%
  select(barcodes, basic_idents) %>%
  column_to_rownames("barcodes")

# Add full idents
pbmc_mock_rd2_annotated <- pbmc_mock_rd2

pbmc_mock_rd2_annotated <- RenameIdents(pbmc_mock_rd2_annotated, annotation_pbmc_rd2[["new_cluster_idents"]])

DimPlot(object = pbmc_mock_rd2_annotated, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

# Add basic idents
seurat_meta <- pbmc_mock_rd2_annotated@meta.data %>%
  rownames_to_column("barcodes")

store_basic_mapping <- store_basic_mapping %>%
  rownames_to_column("barcodes")

joined_meta <- left_join(seurat_meta, store_basic_mapping) %>%
  select(barcodes, basic_idents) %>%
  column_to_rownames("barcodes")

pbmc_mock_rd2_annotated <- AddMetaData(object = pbmc_mock_rd2_annotated, metadata = joined_meta)

DimPlot(object = pbmc_mock_rd2_annotated, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap")

DimPlot(object = pbmc_mock_rd2_annotated, cols = DiscretePalette(36, "polychrome"), label = T, label.size = 6, reduction = "umap", group.by = "basic_idents")

# reorg ident levels
new_ident_levels <- c("T Cell_CD4", "T Cell_CD8", "T Cell_gdT", "T Cell_MAIT", "T Cell_other", "B Cell", "NK Cell", "NK Cell_CD56Br", "Monocyte_CD14",  "Monocyte_CD14_exAM", "Monocyte_CD16","DC", "HSPC", "Platlet")

new_basic_ident_levels <- c("T Cell", "B Cell", "NK Cell", "Monocyte", "DC", "HSPC", "Platlet")


Idents(pbmc_mock_rd2_annotated) <- factor(Idents(pbmc_mock_rd2_annotated), new_ident_levels)

pbmc_mock_rd2_annotated@meta.data$basic_idents <- factor(pbmc_mock_rd2_annotated@meta.data$basic_idents, new_basic_ident_levels)

# Create color palettes
colors_full_idents <- c("dodgerblue", "skyblue", "dodgerblue4", "navy", "skyblue3",
                        "darkgoldenrod2",
                        "palegreen4", "palegreen2",
                        "firebrick1", "firebrick3", "indianred", "tomato4",
                        "magenta", "darkorchid")

colors_basic_idents <- c("dodgerblue",
                         "darkgoldenrod2",
                         "forestgreen",
                         "indianred2",
                         "tomato4",
                         "magenta",
                         "darkorchid")

color_list <- list(colors_full_idents, colors_basic_idents)
names(color_list) <- c("colors_full_idents", "colors_basic_idents")

DimPlot(object = pbmc_mock_rd2_annotated, cols = colors_full_idents, label = T, label.size = 6, reduction = "umap", repel = T)

DimPlot(object = pbmc_mock_rd2_annotated, cols = colors_basic_idents, label = T, label.size = 6, reduction = "umap", group.by = "basic_idents")

# Add exAM Score ------------------------------------------------------------------------------
pbmc_mock_rd2_annotated <- AddModuleScore(object = pbmc_mock_rd2_annotated, features = human_shared_sig, name = "mouse_exAM", search = T)


# Identify Shared DEG by Celltype -------------------------------------------------------------
# DEG Basic Class by Status -----------------------------------------------
Idents(pbmc_mock_rd2_annotated) <- "basic_idents"

t_basic_deg <- FindMarkers(object = pbmc_mock_rd2_annotated, subset.ident = "T Cell", ident.1 = "control", ident.2 = "inhib", group.by = "orig.ident") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2)
t_basic_deg$gene <- rownames(t_basic_deg)
t_deg_list <- t_basic_deg %>%
  pull(gene)
write.csv(t_basic_deg, "07_outputs/t_cell_deg_seurat.csv")

b_basic_deg <- FindMarkers(object = pbmc_mock_rd2_annotated, subset.ident = "B Cell", ident.1 = "control", ident.2 = "inhib", group.by = "orig.ident") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2)
b_basic_deg$gene <- rownames(b_basic_deg)
b_deg_list <- b_basic_deg %>%
  pull(gene)
write.csv(b_basic_deg, "07_outputs/b_cell_deg_seurat.csv")

nk_basic_deg <- FindMarkers(object = pbmc_mock_rd2_annotated, subset.ident = "NK Cell", ident.1 = "control", ident.2 = "inhib", group.by = "orig.ident") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2)
nk_basic_deg$gene <- rownames(nk_basic_deg)
nk_deg_list <- nk_basic_deg %>%
  pull(gene)
write.csv(nk_basic_deg, "07_outputs/nk_cell_deg_seurat.csv")

mono_basic_deg <- FindMarkers(object = pbmc_mock_rd2_annotated, subset.ident = "Monocyte", ident.1 = "control", ident.2 = "inhib", group.by = "orig.ident") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2)
mono_basic_deg$gene <- rownames(mono_basic_deg)
mono_deg_list <- mono_basic_deg %>%
  pull(gene)
write.csv(mono_basic_deg, "07_outputs/mono_cell_deg_seurat.csv")

mono_DC_basic_deg <- FindMarkers(object = pbmc_mock_rd2_annotated, subset.ident = c("Monocyte", "DC"), ident.1 = "control", ident.2 = "inhib", group.by = "orig.ident") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2)
mono_DC_basic_deg$gene <- rownames(mono_DC_basic_deg)
mono_DC_deg_list <- mono_DC_basic_deg %>%
  pull(gene)
write.csv(mono_basic_deg, "07_outputs/mono_DC_cell_deg_seurat.csv")

dc_basic_deg <- FindMarkers(object = pbmc_mock_rd2_annotated, subset.ident = "DC", ident.1 = "control", ident.2 = "inhib", group.by = "orig.ident") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2)
dc_basic_deg$gene <- rownames(dc_basic_deg)
dc_deg_list <- dc_basic_deg %>%
  pull(gene)
write.csv(dc_basic_deg, "07_outputs/dc_cell_deg_seurat.csv")

HSPC_basic_deg <- FindMarkers(object = pbmc_mock_rd2_annotated, subset.ident = "HSPC", ident.1 = "control", ident.2 = "inhib", group.by = "orig.ident") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2)
HSPC_basic_deg$gene <- rownames(HSPC_basic_deg)
HSPC_deg_list <- HSPC_basic_deg %>%
  pull(gene)

platlet_basic_deg <- FindMarkers(object = pbmc_mock_rd2_annotated, subset.ident = "Platlet", ident.1 = "control", ident.2 = "inhib", group.by = "orig.ident") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2)
platlet_basic_deg$gene <- rownames(platlet_basic_deg)
platlet_deg_list <- platlet_basic_deg %>%
  pull(gene)

shared_in_3plus_celltypes <- c("BTG1", "BTG2", "CD69", "JUN", "DDIT4", "DUSP1", "DUSP2", "IER2", "JUN", "JUNB", "MT-ATP8", "MTRNR2L12", "NFKBIA", "RPL26", "RPS29", "TNFAIP3", "ZFP36", "ZFP36L2")

upreg_shared3 <- c("BTG1", "BTG2", "CD69", "JUN", "DDIT4", "DUSP1", "DUSP2", "IER2", "JUN", "JUNB", "MT-ATP8", "NFKBIA", "RPS29", "TNFAIP3", "ZFP36", "ZFP36L2")

pbmc_mock_rd2_annotated <- AddModuleScore(object = pbmc_mock_rd2_annotated, features = list(upreg_shared3), name = "upreg_3")
