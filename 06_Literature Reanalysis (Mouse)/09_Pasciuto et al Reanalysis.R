# Analysis performed using Seurat 3.2.3

# Load R Data files -------------------------------------------------------
load("05_Mendeley Data/hsmzw47kbg-3/Data Resource 1/Supplementary resource 1/WT01.Rdata")
load("05_Mendeley Data/hsmzw47kbg-3/Data Resource 1/Supplementary resource 1/WT02.Rdata")
load("05_Mendeley Data/hsmzw47kbg-3/Data Resource 1/Supplementary resource 1/KO.Rdata")

# Create data list
data_list <- list(WT01.data, WT02.data, KO.data)
names(data_list) <- c("WT1", "WT2", "KO1")

# remove individual data objects
rm(KO.data)
rm(WT01.data)
rm(WT02.data)
gc()

# Merge and create Seurat Object
pasciuto_merged <- Merge_Sparse_Data_All(matrix_list = data_list)

pasciuto_seurat <- CreateSeuratObject(counts = pasciuto_merged, min.cells = 5, min.features = 200)

# Add meta data
pasciuto_seurat@meta.data$genotype[pasciuto_seurat@meta.data$orig.ident == "WT01" | pasciuto_seurat@meta.data$orig.ident == "WT02"] <- "WT"
pasciuto_seurat@meta.data$genotype[pasciuto_seurat@meta.data$orig.ident == "KO"] <- "KO"

pasciuto_seurat <- Add_Mito_Ribo_Seurat(seurat_object = pasciuto_seurat, species = "mouse")

# Analysis Round 1 --------------------------------------------------------
# Normalize RNA data with log normalization
pasciuto_seurat <- NormalizeData(pasciuto_seurat)

# Find and scale variable features
pasciuto_seurat <- FindVariableFeatures(pasciuto_seurat, selection.method = "vst", nfeatures = 1000)

pasciuto_seurat <- ScaleData(pasciuto_seurat, features = VariableFeatures(pasciuto_seurat))
pasciuto_seurat <- RunPCA(object = pasciuto_seurat, npcs = 50)
ElbowPlot(pasciuto_seurat, ndims = 40)

# Cluster and Visualize Data
pasciuto_seurat <- FindNeighbors(object = pasciuto_seurat, dims = 1:20)
pasciuto_seurat <- FindClusters(object = pasciuto_seurat, resolution = 0.3)
pasciuto_seurat <- RunUMAP(object = pasciuto_seurat, dims = 1:20)

FeaturePlot_scCustom(seurat_object = pasciuto_seurat, features = c("P2ry12", "Ms4a7", "Ccr2", "Ms4a4c", "S100a9", "Mrc1"))

# Remove monocytes, macrophages, and neutrophils
pasciuto_seurat_sub <- subset(x = pasciuto_seurat, idents = c(4, 6, 8), invert = T)

# Analysis Round 2 ----------------------------------------------------------------------------
# Normalize RNA data with log normalization
pasciuto_seurat_sub <- NormalizeData(pasciuto_seurat_sub)

# Find and scale variable features
pasciuto_seurat_sub <- FindVariableFeatures(pasciuto_seurat_sub, selection.method = "vst", nfeatures = 1000)

pasciuto_seurat_sub <- ScaleData(pasciuto_seurat_sub, features = VariableFeatures(pasciuto_seurat_sub), vars.to.regress = c("percent_mito", "nCount_RNA"))

pasciuto_seurat_sub <- RunPCA(object = pasciuto_seurat_sub, npcs = 50)
ElbowPlot(pasciuto_seurat_sub, ndims = 40)
pasciuto_seurat_sub <- JackStraw(pasciuto_seurat_sub, num.replicate = 100)
pasciuto_seurat_sub <- ScoreJackStraw(pasciuto_seurat_sub, dims = 1:20)
JackStrawPlot(pasciuto_seurat_sub, dims = 1:20)

# Cluster & Visualize Data
pasciuto_seurat_sub <- FindNeighbors(object = pasciuto_seurat_sub, dims = 1:18)
pasciuto_seurat_sub <- FindClusters(object = pasciuto_seurat_sub, resolution = 0.2)
pasciuto_seurat_sub <- RunUMAP(object = pasciuto_seurat_sub, dims = 1:18)
DimPlot(pasciuto_seurat_sub, cols = DiscretePalette(n = 24), label = T, label.size = 6)

# Add exAM score
pasciuto_seurat_sub <- AddModuleScore(object = pasciuto_seurat_sub, features = shared_sig, name = "exAM_score")

# Save Object ---------------------------------------------------------------------------------
write_rds(pasciuto_seurat_sub, "Pasciuto_Mouse_Microglia_Final.RDS")