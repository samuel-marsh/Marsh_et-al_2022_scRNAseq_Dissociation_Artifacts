# Load Data -----------------------------------------------------------
raw_data_list <- Read_GEO_Delim(data_dir = "/Volumes/marsh_local_2Big/Literature Datasets/Human DataSets/Pasciuto_CD4T_Microglia_2020/Brain_vs_Blood_TCells/GSE146165_RAW/", file_suffix = ".txt.gz", move_genes_rownames = T)

raw_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = c("Brain", "PBMC"))

brain_pbmc_seurat <- CreateSeuratObject(counts = raw_merged, min.cells = 1, min.features = 200)

# Add mito percent
brain_pbmc_seurat <- Add_Mito_Ribo_Seurat(seurat_object = pasciuto_seurat, species = "human")

# QC Filter and Subset
brain_pbmc_seurat <- subset(x = brain_pbmc_seurat, subset = nCount_RNA < 15000 & percent_mito < 15)

# Analysis --------------------------------------------------------
# Normalize RNA data with log normalization
brain_pbmc_seurat <- NormalizeData(brain_pbmc_seurat)
# Find and scale variable features
brain_pbmc_seurat <- FindVariableFeatures(brain_pbmc_seurat, selection.method = "vst", nfeatures = 1000)

brain_pbmc_seurat <- ScaleData(brain_pbmc_seurat, features = VariableFeatures(brain_pbmc_seurat))

# PCA & Clustering
brain_pbmc_seurat <- RunPCA(object = brain_pbmc_seurat, npcs = 50)
ElbowPlot(brain_pbmc_seurat, ndims = 40)
brain_pbmc_seurat <- JackStraw(brain_pbmc_seurat, num.replicate = 100)
brain_pbmc_seurat <- ScoreJackStraw(brain_pbmc_seurat, dims = 1:20)
JackStrawPlot(brain_pbmc_seurat, dims = 1:15)

brain_pbmc_seurat <- FindNeighbors(object = brain_pbmc_seurat, dims = 1:10)
brain_pbmc_seurat <- FindClusters(object = brain_pbmc_seurat, resolution = 0.6)
brain_pbmc_seurat <- RunUMAP(object = brain_pbmc_seurat, dims = 1:10)
DimPlot(brain_pbmc_seurat, cols = DiscretePalette(n = 6), label = T, label.size = 6)
DimPlot(brain_pbmc_seurat, cols = c("navy", "orange"), group.by = "orig.ident", shuffle = T)

# Add additional meta.data column
brain_pbmc_seurat[["tissue"]] <- brain_pbmc_seurat@meta.data$orig.ident

brain_pbmc_seurat@meta.data$tissue <- factor(brain_pbmc_seurat@meta.data$tissue, levels = c("PBMC", "Brain"))

# Convert to human gene symbols (homologous gene names previously confirmed so just using simple gene name capitalization conversion)
human_shared_sig <- lapply(shared_sig, function(x){str_to_upper(x)})

brain_pbmc_seurat <- AddModuleScore(object = brain_pbmc_seurat, features = list(shared_sig), name = "exAM Score", search = TRUE)

# Save Object ---------------------------------------------------------------------------------
write_rds(brain_pbmc_seurat, "Pasciuto_Human_TCell_Final.RDS")