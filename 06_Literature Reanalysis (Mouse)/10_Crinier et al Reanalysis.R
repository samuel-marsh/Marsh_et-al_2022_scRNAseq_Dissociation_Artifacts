# Analysis performed using Seurat 3.2.3

# Load Data & Create Object -------------------------------------------
# Load data
raw_geo_all <- Read_GEO_Delim("~/Downloads/GSE119562_RAW/", file_suffix = ".csv.gz")

# Create sample meta data csv
all_samples <- tibble::tribble(
  ~GEO_ID,      ~Species, ~Tissue, ~Sample, ~Sample_Name,
  "GSM3377673", "human", "spleen", "donor1", "SpD1-Hu",
  "GSM3377674", "human", "spleen", "donor2", "SpD2-Hu",
  "GSM3377675", "human", "spleen", "donor3", "SpD3-Hu",
  "GSM3377676", "human", "blood",  "donor1", "BlD1-Hu",
  "GSM3377677", "human", "blood",  "donor2", "BlD2-Hu",
  "GSM3377678", "human", "blood",  "donor3", "BlD3-Hu",
  "GSM3377679", "mouse", "spleen", "group1", "SpG1-Ms",
  "GSM3377680", "mouse", "spleen", "group2", "SpG2-Ms",
  "GSM3377681", "mouse", "spleen", "group3", "SpG3-Ms",
  "GSM3377682", "mouse", "blood",  "group1", "BlG1-Ms",
  "GSM3377683", "mouse", "blood",  "group2", "BlG2-Ms",
  "GSM3377684", "mouse", "blood",  "group3", "BlG3-Ms"
)
write.csv(all_samples, "crinier_sample_meta_data.csv", row.names = F)

# Pull sample names
mouse_names <- crinier_meta_data %>%
  filter(Species == "mouse") %>%
  pull(Sample_Name)

# Merge data
mouse_raw <- raw_geo_all[7:12]

mouse_merged <- Merge_Sparse_Data_All(matrix_list = mouse_raw, add_cell_ids = mouse_names)

mouse_nk <- CreateSeuratObject(counts = mouse_merged, min.cells = 5, min.features = 200)

# Add relevant metadata
# Add Tissue Information
mouse_nk@meta.data$tissue[mouse_nk@meta.data$orig.ident == "SpG1-Ms" | mouse_nk@meta.data$orig.ident == "SpG2-Ms" | mouse_nk@meta.data$orig.ident == "SpG3-Ms"] <- "Spleen"
mouse_nk@meta.data$tissue[mouse_nk@meta.data$orig.ident == "BlG1-Ms" | mouse_nk@meta.data$orig.ident == "BlG2-Ms" | mouse_nk@meta.data$orig.ident == "BlG3-Ms"] <- "Blood"

mouse_nk <- Add_Mito_Ribo_Seurat(seurat_object = mouse_nk, species = "mouse")

# QC Filter
QC_Plots_Genes(seurat_object = mouse_nk, high_cutoff = 2500, low_cutoff = 200)
QC_Plots_UMIs(seurat_object = mouse_nk, high_cutoff = 7000)
QC_Plots_Mito(seurat_object = mouse_nk, high_cutoff = 10)

# Subset the data by QC metrics
mouse_nk_filtered <- subset(x = mouse_nk, subset = nFeature_RNA < 2500 & nFeature_RNA > 400 & percent_mito < 10 & nCount_RNA < 7000)

# Analysis --------------------------------------------------
# Normalize
mouse_nk_filtered <- NormalizeData(mouse_nk_filtered)
# Find and scale variable features
mouse_nk_filtered <- FindVariableFeatures(mouse_nk_filtered, selection.method = "vst", nfeatures = 1000)

mouse_nk_filtered <- ScaleData(mouse_nk_filtered, features = VariableFeatures(mouse_nk_filtered))
mouse_nk_filtered <- RunPCA(object = mouse_nk_filtered, npcs = 50)
ElbowPlot(mouse_nk_filtered, ndims = 40)
mouse_nk_filtered <- JackStraw(mouse_nk_filtered, num.replicate = 100)
mouse_nk_filtered <- ScoreJackStraw(mouse_nk_filtered, dims = 1:20)
JackStrawPlot(mouse_nk_filtered, dims = 1:20)

# Clustering & DimReduc
mouse_nk_filtered <- FindNeighbors(object = mouse_nk_filtered, dims = 1:19)
mouse_nk_filtered <- FindClusters(object = mouse_nk_filtered, resolution = 0.4)
mouse_nk_filtered <- RunUMAP(object = mouse_nk_filtered, dims = 1:19)

# Add exAM module score
mouse_nk_filtered <- AddModuleScore(object = mouse_nk_filtered, features = shared_sig, name = "exAM_score")

# Save Object ---------------------------------------------------------------------------------
write_rds(mouse_nk_filtered, "Crinier_Mouse_NK_Final.RDS")

