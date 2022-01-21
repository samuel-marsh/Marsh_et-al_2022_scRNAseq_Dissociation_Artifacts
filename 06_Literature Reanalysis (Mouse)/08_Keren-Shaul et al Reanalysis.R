# Analysis performed using Seurat 3.2.3

# Filter Experiment Data of Interest ----------------------------------------------------------
# Create meta data csv with list of plates from only from Trem2KO experiment
# csv file can be found here: https://github.com/samuel-marsh/Marsh_et-al_2022_scRNAseq_Dissociation_Artifacts/blob/master/06_Literature%20Reanalysis%20(Mouse)/Supplemental_CSV_Files/keren_shaul_meta_organized_trem2_exp_only.csv

# Load Data and Create Seurat Object ----------------------------------------------------------
# Read in organized meta data
keren_shaul_meta <- read.csv(file = "keren_shaul_meta_organized_trem2_exp_only.csv", stringsAsFactors = F)

# Create file names for read in
keren_shaul_meta <- keren_shaul_meta %>%
  mutate(file_name = paste0(GEO_ID, "_", Plate_ID))

# Files to import
file_names_import <- keren_shaul_meta %>%
  pull(file_name)

keren_shaul_trem2_raw_data <- Read_GEO_Delim(data_dir = "~/GEO_Downloads/GSE98969_RAW/", file_suffix = ".txt.gz", sample_list = file_names_import, full_names = FALSE, move_genes_rownames = T)

# Create single matrix
plate_ids <- keren_shaul_meta %>%
  pull(Plate_ID)
plate_ids

keren_shaul_raw_merged <- Merge_Sparse_Data_All(matrix_list = keren_shaul_trem2_raw_data, add_cell_ids = plate_ids)

# Remove ERCC Genes
# Calculate ERCC abundances on the raw counts before creating a Seurat object
ERCC.WT.index <- grep(pattern = "^ERCC-", x = rownames(keren_shaul_raw_merged), value = FALSE) # Select row indices and not ERCC names
percent.ERCC.WT <- Matrix::colSums(keren_shaul_raw_merged[ERCC.WT.index, ])/Matrix::colSums(keren_shaul_raw_merged)
# Remove ERCC from count.data
keren_shaul_raw_merged_noERCC <- keren_shaul_raw_merged[-ERCC.WT.index, ]
keren_shaul_seurat <- CreateSeuratObject(counts = keren_shaul_raw_merged_noERCC, min.cells = 5, min.features = 200)

# Add meta data to Seurat Object
keren_shaul_meta <- keren_shaul_meta %>%
  mutate(AD_Trem2_Genotype = paste0(AD_Genotype, "_", Trem2_Genotype))

# Pull object meta
seurat_meta <- keren_shaul_seurat@meta.data

seurat_meta <- seurat_meta %>%
  rownames_to_column("barcodes")

joined_meta <- left_join(seurat_meta, keren_shaul_meta, by = c("orig.ident" = "Plate_ID"))

joined_meta <- joined_meta %>%
  select(-orig.ident, -nCount_RNA, -nFeature_RNA) %>%
  column_to_rownames("barcodes")

keren_shaul_seurat <- AddMetaData(object = keren_shaul_seurat, metadata = joined_meta)

# Add percent mito/ribo
keren_shaul_seurat <- Add_Mito_Ribo_Seurat(seurat_object = keren_shaul_seurat, species = "mouse")

# QC Filter
keren_shaul_seurat <- subset(x = keren_shaul_seurat, subset = nCount_RNA < 7500 & nFeature_RNA < 2500)


# Analysis Round 1 ----------------------------------------------------------------------------
keren_shaul_seurat <- NormalizeData(object = keren_shaul_seurat)
keren_shaul_seurat <- FindVariableFeatures(object = keren_shaul_seurat, selection.method = "vst", nfeatures = 2000)
keren_shaul_seurat <- ScaleData(object = keren_shaul_seurat, features = rownames(keren_shaul_seurat))
keren_shaul_seurat <- RunPCA(object = keren_shaul_seurat, npcs = 75)
keren_shaul_seurat <- JackStraw(object = keren_shaul_seurat, num.replicate = 100)
keren_shaul_seurat <- ScoreJackStraw(object = keren_shaul_seurat, dims = 1:20)
JackStrawPlot(object = keren_shaul_seurat, dims = 1:20)
ElbowPlot(object = keren_shaul_seurat, ndims = 30)
beep(sound = 2)
JackStrawPlot(object = keren_shaul_seurat, dims = 1:20)
keren_shaul_seurat <- FindNeighbors(object = keren_shaul_seurat, dims = 1:15)
keren_shaul_seurat <- FindClusters(object = keren_shaul_seurat, resolution = 0.8)
keren_shaul_seurat <- RunUMAP(object = keren_shaul_seurat, dims = 1:15)
DimPlot(object = keren_shaul_seurat, cols = DiscretePalette(n = 24, palette = "polychrome"), label = T, label.size = 6)


# Filter and subset the data ------------------------------------------------------------------
# Remove contaminating cell types (Macrophages/Monocytes, Astrocytes, Peripheral Immune Cells)
keren_shaul_seurat_filtered <- subset(x = keren_shaul_seurat, idents = c(5, 7, 9), invert = TRUE)
keren_shaul_seurat_filtered <- DietSeurat(object = keren_shaul_seurat_filtered)

# subset to remove Trem2KO samples
keren_shaul_seurat_wt_only <- subset(x = keren_shaul_seurat_filtered, subset = Trem2_Genotype == "WT")

# Analysis Round 2 -----------------------------------------------------------------
keren_shaul_seurat_wt_only <- NormalizeData(object = keren_shaul_seurat_wt_only)
keren_shaul_seurat_wt_only <- FindVariableFeatures(object = keren_shaul_seurat_wt_only, selection.method = "vst", nfeatures = 2000)
keren_shaul_seurat_wt_only <- ScaleData(object = keren_shaul_seurat_wt_only, features = rownames(keren_shaul_seurat_wt_only))
keren_shaul_seurat_wt_only <- RunPCA(object = keren_shaul_seurat_wt_only, npcs = 75)
keren_shaul_seurat_wt_only <- JackStraw(object = keren_shaul_seurat_wt_only, num.replicate = 100)
keren_shaul_seurat_wt_only <- ScoreJackStraw(object = keren_shaul_seurat_wt_only, dims = 1:20)
JackStrawPlot(object = keren_shaul_seurat_wt_only, dims = 1:20)
ElbowPlot(object = keren_shaul_seurat_wt_only, ndims = 30)
beep(sound = 2)
keren_shaul_seurat_wt_only <- FindNeighbors(object = keren_shaul_seurat_wt_only, dims = 1:4)
keren_shaul_seurat_wt_only <- FindClusters(object = keren_shaul_seurat_wt_only, resolution = 1)
keren_shaul_seurat_wt_only <- RunUMAP(object = keren_shaul_seurat_wt_only, dims = 1:4)

# Annotate the data
# Annotation CSV can be found here: https://github.com/samuel-marsh/Marsh_et-al_2022_scRNAseq_Dissociation_Artifacts/blob/master/06_Literature%20Reanalysis%20(Mouse)/Supplemental_CSV_Files/Keren_Shaul_Cluster_Annotation.csv
# Pull data from CSV
wt_only_cluster_annotation <- Pull_Cluster_Annotation(seurat_object = keren_shaul_seurat_wt_only, read_from_file = T, file_name = "Keren_Shaul_Cluster_Annotation.csv")

# rename idents
keren_shaul_seurat_wt_only_renamed <- RenameIdents(object = keren_shaul_seurat_wt_only, wt_only_cluster_annotation[["new_cluster_idents"]])

# relevel the idents
wt_only_new_levels <- c("Homeostatic", "Ifn-Responsive", "AD non-DAM", "DAMa", "DAMb")

Idents(keren_shaul_seurat_wt_only_renamed) <- factor(Idents(keren_shaul_seurat_wt_only_renamed), levels = wt_only_new_levels)

# Add exAM score
keren_shaul_seurat_wt_only <- AddModuleScore(object = keren_shaul_seurat_wt_only, features = shared_sig, name = "exAM_score")

DimPlot(keren_shaul_seurat_wt_only_renamed, label = T, cols = DiscretePalette(n = 24, palette = "polychrome"))

# Save data -----------------------------------------------------------------------------------
write_rds(keren_shaul_seurat_wt_only_renamed, "Kren-Shaul_Microglia_Final.RDS")