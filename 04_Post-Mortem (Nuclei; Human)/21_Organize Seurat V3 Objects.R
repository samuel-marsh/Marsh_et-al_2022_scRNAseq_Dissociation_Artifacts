# Load Packages 
library(tidyverse)
library(Seurat) # Starting with Seurat V3.1.5
library(patchwork)
library(marsh_post.utils)
library(viridis)
library(liger)
library(scCustomize)
library(beepr)

# Load Seurat V3 Objects --------------------------------------------------
marsh_seurat3 <- read_rds("RDS_SeuratV3/marsh_seuratv3_final.RDS")

morabito_seurat3 <- read_rds("RDS_SeuratV3/morabito_seuratv3_final.RDS")

zhou_seurat3 <- read_rds("RDS_SeuratV3/zhou_seuratv3_final.RDS")

leng_ec_seurat3 <- read_rds("RDS_SeuratV3/leng_ec_seuratv3_final.RDS")

leng_sfg_seurat3 <- read_rds("RDS_SeuratV3/leng_sfg_seuratv3_final.RDS")


# Load Annotation ---------------------------------------------------------
# leng SFG 
leng_sfg_final_annotation <- tibble::tribble(
  ~cluster,     ~cell_type,             ~color,
  0L,         "Excitatory",       "dodgerblue",
  1L,    "Oligodendrocyte",           "orange",
  2L,         "Inhibitory",             "navy",
  3L,          "Astrocyte",      "forestgreen",
  4L,    "Oligodendrocyte",           "orange",
  5L,         "Excitatory",       "dodgerblue",
  6L,          "Microglia",             "gold",
  7L,                "OPC",      "darkorange2",
  8L,        "Endothelial",           "orchid",
  9L,          "Astrocyte",      "forestgreen"
)

leng_sfg_colors <- leng_sfg_final_annotation %>% 
  pull(color)

# leng_ec 
leng_ec_final_annotation <- tibble::tribble(
  ~cluster,        ~cell_type,        ~color,
  0L,       "Oligodendrocyte",      "orange",
  1L,             "Microglia",        "gold",
  2L,            "Excitatory",  "dodgerblue",
  3L,                   "OPC", "darkorange2",
  4L,             "Astrocyte", "forestgreen",
  5L,            "Inhibitory",        "navy",
  6L,       "Oligodendrocyte",      "orange",
  7L,            "Inhibitory",        "navy",
  8L,             "Astrocyte", "forestgreen",
  9L,            "Excitatory",  "dodgerblue",
  10L,           "Excitatory",  "dodgerblue",
  11L,           "Excitatory",  "dodgerblue",
  12L,      "Endo/Fibro/Peri",      "orchid"
)

leng_ec_colors <- leng_ec_final_annotation %>% 
  pull(color)

# morabito 
morabito_final_annotation <- tibble::tribble(
  ~cluster,       ~cell_type,        ~color,
  0L,      "Oligodendrocyte",      "orange",
  1L,      "Oligodendrocyte",      "orange",
  2L,           "Excitatory",  "dodgerblue",
  3L,            "Astrocyte", "forestgreen",
  4L,           "Excitatory",  "dodgerblue",
  5L,           "Excitatory",  "dodgerblue",
  6L,           "Inhibitory",        "navy",
  7L,                  "OPC", "darkorange2",
  8L,           "Inhibitory",        "navy",
  9L,      "Oligodendrocyte",      "orange",
  10L,          "Excitatory",  "dodgerblue",
  11L,           "Microglia",        "gold",
  12L,          "Inhibitory",        "navy",
  13L,           "Astrocyte", "forestgreen",
  14L,          "Excitatory",  "dodgerblue",
  15L,          "Inhibitory",        "navy",
  16L,          "Excitatory",  "dodgerblue",
  17L,           "Astrocyte", "forestgreen",
  18L,          "Endo/Fibro",      "orchid",
  19L,     "Oligodendrocyte",      "orange",
  20L,          "Excitatory",  "dodgerblue",
  21L,          "Excitatory",  "dodgerblue",
  22L,     "Oligodendrocyte",      "orange",
  23L,          "Endo/Fibro",      "orchid"
)


morabito_colors <- morabito_final_annotation %>% 
  pull(color)

# Zhou 
zhou_final_annotation <- tibble::tribble(
  ~cluster,  ~cell_type,        ~color,
  0L,     "Oligodendrocyte",      "orange",
  1L,     "Oligodendrocyte",      "orange",
  2L,     "Oligodendrocyte",      "orange",
  3L,          "Excitatory",  "dodgerblue",
  4L,          "Excitatory",  "dodgerblue",
  5L,          "Inhibitory",        "navy",
  6L,           "Astrocyte", "forestgreen",
  7L,          "Excitatory",  "dodgerblue",
  8L,          "Excitatory",  "dodgerblue",
  9L,          "Excitatory",  "dodgerblue",
  10L,         "Excitatory",  "dodgerblue",
  11L,                "OPC", "darkorange2",
  12L,          "Microglia",        "gold",
  13L,          "Astrocyte", "forestgreen",
  14L,         "Excitatory",  "dodgerblue",
  15L,         "Excitatory",  "dodgerblue",
  16L,         "Inhibitory",        "navy",
  17L,    "Endo/Fibro/Peri",      "orchid",
  18L,    "Oligodendrocyte",      "orange",
  19L,        "Excitatory",   "dodgerblue"
)

zhou_colors <- zhou_final_annotation %>% 
  pull(color)

# marsh 
marsh_final_annotation <- tibble::tribble(
  ~cluster, ~cell_type,        ~color,
  0L,    "Oligodendrocyte",      "orange",
  1L,    "Oligodendrocyte",      "orange",
  2L,    "Oligodendrocyte",      "orange",
  3L,         "Excitatory",  "dodgerblue",
  4L,         "Excitatory",  "dodgerblue",
  5L,          "Astrocyte", "forestgreen",
  6L,          "Microglia",        "gold",
  7L,          "Astrocyte", "forestgreen",
  8L,         "Inhibitory",        "navy",
  9L,         "Inhibitory",        "navy",
  10L,               "OPC", "darkorange2",
  11L,       "Endothelial",      "orchid",
  12L,              "PBMC",        "gray",
  13L,        "Fibroblast", "darkorchid3"
)

marsh_colors <- marsh_final_annotation %>% 
  pull(color)

# Add Mito ----------------------------------------------------------------
marsh_seurat3 <- PercentageFeatureSet(marsh_seurat3, pattern = "^MT", col.name = "percent_mito")
morabito_seurat3 <- PercentageFeatureSet(morabito_seurat3, pattern = "^MT", col.name = "percent_mito")
zhou_seurat3 <- PercentageFeatureSet(zhou_seurat3, pattern = "^MT", col.name = "percent_mito")
leng_ec_seurat3 <- PercentageFeatureSet(leng_ec_seurat3, pattern = "^MT", col.name = "percent_mito")
leng_sfg_seurat3 <- PercentageFeatureSet(leng_sfg_seurat3, pattern = "^MT", col.name = "percent_mito")

# Rename Clusters ---------------------------------------------------------
# Combine clusters
marsh_new_ids <- marsh_final_annotation %>% 
  pull(cell_type)
names(marsh_new_ids) <- levels(marsh_seurat3)
marsh_seurat3 <- RenameIdents(marsh_seurat3, marsh_new_ids)

zhou_new_ids <- zhou_final_annotation %>% 
  pull(cell_type)
names(zhou_new_ids) <- levels(zhou_seurat3)
zhou_seurat3 <- RenameIdents(zhou_seurat3, zhou_new_ids)

morabito_new_ids <- morabito_final_annotation %>% 
  pull(cell_type)
names(morabito_new_ids) <- levels(morabito_seurat3)
morabito_seurat3 <- RenameIdents(morabito_seurat3, morabito_new_ids)

leng_ec_new_ids <- leng_ec_final_annotation %>% 
  pull(cell_type)
names(leng_ec_new_ids) <- levels(leng_ec_seurat3)
leng_ec_seurat3 <- RenameIdents(leng_ec_seurat3, leng_ec_new_ids)

leng_sfg_new_ids <- leng_sfg_final_annotation %>% 
  pull(cell_type)
names(leng_sfg_new_ids) <- levels(leng_sfg_seurat3)
leng_sfg_seurat3 <- RenameIdents(leng_sfg_seurat3, leng_sfg_new_ids)

# Relevel New Idents & Plot ------------------------------------------------------
# marsh
QC_Plots_Genes(marsh_seurat3, pt.size = 0)
marsh_levels <- c("Excitatory", "Inhibitory", "Astrocyte", "OPC", "Fibroblast", "Endothelial", "Oligodendrocyte", "Microglia", "PBMC")

Idents(marsh_seurat3) <- factor(Idents(marsh_seurat3), levels= marsh_levels)
QC_Plots_Genes(marsh_seurat3, pt.size = 0)

marsh_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid", "orange", "gold", "gray")
QC_Plots_Genes(marsh_seurat3, pt.size = 0, cols = marsh_renamed_colors)

# Zhou Data
QC_Plots_Genes(zhou_seurat3, pt.size = 0)
zhou_levels <- c("Excitatory", "Inhibitory", "Astrocyte", "OPC", "Endo/Fibro/Peri", "Oligodendrocyte", "Microglia")

Idents(zhou_seurat3) <- factor(Idents(zhou_seurat3), levels= zhou_levels)
QC_Plots_Genes(zhou_seurat3, pt.size = 0)

zhou_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
QC_Plots_Genes(zhou_seurat3, pt.size = 0, cols = zhou_renamed_colors)

# Morabito
QC_Plots_Genes(morabito_seurat3, pt.size = 0)
morabito_levels <- c("Excitatory", "Inhibitory", "Astrocyte", "OPC", "Endo/Fibro", "Oligodendrocyte", "Microglia")

Idents(morabito_seurat3) <- factor(Idents(morabito_seurat3), levels= morabito_levels)
QC_Plots_Genes(morabito_seurat3, pt.size = 0)

morabito_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
QC_Plots_Genes(morabito_seurat3, pt.size = 0, cols = morabito_renamed_colors)

# Leng EC
QC_Plots_Genes(leng_ec_seurat3, pt.size = 0)
leng_ec_levels <- c("Excitatory", "Inhibitory", "Astrocyte", "OPC", "Endo/Fibro/Peri", "Oligodendrocyte", "Microglia")

Idents(leng_ec_seurat3) <- factor(Idents(leng_ec_seurat3), levels= leng_ec_levels)
QC_Plots_Genes(leng_ec_seurat3, pt.size = 0)

leng_ec_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
QC_Plots_Genes(leng_ec_seurat3, pt.size = 0, cols = leng_ec_renamed_colors)

# Leng SFG
QC_Plots_Genes(leng_sfg_seurat3, pt.size = 0)
leng_sfg_levels <- c("Excitatory", "Inhibitory", "Astrocyte", "OPC", "Endothelial", "Oligodendrocyte", "Microglia")

Idents(leng_sfg_seurat3) <- factor(Idents(leng_sfg_seurat3), levels= leng_sfg_levels)
QC_Plots_Genes(leng_sfg_seurat3, pt.size = 0)

leng_sfg_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
QC_Plots_Genes(leng_sfg_seurat3, pt.size = 0, cols = leng_sfg_renamed_colors)

beep(sound = 2)

# Add Subject Meta Data ---------------------------------------------------

# Read all meta
meta_overview <- read.csv("03_Meta Data/Human_Meta_Data.csv", stringsAsFactors = FALSE)

# Marsh meta add 
# Seurat Meta
seurat_meta <- marsh_seurat3@meta.data
sample_ids <- unique(marsh_seurat3@meta.data$orig.ident)

# Filter by sample meta_data
meta_filter <- meta_overview %>% 
  filter(Sample_ID_Analysis %in% sample_ids)

# Move barcodes to column
seurat_meta <- seurat_meta %>% 
  rownames_to_column("barcodes")

joined_meta <- inner_join(seurat_meta, meta_filter, by = c("orig.ident" = "Sample_ID_Analysis"))

# move barcodes back to rownames
joined_meta <- joined_meta %>% 
  column_to_rownames("barcodes") %>% 
  select(-c("orig.ident", "nCount_RNA", "nFeature_RNA"))

# Add back to Seurat
marsh_seurat3 <- AddMetaData(object = marsh_seurat3, metadata = joined_meta)

# zhou meta add 
# Seurat Meta
seurat_meta <- zhou_seurat3@meta.data
sample_ids <- unique(zhou_seurat3@meta.data$orig.ident)

# Filter by sample meta_data
meta_filter <- meta_overview %>% 
  filter(Sample_ID_Analysis %in% sample_ids)

# Move barcodes to column
seurat_meta <- seurat_meta %>% 
  rownames_to_column("barcodes")

joined_meta <- inner_join(seurat_meta, meta_filter, by = c("orig.ident" = "Sample_ID_Analysis"))

# move barcodes back to rownames
joined_meta <- joined_meta %>% 
  column_to_rownames("barcodes") %>% 
  select(-c("orig.ident", "nCount_RNA", "nFeature_RNA"))

# Add back to Seurat
zhou_seurat3 <- AddMetaData(object = zhou_seurat3, metadata = joined_meta)

# morabito meta add 
# Seurat Meta
seurat_meta <- morabito_seurat3@meta.data
sample_ids <- unique(morabito_seurat3@meta.data$orig.ident)

# Filter by sample meta_data
meta_filter <- meta_overview %>% 
  filter(Sample_ID_Analysis %in% sample_ids)

# Move barcodes to column
seurat_meta <- seurat_meta %>% 
  rownames_to_column("barcodes")

joined_meta <- inner_join(seurat_meta, meta_filter, by = c("orig.ident" = "Sample_ID_Analysis"))

# move barcodes back to rownames
joined_meta <- joined_meta %>% 
  column_to_rownames("barcodes") %>% 
  select(-c("orig.ident", "nCount_RNA", "nFeature_RNA"))

# Add back to Seurat
morabito_seurat3 <- AddMetaData(object = morabito_seurat3, metadata = joined_meta)

# leng_ec meta add 
# Seurat Meta
seurat_meta <- leng_ec_seurat3@meta.data
sample_ids <- unique(leng_ec_seurat3@meta.data$orig.ident)

# Filter by sample meta_data
meta_filter <- meta_overview %>% 
  filter(Sample_ID_Analysis %in% sample_ids)

# Move barcodes to column
seurat_meta <- seurat_meta %>% 
  rownames_to_column("barcodes")

joined_meta <- inner_join(seurat_meta, meta_filter, by = c("orig.ident" = "Sample_ID_Analysis"))

# move barcodes back to rownames
joined_meta <- joined_meta %>% 
  column_to_rownames("barcodes") %>% 
  select(-c("orig.ident", "nCount_RNA", "nFeature_RNA"))

# Add back to Seurat
leng_ec_seurat3 <- AddMetaData(object = leng_ec_seurat3, metadata = joined_meta)

# leng_sfg meta add 
# Seurat Meta
seurat_meta <- leng_sfg_seurat3@meta.data
sample_ids <- unique(leng_sfg_seurat3@meta.data$orig.ident)

# Filter by sample meta_data
meta_filter <- meta_overview %>% 
  filter(Sample_ID_Analysis %in% sample_ids)

# Move barcodes to column
seurat_meta <- seurat_meta %>% 
  rownames_to_column("barcodes")

joined_meta <- inner_join(seurat_meta, meta_filter, by = c("orig.ident" = "Sample_ID_Analysis"))

# move barcodes back to rownames
joined_meta <- joined_meta %>% 
  column_to_rownames("barcodes") %>% 
  select(-c("orig.ident", "nCount_RNA", "nFeature_RNA"))

# Add back to Seurat
leng_sfg_seurat3 <- AddMetaData(object = leng_sfg_seurat3, metadata = joined_meta)

# Save renamed & Meta Added Objects ----------------------------------------------------
write_rds(marsh_seurat3, "RDS_SeuratV3/rename_meta_final/marsh_seuratv3_RENAMED.RDS")
write_rds(zhou_seurat3, "RDS_SeuratV3/rename_meta_final/zhou_seuratv3_RENAMED.RDS")
write_rds(morabito_seurat3, "RDS_SeuratV3/rename_meta_final/morabito_seuratv3_RENAMED.RDS")
write_rds(leng_ec_seurat3, "RDS_SeuratV3/rename_meta_final/leng_ec_seuratv3_RENAMED.RDS")
write_rds(leng_sfg_seurat3, "RDS_SeuratV3/rename_meta_final/leng_sfg_seuratv3_RENAMED.RDS")