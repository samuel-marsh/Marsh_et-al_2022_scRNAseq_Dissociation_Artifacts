# 1.0 Load packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# 2.0 Load V3 All Datasets ----------------------------------------------------
marsh_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/marsh_seuratv3_RENAMED.RDS")
zhou_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/zhou_seuratv3_RENAMED.RDS")
morabito_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/morabito_seuratv3_RENAMED.RDS")
leng_ec_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/leng_ec_seuratv3_RENAMED.RDS")
leng_sfg_seurat3 <- read_rds("RDS_SeuratV3/rename_meta_final/leng_sfg_seuratv3_RENAMED_b.RDS")

beep(sound = 2)

# 3.0 Load Colors -------------------------------------------------------------
marsh_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid", "orange", "gold", "gray")
zhou_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
morabito_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
leng_ec_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")
leng_sfg_renamed_colors <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "orchid", "orange", "gold")

# 4.0 Load Gene Lists ---------------------------------------------------------
# 4.1 Load Mouse DEG Lists ------------------------------------------------
# Add Mouse list from SI Table XX
mouse_myeloid_list <- read.csv("Mouse_Data_Import/mouse_DEG_list_myeloid.csv", stringsAsFactors = FALSE, header = TRUE)

# Pull myeloid list
mouse_myeloid_list <- pull(mouse_myeloid_list, myeloid)

# Human homolog of Hist2h2aa1 is HIST2H2AA4
mouse_myeloid_list <- gsub("Hist2h2aa1", "HIST2H2AA4", mouse_myeloid_list)
print(mouse_myeloid_list)

# Convert to uppercase for human genes
mouse_myeloid_list_HUMAN <- mouse_myeloid_list %>% 
  str_to_upper()

# Pull all CNS list 
mouse_all_cns_list <- read.csv("Mouse_Data_Import/mouse_DEG_list_all_cns.csv", stringsAsFactors = FALSE, header = TRUE)

mouse_all_cns_list <- pull(mouse_all_cns_list, all_cns)
# Human homolog of 1500015O10Rik is C2orf40
mouse_all_cns_list <- gsub("1500015O10Rik", "C2orf40", mouse_all_cns_list)
print(mouse_all_cns_list)

mouse_all_cns_list_HUMAN <- mouse_all_cns_list %>% 
  str_to_upper() %>% 
  gsub("C2ORF40", "C2orf40", .)
print(mouse_all_cns_list_HUMAN)

# Create Combined Mouse Myeloid and All CNS list (No duplicates)
mouse_combined_list_HUMAN <- union(mouse_myeloid_list_HUMAN, mouse_all_cns_list_HUMAN)

# Remove mouse lists
rm(mouse_myeloid_list)
rm(mouse_all_cns_list)

# 4.2 Load LIGER Factor Lists ----------------------------------------------------------
liger_mg_list <- read_rds("subcluster_factor_gene_list/new_thresholded_lists/liger_micro_cutoff_gene_list.RDS")

liger_astro_list <- read_rds("subcluster_factor_gene_list/new_thresholded_lists/liger_astro_cutoff_gene_list.RDS")

liger_mg_list <- list(as.character(liger_mg_list[["micro_factor_cutoff"]]))
liger_astro_list <- list(as.character(liger_astro_list[["astro_factor_cutoff"]]))

beep(sound = 2)

# 5.0 Add Module Scores ----------------------------------------------------
# 5.1 Marsh Dataset -------------------------------------------------------
# Mouse Lists
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# 5.2 zhou Dataset -------------------------------------------------------
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# 5.3 Morabito Dataset -------------------------------------------------------
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# 5.4 leng ec Dataset -------------------------------------------------------
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# 5.5 leng sfg Dataset -------------------------------------------------------
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(mouse_myeloid_list_HUMAN), name = "mouse_myeloid", search = TRUE)
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(mouse_all_cns_list_HUMAN), name = "mouse_all_cns", search = TRUE)
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(mouse_combined_list_HUMAN), name = "mouse_combined", search = TRUE)

# 5.6 Add LIGER Factor Scores -------------------------------------------------------
# micro score
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(liger_mg_list), search = TRUE, name = "liger_mg_factor")
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(liger_mg_list), search = TRUE, name = "liger_mg_factor")
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(liger_mg_list), search = TRUE, name = "liger_mg_factor")
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(liger_mg_list), search = TRUE, name = "liger_mg_factor")
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(liger_mg_list), search = TRUE, name = "liger_mg_factor")

# astro score
marsh_seurat3 <- AddModuleScore(object = marsh_seurat3, features = list(liger_astro_list), search = TRUE, name = "liger_astro_factor")
zhou_seurat3 <- AddModuleScore(object = zhou_seurat3, features = list(liger_astro_list), search = TRUE, name = "liger_astro_factor")
morabito_seurat3 <- AddModuleScore(object = morabito_seurat3, features = list(liger_astro_list), search = TRUE, name = "liger_astro_factor")
leng_ec_seurat3 <- AddModuleScore(object = leng_ec_seurat3, features = list(liger_astro_list), search = TRUE, name = "liger_astro_factor")
leng_sfg_seurat3 <- AddModuleScore(object = leng_sfg_seurat3, features = list(liger_astro_list), search = TRUE, name = "liger_astro_factor")

# 6.0 Save with new scores added ----------------------------------------------
write_rds(marsh_seurat3, "RDS_SeuratV3/rename_scored/marsh_seuratv3_RENAMED_SCORED.RDS")
write_rds(zhou_seurat3, "RDS_SeuratV3/rename_scored/zhou_seuratv3_RENAMED_SCORED.RDS")
write_rds(morabito_seurat3, "RDS_SeuratV3/rename_scored/morabito_seuratv3_RENAMED_SCORED.RDS")
write_rds(leng_ec_seurat3, "RDS_SeuratV3/rename_scored/leng_ec_seuratv3_RENAMED_SCORED.RDS")
write_rds(leng_sfg_seurat3, "RDS_SeuratV3/rename_scored/leng_sfg_seuratv3_RENAMED_SCORED.RDS")

beep(sound = 2)