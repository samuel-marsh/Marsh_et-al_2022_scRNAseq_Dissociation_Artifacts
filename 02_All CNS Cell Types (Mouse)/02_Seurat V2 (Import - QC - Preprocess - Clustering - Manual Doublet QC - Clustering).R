# Seurat V2 Analysis -------------------------------------------------------------------
# Originally performed between 2019-08-28 and 2019-10-22 using packages as installed on or prior to start date
    # Originally performed using R 3.5.3
# Script validated using packrat and versions package with date 2019-08-28 (see prior script for setup)
    # Validated with R 3.6.1 using the same setup
    # Seurat 2.3.4 installation as described in previous script

# Load data
raw_data <- Read10X("~/Dropbox (BCH)/untitled folder/RAW_DATA_EXP_17/Seurat 2 Compat/")
# Modified Cell Ranger 3.0 output to rename the features file into genes for reading into Seurat 2.3.4

# Create Seurat object 
all_cns_cells <- CreateSeuratObject(raw.data = raw_data, min.cells = 10, min.genes = 200, project = "all_cns_cells", names.delim = "-", names.field = 2)

# Meta data
meta_data <- data.frame(all_cns_cells@meta.data)
meta_data_tb <- meta_data %>% 
  rownames_to_column(var = "cell_name") %>% 
  as_tibble()

# Create string of which samples are INHIB+
INHIB <- c("3", "4")

# Create new column based on string
meta_data_tb <- mutate(meta_data_tb, Treatment = ifelse(orig.ident %in% INHIB, "INHIB", "NO_INHIB"))

# Convert back to data frame
meta_data_mod <- meta_data_tb %>% 
  select("cell_name", "Treatment") %>% 
  data.frame() %>% 
  column_to_rownames(var = "cell_name")

# Remove unesscessary meta_data objects
rm(INHIB)
rm(meta_data)
rm(meta_data_tb)

# Add new metadata to seurat object
all_cns_cells <- AddMetaData(all_cns_cells,
                             metadata = meta_data_mod,
                             col.name = Treatment)

mito_genes <- grep(pattern = "^mt-", x = rownames(x = raw_data), value = TRUE)
percent_mito <- Matrix::colSums(all_cns_cells@raw.data[mito_genes, ]) / Matrix::colSums(all_cns_cells@raw.data)

all_cns_cells <- AddMetaData(object = all_cns_cells, metadata = percent_mito, col.name = "percent_mito")

# QC thresholds
all_cns_cells <- FilterCells(object = all_cns_cells, subset.names = c("nGene", "percent_mito", "nUMI"), low.thresholds = c(200, -Inf, -Inf), high.thresholds = c(4000, 0.25, 12000))

# normalize & var genes
all_cns_cells <- NormalizeData(object = all_cns_cells, normalization.method = "LogNormalize", scale.factor = 1e4, display.progress = TRUE)
all_cns_cells <- FindVariableGenes(object = all_cns_cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = TRUE)

# Scale the dataset & PCA
all_cns_cells <- ScaleData(object = all_cns_cells, vars.to.regress = c("nUMI", "percent_mito"))
all_cns_cells <- RunPCA(object = all_cns_cells, pc.genes = all_cns_cells@var.genes, pcs.compute = 75)
all_cns_cells <- ProjectPCA(object = all_cns_cells, do.print = FALSE)

# Cluster and Plot
all_cns_cells <- FindClusters(object = all_cns_cells, reduction.type = "pca", dims.use = 1:40, resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
all_cns_cells <- RunTSNE(object = all_cns_cells, dims.use = 1:40)
TSNEPlot(object = all_cns_cells, do.label = TRUE)

# Visualize low wuality clusters 2nd way
VlnPlot(all_cns_cells, features.plot = "nGene", x.lab.rot = TRUE, do.return = TRUE) + 
  geom_hline(yintercept = 600, linetype = "dashed", color = "red")

# Subset out the low quality clusters
all_cns_cells <- SubsetData(all_cns_cells, ident.remove = c(7, 10, 17), do.clean = TRUE)


# Reanalyze without the low quality cells ---------------------------------
all_cns_cells <- NormalizeData(object = all_cns_cells, normalization.method = "LogNormalize", scale.factor = 1e4, display.progress = TRUE)

# Find Variable genes in the data
all_cns_cells <- FindVariableGenes(object = all_cns_cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = TRUE)

all_cns_cells <- ScaleData(object = all_cns_cells, vars.to.regress = c("nUMI", "percent_mito"))

all_cns_cells <- RunPCA(object = all_cns_cells, pc.genes = all_cns_cells@var.genes, pcs.compute = 75)

all_cns_cells <- ProjectPCA(object = all_cns_cells, do.print = FALSE)

all_cns_cells <- FindClusters(object = all_cns_cells, reduction.type = "pca", dims.use = 1:40, resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

all_cns_cells <- RunTSNE(object = all_cns_cells, dims.use = 1:40)

TSNEPlot(object = all_cns_cells, do.label = TRUE)


# Examination of gene expression revealed likely microglial contamination/doublets in other clusters.  Remove contaminating cells before final analyses

# Manual filtering --------------------------------------------------------
# Astrocytes
astroytes_subset <- SubsetData(all_cns_cells, ident.use = c(2, 4), do.clean = TRUE)
astroytes_subset <- NormalizeData(object = astroytes_subset, normalization.method = "LogNormalize", scale.factor = 1e4, display.progress = TRUE)
astroytes_subset <- FindVariableGenes(object = astroytes_subset, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = TRUE)
astroytes_subset <- ScaleData(object = astroytes_subset, vars.to.regress = c("nUMI", "percent_mito"))
astroytes_subset <- RunPCA(object = astroytes_subset, pc.genes = astroytes_subset@var.genes, pcs.compute = 75)
astroytes_subset <- ProjectPCA(object = astroytes_subset, do.print = FALSE)
astroytes_subset <- FindClusters(object = astroytes_subset, reduction.type = "pca", dims.use = 1:20, resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
astroytes_subset <- RunTSNE(object = astroytes_subset, dims.use = 1:20)
TSNEPlot(object = astroytes_subset, do.label = TRUE)
beep(sound = 2)

astroytes_filtered <- FilterCells(astroytes_subset, subset.names = "Tmem119", high.thresholds = 0.1)
astroytes_filtered <- SubsetData(astroytes_filtered, ident.remove = 2, do.clean = TRUE)

# Endothelial/Pericytes
endo_subset <- SubsetData(all_cns_cells, ident.use = c(1, 8), do.clean = TRUE)
endo_subset <- NormalizeData(object = endo_subset, normalization.method = "LogNormalize", scale.factor = 1e4, display.progress = TRUE)
endo_subset <- FindVariableGenes(object = endo_subset, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = TRUE)
endo_subset <- ScaleData(object = endo_subset, vars.to.regress = c("nUMI", "percent_mito"))
endo_subset <- RunPCA(object = endo_subset, pc.genes = endo_subset@var.genes, pcs.compute = 75)
endo_subset <- ProjectPCA(object = endo_subset, do.print = FALSE)
endo_subset <- FindClusters(object = endo_subset, reduction.type = "pca", dims.use = 1:10, resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
endo_subset <- RunTSNE(object = endo_subset, dims.use = 1:10)
TSNEPlot(object = endo_subset, do.label = TRUE)
beep(sound = 2)

endo_filtered <- FilterCells(endo_subset, subset.names = "Tmem119", high.thresholds = 0.1)
endo_filtered <- SubsetData(endo_filtered, ident.remove = 4, do.clean = TRUE)

# Remove endo and astrocyte unfiltered clusters from original object
all_cns_cells <- SubsetData(all_cns_cells, ident.remove = c(2, 4, 1, 8), do.clean = TRUE)
all_cns_cells <- MergeSeurat(all_cns_cells, astroytes_filtered)
all_cns_cells <- MergeSeurat(all_cns_cells, endo_filtered)

# Reanalyze dataset
all_cns_cells <- NormalizeData(object = all_cns_cells, normalization.method = "LogNormalize", scale.factor = 1e4, display.progress = TRUE)

all_cns_cells <- FindVariableGenes(object = all_cns_cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = TRUE)
length(x = all_cns_cells@var.genes)

all_cns_cells <- ScaleData(object = all_cns_cells, vars.to.regress = c("nUMI", "percent_mito"))

all_cns_cells <- RunPCA(object = all_cns_cells, pc.genes = all_cns_cells@var.genes, pcs.compute = 75)

all_cns_cells <- ProjectPCA(object = all_cns_cells, do.print = FALSE)

PCElbowPlot(object = all_cns_cells, num.pc = 75)
beep(sound = 2)

all_cns_cells <- FindClusters(object = all_cns_cells, reduction.type = "pca", dims.use = 1:40, resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

all_cns_cells <- RunTSNE(object = all_cns_cells, dims.use = 1:40)

# Filter
all_cns_cells <- SubsetData(all_cns_cells, ident.remove = 17, do.clean = TRUE)
all_cns_cells <- NormalizeData(object = all_cns_cells, normalization.method = "LogNormalize", scale.factor = 1e4, display.progress = TRUE)
all_cns_cells <- FindVariableGenes(object = all_cns_cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = TRUE)
all_cns_cells <- ScaleData(object = all_cns_cells, vars.to.regress = c("nUMI", "percent_mito"))
all_cns_cells <- RunPCA(object = all_cns_cells, pc.genes = all_cns_cells@var.genes, pcs.compute = 75)
all_cns_cells <- ProjectPCA(object = all_cns_cells, do.print = FALSE)
PCElbowPlot(object = all_cns_cells, num.pc = 75)
beep(sound = 2)

all_cns_cells <- FindClusters(object = all_cns_cells, reduction.type = "pca", dims.use = 1:40, resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
all_cns_cells <- RunTSNE(object = all_cns_cells, dims.use = 1:40)
TSNEPlot(object = all_cns_cells, do.label = TRUE, pt.size = 1)
beep(sound = 2)


# Some microglia/doublets in the neuron cluster
neuron_subset<- SubsetData(all_cns_cells, ident.use = c(7))
TSNEPlot(neuron_subset, do.label = TRUE)
neuron_filtered <- FilterCells(neuron_subset, subset.names = "Tmem119", high.thresholds = 0.1)
neuron_filtered <- FilterCells(neuron_filtered, subset.names = "P2ry12", high.thresholds = 0.1)

neuron_filtered <- SubsetData(neuron_filtered, ident.use = c(7), do.clean = TRUE)

# Remove unfiltered neurons
all_cns_cells <- SubsetData(all_cns_cells, ident.remove = c(7), do.clean = TRUE)
# merge filtered neurons
all_cns_cells <- MergeSeurat(all_cns_cells, neuron_filtered, do.normalize = FALSE)

# Final Analysis & Clustering ---------------------------------------------
all_cns_cells <- NormalizeData(object = all_cns_cells, normalization.method = "LogNormalize", scale.factor = 1e4, display.progress = TRUE)
all_cns_cells <- FindVariableGenes(object = all_cns_cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = TRUE)
length(x = all_cns_cells@var.genes)
all_cns_cells <- ScaleData(object = all_cns_cells, vars.to.regress = c("nUMI", "percent_mito"))
all_cns_cells <- RunPCA(object = all_cns_cells, pc.genes = all_cns_cells@var.genes, pcs.compute = 75)
all_cns_cells <- ProjectPCA(object = all_cns_cells, do.print = FALSE)
PCElbowPlot(object = all_cns_cells, num.pc = 75); beep(sound = 2)
all_cns_cells <- FindClusters(object = all_cns_cells, reduction.type = "pca", dims.use = 1:40, resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
all_cns_cells <- RunTSNE(object = all_cns_cells, dims.use = 1:40)

TSNEPlot(object = all_cns_cells, do.label = TRUE, pt.size = 1)

# Annotate clusters with cell type names
cluster_numbers <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
clusters_annotated <- c("Microglia", "Endothelial", "Astrocyte 1", "Oligodendrocytes", "Astrocyte 2", "Epithelial/CP", "Microglia (exAM)", "Neurons", "Pericytes", "Mono/Mac", "Ependymal Cells", "Neural Progenitors", "Fibroblasts", "OPC", "Bergmann Glia", "Erythrocytes", "OEC")

all_cns_cells@ident <- plyr::mapvalues(x = all_cns_cells@ident, from = cluster_numbers, to = clusters_annotated)

# Reorder the clusters for plotting
new_cluster_order <- c("Microglia", "Microglia (exAM)", "Mono/Mac", "Astrocyte 1", "Astrocyte 2", "Bergmann Glia", "Oligodendrocytes", "OPC",  "Endothelial", "Pericytes", "Epithelial/CP", "Ependymal Cells", "Fibroblasts", "Neurons", "Neural Progenitors", "OEC", "Erythrocytes")

all_cns_cells@ident <- factor(all_cns_cells@ident, levels = new_cluster_order)

