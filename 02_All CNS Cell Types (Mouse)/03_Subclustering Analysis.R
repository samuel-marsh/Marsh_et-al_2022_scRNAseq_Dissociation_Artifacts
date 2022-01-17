# Subclustering Analysis was performed using following packages:
    # R 3.6.1
    # Seurat 3.1.5 

# Analysis was performed using saved .RDS Seurat Object as starting point.
    # Original object was created as described in preceding code using Seurat 2.3.4
    # Object was read into R and updated to Seurat V3 using `UpdateSeuratObject()` before conducting subclustering analysis below.`


# Subset Decisions --------------------------------------------------------
# Subset cell types with >= 80 cells per replicate 
# Dimplot of original object
all_cns_cells <- UpdateSeuratObject(all_cns_cells)

DimPlot(all_cns_cells, reduction = "tsne", label = TRUE)

all_stats <- Cluster_Stats_All_Samples(all_cns_cells)

# Create Subset Objects
myeloid <- subset(x = all_cns_cells, idents = c("Microglia", "Microglia (exAM)", "Mono/Mac"))
DimPlot(myeloid)

endo <- subset(x = all_cns_cells, idents = c("Endothelial", "Pericytes"))
DimPlot(endo)

astro <- subset(x = all_cns_cells, idents = c("Astrocyte 1", "Astrocyte 2", "Bergmann Glia"))
DimPlot(astro)

epi <- subset(x = all_cns_cells, idents = c("Epithelial/CP"))
DimPlot(epi)

oligo <- subset(x = all_cns_cells, idents = c("Oligodendrocytes", "OPC"))
DimPlot(oligo)

neuron <- subset(x = all_cns_cells, idents = c("Neurons", "Neural Progenitors"))
DimPlot(neuron)

# Myeloid Subclustering ---------------------------------------------------
#  Normalize Data 
myeloid <- NormalizeData(myeloid, normalization.method = "LogNormalize", scale.factor = 10000)
myeloid <- FindVariableFeatures(myeloid, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
# Scale the data
all_genes <- rownames(myeloid)
myeloid <- ScaleData(myeloid, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))
myeloid <- RunPCA(myeloid, features = VariableFeatures(object = myeloid), npcs = 50)
myeloid <- JackStraw(myeloid)
myeloid <- ScoreJackStraw(myeloid, dims = 1:20)
beep(sound = 2)
# Clustering
myeloid <- FindNeighbors(myeloid, dims = 1:15)
myeloid <- FindClusters(myeloid, resolution = 0.2)
myeloid <- RunTSNE(myeloid, dims = 1:15)

DimPlot(myeloid, label = TRUE, reduction = "tsne")
p_myeloid <- DimPlot(myeloid, label = FALSE, reduction = "tsne", group.by = "Treatment", cols = c("navy", "orange"), pt.size = 4, shuffle = TRUE)
p_myeloid

# Astro Subclustering ---------------------------------------------------
# Normalize Data 
astro <- NormalizeData(astro, normalization.method = "LogNormalize", scale.factor = 10000)
astro <- FindVariableFeatures(astro, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.75, Inf))
# Scale the data
all_genes <- rownames(astro)
astro <- ScaleData(astro, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))
astro <- RunPCA(astro, features = VariableFeatures(object = astro), npcs = 50)
astro <- JackStraw(astro)
astro <- ScoreJackStraw(astro, dims = 1:20)
beep(sound = 2)
# Clustering
astro <- FindNeighbors(astro, dims = 1:7)
astro <- FindClusters(astro, resolution = 0.2)
astro <- RunTSNE(astro, dims = 1:7)

DimPlot(astro, label = TRUE, reduction = "tsne")
p_astro <- DimPlot(astro, label = FALSE, reduction = "tsne", group.by = "Treatment", cols = c("navy", "orange"), pt.size = 4, shuffle = TRUE)
p_astro

# Endo Subclustering ---------------------------------------------------
# Normalize Data 
endo <- NormalizeData(endo, normalization.method = "LogNormalize", scale.factor = 10000)
endo <- FindVariableFeatures(endo, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.75, Inf))
# Scale the data
all_genes <- rownames(endo)
endo <- ScaleData(endo, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))
endo <- RunPCA(endo, features = VariableFeatures(object = endo), npcs = 50)
endo <- JackStraw(endo)
endo <- ScoreJackStraw(endo, dims = 1:20)
beep(sound = 2)
# Clustering
endo <- FindNeighbors(endo, dims = 1:8)
endo <- FindClusters(endo, resolution = 0.2)
endo <- RunTSNE(endo, dims = 1:8)

DimPlot(endo, label = TRUE, reduction = "tsne")
p_endo <- DimPlot(endo, label = FALSE, reduction = "tsne", group.by = "Treatment", cols = c("navy", "orange"), pt.size = 4, shuffle = TRUE)
p_endo

# Oligo Subclustering ---------------------------------------------------
# Normalize Data 
oligo <- NormalizeData(oligo, normalization.method = "LogNormalize", scale.factor = 10000)
oligo <- FindVariableFeatures(oligo, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))
# Scale the data
all_genes <- rownames(oligo)
oligo <- ScaleData(oligo, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))
oligo <- RunPCA(oligo, features = VariableFeatures(object = oligo), npcs = 50)
oligo <- JackStraw(oligo)
oligo <- ScoreJackStraw(oligo, dims = 1:20)
beep(sound = 2)
# Clustering
oligo <- FindNeighbors(oligo, dims = 1:13)
oligo <- FindClusters(oligo, resolution = 0.6)
oligo <- RunTSNE(oligo, dims = 1:13)

DimPlot(oligo, label = TRUE, reduction = "tsne")
p_oligo <- DimPlot(oligo, label = FALSE, reduction = "tsne", group.by = "Treatment", cols = c("navy", "orange"), pt.size = 4, shuffle = TRUE)
p_oligo

# Epi Subclustering ---------------------------------------------------
# Normalize Data 
epi <- NormalizeData(epi, normalization.method = "LogNormalize", scale.factor = 10000)
epi <- FindVariableFeatures(epi, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))
# Scale the data
all_genes <- rownames(epi)
epi <- ScaleData(epi, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))
epi <- RunPCA(epi, features = VariableFeatures(object = epi), npcs = 50)
epi <- JackStraw(epi)
epi <- ScoreJackStraw(epi, dims = 1:20)
beep(sound = 2)
# Clustering
epi <- FindNeighbors(epi, dims = 1:3)
epi <- FindClusters(epi, resolution = 0.4)
epi <- RunTSNE(epi, dims = 1:3)

DimPlot(epi, label = TRUE, reduction = "tsne")
p_epi <- DimPlot(epi, label = FALSE, reduction = "tsne", group.by = "Treatment", cols = c("navy", "orange"), pt.size = 4, shuffle = TRUE)
p_epi

# Neuron Subclustering ---------------------------------------------------
# Normalize Data 
neuron <- NormalizeData(neuron, normalization.method = "LogNormalize", scale.factor = 10000)
neuron <- FindVariableFeatures(neuron, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(1, Inf))
# Scale the data
all_genes <- rownames(neuron)
neuron <- ScaleData(neuron, features = all_genes, vars.to.regress = c("nCount_RNA", "percent_mito"))
neuron <- RunPCA(neuron, features = VariableFeatures(object = neuron), npcs = 50)
neuron <- JackStraw(neuron)
neuron <- ScoreJackStraw(neuron, dims = 1:20)
beep(sound = 2)
# Clustering
neuron <- FindNeighbors(neuron, dims = 1:3)
neuron <- FindClusters(neuron, resolution = 0.4)
neuron <- RunTSNE(neuron, dims = 1:3)

DimPlot(neuron, label = TRUE, reduction = "tsne")
p_neuron <- DimPlot(neuron, label = FALSE, reduction = "tsne", group.by = "Treatment", cols = c("navy", "orange"), pt.size = 4, shuffle = TRUE)
p_neuron

# Finalized ---------------------------------------------------------------
# Plot by Cell Type and Treatment
wrap_plots(p_myeloid, p_astro, p_oligo, p_neuron, p_endo, p_epi)
