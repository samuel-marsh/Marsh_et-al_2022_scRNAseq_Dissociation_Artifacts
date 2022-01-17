# Microglia ---------------------------------------------------------------
micro_sub2 <- subsetLiger(object = marsh_post_liger_cleaned, clusters.use = micro_clu2, remove.missing = FALSE)
micro_sub2 <- normalize(micro_sub2)
micro_sub2 <- selectGenes(micro_sub2, do.plot = T, num.genes = 200)
micro_sub2 <- scaleNotCenter(micro_sub2)
micro_sub2 <- online_iNMF(micro_sub2, k = 25, lambda = 7, miniBatch_size = 300)
micro_sub2 <- quantile_norm(micro_sub2)
micro_sub2 <- clusterLouvainJaccard(micro_sub2,resolution = 1.5)
micro_sub2 <- runUMAP(micro_sub2, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(micro_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

umap_dim[[1]] + 
  scale_color_manual(values = c("black", "NA", "NA")) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning_rd2/micro_sub_factors_rd2.pdf")
plotFactors(liger_object = micro_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = micro_sub2, "OPALIN", plot.by = "none")
plotGene_keep_scale(object = micro_sub2, "ST18", plot.by = "none", zero.color = "lightgray")
plotGene_keep_scale(object = micro_sub2, "PTPRC", plot.by = "none")

seurat_sub2 <- ligerToSeurat(object = micro_sub2)
clu7 <- prestowrapper(object = seurat_sub2, ident.1 = 7)

micro_doubl2 <- c("12")
doublet_micro2 <- data.frame(micro_sub2@clusters, stringsAsFactors = FALSE)
names(doublet_micro2)[1] <- "cluster"
doublet_micro2 <- doublet_micro2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% micro_doubl2) %>% 
  pull(barcodes)

# Oligo ---------------------------------------------------------------
oligo_sub2 <- subsetLiger(object = marsh_post_liger_cleaned, clusters.use = oligo_clu2, remove.missing = FALSE)
oligo_sub2 <- normalize(oligo_sub2)
oligo_sub2 <- selectGenes(oligo_sub2, do.plot = T, num.genes = 400)
oligo_sub2 <- scaleNotCenter(oligo_sub2)
oligo_sub2 <- online_iNMF(oligo_sub2, k = 25, lambda = 5, miniBatch_size = 5000)
oligo_sub2 <- quantile_norm(oligo_sub2, knn_k = 15)
oligo_sub2 <- clusterLouvainJaccard(oligo_sub2,resolution = 1)
oligo_sub2 <- runUMAP(oligo_sub2)

umap_dim <- plotByDatasetAndCluster(oligo_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("plots/cleaning_rd2/oligo_sub_factors_rd2_NEW.pdf")
plotFactors(liger_object = oligo_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = oligo_sub2, "COL18A1", plot.by = "none")
plotGene_keep_scale(object = oligo_sub2, "OPALIN", plot.by = "none")
plotGene_keep_scale(object = oligo_sub2, "SAMMSON", plot.by = "none")
plotGene_keep_scale(object = oligo_sub2, "FLT1", plot.by = "none")
plotFeature(object = oligo_sub2, feature = "nGene")

sub9 <- subsetLiger(oligo_sub2, clusters.use = 9)
plotByDatasetAndCluster(sub9, return.plots = TRUE, do.legend = TRUE, text.size = 6)
plotGene_keep_scale(object = sub9, "PTPRC", plot.by = "none")

seurat_sub2 <- ligerToSeurat(object = oligo_sub2)
clu7 <- prestowrapper(object = seurat_sub2, ident.1 = 7)
beep(sound = 2)

# Excit---------------------------------------------------------------
excit_sub2 <- subsetLiger(object = marsh_post_liger_cleaned, clusters.use = excit_clu2, remove.missing = FALSE)
excit_sub2 <- normalize(excit_sub2)
excit_sub2 <- selectGenes(excit_sub2, do.plot = T, num.genes = 300)
excit_sub2 <- scaleNotCenter(excit_sub2)
excit_sub2 <- online_iNMF(excit_sub2, k = 15, lambda = 5, miniBatch_size = 2000)
excit_sub2 <- quantile_norm(excit_sub2)
excit_sub2 <- clusterLouvainJaccard(excit_sub2,resolution = 0.5)
excit_sub2 <- runUMAP(excit_sub2)

umap_dim <- plotByDatasetAndCluster(excit_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning_rd2/excit_sub_factors_rd2_NEW.pdf")
plotFactors(liger_object = excit_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = excit_sub2, "ST18", plot.by = "none")
plotGene_keep_scale(object = excit_sub2, "GFAP", plot.by = "none")
plotGene_keep_scale(object = excit_sub2, "PDGFRA", plot.by = "none")
plotFeature(object = excit_sub2, feature = "nGene")

seurat_sub2 <- ligerToSeurat(object = excit_sub2)
clu13 <- prestowrapper(seurat_sub2, ident.1 = 13)

# inhib ---------------------------------------------------------------
inhib_sub2 <- subsetLiger(object = marsh_post_liger_cleaned, clusters.use = inhib_clu2, remove.missing = FALSE)
inhib_sub2 <- normalize(inhib_sub2)
inhib_sub2 <- selectGenes(inhib_sub2, do.plot = T, num.genes = 350)
inhib_sub2 <- scaleNotCenter(inhib_sub2)
inhib_sub2 <- online_iNMF(inhib_sub2, k = 15, lambda = 5, miniBatch_size = 1000)
inhib_sub2 <- quantile_norm(inhib_sub2, knn_k = 15)
inhib_sub2 <- clusterLouvainJaccard(inhib_sub2,resolution = 2.5)
inhib_sub2 <- runUMAP(inhib_sub2, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(inhib_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 36, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("plots/cleaning_rd2/inhib_sub_factors_rd2_NEW.pdf")
plotFactors(liger_object = inhib_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = inhib_sub2, "GAD1", plot.by = "none")
plotGene_keep_scale(object = inhib_sub2, "GAD2", plot.by = "none")
plotGene_keep_scale(object = inhib_sub2, "KLK6", plot.by = "none")

inhib_doubl2 <- c("25", "8", "19")
doublet_inhib2 <- data.frame(inhib_sub2@clusters, stringsAsFactors = FALSE)
names(doublet_inhib2)[1] <- "cluster"
doublet_inhib2 <- doublet_inhib2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% inhib_doubl2) %>% 
  pull(barcodes)

# Astro ---------------------------------------------------------------
astro_sub2 <- subsetLiger(object = marsh_post_liger_cleaned, clusters.use = astro_clu2, remove.missing = FALSE)
astro_sub2 <- normalize(astro_sub2)
astro_sub2 <- selectGenes(astro_sub2, do.plot = T, num.genes = 200)
astro_sub2 <- scaleNotCenter(astro_sub2)
astro_sub2 <- online_iNMF(astro_sub2, k = 15, lambda = 5, miniBatch_size = 1000)
astro_sub2 <- quantile_norm(astro_sub2, knn_k = 15)
astro_sub2 <- clusterLouvainJaccard(astro_sub2,resolution = 0.7)
astro_sub2 <- runUMAP(astro_sub2, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(astro_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning_rd2/astro_sub_factors_rd2_NEW.pdf")
plotFactors(liger_object = astro_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = astro_sub2, "FGFR3", plot.by = "none")
plotGene_keep_scale(object = astro_sub2, "GJA1", plot.by = "none")
plotGene_keep_scale(object = astro_sub2, "ST18", plot.by = "none")

seurat_sub2 <- ligerToSeurat(astro_sub2)
clu9 <- prestowrapper(seurat_sub2, ident.1 = 9)

astro_doubl2 <- c("25", "8", "19")
doublet_astro2 <- data.frame(astro_sub2@clusters, stringsAsFactors = FALSE)
names(doublet_astro2)[1] <- "cluster"
doublet_astro2 <- doublet_astro2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_doubl2) %>% 
  pull(barcodes)

# OPC ---------------------------------------------------------------
opc_sub2 <- subsetLiger(object = marsh_post_liger_cleaned, clusters.use = opc_clu2, remove.missing = FALSE)
opc_sub2 <- normalize(opc_sub2)
opc_sub2 <- selectGenes(opc_sub2, do.plot = T, num.genes = 200)
opc_sub2 <- scaleNotCenter(opc_sub2)
opc_sub2 <- online_iNMF(opc_sub2, k = 15, lambda = 5, miniBatch_size = 400)
opc_sub2 <- quantile_norm(opc_sub2, knn_k = 15)
opc_sub2 <- clusterLouvainJaccard(opc_sub2,resolution = 0.5)
opc_sub2 <- runUMAP(opc_sub2, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(opc_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning_rd2/opc_sub_factors_rd2.pdf")
plotFactors(liger_object = opc_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = opc_sub2, "TNR", plot.by = "none")
plotGene_keep_scale(object = opc_sub2, "PTPRC", plot.by = "none")
plotGene_keep_scale(object = opc_sub2, "SYT1", plot.by = "none")

# Endo ---------------------------------------------------------------
endo_sub2 <- subsetLiger(object = marsh_post_liger_cleaned, clusters.use = endo_clu2, remove.missing = FALSE)
endo_sub2 <- normalize(endo_sub2)
endo_sub2 <- selectGenes(endo_sub2, do.plot = T, num.genes = 500)
endo_sub2 <- scaleNotCenter(endo_sub2)
endo_sub2 <- online_iNMF(endo_sub2, k = 25, lambda = 5, miniBatch_size = 50)
endo_sub2 <- quantile_norm(endo_sub2)
endo_sub2 <- clusterLouvainJaccard(endo_sub2,resolution = 0.9)
endo_sub2 <- runUMAP(endo_sub2)

umap_dim <- plotByDatasetAndCluster(endo_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning_rd2/endo_sub_factors_rd2.pdf")
plotFactors(liger_object = endo_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = endo_sub2, "FLT1", plot.by = "none")
plotGene_keep_scale(object = endo_sub2, "COL1A2", plot.by = "none")
plotGene_keep_scale(object = endo_sub2, "OPALIN", plot.by = "none")

# Elim Doublets
endo_doubl2 <- c("3")
doublet_endo2 <- data.frame(endo_sub2@clusters, stringsAsFactors = FALSE)
names(doublet_endo2)[1] <- "cluster"
doublet_endo2 <- doublet_endo2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% endo_doubl2) %>% 
  pull(barcodes)

# Fibro ---------------------------------------------------------------
fibro_sub2 <- subsetLiger(object = post_liger_cleaned, clusters.use = fibro_clu2, remove.missing = FALSE)
fibro_sub2 <- normalize(fibro_sub2)
fibro_sub2 <- selectGenes(fibro_sub2, do.plot = T, num.genes = 300)
fibro_sub2 <- scaleNotCenter(fibro_sub2)
#fibro_sub2 <- online_iNMF(fibro_sub2, k = 25, lambda = 5, miniBatch_size = 75)
fibro_sub2 <- optimizeALS(fibro_sub2, k = 5, lambda = 5, nrep = 1)
fibro_sub2 <- quantile_norm(fibro_sub2)
fibro_sub2 <- clusterLouvainJaccard(fibro_sub2,resolution = 2)
fibro_sub2 <- runUMAP(fibro_sub2)

umap_dim <- plotByDatasetAndCluster(fibro_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning_rd2/fibro_sub_factors_rd2.pdf")
plotFactors(liger_object = fibro_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = fibro_sub2, "DCN", plot.by = "none")
plotGene_keep_scale(object = fibro_sub2, "FLT1", plot.by = "none")
plotGene_keep_scale(object = fibro_sub2, "ST18", plot.by = "none")

# Elim Doublets
fibro_doubl2 <- c("0")
doublet_fibro2 <- data.frame(fibro_sub2@clusters, stringsAsFactors = FALSE)
names(doublet_fibro2)[1] <- "cluster"
doublet_fibro2 <- doublet_fibro2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% fibro_doubl2) %>% 
  pull(barcodes)

# Immune Peripheral ---------------------------------------------------------------
immune_sub2 <- subsetLiger(object = post_liger_cleaned, clusters.use = immune_clu2, remove.missing = FALSE)
immune_sub2 <- normalize(immune_sub2)
immune_sub2 <- selectGenes(immune_sub2, do.plot = T, num.genes = 200)
immune_sub2 <- scaleNotCenter(immune_sub2)
#immune_sub2 <- online_iNMF(immune_sub2, k = 25, lambda = 5, miniBatch_size = 25)
immune_sub2 <- optimizeALS(immune_sub2, k = 10, lambda = 5, nrep = 1) #same error as fibro cluster original
immune_sub2 <- quantile_norm(immune_sub2, knn_k = 10) 
immune_sub2 <- clusterLouvainJaccard(immune_sub2,resolution = 2)
immune_sub2 <- runUMAP(immune_sub2)

umap_dim <- plotByDatasetAndCluster(immune_sub2, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning_rd2/immune_sub_factors_rd2.pdf")
plotFactors(liger_object = immune_sub2, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = immune_sub2, "PTPRC", plot.by = "none")
plotGene_keep_scale(object = immune_sub2, "LYZ", plot.by = "none")
plotGene_keep_scale(object = immune_sub2, "CD3E", plot.by = "none")
plotGene_keep_scale(object = immune_sub2, "CX3CR1", plot.by = "none")

# Elim Doublets
immune_doubl2 <- c("0")
doublet_immune2 <- data.frame(immune_sub2@clusters, stringsAsFactors = FALSE)
names(doublet_immune2)[1] <- "cluster"
doublet_immune2 <- doublet_immune2 %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% immune_doubl2) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
doublet_rd2_new <- data.frame(marsh_post_liger_cleaned@clusters, stringsAsFactors = FALSE)
names(doublet_rd2_new)[1] <- "cluster"
doublet_rd2_new <- doublet_rd2_new %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% doublet_clu2) %>% 
  pull(barcodes)

# Compile doublets into single vector
doublets_complete2 <- c(doublet_astro2, doublet_endo2, doublet_inhib2, doublet_rd2_new, doublet_fibro2, doublet_immune2)

# Barcodes already in cell.data slot from previous round

# Remove doublet cells from Round02 object and save as cleaned2 object
marsh_post_liger_cleaned2 <- subsetLiger(object = marsh_post_liger_cleaned, cells.use = setdiff(unique(marsh_post_liger_cleaned@cell.data$barcodes), doublets_complete2), remove.missing = FALSE)
