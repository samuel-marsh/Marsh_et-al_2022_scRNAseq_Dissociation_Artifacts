# Doublet cleaning 

# Microglia ---------------------------------------------------------------
micro_sub <- subsetLiger(object = marsh_post_liger, clusters.use = micro_clu, remove.missing = FALSE)
micro_sub <- normalize(micro_sub)
micro_sub <- selectGenes(micro_sub, do.plot = T, num.genes = 400)
micro_sub <- scaleNotCenter(micro_sub)
micro_sub <- online_iNMF(micro_sub, k = 25, lambda = 7, miniBatch_size = 200)
micro_sub <- quantile_norm(micro_sub)
micro_sub <- clusterLouvainJaccard(micro_sub,resolution = 2.5)
micro_sub <- runUMAP(micro_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(micro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("plots/cleaning/micro_sub_factors_NEW.pdf")
plotFactors(liger_object = micro_sub, num.genes = 8, plot.tsne = T)
dev.off()

# Remove mito clustering factor to sep neurons better
micro_sub <- quantile_norm(micro_sub, dims.use = setdiff(1:25, 3))
micro_sub <- clusterLouvainJaccard(micro_sub,resolution = 2.5)
micro_sub <- runUMAP(micro_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(micro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("plots/cleaning/micro_sub_factors_NEW2.pdf")
plotFactors(liger_object = micro_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = micro_sub, "CX3CR1", plot.by = "none")
plotGene_keep_scale(object = micro_sub, "PTPRC", plot.by = "none", zero.color = "lightgray")
plotGene_keep_scale(object = micro_sub, "SYT1", plot.by = "none")

seurat_sub <- ligerToSeurat(object = micro_sub)
clu7 <- prestowrapper(object = seurat_sub, ident.1 = 7)

# Elim Doublets
# ADD NOTES
micro_doubl <- c("8", "14", "17", "1", "21", "19", "15")
doublet_micro <- data.frame(micro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_micro)[1] <- "cluster"
doublet_micro <- doublet_micro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% micro_doubl) %>% 
  pull(barcodes)

# Oligo ---------------------------------------------------------------
oligo_sub <- subsetLiger(object = marsh_post_liger, clusters.use = oligo_clu, remove.missing = FALSE)
oligo_sub <- normalize(oligo_sub)
oligo_sub <- selectGenes(oligo_sub, do.plot = T, num.genes = 400)
oligo_sub <- scaleNotCenter(oligo_sub)
oligo_sub <- online_iNMF(oligo_sub, k = 25, lambda = 5, miniBatch_size = 5000)
oligo_sub <- quantile_norm(oligo_sub, knn_k = 15)
oligo_sub <- clusterLouvainJaccard(oligo_sub,resolution = 0.5)
oligo_sub <- runUMAP(oligo_sub)

umap_dim <- plotByDatasetAndCluster(oligo_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("plots/cleaning/oligo_sub_factors.pdf")
plotFactors(liger_object = oligo_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = oligo_sub, "OPALIN", plot.by = "none")
plotGene_keep_scale(object = oligo_sub, "COL18A1", plot.by = "none")
plotGene_keep_scale(object = oligo_sub, "SYT1", plot.by = "none")
plotFeature(object = oligo_sub, feature = "nGene")

sub17 <- subsetLiger(oligo_sub, clusters.use = 17)
plotByDatasetAndCluster(sub17, return.plots = TRUE, do.legend = TRUE, text.size = 6)

seurat_sub <- ligerToSeurat(object = oligo_sub)
clu8 <- prestowrapper(object = seurat_sub, ident.1 = 8, ident.2 = 0)
clu12 <- prestowrapper(object = seurat_sub, ident.1 = 12)
sub12 <- subsetLiger(oligo_sub, clusters.use = 12)
plotByDatasetAndCluster(sub12, return.plots = TRUE, do.legend = TRUE, text.size = 6)

clu14 <- prestowrapper(object = seurat_sub, ident.1 = 14)

clu18 <- prestowrapper(object = seurat_sub, ident.1 = 18)
beep(sound = 2)

# Elim Doublets
oligo_doubl <- c("7", "18", "11", "17")
doublet_oligo <- data.frame(oligo_sub@clusters, stringsAsFactors = FALSE)
names(doublet_oligo)[1] <- "cluster"
doublet_oligo <- doublet_oligo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% oligo_doubl) %>% 
  pull(barcodes)

# Excit---------------------------------------------------------------
excit_sub <- subsetLiger(object = marsh_post_liger, clusters.use = excit_clu, remove.missing = FALSE)
excit_sub <- normalize(excit_sub)
excit_sub <- selectGenes(excit_sub, do.plot = T, num.genes = 400)
excit_sub <- scaleNotCenter(excit_sub)
excit_sub <- online_iNMF(excit_sub, k = 25, lambda = 5, miniBatch_size = 2000)
excit_sub <- quantile_norm(excit_sub)
excit_sub <- clusterLouvainJaccard(excit_sub,resolution = 1)
excit_sub <- runUMAP(excit_sub)

umap_dim <- plotByDatasetAndCluster(excit_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "polychrome")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("plots/cleaning/excit_sub_factors.pdf")
plotFactors(liger_object = excit_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = excit_sub, "SLC17A7", plot.by = "none")
plotGene_keep_scale(object = excit_sub, "RBFOX3", plot.by = "none")
plotGene_keep_scale(object = excit_sub, "FEZF2", plot.by = "none")
plotFeature(object = excit_sub, feature = "nGene")

seurat_sub <- ligerToSeurat(object = excit_sub)
clu13 <- prestowrapper(seurat_sub, ident.1 = 13)
clu12 <- prestowrapper(seurat_sub, ident.1 = 12)

# Elim Doublets
excit_doubl <- c("11", "18", "3", "16", "19")
doublet_excit <- data.frame(excit_sub@clusters, stringsAsFactors = FALSE)
names(doublet_excit)[1] <- "cluster"
doublet_excit <- doublet_excit %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% excit_doubl) %>% 
  pull(barcodes)

# inhib ---------------------------------------------------------------
inhib_sub <- subsetLiger(object = marsh_post_liger, clusters.use = inhib_clu, remove.missing = FALSE)
inhib_sub <- normalize(inhib_sub)
inhib_sub <- selectGenes(inhib_sub, do.plot = T, num.genes = 500)
inhib_sub <- scaleNotCenter(inhib_sub)
inhib_sub <- online_iNMF(inhib_sub, k = 25, lambda = 5, miniBatch_size = 1000)
inhib_sub <- quantile_norm(inhib_sub, knn_k = 15)
inhib_sub <- clusterLouvainJaccard(inhib_sub,resolution = 0.8)
inhib_sub <- runUMAP(inhib_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(inhib_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

beep(sound = 2)

pdf("plots/cleaning/inhib_sub_factors.pdf")
plotFactors(liger_object = inhib_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = inhib_sub, "GAD1", plot.by = "none")
plotGene_keep_scale(object = inhib_sub, "HOMER3", plot.by = "none")
plotGene_keep_scale(object = inhib_sub, "SYT1", plot.by = "none")

# Elim Doublets
inhib_doubl <- c("15")
doublet_inhib <- data.frame(inhib_sub@clusters, stringsAsFactors = FALSE)
names(doublet_inhib)[1] <- "cluster"
doublet_inhib <- doublet_inhib %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% inhib_doubl) %>% 
  pull(barcodes)

# Astro ---------------------------------------------------------------
astro_sub <- subsetLiger(object = marsh_post_liger, clusters.use = astro_clu, remove.missing = FALSE)
astro_sub <- normalize(astro_sub)
astro_sub <- selectGenes(astro_sub, do.plot = T, num.genes = 400)
astro_sub <- scaleNotCenter(astro_sub)
astro_sub <- online_iNMF(astro_sub, k = 25, lambda = 5, miniBatch_size = 1000)
astro_sub <- quantile_norm(astro_sub, knn_k = 15)
astro_sub <- clusterLouvainJaccard(astro_sub,resolution = 0.7)
astro_sub <- runUMAP(astro_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(astro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning/astro_sub_factors.pdf")
plotFactors(liger_object = astro_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = astro_sub, "GFAP", plot.by = "none")
plotGene_keep_scale(object = astro_sub, "VIM", plot.by = "none")
plotGene_keep_scale(object = astro_sub, "SYT1", plot.by = "none")
plotFeature(astro_sub, feature = "nGene", by.dataset = F)

seurat_sub <- ligerToSeurat(astro_sub)
clu9 <- prestowrapper(seurat_sub, ident.1 = 9)

# Elim Doublets
astro_doubl <- c("6", "10")
doublet_astro <- data.frame(astro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_astro)[1] <- "cluster"
doublet_astro <- doublet_astro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% astro_doubl) %>% 
  pull(barcodes)

# OPC ---------------------------------------------------------------
opc_sub <- subsetLiger(object = marsh_post_liger, clusters.use = opc_clu, remove.missing = FALSE)
opc_sub <- normalize(opc_sub)
opc_sub <- selectGenes(opc_sub, do.plot = T, num.genes = 500)
opc_sub <- scaleNotCenter(opc_sub)
opc_sub <- online_iNMF(opc_sub, k = 25, lambda = 5, miniBatch_size = 400)
opc_sub <- quantile_norm(opc_sub, knn_k = 15)
opc_sub <- clusterLouvainJaccard(opc_sub,resolution = 0.5)
opc_sub <- runUMAP(opc_sub, min_dist = 0.3, n_neighbors = 15)

umap_dim <- plotByDatasetAndCluster(opc_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning/opc_sub_factors.pdf")
plotFactors(liger_object = opc_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = opc_sub, "TNR", plot.by = "none")
plotGene_keep_scale(object = opc_sub, "PTPRC", plot.by = "none")
plotGene_keep_scale(object = opc_sub, "OPALIN", plot.by = "none")

# Elim Doublets
opc_doubl <- c("3", "5", "4")
doublet_opc <- data.frame(opc_sub@clusters, stringsAsFactors = FALSE)
names(doublet_opc)[1] <- "cluster"
doublet_opc <- doublet_opc %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% opc_doubl) %>% 
  pull(barcodes)

# Endo ---------------------------------------------------------------
endo_sub <- subsetLiger(object = marsh_post_liger, clusters.use = endo_clu, remove.missing = FALSE)
endo_sub <- normalize(endo_sub)
endo_sub <- selectGenes(endo_sub, do.plot = T, num.genes = 500)
endo_sub <- scaleNotCenter(endo_sub)
endo_sub <- online_iNMF(endo_sub, k = 25, lambda = 5, miniBatch_size = 50)
endo_sub <- quantile_norm(endo_sub)
endo_sub <- clusterLouvainJaccard(endo_sub,resolution = 0.9)
endo_sub <- runUMAP(endo_sub)

umap_dim <- plotByDatasetAndCluster(endo_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning/endo_sub_factors.pdf")
plotFactors(liger_object = endo_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = endo_sub, "FLT1", plot.by = "none")
plotGene_keep_scale(object = endo_sub, "DCN", plot.by = "none")
plotGene_keep_scale(object = endo_sub, "CLDN5", plot.by = "none")

# Elim Doublets
endo_doubl <- c("2", "5", "4")
doublet_endo <- data.frame(endo_sub@clusters, stringsAsFactors = FALSE)
names(doublet_endo)[1] <- "cluster"
doublet_endo <- doublet_endo %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% endo_doubl) %>% 
  pull(barcodes)

# Fibro ---------------------------------------------------------------
fibro_sub <- subsetLiger(object = marsh_post_liger, clusters.use = fibro_clu, remove.missing = FALSE)
fibro_sub <- normalize(fibro_sub)
fibro_sub <- selectGenes(fibro_sub, do.plot = T, num.genes = 300)
fibro_sub <- scaleNotCenter(fibro_sub)
fibro_sub <- online_iNMF(fibro_sub, k = 25, lambda = 5, miniBatch_size = 75)
fibro_sub <- quantile_norm(fibro_sub)
fibro_sub <- clusterLouvainJaccard(fibro_sub,resolution = 0.9)
fibro_sub <- runUMAP(fibro_sub)

umap_dim <- plotByDatasetAndCluster(fibro_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning/fibro_sub_factors.pdf")
plotFactors(liger_object = fibro_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = fibro_sub, "COL1A2", plot.by = "none")
plotGene_keep_scale(object = fibro_sub, "CLDN5", plot.by = "none")
plotGene_keep_scale(object = fibro_sub, "PLP1", plot.by = "none")

# Elim Doublets
fibro_doubl <- c("1", "3", "4", "6", "0", "2")
doublet_fibro <- data.frame(fibro_sub@clusters, stringsAsFactors = FALSE)
names(doublet_fibro)[1] <- "cluster"
doublet_fibro <- doublet_fibro %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% fibro_doubl) %>% 
  pull(barcodes)

# Immune Peripheral ---------------------------------------------------------------
immune_sub <- subsetLiger(object = marsh_post_liger, clusters.use = immune_clu, remove.missing = FALSE)
immune_sub <- normalize(immune_sub)
immune_sub <- selectGenes(immune_sub, do.plot = T, num.genes = 400)
immune_sub <- scaleNotCenter(immune_sub)
immune_sub <- online_iNMF(immune_sub, k = 25, lambda = 5, miniBatch_size = 25)
immune_sub <- quantile_norm(immune_sub)
immune_sub <- clusterLouvainJaccard(immune_sub,resolution = 1)
immune_sub <- runUMAP(immune_sub)

umap_dim <- plotByDatasetAndCluster(immune_sub, return.plots = TRUE, do.legend = TRUE, text.size = 6)
umap_palette <- DiscretePalette(n = 26, palette = "alphabet")
umap_dim <- umap_dim[[2]] + 
  scale_color_manual(values = umap_palette) +
  theme_classic() + 
  theme(legend.position = "right") +
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
umap_dim

pdf("plots/cleaning/immune_sub_factors.pdf")
plotFactors(liger_object = immune_sub, num.genes = 8, plot.tsne = T)
dev.off()

plotGene_keep_scale(object = immune_sub, "PTPRC", plot.by = "none")
plotGene_keep_scale(object = immune_sub, "LYZ", plot.by = "none")
plotGene_keep_scale(object = immune_sub, "CD3E", plot.by = "none")
plotGene_keep_scale(object = immune_sub, "PLP1", plot.by = "none")

# Elim Doublets
immune_doubl <- c("2", "3", "0")
doublet_immune <- data.frame(immune_sub@clusters, stringsAsFactors = FALSE)
names(doublet_immune)[1] <- "cluster"
doublet_immune <- doublet_immune %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% immune_doubl) %>% 
  pull(barcodes)

# Unknown ---------------------------------------------------------------
# Most likely doublets
# Elim Doublets
unknown_doubl <- unknown_clu
doublet_unknown <- data.frame(marsh_post_liger@clusters, stringsAsFactors = FALSE)
names(doublet_unknown)[1] <- "cluster"
doublet_unknown <- doublet_unknown %>% 
  rownames_to_column(var = "barcodes") %>% 
  filter(cluster %in% unknown_doubl) %>% 
  pull(barcodes)

# Compile Doublets and Clean ----------------------------------------------
# Compile doublet list from all clusters
doublets_complete <- c(doublet_astro, doublet_endo, doublet_excit, doublet_fibro, doublet_immune, doublet_inhib, doublet_micro, doublet_oligo, doublet_opc, doublet_unknown)

# Add barcodes to cell.data slot
marsh_post_liger@cell.data$barcodes <- rownames(marsh_post_liger@cell.data)

# Remove doublet cells from original liger object and save as cleaned object
marsh_post_liger_cleaned <- subsetLiger(object = marsh_post_liger, cells.use = setdiff(unique(marsh_post_liger@cell.data$barcodes), doublets_complete), remove.missing = FALSE)

write_rds(marsh_post_liger_cleaned, "RDS_Objects/marsh_post_liger_round01_cleaned.RDS")