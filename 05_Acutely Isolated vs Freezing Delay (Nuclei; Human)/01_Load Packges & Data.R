# Load Packages -----------------------------------------------------------
library(Seurat) #v2.3.4
library(liger)
library(NNLM)
library(stringr)
library(presto)
library(gridExtra)

# Load Data ---------------------------------------------------------------
work.dir <- '/home/tkamath/smarsh_mgs/SMarsh_MG/'
lib.list <- c('pEXP19sHSrCTXt6iSAMPLE03d20190601dapi','pEXP20sHSrCTXt6iSAMPLE01d20190601O1',
              'pEXP20sHSrCTXt0iSAMPLE01d20190601L1','pEXP20sHSrCTXt0iSAMPLE03d20190601N1')
lib.names <- c('T6C','T6A','T0A','T0C')
dge.full <- lapply(c(1:length(lib.list)), function(x){
  dge.1 <- readMM(paste0(work.dir,lib.list[x],'/filtered_feature_bc_matrix/matrix.mtx.gz') )
  bx2 <- read.table(paste0(work.dir,lib.list[x],'/filtered_feature_bc_matrix/barcodes.tsv.gz'))
  genes2 <- read.table(paste0(work.dir,lib.list[x],'/filtered_feature_bc_matrix/features.tsv.gz'))
  rownames(dge.1) <- genes2$V2
  bx.use <- as.vector(unlist(lapply(bx2$V1,function(y){paste0(lib.names[x],"_",y)})))
  colnames(dge.1) <- bx.use
  gc()
  print(paste0('Returning ',as.character(x)))
  return(dge.1)
})
names(dge.full) <- lib.names

# Liger Round 01 ----------------------------------------------------------
a.fresh.use = createLiger(dge.full)
a.fresh.use = liger::normalize(a.fresh.use)
a.fresh.use = selectGenes(a.fresh.use, do.plot = T,num.genes = 900)
a.fresh.use = scaleNotCenter(a.fresh.use)
a.fresh.use = online_iNMF(a.fresh.use,k = 40,lambda = 10,miniBatch_size = 20000)
a.fresh.use = quantile_norm(a.fresh.use)
a.fresh.use = clusterLouvainJaccard(a.fresh.use,resolution = 0.5)
a.fresh.use = runUMAP(a.fresh.use)

# plotting
pdf('factors_new.pdf')
plotFactors(a.fresh.use, num.genes = 8, plot.tsne = T)
dev.off()

marker.genes <- c('OLIG1','CX3CR1','AQP4','TNR','THEMIS','SLC17A7','GAD2','CLDN5')
pdf('markergenes.pdf')
lapply(marker.genes,function(x){print(plotGene(a.fresh.use, gene = x, by.dataset = F))})
dev.off()

p1 <- plotByDatasetAndCluster(a.fresh.use,return.plots = T,do.legend = T)
pdf('tsne.pdf')
p1[[1]]
p1[[2]]
dev.off()

# Convert to Seurat and Annotate
a.fresh.seurat <- ligerToSeurat(a.fresh.use)
fresh.markers <- prestowrapper(a.fresh.seurat, log.FC = log(1.2),all.clusters = T)

# 0 - Oligodendrocyte
# 1 - Oligodendrocyte
# 7 - Oligodendrocyte
# 11 - Oligo?
# 9 - OPC
# 14 - OPC
# 4 - TP73+/tumor
# 12 - tumor+?
# 3 - Microglia
# 5 - Microglia
# 6 - Microglia
# 18 - Microglia
# 8 - Ex Neuron
# 15 - Ex Neuron
# 13 - Interneuron
# 21 - Interneuron
# 10 - PBMCs
# 2 - Tumor TP73?/ astrocyte
# 16 - Astrocyte
# 17 - REMOVE!
# 19 - Endothelialcells
# 23- Fibroblast
# 20 - Dividing (?)
# 22 - NEURON?
# 24 - MOBP
# 25 - MOBP

# Clean object and add annotation
a.fresh.clean <- subsetLiger(a.fresh.use,clusters.use = setdiff(levels(a.fresh.use@clusters),c(17)),
                             remove.missing = F)
levels(a.fresh.clean@clusters) <- c('Oligodendrocyte','Oligodendrocyte','Astrocyte','Microglia',
                                    'Tumor_TP73','Microglia','Microglia','Oligodendrocyte',
                                    'Ex_Neuron','OPC','PBMC','Oligodendrocyte',
                                    'Tumor_TP73','Interneuron','OPC','Ex_Neuron',
                                    'Astrocyte','Microglia','Endothelial',
                                    'Dividing_cells','Interneuron','Interneuron', 'Fibroblast',
                                    'Oligodendrocyte_MOBP','Oligodendrocyte_MOBP')

# Add mito and subset
a.fresh.clean <- AddMito(a.fresh.clean,species = 'human')
plotFeature(a.fresh.clean,feature = 'percent.mito',by.dataset = F)

idx.use <- rownames(a.fresh.clean@cell.data)[which(a.fresh.clean@cell.data$percent.mito <= 0.1)]
a.fresh.clean2 <- subsetLiger(a.fresh.clean, cells.use = idx.use,remove.missing = F)

plotByDatasetAndCluster(a.fresh.clean2)

# Add additional meta data
a.fresh.clean@cell.data$individual <- as.vector(unlist(lapply(a.fresh.clean@cell.data$dataset,function(x){
  if(grepl('C',x))
    return('C')
  else
    return('A')
})))

################################### 
##### Subcluster MGs ##############
###################################
a.mg <- subsetLiger(a.fresh.use, clusters.use = c(3,5,6,18))
a.mg = liger::normalize(a.mg)
a.mg = selectGenes(a.mg, do.plot = T,num.genes = 600)
a.mg = scaleNotCenter(a.mg)
a.mg <- optimizeALS(a.mg,k = 20,lambda = 5)
a.mg = quantile_norm(a.mg)
a.mg = clusterLouvainJaccard(a.mg,resolution = 0.5)
a.mg = runUMAP(a.mg)

# plotting
p1 <- plotByDatasetAndCluster(a.mg,return.plots = T,do.legend = T)
pdf('tsne_mg.pdf')
p1[[1]]
p1[[2]]
dev.off()

# Clean Mg object
# Remove clusters: 5,8,9
a.mgcleaned <- subsetLiger(a.mg, clusters.use = setdiff(unique(a.mg@clusters),c(5,8,9) ))
a.mgcleaned = quantile_norm(a.mgcleaned)
a.mgcleaned = clusterLouvainJaccard(a.mgcleaned,resolution = 0.6)
a.mgcleaned = runUMAP(a.mgcleaned)
p1 <- plotByDatasetAndCluster(a.mgcleaned,return.plots = T,do.legend = T)
pdf('tsne_mgcleaned.pdf')
p1[[1]]
p1[[2]]
dev.off()

# plotting factors
pdf('/home/tkamath/smarsh_mgs/MG/factorplot_mgcleaned.pdf')
plotFactors(a.mgcleaned, plot.tsne = T)
dev.off()

# Add mito %
a.mgcleaned <- AddMito(a.mgcleaned, species = 'human')

# Save final object
saveRDS(a.mgcleaned,'MG/a.mgcleaned.rds')

################################
########Astrocytes##############
################################
a.astro <- subsetLiger(a.fresh.use, clusters.use = c(2,16))
a.astro = liger::normalize(a.astro)
a.astro = selectGenes(a.astro, do.plot = T,num.genes = 700)
a.astro = scaleNotCenter(a.astro)
a.astro <- optimizeALS(a.astro,k = 20,lambda = 5)
a.astro = quantile_norm(a.astro)
a.astro = clusterLouvainJaccard(a.astro,resolution = 0.4)
a.astro = runUMAP(a.astro)

# plotting
p1 <- plotByDatasetAndCluster(a.astro,return.plots = T,do.legend = T)
pdf('tsne_astro.pdf')
p1[[1]]
p1[[2]]
dev.off()

pdf('factor_astro.pdf')
plotFactors(a.astro, plot.tsne = T)
dev.off()

# Clean astros
# Remove 3,6,..
a.astrocleaned <- subsetLiger(a.astro, clusters.use = setdiff(levels(a.astro@clusters),c(3,6)))
a.astrocleaned <- optimizeALS(a.astrocleaned,k = 15,lambda = 10)
a.astrocleaned = quantile_norm(a.astrocleaned)
a.astrocleaned = clusterLouvainJaccard(a.astrocleaned,resolution = 0.3)
a.astrocleaned = runUMAP(a.astrocleaned)
p1 <- plotByDatasetAndCluster(a.astrocleaned,return.plots = T,do.legend = T)
pdf('tsne_astro_cleaned.pdf')
p1[[1]]
p1[[2]]
dev.off()

# Add mito % and subset
a.astrocleaned <- AddMito(a.astrocleaned,species = 'human')
plotFeature(a.astrocleaned,feature = 'percent.mito',by.dataset = F)

idx.use <- rownames(a.astrocleaned@cell.data)[which(a.astrocleaned@cell.data$percent.mito <= 0.1)]
a.astrocleaned <- subsetLiger(a.astrocleaned,cells.use = idx.use)

# save final object
saveRDS(a.astrocleaned,'a.astrocleaned.rds')