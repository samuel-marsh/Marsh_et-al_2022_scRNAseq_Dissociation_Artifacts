# 1.0 Load Packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(marsh.utils)
library(viridis)
library(liger)
library(beepr)
library(scCustomize)

# 2.0 Load Liger Subsetted Objects --------------------------------------------
micro_merged_liger_cleaned <- read_rds("RDS_subset_merge_liger/micro_subset_cleaned.RDS")
astro_merged_liger_cleaned <- read_rds("RDS_subset_merge_liger/astro_subset_cleaned3_C.RDS")
oligo_merged_liger_cleaned <- read_rds("RDS_subset_merge_liger/oligo_merge_cleaned_b.RDS")
opc_merged_liger_cleaned2 <- read_rds("RDS_subset_merge_liger/opc_subset_cleaned2.RDS")
excit_merged_liger_cleaned2 <- read_rds("RDS_subset_merge_liger/excit_subset_cleaned2.RDS")
inhib_merged_liger_cleaned3 <- read_rds("RDS_subset_merge_liger/inhib_subset_cleaned3.RDS")

# Factor Numbers by cell type
microglia - 6
astrocyte - 15
Oligo - 18
OPC - 8
Excit - 11
Inhib - 18

# 3.0 Micro Factor ---------------------------------------------------------
# plot factor
micro_plots2 <- plotGeneLoadings(micro_merged_liger_cleaned, dataset1 = "Marsh", dataset2 = "Zhou", return.plots = T, pt.size = 0.5, num.genes.show = 25)
micro_plots2[[6]]

# pull factor and overall gene list
microW <- t(micro_merged_liger_cleaned@W)
microW <- data.frame(microW)
microW <- microW %>% 
  rownames_to_column(var = "gene")

microW <- microW %>% 
  select(gene, X6)
microW <- microW %>% 
  arrange(desc(X6))
genes <- microW$gene
microW$gene <- factor(microW$gene, levels = genes)

# Plot factor and select threshold cutoff
ggplot(data = microW, aes(x = gene, y = X6)) +
  geom_point() +
  geom_hline(yintercept = 0.67)

# Extract gene list
micro_genes_above_cutoff <- sum(microW$X6 > 0.67)

micro_factor_cutoff <- data.frame(top_genes_by_factor(liger_object = micro_merged_liger_cleaned, liger_factor = 6, num_genes = micro_genes_above_cutoff))

colnames(micro_factor_cutoff) <- "micro_factor_cutoff"

# save factor list
write_rds(micro_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_micro_cutoff_gene_list.RDS")
write.csv(micro_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_micro_cutoff_gene_list.csv")

#4.0 Astrpcyte Factor ---------------------------------------------------------
# Plot factor
astro_plots2 <- plotGeneLoadings(astro_merged_liger_cleaned, dataset1 = "Marsh", dataset2 = "Zhou", return.plots = T, pt.size = 0.5, num.genes.show = 25)
astro_plots2[[15]]

# pull factor and genes
astroW <- t(astro_merged_liger_cleaned@W)
astroW <- data.frame(astroW)
astroW <- astroW %>% 
  rownames_to_column(var = "gene")

astroW <- astroW %>% 
  select(gene, X15)

astroW <- astroW %>% 
  arrange(desc(X15))

genes <- astroW$gene

astroW$gene <- factor(astroW$gene, levels = genes)

# plot and select threshold
ggplot(data = astroW, aes(x = gene, y = X15)) +
  geom_point() +
  geom_hline(yintercept = 1.7)

# Extract gene list
astro_genes_above_cutoff <- sum(astroW$X15 > 1.7)

astro_factor_cutoff <- data.frame(top_genes_by_factor(liger_object = astro_merged_liger_cleaned, liger_factor = 15, num_genes = astro_genes_above_cutoff))

colnames(astro_factor_cutoff) <- "astro_factor_cutoff"

# save factor list
write_rds(astro_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_astro_cutoff_gene_list.RDS")
write.csv(astro_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_astro_cutoff_gene_list.csv")

# 5.0 Oligodendrocyte Factor ---------------------------------------------------------
# Plot Factor
oligo_plots2 <- plotGeneLoadings(oligo_merged_liger_cleaned, dataset1 = "Marsh", dataset2 = "Zhou", return.plots = T, pt.size = 0.5, num.genes.show = 25)
oligo_plots2[[18]]

# pull factor and gene list
oligoW <- t(oligo_merged_liger_cleaned@W)
oligoW <- data.frame(oligoW)
oligoW <- oligoW %>% 
  rownames_to_column(var = "gene")

oligoW <- oligoW %>% 
  select(gene, X18)

oligoW <- oligoW %>% 
  arrange(desc(X18))

genes <- oligoW$gene

oligoW$gene <- factor(oligoW$gene, levels = genes)

# Plot and threshold factor
ggplot(data = oligoW, aes(x = gene, y = X18)) +
  geom_point() +
  geom_hline(yintercept = 0.2)

# Extract gene list
oligo_genes_above_cutoff <- sum(oligoW$X18 > 0.2)

oligo_factor_cutoff <- data.frame(top_genes_by_factor(liger_object = oligo_merged_liger_cleaned, liger_factor = 18, num_genes = oligo_genes_above_cutoff))

colnames(oligo_factor_cutoff) <- "oligo_factor_cutoff"

# save factor list
write_rds(oligo_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_oligo_cutoff_gene_list.RDS")
write.csv(oligo_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_oligo_cutoff_gene_list.csv")

# 6.0 OPC Factor ---------------------------------------------------------
# Plot factor
opc_plots2 <- plotGeneLoadings(opc_merged_liger_cleaned2, dataset1 = "Marsh", dataset2 = "Zhou", return.plots = T, pt.size = 0.5, num.genes.show = 25)
opc_plots2[[9]]

# pull factor and gene list
opcW <- t(opc_merged_liger_cleaned2@W)
opcW <- data.frame(opcW)
opcW <- opcW %>% 
  rownames_to_column(var = "gene")

opcW <- opcW %>% 
  select(gene, X9)

opcW <- opcW %>% 
  arrange(desc(X9))

genes <- opcW$gene

opcW$gene <- factor(opcW$gene, levels = genes)

# Plot and threshold factor
ggplot(data = opcW, aes(x = gene, y = X9)) +
  geom_point() +
  geom_hline(yintercept = 2.4)

# Extract gene list
opc_genes_above_cutoff <- sum(opcW$X9 > 2.4)

opc_factor_cutoff <- data.frame(top_genes_by_factor(liger_object = opc_merged_liger_cleaned2, liger_factor = 9, num_genes = opc_genes_above_cutoff))

colnames(opc_factor_cutoff) <- "opc_factor_cutoff"

# save factor
write_rds(opc_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_opc_cutoff_gene_list.RDS")
write.csv(opc_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_opc_cutoff_gene_list.csv")

# 7.0 Excitatory Neuron Factor ---------------------------------------------------------
excit_plots2 <- plotGeneLoadings(excit_merged_liger_cleaned2, dataset1 = "Marsh", dataset2 = "Zhou", return.plots = T, pt.size = 0.5, num.genes.show = 25)
excit_plots2[[11]]

# Extract factor and genes
excitW <- t(excit_merged_liger_cleaned2@W)
excitW <- data.frame(excitW)
excitW <- excitW %>% 
  rownames_to_column(var = "gene")

excitW <- excitW %>% 
  select(gene, X11)

excitW <- excitW %>% 
  arrange(desc(X11))

genes <- excitW$gene

excitW$gene <- factor(excitW$gene, levels = genes)

# Plot and threshold factor
ggplot(data = excitW, aes(x = gene, y = X11)) +
  geom_point() +
  geom_hline(yintercept = 1.5)

# Extract gene list
excit_genes_above_cutoff <- sum(excitW$X11 > 1.5)

excit_factor_cutoff <- data.frame(top_genes_by_factor(liger_object = excit_merged_liger_cleaned2, liger_factor = 11, num_genes = excit_genes_above_cutoff))

colnames(excit_factor_cutoff) <- "excit_factor_cutoff"

# save factor
write_rds(excit_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_excit_cutoff_gene_list.RDS")
write.csv(excit_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_excit_cutoff_gene_list.csv")

# 8.0 Inhibitory Neuron Factor ---------------------------------------------------------
# Plot factor
inhib_plots2 <- plotGeneLoadings(inhib_merged_liger_cleaned3, dataset1 = "Marsh", dataset2 = "Zhou", return.plots = T, pt.size = 0.5, num.genes.show = 25)
inhib_plots2[[18]]

# Pull factor and genes
inhibW <- t(inhib_merged_liger_cleaned3@W)
inhibW <- data.frame(inhibW)
inhibW <- inhibW %>% 
  rownames_to_column(var = "gene")

inhibW <- inhibW %>% 
  select(gene, X18)

inhibW <- inhibW %>% 
  arrange(desc(X18))

genes <- inhibW$gene

inhibW$gene <- factor(inhibW$gene, levels = genes)

# plot and threshold factor
ggplot(data = inhibW, aes(x = gene, y = X18)) +
  geom_point() +
  geom_hline(yintercept = 0.58)

# Extract gene list
inhib_genes_above_cutoff <- sum(inhibW$X18 > 0.58)

inhib_factor_cutoff <- data.frame(top_genes_by_factor(liger_object = inhib_merged_liger_cleaned3, liger_factor = 18, num_genes = inhib_genes_above_cutoff))

colnames(inhib_factor_cutoff) <- "inhib_factor_cutoff"

# save factor list
write_rds(inhib_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_inhib_cutoff_gene_list.RDS")
write.csv(inhib_factor_cutoff, "subcluster_factor_gene_list/new_thresholded_lists/liger_inhib_cutoff_gene_list.csv")