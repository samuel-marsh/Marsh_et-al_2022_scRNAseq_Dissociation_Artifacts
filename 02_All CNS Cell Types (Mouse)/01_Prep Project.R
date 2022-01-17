# Start new RStudio Project File
# Isolate libraries for package versioning

# Load packrat & start
library(packrat)
init()

# Install versions package
install.packages("versions")

# Load versions and install properly dataed packages
library(versions)
install.dates("tidyverse", "2019-08-28")

install.dates("BiocManager", "2019-08-28")

# Install Seurat dependecies
install.dates(c("cowplot", "ROCR", "mixtools", "lars", "ica", "tsne", "Rtsne", "fpc", "ape", "pbapply", "igraph", "RANN", "irlba", "gplots", "dtw", "SDMTools", "plotly", "Hmisc", "ggridges", "metap", "lmtest", "fitdistrplus", "png", "doSNOW", "reticulate", "foreach", "hdf5r", "RcppEigen", "RcppProgress"), "2019-08-28")

# Install Seurat from CRAN archive
#wget https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_2.3.4.tar.gz
install.packages("Seurat_2.3.4.tar.gz", repos = NULL, type = "source")