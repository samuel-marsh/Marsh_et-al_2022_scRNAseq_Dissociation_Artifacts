# clusterLouvainJaccard function
    # Now deprecated in LIGER

clusterLouvainJaccard = function(object, resolution = 0.1, k.param=30, n.iter = 10, n.start = 10,
                                 print.output = F, ...)
{
  if (!require("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  temp.seurat = CreateSeuratObject(t(Reduce(rbind,object@scale.data)))
  temp.seurat@scale.data = t(Reduce(rbind,object@scale.data))
  rownames(object@H.norm)=colnames(temp.seurat@scale.data)
  temp.seurat@dr$NMF=new(Class="dim.reduction",cell.embeddings=object@H.norm,key="NMF")
  temp.seurat <- FindClusters(object = temp.seurat, reduction.type = "NMF",
                              dims.use = 1:ncol(object@H.norm),force.recalc=T,
                              save.SNN = T,resolution=resolution,k.param=k.param,
                              n.iter = n.iter, n.start = n.start, print.output = print.output, ...)
  object@clusters = temp.seurat@ident
  return(object)
}



# Custom Plotting by Sample LIGER
DimPlot_All_Samples_Liger <- function(liger_object, cell_data_column = "dataset", color = "black", num_col = NULL, label_size = 14, ...){
  # Pull dataset info
  column_list <- as.character(unique(liger_object@cell.data[, "dataset"]))
  # Pull reduc coordinates
  reduc_coordinates <- data.frame(liger_object@tsne.coords)
  colnames(reduc_coordinates) <- c("dr1", "dr2")
  x_axis <- c(min(reduc_coordinates[, 1]),
              max(reduc_coordinates[, 1]))
  y_axis <- c(min(reduc_coordinates[, 2]),
              max(reduc_coordinates[, 2]))
  # Plotting Function
  Dim_by_Sample_Liger <- function(liger_object, column_list){
    # Pull Cell Barcodes
    cell_barcodes <- row.names(liger_object@cell.data)[which(liger_object@cell.data[, cell_data_column] == column_list)]
    
    # Subset liger
    temp_liger <- suppressWarnings(subsetLiger(object = liger_object, cells.use = cell_barcodes, remove.missing = FALSE))
    
    # Plot
    plot <- plotByDatasetAndCluster(temp_liger, return.plots = TRUE, do.legend = TRUE, text.size = 6)
    plot <- plot[[1]] +
      scale_color_manual(values = "black") +
      ggtitle(paste(column_list)) +
      theme(plot.title = element_text(hjust = 0.5, size = label_size),
            legend.position = "none") +
      xlim(x_axis) +
      ylim(y_axis)
  }
  # Iterate and create plot list
  all_plots <- map(column_list, ~Dim_by_Sample_Liger(liger_object, .x))
  # Wrap Plots into single output
  wrap_plots(all_plots, ncol = num_col)
}


