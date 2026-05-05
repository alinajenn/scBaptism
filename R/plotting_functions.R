# Viz file

#' plot_multiple_tSNE
#'
#'
#' @param sce_query SingleCellExperiment object with the annotations to plot in the colData
#' @param labels_vector Vector with names of the colData
#'
#' @importFrom cowplot plot_grid
#' @importFrom scater plotTSNE
#'
#' @returns plot with all selected annotations plotted with tSNE dimensionality reduction
#'
#' @export
#'
plot_multiple_tSNE <- function(sce_query, labels_vector) {

  #plot all selected columns
  plot <- lapply(labels_vector, function(arg){scater::plotTSNE(sce_query, color_by = arg)})

  #display all plots in one graphic
  plot <- cowplot::plot_grid(plotlist = plot, nrow = round(sqrt(length(vec)), digits = 0))

  return(plot)
}

#confusion matrix

#alluvial plot
