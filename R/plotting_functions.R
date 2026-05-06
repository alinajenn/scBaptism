# file for all the plots
#
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

plot_confusion_matrix <- function(sce_query, first_tool, second_tool) {

  #turn all the labels into a table
  anno_table <- table(sce_query[[first_tool]], sce_query[[second_tool]])


  #delete all the row with value 0 (good idea??? probably not for the confusion matrix)

  #plot the confusion matrix
  plot <- anno_tabledf |>
          ggplot(aes(x = Var2, y = Var1)) +
          geom_tile(aes(fill = Freq)) +
          labs(x = first_tool,
               y = second_tool,
               title = "Comparing annotations",
               fill = "Number of cells")

  return(plot)

}


  #not the nicest color scale

#alluvial plot


  #anno_table <- table(sce_query[[first_tool]], sce_query[[second_tool]])
