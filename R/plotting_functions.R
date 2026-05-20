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
#' @examples
#' #load example SCE from iUSEiSEE package
#'
#' sce_annotated <- readRDS(
#'  file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' #plot two of the exisiting annotations
#'
#' plot <- plot_multiple_tSNE(sce_query = sce_annotated,
#'                           labels_vector = c("labels_main", "labels_fine"))
#'
#' plot
#'
#'
plot_multiple_tSNE <- function(sce_query, labels_vector) {

  #plot all selected columns
  plot <- lapply(labels_vector, function(arg){scater::plotTSNE(sce_query, color_by = arg)})

  #display all plots in one graphic
  plot <- cowplot::plot_grid(plotlist = plot, nrow = round(sqrt(length(labels_vector)), digits = 0))

  return(plot)
}



#' plot_confusion_matrix
#'
#'
#' @param sce_query SingleCellExperiment object with the annotations to plot in the colData
#' @param first_tool name of the first annotation column in the colData of sce_query
#' @param second_tool name of the second annotation column in the colData of sce_query
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient labs theme_minimal
#'
#' @returns confusion matrix comparing the annotation of two selected annotation columns
#'
#' @export
#'
#' @examples
#' #load example SCE from iUSEiSEE package
#'
#' sce_annotated <- readRDS(
#'  file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' #plot two of the exisiting annotations
#'
#' plot <- plot_confusion_matrix(sce_query = sce_annotated,
#'                               first_tool = "labels_main",
#'                               second_tool = "labels_fine")
#'
#' plot
#'
plot_confusion_matrix <- function(sce_query, first_tool, second_tool) {

  #turn all the labels into a table

  anno_table <- table(sce_query[[first_tool]], sce_query[[second_tool]])

  #plot the confusion matrix
  plot <- as.data.frame(anno_table) |>
          ggplot2::ggplot(ggplot2::aes(x = Var2, y = Var1)) +
          ggplot2::geom_tile(ggplot2::aes(fill = Freq)) +
          ggplot2::scale_fill_gradient(low = "white", high = "red") +
          ggplot2::labs(x = first_tool,
                       y = second_tool,
                       title = "Comparing annotations",
                       fill = "Number of cells") +
          ggplot2::theme_minimal()


return(plot)

}



#' plot_alluvial
#'
#'
#' @param sce_query SingleCellExperiment object with the annotations to plot in the colData
#' @param first_tool name of the first annotation column in the colData of sce_query
#' @param second_tool name of the second annotation column in the colData of sce_query
#' @param threshold number defining below which percentage of occurrence cells get filtered out (default 0.02, so 2%)
#'
#' @importFrom ggplot2 ggplot aes labs theme element_blank after_stat geom_text
#' @importFrom dplyr count filter mutate
#' @importFrom magrittr %>%
#' @importFrom ggalluvial geom_alluvium StatStratum
#' @importFrom cowplot theme_minimal_hgrid
#'
#' @returns alluvial plot comparing the annotation of two selected annotation columns
#'
#' @export
#'
#' @examples
#' #load example SCE from iUSEiSEE package
#'
#' sce_annotated <- readRDS(
#'  file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' #plot two of the exisiting annotations
#'
#' plot <- plot_alluvial(sce_query = sce_annotated,
#'                       first_tool = "labels_main",
#'                       second_tool = "labels_fine",
#'                       threshold = 0.02)
#'
#' plot
#'
plot_alluvial <- function(sce_query, first_tool, second_tool, threshold = 0.02) {

  #get the data frame from the sce_query
  df_compare <- data.frame(
    first_anno = sce_query[[first_tool]],
    second_anno = sce_query[[second_tool]]
    )

  filtered_first <- df_compare %>%
    dplyr::count(first_anno) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::filter(prop < threshold) %>%
    dplyr::pull(first_anno)

  filtered_second <- df_compare %>%
    dplyr::count(second_anno) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::filter(prop < threshold) %>%
    dplyr::pull(second_anno)

  df_compare2 <- df_compare %>%
    dplyr::mutate(
      first_anno = ifelse(first_anno %in% filtered_first, "other", first_anno),
      second_anno = ifelse(second_anno %in% filtered_second, "other", second_anno)
    )


  df_agg <- df_compare2 %>%
    dplyr::count(first_anno, second_anno)


  plot <- ggplot2::ggplot(df_agg,
                  ggplot2::aes(axis1 = first_anno, axis2 = second_anno, y = n)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = first_anno)) +
    ggalluvial::geom_stratum() +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(stratum)),
                       size = 2.5) +
    cowplot::theme_minimal_hgrid() +
    ggplot2::theme(axis.title.x= ggplot2::element_blank(),
                  axis.text.x= ggplot2::element_blank(),
                  axis.ticks.x= ggplot2::element_blank()) +
    ggplot2::labs(
         y = "Number of cells",
         title = "Comparing annotations",
       fill = "Cell type")

  return(plot)
}

