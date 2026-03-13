#' run_CelliDmk
#'
#' @param sce_query SCE to be annotated
#' @param markers_list List of marker genes
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#' @param verbose display message after annotation is finished
#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' export
#'
#'@importFrom CelliD RunCellHGT
#'@importFrom CelliD RunMCA
#'@importFrom SingeCellExperiment colData
#'
#' @examples
#'
#'
#' library(iUSEiSEE)
#' library(dplyr)
#'
#' load SCE from iUSEiSEE
#'
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' run the annotation
#' sce_annotated <- run_tool(sce_annotated, markers_lists)
#'
#' plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_CelliDmk_labels")
#'
#'
#'@family hybrid family
run_CelliDmk <- function(sce_query,
                         markers_list,
                         return_extra_info = FALSE,
                         verbose = FALSE,
                         ...)

{

  # checks ----------------------------------------------------------------


  # transformation --------------------------------------------------------
  #sce input is fine

  #input markers: list of vectors (characters)


  # running annotation-----------------------------------------------------

  #dimansionality reduction MCA is neccessary for running CelliD
  sce_query <- CelliD::RunMCA(sce_query)


  # run marker based CelliS version

  result_marker <- CelliD::RunCellHGT(sce_query,
                              pathways = markers_list,
                              dims = 1:50,
                              n.features = 200
                              )


  # return input SCE with new annotation-----------------------------------

  #for each cell, assess the signature with the lowest corrected p-value (max -log10 corrected p-value)
  marker_prediction <- rownames(result_marker)[apply(result_marker, 2, which.max)]

  # for each cell, evaluate if the lowest p-value is significant
  marker_prediction_signif <- ifelse(apply(result_marker, 2, max)>2, yes = marker_prediction, "unassigned")

  # add new annotation to SCE query
  SummarizedExperiment::colData(sce_query)$scb_CelliDmk_labels <- marker_prediction_signif



  #message("CelliD annotation done")
  if(verbose) message("CelliD annotation done")

  return(sce_query)

}
