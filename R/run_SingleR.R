#' run_SinlgeR
#'
#' @param sce_query SCE object to be annotated
#' @param reference SCE object that acts as a reference
#' @param ref_labels List of gene labels or column from the references colData
#'
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom SingleR SingleR
#' @importFrom SingleCellExperiment colData
#'
#' @examples
#'
#' library(iUSEiSEE)
#' library(dplyr)
#'
#'
#' # load SCE from iuseisee
#'
#' sce_annotated <- readRDS(
#'   file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE")
#' )
#'
#' #using the labels_main of the example to run another annotation with SingleR
#' sce_query <- run_SingleR(sce_query = sce_annotated, reference = sce_annotated, ref_labels = sce_annotated$labels_main)
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "labels_main")
#'
#'
#'@family reference-based family
run_SingleR <- function(sce_query, #SummarizedExperiment
                        reference, #SummarizedExperiment
                        ref_labels,#List or column of your SCE etc
                        return_extra_info = FALSE,
                        verbose = FALSE,
                      ...)

{

  # checks ----------------------------------------------------------------


  # transformation --------------------------------------------------------
  #ref needs to be normalized and log-transformed
  #query: sce is accepted as input

  #labels column from sce or list, no transformation needed

  # running annotation-----------------------------------------------------
  SingleR_res <- SingleR::SingleR(test = sce_query,
                                ref = reference,
                                labels = ref_labels,
                                ...
                                )

  # return input SCE with new annotation-----------------------------------


  SummarizedExperiment::colData(sce_query)$scb_SingleR_labels <- SingleR_res$labels


  #additional data SingleR provides along the annotation
  if (return_extra_info){
    SummarizedExperiment::colData(sce_query)$scb_SingleR_delta.next <- SingleR_res$delta.next
    SummarizedExperiment::colData(sce_query)$scb_SingleR_scores <- SingleR_res$scores
    SummarizedExperiment::colData(sce_query)$scb_SingleR_prunded.labels <- SingleR_res$prunded.lables
  }




  #sce_out add to this sce
  if(verbose) message("SingleR annotation done")
  return(sce_query)



}
