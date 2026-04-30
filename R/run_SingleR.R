#' run_SinlgeR
#'
#' @param sce_query SCE object to be annotated
#' @param reference SCE object that acts as a reference
#' @param ref_labs List of gene labels or column from the references colData
#' @param return_extra_info if TRUE, adds additional metadata from the annotation (delta.next, scores, prunded.labels)
#' @param verbose display message after annotation is finished
#' @param clusters vector or factor of cluster identities for each cell in query
#' @param restrict character vector of gene names used for marker selection, default is all genes
#' @param check.missing.test  binary, whether rows of query should be tested for missing values
#' @param check.missing.ref  binary, whether rows of reference should be tested for missing values
#' @param num.threads Integer specifying number of threads to use for index building & classification
#' @param BPPARAM BiocParallelParam specifying how parallelization should be performed
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom SingleR SingleR
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#'
#' library(iUSEiSEE)
#' library(dplyr)
#'
#'
#' # load SCE from iUSeiSEE
#'
#' sce_annotated <- readRDS(
#'   file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' #using the labels_main of the example to run another annotation with SingleR
#' sce_annotated <- run_SingleR(sce_query = sce_annotated,
#'                           reference = sce_annotated,
#'                           ref_labs = "labels_main")
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_SingleR_labels")
#'
#'
#'@family reference-based family
run_SingleR <- function(sce_query,
                        reference,
                        ref_labs,
                        clusters = NULL,
                        restrict = NULL,
                        check.missing.test = FALSE,
                        check.missing.ref = FALSE,
                        num.threads = bpnworkers(BPPARAM),
                        BPPARAM = SerialParam(),
                        return_extra_info = FALSE,
                        verbose = FALSE
                        )

{

  # checks ----------------------------------------------------------------


  # transformation --------------------------------------------------------

  #query & reference: SCE is accepted as input

  #extract labels column from reference SCE

  labels_col <- SummarizedExperiment::colData(reference)[[ref_labs]]

  # running annotation-----------------------------------------------------
  SingleR_res <- SingleR::SingleR(test = sce_query,
                                ref = reference,
                                labels = labels_col,
                                clusters = clusters,
                                restrict = restrict,
                                check.missing.test = check.missing.test,
                                check.missing.ref = check.missing.ref,
                                num.threads = num.threads,
                                BPPARAM = BPPARAM
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
