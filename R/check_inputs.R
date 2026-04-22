#' check_inputs
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE with annotations
#' @param ref_labs existing annotation in the reference
#'
#' @returns
#'
#' export
#'
#' @examples
#'
#' library(iUSEiSEE)
#' library(dplyr)
#'
#' # load SCE from iUSEiSEE
#'
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' # run the annotation
#' sce_annotated <- run_CelliDref(sce_annotated, markers_lists)
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_CelliDref_labels")
#'
#'
###############################documentation missing for all


  # sce_query ----------------------------------------------------------

.check_query <- function(sce_query){

  #check if query is provided
  if(is.null(sce_query)) {
    stop("Please provide a query object.")
  }

  if(!(class(sce_query) == "SingleCellExperiment")) {
    stop("Query input is not a SingleCellExperiment object.")
  }

  if(!("counts" %in% names(SummarizedExperiment::assays(sce_query)))) {
    stop("Query: no assay called 'counts' found. Rename your assay or add the counts")
  }

  if(!("logcounts" %in% names(SummarizedExperiment::assays(sce_query)))) {
    stop("Query: no assay called 'logcounts' found. Rename your assay or calculate the logcounts.")
  }

}

  # anno_methods -------------------------------------------------------

.check_anno_methods <- function(anno_methods){
  #is it provided?
  if(is.null(anno_methods)) {
    stop("Please provide a character vector with the annotation methods you want to use")
  }

  #is it characters?

  if(!is.character(anno_methods)){
    stop("anno_methods needs to be a characters or a vector of characters")
  }

  #are all the provided methods valid?

  anno_methods_input_length <- length(anno_methods)

  valid_anno_methods <- c("SingleR", "Seurat", "clustifyr", "scPred", "scClassify", "scmap", "CelliDref", "CIA", "SCINA", "CelliDmk")

  anno_methods <- anno_methods[anno_methods %in% valid_anno_methods]


  #any methods left after discarding non-valid ones?
  if(is.null(anno_methods)) {
    stop("No valid annotation methods provided. Please select one or more of: SingleR, Seurat, clustifyr, scPred, scClassify, scmap, CelliDref, CIA, SCINA, CelliDmk")
  }

  #send message to user if some methods were deleted due to invalidity

  if((length(anno_methods)) < anno_methods_input_length & (length(anno_methods > 0))){
    message("Some annotation methods provided are not valid options, proceeding with only the valid methods")
  }

  if(length(anno_methods) == 0) {
    stop("No valid annotations methods specified")
  }

  return(anno_methods)
}

  # reference ----------------------------------------------------------


.check_reference <- function(reference, ref_labs){
  if(!(class(reference) == "SingleCellExperiment") & !is.null(reference)) {
    stop("Reference input is not a SingleCellExperiment object.")
  }


  if(!is.null(reference)) {
    if(!("counts" %in% names(SummarizedExperiment::assays(reference)))) {
     stop("reference: no assay called 'counts' found. Rename your assay or add the counts")
    }
  }

  if(!is.null(reference)) {
   if(!("logcounts" %in% names(SummarizedExperiment::assays(reference)))) {
      stop("reference: no assay called 'logcounts' found. Rename your assay or calculate the logcounts.")
    }
  }

  #check: if there are ref_labs, is the reference also provided


  if(!(is.null(ref_labs)) & is.null(reference)) {
    stop("Please provide the reference SingleCellExperiment object")
  }


  #check if ref_labs is character of length 1

  if(!is.character(ref_labs) & !is.null(ref_labs)) {
    stop("ref_labs needs to be a character")

    if(length(ref_labs) != 1) {
      stop("ref_labs needs to be a character of length 1")
    }
  }

  #check if ref_labs is colData in reference

  if(!is.null(reference)) {
   if(is.null(SingleCellExperiment::colData(reference)[[ref_labs]])){
      stop("ref_labs needs to be the name of colData column in reference")
   }

  }
}

  # markers_list --------------------------------------------------------


.check_markers_list <- function(markers_list){

  # check whether an provided object is a list
  if(!is.list(markers_list) & !is.null(markers_list)) {
    stop("markers_list is not a list")
  }

}



