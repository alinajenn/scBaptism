# file for all the functions checking inputs provided by the user




  # sce_query ----------------------------------------------------------


#' basic checks whether sce_query is is provided, is an SCE and contains the neccessary assays
#'
#' @param sce_query SingleCellExperiment object with assays counts and logcounts
#'
#' @returns no return, execution of package is stopped if sce_query is not a suitable input
#'
#' @noRd
.check_query <- function(sce_query){

  #check if query is provided
  if(is.null(sce_query)) {
    stop("Please provide a query object.")
  }

  if(!(is(sce_query, "SingleCellExperiment"))) {
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

#' checks if anno_methods is  the correct format (list of strings), if all strings are valid methods
#' and discards any invalid methods
#'
#' @param anno_methods list of strings
#'
#' @returns anno_methods only containing valid methods
#'
#' @noRd
.check_anno_methods <- function(anno_methods){
  #is it provided?
  if(is.null(anno_methods)) {
    stop("Please provide a character vector with the annotation methods you want to use")
  }

  #is it characters?

  if(!is.character(anno_methods)){
    stop("anno_methods needs to be a character or a vector of characters")
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

#' basic checks whether reference is a SCE and if it has annotation data in ref_labs colData
#'
#' @param reference SCE with counts, logcounts and colData that contains annotation
#' @param ref_labs string containing the name of annotation colData of reference
#'
#' @returns no return, execution of package is stopped if reference or ref_labs are not a reasonable inputs
#'
#' @noRd
.check_reference <- function(reference, ref_labs){
  if(!(is(reference, "SingleCellExperiment")) & !is.null(reference)) {
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
   if(!(ref_labs %in% names(colData(reference)))){
      stop("ref_labs needs to be the name of colData column in reference")
   }

  }
}





  # markers_list --------------------------------------------------------


#' basic checks whether markers_list is actually a nested list and contains at least one item
#'
#' @param markers_list nested list of characters
#'
#' @returns no return, execution of package is stopped if markers_list is not a suitable input
#'
#' @noRd
.check_markers_list <- function(markers_list){

  # check whether an provided object is a nested list
  if(!is.list(markers_list) & !is.null(markers_list)) {
    stop("markers_list is not a list")

    if (!all(vapply(markers_list, is.list, logical(1)))) {
      stop("markers_list is not a nested list")
    }
  }

  if(!is.null(markers_list)) {
    if(length(markers_list) < 1){
      stop("markers_list does not contain anything")
    }
  }

}



