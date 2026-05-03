
#' checks whether the input parameters provided by the user match with the annotation methods they selected
#' marker based tools need a markers_list, reference based tools need a reference with an annotation column
#'
#' @param reference SingleCellExperiment object with annotation column
#' @param markers_list List of marker genes
#' @param anno_methods vector with names of annotation methods the user wants to perform
#'
#' @returns anno_methods_final, a list with only the methods that can be run with the inputs provided
#'
#' @noRd

.match_input_methods <- function(anno_methods, reference, markers_list)

  {

  anno_methods_final <- anno_methods

  #all valid methods, sorted by reference needed vs markers needed

  anno_referencebased <- c("SingleR", "Seurat", "clustifyr", "scPred", "scClassify", "scmap", "CelliDref")
  anno_markerbased <- c("CIA", "SCINA", "CelliDmk")

  #check if a reference is provided, when user selects reference-based tools

  if(is.null(reference) & any(anno_methods_final %in% anno_referencebased)) {
    warning("No reference provided. Annotation will be performed with selected marker-based methods only")
    anno_methods_final <- setdiff(anno_methods_final, anno_referencebased)
    anno_methods_excluded <- setdiff(anno_referencebased, anno_methods_final)
    if(length(anno_methods_final) == 0) {
      stop("Please provide a reference")
    }
  }

  #check if markers list is provided, when user selects marker-based tools
  if(is.null(markers_list) & any(anno_methods_final %in% anno_markerbased)) {
    warning("No markers list  provided. Annotations cannot be provided without")

    anno_methods_final <- setdiff(anno_methods_final, anno_markerbased)
    if(length(anno_methods_final) == 0) {
      stop("Please provide a list of markers")
    }
  }

  return(anno_methods_final)
}
