#' run_scBaptism
#'
#' @param sce_query SCE to be annotated
#' @param anno_methods vector with all annotation methods to be used for annotation, default is all methods
#' @param markers_list list of marker genes
#' @param reference SingleCellExperiment object with annotation information column
#' @param ref_labs name of column with annotation information
#' @param return_extra_info if TRUE, adds additional metadata from the annotation when available
#' @param verbose display message after each annotation is finished
#'
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'


run_scBaptism <- function(sce_query = NULL,
                          anno_methods = NULL,
                          markers_list = NULL,
                          reference = NULL,
                          ref_labs = NULL,
                          verbose = FALSE,
                          return_extra_info = FALSE,
                          ...)

  {

  # checks -------------------------------------------------------------

  # checks of input parameters

  .check_query(sce_query)
  anno_methods <- .check_anno_methods(anno_methods)
  .check_markers_list(markers_list)
  .check_reference(reference, ref_labs)

  #check is selected tools match provided inputs

  anno_methods_final <- anno_methods

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
    warning("No markers list provided. and you wanted any marker based methods")

    anno_methods_final <- setdiff(anno_methods_final, anno_markerbased)
    if(length(anno_methods_final) == 0) {
      stop("Please provide a list of markers")
    }
  }



  # run selected tools -------------------------------------------------
  for (item in anno_methods_final) {


  sce_query <- switch(item,
  'CIA' = sce_query <- run_CIA(sce_query = sce_query, markers_list = markers_list),
  'SCINA' = sce_query <- run_SCINA(sce_query = sce_query, markers_list = markers_list),
  'CelliDmk' = sce_query <- run_CelliDmk(sce_query = sce_query, markers_list = markers_list),
  'SingleR' = sce_query <- run_SingleR(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'Seurat' = sce_query <- run_Seurat(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'clustifyr' = sce_query <- run_clustifyr(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'scPred' = sce_query <- run_scPred(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'scClassify' = sce_query <- run_scClassify(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'scmap' = sce_query <- run_scmap(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'CelliDref' = sce_query <- run_CelliDref(sce_query = sce_query, reference = reference, ref_labs = ref_labs)
  )

  }

  return(sce_query)
}

