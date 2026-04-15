#' run_scBaptism
#'
#' @param sce_query SCE to be annotated
#' @param anno_methods vector with all annotation methods to be used for annotation, default is all methods
#' @param markers_list list of marker genes
#' @param reference SingleCellExperiment object with annotation information column
#' @param ref_labs name of column with annotation information
#' @param name description
#' @param return_extra_info if TRUE, adds additional metadata from the annotation when available
#' @param verbose display message after each annotation is finished
#'
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom
#'
#' examples
#'
#'
#'


run_scBaptism <- function(sce_query,
                          anno_methods,
                          markers_list,
                          reference,
                          ref_labs,
                          log_name_r,
                          log_name_q,
                          verbose = FALSE,
                          return_extra_info = FALSE,
                          ...)

  {

  # checks -------------------------------------------------------------
  #log counts?
  #do selected tools match provided inputs?

  # run selected tools -------------------------------------------------
  for (item in anno_methods) {


  sce_query <- switch(item,
  'CIA' = sce_query <- run_CIA(sce_query = sce_query, markers_list = markers_list),
  'SCINA' = sce_query <- run_SCINA(sce_query = sce_query, markers_list = markers_list),
  'CelliDmk' = sce_query <- run_CelliDmk(sce_query = sce_query, markers_list = markers_list),
  'SingleR' = sce_query <- run_SingleR(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'Seurat' = sce_query <- run_Seurat(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'clustifyr' = sce_query <- run_clustifyr(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'scPred' = sce_query <- run_scPred(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'scClassify' = sce_query <- run_scClassify(sce_query = sce_query, log_name_r = log_name_r, log_name_q = log_name_q, reference = reference, ref_labs = ref_labs),
  'scmap' = sce_query <- run_scmap(sce_query = sce_query, reference = reference, ref_labs = ref_labs),
  'CelliDref' = sce_query <- run_CelliDref(sce_query = sce_query, reference = reference, ref_labs = ref_labs)
  )

  }

  return(sce_query)
}
