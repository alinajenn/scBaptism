#' run_scBaptism
#'
#' @param sce_query SCE to be annotated
#' @param anno_methods vector with all annotation methods to be used for annotation, default is all methods
#' @param markers_list list of marker genes
#' @param reference SingleCellExperiment object with annotation information column
#' @param ref_labs name of column with annotation information
#' @param return_extra_info if TRUE, adds additional metadata from the annotation when available
#' @param verbose display message after each annotation is finished
#' @param ... further arguments passed to other methods
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @examples
#' library(iUSEiSEE)
#' library(Seurat)
#' library(dplyr)
#'
#'
#' #load example SCE from iUSEiSEE package (also used as reference)
#'
#' sce_annotated <- readRDS(
#'  file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' #get the markers list using Seurat
#' myseu <- Seurat::as.Seurat(sce_annotated)
#' myseu <- Seurat::ScaleData(myseu)
#' Seurat::Idents(myseu) <- "labels_main"
#'
#' seu_all_markers <- Seurat::FindAllMarkers(myseu, test.use = "wilcox", only.pos = TRUE,
#'                                          min.pct = 0.25, logfc.threshold = 0.25)
#'
#' top_k_markers <- 50
#'
#' markers_lists <- seu_all_markers %>%
#'  dplyr::group_by(cluster) %>%
#'  dplyr::top_n(n = top_k_markers, wt = avg_log2FC)
#' markers_lists <- split(markers_lists$gene, markers_lists$cluster)
#'
#' #run scBaptism to annotate with SingleR and CIA
#' anno_result <- run_scBaptism(sce_query = sce_annotated,
#'                              anno_methods = c("SCINA", "CIA"),
#'                              markers_list = markers_lists,
#'                              reference = sce_annotated,
#'                              ref_labs = "labels_main")
#'


run_scBaptism <- function(sce_query = NULL,
                          anno_methods = NULL,
                          markers_list = NULL,
                          reference = NULL,
                          ref_labs = NULL,
                          verbose = FALSE,
                          return_extra_info = FALSE,
                          ...){

  # checks -------------------------------------------------------------

  # checks of each input parameter

  .check_query(sce_query)
  anno_methods <- .check_anno_methods(anno_methods)
  .check_markers_list(markers_list)
  .check_reference(reference, ref_labs)


  #check if selected tools match the  provided inputs

  anno_methods_final <- .match_input_methods(anno_methods, reference, markers_list)

  # prepare parameter inputs ---------------------------------------------------

  # get all arguments users might have passed through the ...

  dots <- list(...)


  # get the extra arguments sorted into lists to pass them to the correct function



  args_CIA <- dots[names(dots) %in% names(formals(run_CIA))]
  args_SCINA <- dots[names(dots) %in% names(formals(run_SCINA))]
  args_CelliDmk <- dots[names(dots) %in% names(formals(run_CelliDmk))]
  args_SingleR <- dots[names(dots) %in% names(formals(run_SingleR))]
  args_Seurat <- dots[names(dots) %in% names(formals(run_Seurat))]
  args_clustifyr <- dots[names(dots) %in% names(formals(run_clustifyr))]
  args_scPred <- dots[names(dots) %in% names(formals(run_scPred))]
  args_scClassify <- dots[names(dots) %in% names(formals(run_scClassify))]
  args_scmap <- dots[names(dots) %in% names(formals(run_scmap))]
  args_CelliDref <- dots[names(dots) %in% names(formals(run_CelliDref))]



  # run selected tools ------------------------------------------------------
  for (item in anno_methods_final) {

    sce_query <- switch(EXPR = item,
    'CIA' = do.call(run_CIA,
                           c(list(sce_query = sce_query,
                                  markers_list = markers_list,
                                  verbose = verbose),
                            args_CIA)),

    'SCINA' =  do.call(run_SCINA,
                               c(list(sce_query = sce_query,
                                      markers_list = markers_list,
                                      verbose = verbose,
                                      return_extra_info = return_extra_info),
                               args_SCINA)),

    'CelliDmk' =  do.call(run_CelliDmk,
                                      c(list(sce_query = sce_query,
                                             markers_list = markers_list,
                                             verbose = verbose,
                                             return_extra_info = return_extra_info),
                                        args_CelliDmk)),

    'SingleR' =  do.call(run_SingleR,
                                     c(list(sce_query = sce_query,
                                            reference = reference,
                                            ref_labs = ref_labs,
                                            verbose = verbose,
                                            return_extra_info = return_extra_info),
                                       args_SingleR)),

    'Seurat' = do.call(run_Seurat,
                                    c(list(sce_query = sce_query,
                                           reference = reference,
                                           ref_labs = ref_labs,
                                           verbose = verbose,
                                           return_extra_info = return_extra_info),
                                      args_Seurat)),

    'clustifyr' =  do.call(run_clustifyr,
                                       c(list(sce_query = sce_query,
                                              reference = reference,
                                              ref_labs = ref_labs,
                                              verbose = verbose),
                                         args_clustifyr)),

    'scPred' =  do.call(run_scPred,
                                    c(list(sce_query = sce_query,
                                           reference = reference,
                                           ref_labs = ref_labs,
                                           verbose = verbose,
                                           return_extra_info = return_extra_info),
                                      args_scPred)),

    'scClassify' =  do.call(run_scClassify,
                                        c(list(sce_query = sce_query,
                                               reference = reference,
                                               ref_labs = ref_labs,
                                               verbose = verbose,
                                               return_extra_info = return_extra_info),
                                          args_scClassify)),

    'scmap' =  do.call(run_scmap,
                                   c(list(sce_query = sce_query,
                                          reference = reference,
                                          ref_labs = ref_labs,
                                          verbose = verbose,
                                          return_extra_info = return_extra_info),
                                     args_scmap)),

    'CelliDref' =  do.call(run_CelliDref,
                                       c(list(sce_query = sce_query,
                                              reference = reference,
                                              ref_labs = ref_labs,
                                              verbose = verbose,
                                              return_extra_info = return_extra_info),
                                         args_CelliDref))
    )

  }

  # return result --------------------------------------------------------------------

  return(sce_query)
}

