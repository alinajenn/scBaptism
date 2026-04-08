#' run_scClassify
#'
#' @param sce_query SCE to be annotated
#' @param log_name_q string with name of log-normalized query data
#' @param reference SCE object that acts as a reference
#' @param log_name_r string with name of log-normalized reference data
#' @param ref_labels Column from the references colData
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#' @param verbose display message after annotation is finished
#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#'@importFrom scClassify train_scClassify
#'@importFrom scClassify predict_scClassify
#'@importFrom SummarizedExperiment assay
#'@importFrom SummarizedExperiment colData
#'
#'
#' @examples
#'
#'library(iUSEiSEE)
#'library(dplyr)
#'# load SCE from iUSEiSEE
#'
#'sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#'#run the annotation
#'#sce_annotated <- run_scClassify(sce_query = sce_annotated,
#'                                 log_name_q = "logcounts",
#'                                 reference = sce_annotated,
#'                                 log_name_r = "logcounts",
#'                                 ref_labs = "labels_main")
#'
#'# plot the existing annotation with scater(t-SNE)
#'scater::plotTSNE(sce_annotated, color_by = "scb_scClassify_labels")
#'
#'
#'@family reference family
run_scClassify <- function(sce_query,
                           log_name_q, #String name of the log normalized matrices (query)
                           reference,
                           log_name_r,
                           ref_labs,
                           return_extra_info = FALSE,
                           verbose = FALSE,
                           ...)

{

  # transformation ----------------------------------------------------------------

  #take logcount matrices from reference & query

  log_query <- SummarizedExperiment::assay(sce_query, log_name_q)
  log_ref <- SummarizedExperiment::assay(reference, log_name_r)

  #get the labels of the reference

  ref_labels <- SummarizedExperiment::colData(reference)[[ref_labs]]

  # checks --------------------------------------------------------

  #check wether the matrices are dcG matrices and if not transform
  #scClassify requires dcG matrices as input

 if (!is(log_query, "dgCMatrix")) {

   log_query <- as(log_query, "dgCMatrix")
   print ("Transforming query into dgCMatrix")

 } else {
   print ("query is already dcGMatrix, proceeding directly")
 }

  if (!is(log_ref, "dgCMatrix")) {

    log_ref <- as(log_ref, "dgCMatrix")
    print ("Transforming reference into dgC matrix")

  } else {
    print ("reference is already dcGMatrix, proceeding directly")
  }



  #train our own model ----------------------------------------------------------

  classifier <- scClassify::train_scClassify(exprsMat_train = log_ref,
                                 cellTypes_train = ref_labels,
                                 selectFeatures = c("limma"), #different kinds of DE gene analysis methods
                                 returnList = FALSE)

  # running annotation-----------------------------------------------------
  anno_res <- scClassify::predict_scClassify(exprsMat_test = log_query,
                                 trainRes = classifier,
                                 cellTypes_test = NULL,
                                 algorithm = "WKNN",
                                 features = c("limma"),
                                 similarity = c("pearson"),
                                 prob_threshold = 0.7,
                                 verbose = FALSE)


  # return input SCE with new annotation-----------------------------------

  SummarizedExperiment::colData(sce_query)$scb_scClassify_labels <- anno_res$pearson_WKNN_limma$predRes


  #message("xx annotation done")
  if(verbose) message("scClassify annotation done")

  return(sce_query)

}
