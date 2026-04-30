#' run_scClassify
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE object that acts as a reference
#' @param ref_labs Column from the references colData
#' @param selectFeatures string vector deciding which methods to use for feature selection during training. Defaults to "limma", other options:"DV", "DD", "chisq", "BI", "Cepo"
#' @param cellTypes_test A list or a vector indicates cell types of the query datasets (Optional).
#' @param k An integer indicates the number of neighbors
#' @param prob_threshold A numeric indicates the probability threshold for KNN/WKNN/DWKNN.
#' @param cor_threshold_static A numeric indicates the static correlation threshold.
#' @param cor_threshold_high A numeric indicates the highest correlation threshold
#' @param scClassify_features A vector indicates the gene selection method, set as "limma" by default. This should be one or more of "limma", "DV", "DD", "chisq", "BI".
#' @param algorithm A vector indicates the KNN method that are used, set as "WKNN" by default. This should be one or more of "WKNN", "KNN", "DWKNN".
#' @param similarityvA vector indicates the similarity measure that are used, set as "pearson" by default. This should be one or more of "pearson", "spearman", "cosine", "jaccard", "kendall", "binomial", "weighted_rank","manhattan"
#' @param cutoff_methodvA vector indicates the method to cutoff the correlation distribution. Set as "dynamic" by default.
#' @param weighted_ensemble A logical input indicates in ensemble learning, whether the results is combined by a weighted score for each base classifier.
#' @param weights A vector indicates the weights for ensemble
#' @param parallel A logical input indicates whether running in paralllel or not
#' @param BPPARAM A BiocParallelParam class object from the BiocParallel package is used. Default is SerialParam().
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#' @param verbose display message after annotation is finished
#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @importFrom methods as is
#' @importFrom scClassify predict_scClassify train_scClassify
#' @importFrom SummarizedExperiment assay assays colData
#'
#' @export
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
#' sce_annotated <- run_scClassify(sce_query = sce_annotated,
#'                                 reference = sce_annotated,
#'                                 ref_labs = "labels_main")
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_scClassify_labels")
#'
#'
#'@family reference family
run_scClassify <- function(sce_query,
                           reference,
                           ref_labs,
                           cellTypes_test = NULL,
                           k = 10,
                           prob_threshold = 0.7,
                           cor_threshold_static = 0.5,
                           cor_threshold_high = 0.7,
                           scClassify_features = "limma",
                           algorithm = "WKNN",
                           similarity = "pearson",
                           cutoff_method = c("dynamic", "static"),
                           weighted_ensemble = FALSE,
                           weights = NULL,
                           parallel = FALSE,
                           BPPARAM = BiocParallel::SerialParam(),
                           return_extra_info = FALSE,
                           verbose = FALSE
                           )

{

  # checks -----------------------------------------------------------------------

  #are the logcounts called logcounts?

  if(!("logcounts" %in% names(assays(sce_query)))) {
    stop("no assay called 'logcounts' found. Rename your assay or calculate the logcounts.")
  }

  # transformation ----------------------------------------------------------------

  #take logcounts matrices from reference & query

  log_query <- SummarizedExperiment::assay(sce_query, "logcounts")
  log_ref <- SummarizedExperiment::assay(reference, "logcounts")

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
                                 selectFeatures = scClassify_features, #different kinds of DE gene analysis methods
                                 returnList = FALSE)

  # running annotation-----------------------------------------------------
  anno_res <- scClassify::predict_scClassify(exprsMat_test = log_query,
                                 trainRes = classifier,
                                 cellTypes_test = cellTypes_test,
                                 k = k,
                                 prob_threshold = prob_threshold,
                                 cor_threshold_static = cor_threshold_static,
                                 cor_threshold_high = cor_threshold_high,
                                 algorithm = algorithm,
                                 features = scClassify_features, #needs to be same as used for train_scClassify
                                 similarity = similarity,
                                 cutoff_method = cutoff_method,
                                 weighted_ensemble = weighted_ensemble,
                                 weights = weights,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 verbose = verbose
                                 )


  # return input SCE with new annotation-----------------------------------

  SummarizedExperiment::colData(sce_query)$scb_scClassify_labels <- anno_res$pearson_WKNN_limma$predRes


  #message
  if(verbose) message("scClassify annotation done")

  return(sce_query)

}
