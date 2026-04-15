#' run_clustifyr
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE object that acts as a reference
#' @param ref_labs List of gene labels or column from the references colData
#' @param verbose display message after annotation is finished
#'
#'
#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' export
#'
#' @importFrom clustifyr clustify
#' @importFrom clustifyr object_ref
#' @importFrom SummarizedExperiment colData
#'
#'
#' @examples
#'
#' library(iUSEiSEE)
#' library(Seurat)
#' library(dplyr)
#'
#' #load SCE from iUSEiSEE
#'
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' #run the annotation
#' sce_annotated <- run_clustifyr(sce_query = sce_annotated,
#'                                reference = sce_annotated,
#'                                ref_labs = "labels_main")
#'
#' #plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_clustifyr_labels")
#'
#'
#'@family reference-based family
#'
run_clustifyr <- function(sce_query,
                          reference,
                          ref_labs, #string of name of annotation column of SCE
                          verbose = FALSE,
                          ...)

{

  # checks ----------------------------------------------------------------


  # transformation --------------------------------------------------------

  #input is fine as SCE

  #bulid the reference with object_ref() function from clustifyr

  transf_ref <- clustifyr::object_ref(
    input = reference,
    cluster_col = ref_labs # they call it cluster_col, but ask for cell identities
  )


  # running annotation-----------------------------------------------------

  sce_query <- clustifyr::clustify(
    input = sce_query, # an SCE object
    ref_mat = transf_ref, # matrix of RNA-seq expression data for each cell type
    cluster_col = ref_labs, # name of column in meta.data containing cell clusters
    obj_out = TRUE, # output SCE object with cell type inserted as "type" column
    ...)


  # return input SCE with new annotation-----------------------------------


  #gets added automatically to the sce_query, but need to rename column

  names(SummarizedExperiment::colData(sce_query))[which(names(SummarizedExperiment::colData(sce_query))=="type")]="scb_clustifyr_labels"

  #rename column for correlation coefficient

  names(SummarizedExperiment::colData(sce_query))[which(names(SummarizedExperiment::colData(sce_query))=="r")]="scb_clustifyr_corr_coef"

  #message
  if(verbose) message("clustifyr annotation done")

  return(sce_query)

}
