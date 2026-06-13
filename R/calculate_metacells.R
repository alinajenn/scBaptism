#' calculate_metacells
#'
#' @param sce_query SingleCellExperiment object with logcounts to calculate metacells on
#' @param annotation_col name of annotation column in sce_query, set to NULL if no annotation should be transferred to metacell level
#' @param gamma set the graining level for the metacells, aka proportion of single cells to meta cells
#' @param n.pc number of PCAs to be used during the calculation of the metacell identities, default is 10
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SuperCell SCimplify supercell_GE supercell_2_sce
#' @importFrom scran buildSNNGraph
#' @importFrom igraph cluster_walktrap
#'
#' @returns SCE with logcounts on metacell level, also clusters, PCA & tSNE
#'
#' @export
#'
#' @examples
#' # create a SCE on metacell level without prior annotation
#' #load example SCE from iUSEiSEE package
#'
#' sce_annotated <- readRDS(
#'  file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' sce_metacells <- calculate_metacells(sce_query = sce_annotated,
#'                                      annotation_col = NULL)
#'
#' sce_metacells
#'
#'
calculate_metacells <- function(sce_query,
                                annotation_col = NULL,
                                gamma = 10, #SCimplify graining level
                                n.pc = 10
                                ){

  #get the log normalized expression matrix

  query_matrix <- SummarizedExperiment::assay(sce_query, "logcounts")

  #run the metacell calculations, metacell identity for each single cell

  #with or without existing annotation

  if (!is.null(annotation_col)) {
    mc_identity <- SuperCell::SCimplify(X = query_matrix,
                                       cell.annotation = sce_query[[annotation_col]],
                                       gamma = gamma,
                                       n.pc = n.pc)
  } else {
    mc_identity <- SuperCell::SCimplify(X = query_matrix,
                                        gamma = gamma,
                                        n.pc = n.pc)
  }

  #get average gene expression for each metacell, in matrix format

  mc_matrix <- SuperCell::supercell_GE(ge = query_matrix,
                                       groups = mc_identity$membership,
                                       mode = "average")


  sce_metacells <- SuperCell::supercell_2_sce(SC.GE = mc_matrix,
                                              SC = mc_identity)

  #add the annotation to SCE, if provided

  if (!is.null(annotation_col)) {
    SummarizedExperiment::colData(sce_metacells)$annotation <- mc_identity$SC.cell.annotation.
  }

  #calculate clusters (needed for running clustifyr)
  graph <- scran::buildSNNGraph(sce_metacells, k=10, use.dimred = 'PCA')
  clust <- igraph::cluster_walktrap(graph)$membership
  SingleCellExperiment::colLabels(sce_metacells) <- factor(clust)

  names(SummarizedExperiment::colData(sce_metacells))[which(names(SummarizedExperiment::colData(sce_metacells))=="label")]="clusters"


  #calcuate tSNE (useful for plotting)
  sce_metacells <- scater::runTSNE(sce_metacells)

  return(sce_metacells)

}
