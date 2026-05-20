#' calculate_metacells
#'
#' @param sce_query SingleCellExperiment object with logcounts to calculate metacells on
#' @param gamma set the graining level for the metacells, aka proportion of single cells to meta cells
#'
#' @importFrom SummarizedExperiment assay ncol colData
#' @importFrom scater plotTSNE runPCA
#' @importFrom SuperCell SCimplify supercell_GE
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scran buildSNNGraph
#' @importFrom igraph cluster_walktrap
#'
#' @returns SCE with logcounts on metacell level, also clusters, PCA & tSNE
#'
#' @export
#'
#' @examples
#' #load example SCE from iUSEiSEE package
#'
#' sce_annotated <- readRDS(
#'  file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' sce_metacells <- calculate_metacells(sce_query = sce_annotated,
#'                                      gamma = 10)
#'
#' sce_metacells
#'
#'
calculate_metacells <- function(sce_query,
                                gamma = 10 #SCimplify graining level
                                ){

  #get the log normalized expression matrix

  query_matrix <- SummarizedExperiment::assay(sce_query, "logcounts")

  #run the metacell calculations, metacell identity for each single cell

  mc_identity <- SuperCell::SCimplify(X = query_matrix,
                                      gamma = gamma)

  #get average gene expression for each metacell, in matrix format

  mc_matrix <- SuperCell::supercell_GE(ge = query_matrix,
                                       groups = mc_identity$membership,
                                       mode = "average")


  #create new SCE with the metacell expression matrix
  sce_metacells <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = mc_matrix))

  #add cellnames (needed for running scClassify)
  colnames(sce_metacells) <- paste0("cell_", seq_len(SummarizedExperiment::ncol(sce_metacells)))

  #run PCA
  sce_metacells <- scater::runPCA(sce_metacells)

  #calculate clusters (needed for running clustifyr)
  graph <- scran::buildSNNGraph(sce_metacells, k=10, use.dimred = 'PCA')
  clust <- igraph::cluster_walktrap(graph)$membership
  SingleCellExperiment::colLabels(sce_metacells) <- factor(clust)

  names(SummarizedExperiment::colData(sce_metacells))[which(names(SummarizedExperiment::colData(sce_metacells))=="label")]="clusters"


  #calcuate tSNE (useful for plotting)
  sce_metacells <- scater::runTSNE(sce_metacells)

  return(sce_metacells)

}
