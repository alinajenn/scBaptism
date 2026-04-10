#' run_scmap
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE object that acts as a reference
#' @param ref_labs Column from the references colData
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#' @param verbose display message after annotation is finished
#'
#'
#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom scmap scmapCluster
#' @importFrom scmap selectFeatures
#' @importFrom scmap indexCell
#' @importFrom scmap indexCluster
#' @importFrom scmap scmapCell2Cluster
#' @importFrom scmap scmapCell
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @examples
#'
#' library(iUSEiSEE)
#' library(Seurat)
#' library(dplyr)
#'
#' #load SCE from iuseisee
#'
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' #run the annotation
#' #sce_annotated <- run_scmap(sce_query= sce_annotated, reference = sce_annotated, ref_labs = "labels_main)
#'
#' #plot the new annotations with scater(t-SNE)
#' #cell wise annotation
#' scater::plotTSNE(sce_annotated, color_by = "scb_scmap_labels")
#'
#' #cluster wise annotation
#'scater::plotTSNE(sce_annotated, color_by = "scb_scmap_clusterlabels")
#'
#'@family classical machine learning family
run_scmap <- function(sce_query,
                      reference, #SCE
                      ref_labs,
                      return_extra_info = FALSE,
                      verbose = FALSE,
                      ...)

{

  # checks ----------------------------------------------------------------


  # transformation --------------------------------------------------------

  #input can be SCE for reference & query



  #prepare reference
  #use gene names as feature symbols
  SummarizedExperiment::rowData(reference)$feature_symbol <- rownames(reference)

  #if needed, remove duplicated features
  reference[!duplicated(SummarizedExperiment::rownames(reference)), ]

  #select features & perform clustering
  reference <- scmap::selectFeatures(reference, suppress_plot = TRUE)

  reference <- scmap::indexCluster(reference, cluster_col = ref_labs)

  #prepare query (features as rownames)
  SummarizedExperiment::rowData(sce_query)$feature_symbol <- rownames(sce_query)



  # running annotations-----------------------------------------------------

  #cluster-wise
  scmapCluster_results <- scmap::scmapCluster(
    projection = sce_query,
    index_list = list(
      refinfo = metadata(reference)$scmap_cluster_index
    )
  )


  #cell-wise

  #index reference cell wise

  reference <- scmap::indexCell(reference)

  scmapCell_results <- scmap::scmapCell(
    sce_query,
    list(
      refinfo = metadata(reference)$scmap_cell_index
    )
  )

  scmapCell_clusters <- scmap::scmapCell2Cluster(
    scmapCell_results,
    list(
      as.character(SummarizedExperiment::colData(reference)[[ref_labs]])
    )
  )


  # return input SCE with new annotation-----------------------------------

  #cluster-wise
  SummarizedExperiment::colData(sce_query)[["scb_scmap_clusterlabels"]] <- scmapCluster_results$scmap_cluster_labs[,'refinfo']

  #cell-wise
  SummarizedExperiment::colData(sce_query)[["scb_scmap_labels"]] <- scmapCell_clusters$scmap_cluster_labs[,'refinfo']


#  if(return_extra_info){
 #   SummarizedExperiment::colData(sce_query)$scb_scmap_cluster_siml <- scmapCluster_results$scmap_cluster_siml
  #  SummarizedExperiment::colData(sce_query)$scb_scmap_combined_labels <- scmapCluster_results$combined_labs

   # SummarizedExperiment::colData(sce_query)$scb_scmap_combined_labels <- scmapCell_clusters$combined_labs
    #SummarizedExperiment::colData(sce_query)$scb_scmap_combined_labels <- scmapCell_clusters$combined_labs
  #}

  #message
  if(verbose) message("scmap annotation done")

  return(sce_query)

}
