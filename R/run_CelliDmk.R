#' run_CelliDmk
#'
#' @param sce_query SCE to be annotated
#' @param markers_list List of marker genes
#' @param reduction name of the MCA reduction
#' @param n.features integer of top n features to consider for hypergeometric test
#' @param CelliD_features vector of features to calculate the gene ranking by default will take everything in the selected mca reduction.
#' @param CelliD_dims MCA dimensions to use to compute n.features top genes.
#' @param minSize minimum number of overlapping genes in geneset and
#' @param log.trans if TRUE tranform the pvalue matrix with -log10 and convert it to sparse matrix
#' @param p.adjust if TRUE apply Benjamini Hochberg correction to p-value
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#' @param verbose display message after annotation is finished

#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#'@importFrom CelliD RunCellHGT
#'@importFrom CelliD RunMCA
#'@importFrom SingleCellExperiment colData
#'
#' @examples
#'
#'
#' library(iUSEiSEE)
#' library(dplyr)
#' library(Seurat)
#'
#' # load SCE from iUSEiSEE
#'
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
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
#'
#' # run the annotation
#' sce_annotated <- run_CelliDmk(sce_query = sce_annotated,
#'                               markers_list = markers_lists)
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_CelliDmk_labels")
#'
#'
#'@family hybrid family
run_CelliDmk <- function(sce_query,
                         markers_list,
                         reduction = "MCA",
                         n.features = 200,
                         CelliD_features = NULL,
                         CelliD_dims = seq(50),
                         minSize = 10,
                         log.trans = TRUE,
                         p.adjust = TRUE,
                         return_extra_info = FALSE,
                         verbose = FALSE
                         )

{

  # checks ----------------------------------------------------------------


  # transformation --------------------------------------------------------
  #sce input is fine

  #input markers: list of vectors (characters)


  # running annotation-----------------------------------------------------

  #dimansionality reduction MCA is neccessary for running CelliD
  sce_query <- CelliD::RunMCA(sce_query)


  # run marker based CelliD version

  result_marker <- CelliD::RunCellHGT(sce_query,
                              pathways = markers_list,
                              reduction = reduction,
                              n.features = n.features,
                              features = CelliD_features,
                              dims = CelliD_dims,
                              minSize = minSize,
                              log.trans = log.trans,
                              p.adjust = p.adjust
                              )


  # return input SCE with new annotation-----------------------------------

  #for each cell, assess the signature with the lowest corrected p-value (max -log10 corrected p-value)
  marker_prediction <- rownames(result_marker)[apply(result_marker, 2, which.max)]

  # for each cell, evaluate if the lowest p-value is significant
  marker_prediction_signif <- ifelse(apply(result_marker, 2, max)>2, yes = marker_prediction, "unassigned")

  # add new annotation to SCE query
  SummarizedExperiment::colData(sce_query)$scb_CelliDmk_labels <- marker_prediction_signif



  #message("CelliD annotation done")
  if(verbose) message("CelliD (marker based) annotation done")

  return(sce_query)

}
