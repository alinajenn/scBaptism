#' run_Seurat
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE with the reference for the annotation
#' @param ref_labs of the column of the reference, that contains annotation information
#' @param n_pcs number of PCs to compute on reference
#' @param query.assay name of the Assay to use from query
#' @param weight.reduction Dimensional reduction to use for the weighting anchors. Options are:pcaproject, lsiproject, pca, cca, custom DimReduc
#' @param l2.norm Perform L2 normalization on the cell embeddings after dimensional reduction
#' @param Seurat_dims Set of dimensions to use in the anchor weighting procedure. If NULL, the same dimensions that were used to find anchors will be used for weighting.
#' @param k.weight Number of neighbors to consider when weighting anchors
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param n.trees More trees gives higher precision when using annoy approximate nearest neighbor search
#' @param slot Slot to store the imputed data. Must be either "data" (default) or "counts"
#' @param prediction.assay Return an Assay object with the prediction scores for each class stored in the data slot.
#' @param store.weights Optionally store the weights matrix used for predictions in the returned query object.
#' @param verbose display message after annotation is finished
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#'
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom Seurat as.Seurat FindVariableFeatures TransferData ScaleData Idents
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment colData
#'
#'
#' @examples
#'
#' #loading example dataset
#' library(iUSEiSEE)
#'
#' sce_annotated <-
#'   readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#'
#' sce_annotated <- run_Seurat(sce_query = sce_annotated,
#'                             reference = sce_annotated,
#'                             ref_labs = "labels_main")
#'
#' #tSNE plot of the result
#' scater::plotTSNE(sce_annotated, color_by = "scb_Seurat_labels")
#'
#'
#'@family reference-based family
run_Seurat <- function(sce_query,
                       reference,
                       ref_labs,
                       query.assay = NULL,
                       weight.reduction = "pcaproject",
                       l2.norm = FALSE,
                       Seurat_dims = NULL,
                       k.weight = 50,
                       sd.weight = 1,
                       eps = 0,
                       n.trees = 50,
                       slot = "data",
                       prediction.assay = FALSE,
                       store.weights = TRUE,
                       n_pcs = 30,
                       return_extra_info = FALSE,
                       verbose = FALSE
                       )

{

  #transform input into Seurat objects-------------------------------------
  seurat_query <- Seurat::as.Seurat(sce_query)
  seurat_ref <- Seurat::as.Seurat(reference)

  # scale data
  seurat_query <- Seurat::ScaleData(seurat_query)
  seurat_ref <- Seurat::ScaleData(seurat_ref)

  # set Idents (aka existing annotation from the reference)
  Seurat::Idents(seurat_ref) <- colData(reference)[[ref_labs]]

  # run annotation---------------------------------------------------------
  seurat_query <- Seurat::FindVariableFeatures(object = seurat_query)
  seurat_ref <- Seurat::FindVariableFeatures(object = seurat_ref)

  anchors <- Seurat::FindTransferAnchors(reference = seurat_ref,
                                         query = seurat_query,
                                         dims = 1:n_pcs)

  seurat_res <- Seurat::TransferData(anchorset = anchors,
                                     refdata = Seurat::Idents(seurat_ref),
                                     reference = seurat_ref,
                                     query = seurat_query,
                                     query.assay,
                                     dims = Seurat_dims,
                                     weight.reduction = weight.reduction,
                                     l2.norm = l2.norm,
                                     k.weight = k.weight,
                                     sd.weight = sd.weight,
                                     eps = eps,
                                     n.trees = n.trees,
                                     slot = slot,
                                     prediction.assay = prediction.assay,
                                     store.weights = store.weights
                                     )


  # return input SCE with new annotation-----------------------------------
  SummarizedExperiment::colData(sce_query)$scb_Seurat_labels <- seurat_res$predicted.id


  if (return_extra_info){
    SummarizedExperiment::colData(sce_query)$scb_Seurat_score_max <- seurat_res$prediction.score.max
  }

  if(verbose) message("Seurat annotation done")

  return(sce_query)

}
