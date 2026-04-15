#' run_Seurat
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE with the reference for the annotation
#' @param ref_labs of the column of the reference, that contains annotation information
#' @param n_pcs number of PCs to compute on reference
#' @param verbose display message after annotation is finished
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#'
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'

#' @importFrom Seurat as.Seurat
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat FindTransferAnchors
#' @importFrom Seurat TransferData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat Idents
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment colData
#'
#'
#' @examples
#'
#' #loading example dataset
#' library(iUSEiSEE)
#'
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#'
#' sce_annotated <- run_Seurat(sce_query = sce_annotated, reference = sce_annotated, ref_labs = "labels_main")
#'
#' #tSNE plot of the result
#' scater::plotTSNE(sce_annotated, color_by = "scb_Seurat_labels")
#'
#'
#'@family reference-based family
run_Seurat <- function(sce_query, #SCE
                       reference, #SCE with labels
                       ref_labs,# string name of annotation column of SCE
                       n_pcs = 30, # 30 is the default from Seurat
                       return_extra_info = FALSE,
                       verbose = FALSE,
                       ...)

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
                                     dims = 1:n_pcs)


  # return input SCE with new annotation-----------------------------------
  SummarizedExperiment::colData(sce_query)$scb_Seurat_labels <- seurat_res$predicted.id


  if (return_extra_info){
    SummarizedExperiment::colData(sce_query)$scb_Seurat_score_max <- seurat_res$prediction.score.max
  }
  #also returns predictions.score.celltype, worth adding them? (and how? possibly AddMetaData())

  if(verbose) message("Seurat annotation done")

  return(sce_query)

}
