#' run_CelliDref
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE with annotations
#' @param ref_lab existing annotation in the reference
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#' @param verbose display message after annotation is finished
#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' export
#'
#'@importFrom CelliD RunCellHGT
#'@importFrom SummarizedExperiment colData
#'@importFrom CelliD RunMCA
#'
#'
#' @examples
#'
#'
#'library(iUSEiSEE)
#'library(dplyr)
#'
#'# load SCE from iUSEiSEE
#'
#'sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#'#run the annotation
#'#sce_annotated <- run_tool(sce_annotated, markers_lists)
#'
#'# plot the existing annotation with scater(t-SNE)
#'scater::plotTSNE(sce_annotated, color_by = "scb_CelliDref_labels")
#'
#'
#'@family hybrid family
run_CelliDref <- function(sce_query,
                       reference,
                       ref_lab,
                       return_extra_info = FALSE,
                       verbose = FALSE,
                       ...)

{

  # checks ----------------------------------------------------------------


  # transformation --------------------------------------------------------
  #sce input would be fine, but some of the later functions depend on Seurat object

  #reference:

  seurat_query <- Seurat::as.Seurat(sce_query)
  seurat_ref <- Seurat::as.Seurat(reference)

  # running annotation-----------------------------------------------------

  #Restricting to protein-coding genes: Skipped for now
  #Normalize & Scale Data with Seurat functions? Skipped for now
  #RUN MCA: necessary for using the tool


  seurat_query <- Seurat::NormalizeData(seurat_query)
  seurat_query <- Seurat::FindVariableFeatures(seurat_query)
  seurat_query <- Seurat::ScaleData(seurat_query)
  seurat_query <- CelliD::RunMCA(seurat_query)

  #prepare

  #runs other dim reds (that we already have)


  #prepare ref
  #extract per cell gene signatures

  seurat_ref <- Seurat::NormalizeData(seurat_ref)
  seurat_ref <- Seurat::ScaleData(seurat_ref, features = rownames(seurat_ref))
  seurat_ref <- CelliD::RunMCA(seurat_ref)

  ref_cell_gs <- CelliD::GetCellGeneSet(seurat_ref, dims = 1:50, n.features = 200)


  #1 run reference version

  result_ref <- CelliD::RunCellHGT(seurat_query,
                                      pathways = ref_cell_gs,
                                      dims = 1:50,
                                      n.features = 200
  )



  # return input SCE with new annotation-----------------------------------



  #for each cell, assess the signature with the lowest corrected p-value (max -log10 corrected p-value)
  ref_match <- rownames(result_ref)[apply(result_ref, 2, which.max)]

  ref_prediction <- seurat_ref[[]][, ref_lab][ref_match]


  # for each cell, evaluate if the lowest p-value is significant
  ref_prediction_signif <- ifelse(apply(result_ref, 2, max)>2, yes = ref_prediction, "unassigned")

  # add new annotation to SCE query
  SummarizedExperiment::colData(sce_query)$scb_CelliDref_labels <- ref_prediction_signif


  #message("CelliD annotation done")
  if(verbose) message("CelliD annotation done")

  return(sce_query)

}
