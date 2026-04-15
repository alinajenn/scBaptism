#' run_scPred
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE object that acts as a reference
#' @param ref_labs Column from the references colData
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#' @param verbose display message after annotation is finished
#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' export
#'
#'@importFrom scPred scPredict
#'@importFrom scPred getFeatureSpace
#'@importFrom scPred trainModel
#'@importFrom scPred get_scpred
#'@importFrom SeuratObject RenameAssays
#'@importFrom Seurat as.Seurat
#'@importFrom Seurat NormalizeData
#'@importFrom Seurat FindVariableFeatures
#'@importFrom Seurat ScaleData
#'@importFrom Seurat RunPCA
#'@importFrom Seurat RunUMAP
#'@importFrom magrittr %>%
#'@importFrom SummarizedExperiment colData
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
#'#sce_annotated <- run_scPred(sce_query = sce_annotated, reference = sce_annotated, ref_labs = "labels_main")
#'
#'# plot the existing annotation with scater(t-SNE)
#'scater::plotTSNE(sce_annotated, color_by = "scb_scPred_labels")
#'
#'
#'@family classical machine learning family
run_scPred <- function(sce_query,
                       reference,
                       ref_labs,
                       threshold = 0.55,
                       max.iter.harmony = 20,
                       recompute_alignment = TRUE,
                       seed = 66,
                       return_extra_info = FALSE,
                       verbose = FALSE,
                       ...)

{

  # checks ----------------------------------------------------------------

  # this tool offers the option to reclassify, give as option?
  # could be done via the plot output

  # transformation --------------------------------------------------------

  # Ref & Query should both be Seurat objects

  seurat_query <- Seurat::as.Seurat(sce_query)
  seurat_ref <- Seurat::as.Seurat(reference)

  # Assay needs to be renamed to data (from originalexp), needed in later function

  seurat_query <- SeuratObject::RenameAssays(seurat_query, assay.name = 'originalexp', new.assay.name = 'data')
  seurat_ref <- SeuratObject::RenameAssays(seurat_ref, assay.name = 'originalexp', new.assay.name = 'data')

  # normalize, dimred for ref
  seurat_ref <- seurat_ref %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA() %>%
    Seurat::RunUMAP(dims = 1:30)

  #train classifier

  seurat_ref <- scPred::getFeatureSpace(seurat_ref, ref_labs)

  seurat_ref <- scPred::trainModel(seurat_ref)

  scPred::get_scpred(seurat_ref)

  # running annotation-----------------------------------------------------

  #normalize query

  seurat_query <- Seurat::NormalizeData(seurat_query)

  #run annotation

  seurat_query <- scPred::scPredict(new = seurat_query,
                                    reference = seurat_ref,
                                    threshold = 0.55,
                                    max.iter.harmony = 20,
                                    recompute_alignment = TRUE,
                                    seed = 66
                                    )

  #include loop around re-annotation?



  # return input SCE with new annotation-----------------------------------

  SummarizedExperiment::colData(sce_query)$scb_scPred_labels <- seurat_query$scpred_prediction

  if(return_extra_info){
    SummarizedExperiment::colData(sce_query)$scb_scPred_max_score <- seurat_query$scpred_max
    SummarizedExperiment::colData(sce_query)$scb_scPred_no_rejection <- seurat_query$scpred_no_rejection
  }


  #message("scPred annotation done")
  if(verbose) message("scPred annotation done")

  return(sce_query)

}
