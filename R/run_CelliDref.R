#' run_CelliDref
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE with annotations
#' @param ref_labs existing annotation in the reference
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
#' @importFrom CelliD RunCellHGT GetCellGeneSet RunMCA
#' @importFrom SummarizedExperiment colData
#' @importFrom Seurat as.Seurat ScaleData FindVariableFeatures NormalizeData
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
#' sce_annotated <- run_CelliDref(sce_query = sce_annotated,
#'                                reference = sce_annotated,
#'                                ref_labs = "labels_main")
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_CelliDref_labels")
#'
#'
#'@family hybrid family
run_CelliDref <- function(sce_query,
                       reference,
                       ref_labs,
                       reduction = "mca",
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
  #sce input would be fine, but some of the later functions depend on Seurat object

  #transform query and reference into Seurat objects

  seurat_query <- Seurat::as.Seurat(sce_query)
  seurat_ref <- Seurat::as.Seurat(reference)

  # running annotation-----------------------------------------------------

  # prepare query

  seurat_query <- Seurat::NormalizeData(seurat_query)
  seurat_query <- Seurat::FindVariableFeatures(seurat_query)
  seurat_query <- Seurat::ScaleData(seurat_query)
  seurat_query <- CelliD::RunMCA(seurat_query)


  #prepare ref

  seurat_ref <- Seurat::NormalizeData(seurat_ref)
  seurat_ref <- Seurat::ScaleData(seurat_ref, features = rownames(seurat_ref))
  seurat_ref <- CelliD::RunMCA(seurat_ref)

  #extract per cell gene signatures

  ref_cell_gs <- CelliD::GetCellGeneSet(seurat_ref, dims = 1:50, n.features = 200)


  # run reference based annotation of CelliD

  result_ref <- CelliD::RunCellHGT(X = seurat_query,
                                  pathways = ref_cell_gs,
                                  reduction = reduction,
                                  n.features = n.features,
                                  features = CelliD_features,
                                  dims = CelliD_dims,
                                  minSize = minSize,
                                  log.trans = log.trans,
                                  p.adjust = p.adjust
                                  )



  # return input SCE with new annotation-----------------------------------

  # for each cell, assess the signature with the lowest corrected p-value (max -log10 corrected p-value)
  ref_match <- rownames(result_ref)[apply(result_ref, 2, which.max)]

  ref_prediction <- seurat_ref[[]][ref_match, ref_labs]


  # for each cell, evaluate if the lowest p-value is significant
  ref_prediction_signif <- ifelse(apply(result_ref, 2, max)>2, yes = ref_prediction, "unassigned")

  # add new annotation to SCE query
  SummarizedExperiment::colData(sce_query)$scb_CelliDref_labels <- ref_prediction_signif


  # message("CelliD annotation done")
  if(verbose) message("CelliD annotation done")

  return(sce_query)

}
