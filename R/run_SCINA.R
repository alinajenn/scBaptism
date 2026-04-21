#' run_SCINA
#'
#' @param sce_query SCE to be annotated
#' @param markers_list List of marker genes
#' @param return_extra_info if TRUE, adds additional metadata from the annotation
#' @param verbose display message after annotation is finished
#'
#'
#' @returns sce_query: a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom SCINA SCINA
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' library(iUSEiSEE)
#' library(Seurat)
#' library(dplyr)
#'
#'
#' #load example SCE from iUSEiSEE package
#'
#' sce_annotated <- readRDS(
#'  file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
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
#' # run the annotation
#' run_SCINA(sce_query = sce_annotated, markers_list =  markers_lists)
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_SCINA_labels")
#'
#' @family marker family


run_SCINA <- function(sce_query,
                      markers_list,
                      return_extra_info = FALSE,
                      verbose = FALSE,
                      max_iter = 1000,
                      convergence_n = 12,
                      convergence_rate = 0.999,
                      sensitivity_cutoff = 0.9,
                      ...) {

  # checks ----------------------------------------------------------------

  if(!hasArg(sce_query)) {
    stop("please provide a query")
  }

  #do logcounts exist and are called logcounts?

  if(!("logcounts" %in% names(assays(sce_query)))) {
    stop("no assay called 'logcounts' found. Rename your assay or calculate the logcounts.")
  }

  # transformation --------------------------------------------------------

  #select log-normalized counts assay from input SCE
   transf_query <- SummarizedExperiment::assay(sce_query, "logcounts")

  # markers list needs to be provided by user

  # run annotation --------------------------------------------------------
  anno_res <- SCINA::SCINA(exp = transf_query,
                           signatures = markers_list,
                           max_iter,
                           convergence_n,
                           convergence_rate,
                           sensitivity_cutoff,
                           ...)

  # return annotation ------------------------------------------------------

  #add SCINA annotation & probabilities column to SCE
  SummarizedExperiment::colData(sce_query)$scb_SCINA_labels <- anno_res$cell_labels

  #probablities are calculated by SCINA
  if(return_extra_info) {
    SummarizedExperiment::colData(sce_query)$scb_SCINA_prob <- anno_res$probabilities}

  if(verbose) message("SCINA annotation done")
  return(sce_query)

}


