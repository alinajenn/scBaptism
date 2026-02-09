#' run_SCINA
#'
#' @param sce_query
#' @param markers_list
#'
#'
#' @return sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom SCINA SCINA
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' # load libraries
#'
#'library(iUSEiSEE)
#'library(Seurat)
#'library(dplyr)
#'
#'
#'# load example SCE from iuseisee package
#'
#'sce_annotated <- readRDS(
#'  file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#'# get the markers list
#'myseu <- Seurat::as.Seurat(sce_annotated)
#'myseu <- Seurat::ScaleData(myseu)
#'Seurat::Idents(myseu) <- "labels_main"
#'
#'seu_all_markers <- Seurat::FindAllMarkers(myseu, test.use = "wilcox", only.pos = TRUE,
#'                                          min.pct = 0.25, logfc.threshold = 0.25)
#'
#'top_k_markers <- 50
#'
#'markers_lists <- seu_all_markers %>%
#'  group_by(cluster) %>%
#'  dplyr::top_n(n = top_k_markers, wt = avg_log2FC)
#'markers_lists <- split(markers_lists$gene, markers_lists$cluster)
#'
#'# run the annotation
#'run_SCINA(sce_query = sce_annotated, markers_list =  markers_lists)
#'
#' @family marker family


run_SCINA <- function(sce_query,
                      markers_list,
                      max_iter = 1000,
                      convergence_n = 12,
                      convergence_rate = 0.999,
                      sensitivity_cutoff = 0.9,
                      ...) {



  # 1) TRANSFROM: select logcounts from the SCE
   transf_query <- SummarizedExperiment::assay(sce_query, "logcounts")

  # makers list provided by user

  # 2) RUN
  anno_res <- SCINA::SCINA(exp = transf_query,
                           signatures = markers_list,
                           max_iter,
                           convergence_n,
                           convergence_rate,
                           sensitivity_cutoff,
                           ...)

  # 3) RETURN: add SCINA annotation & probabilities column to SCE
  colData(sce_query)$scb_SCINA_res <- anno_res$cell_labels

  #probablities are calculated by SCINA
  colData(sce_query)$scb_SCINA_prob <- anno_res$probabilities

  return(sce_query)

}


