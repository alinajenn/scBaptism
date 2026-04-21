#' run_CIA
#'
#' @param sce_query SCE to be annotated
#' @param markers_list List of marker genes
#' @param similarity_threshold threshold whether the highest score is significantly higher than others
#' @param column_name name of the column that contains the final annotation
#' @param n_cpus number of cpu cores used
#' @param verbose display message after annotation is finished
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom CIA CIA_classify
#'
#' @examples
#'
#' library(iUSEiSEE)
#' library(Seurat)
#' library(dplyr)
#'
#' # load SCE from iuseisee
#'
#' sce_annotated <- readRDS(
#'   file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE")
#' )
#'
#' # get the markers list using Seurat
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
#' sce_annotated <- run_CIA(sce_query = sce_annotated, markers_list = markers_lists)
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_CIA_labels")
#'
#'
#'@family marker family
run_CIA <- function(sce_query,
                    markers_list, #markers input (list)
                    similarity_threshold = 0, #is highest score significantly higher? If not then cell=unassigned
                    n_cpus = 1, # number of cpu cores, default
                    verbose = FALSE)

{

  # checks ----------------------------------------------------



  # transformation --------------------------------------------

  #input as SCE is accepted

  # running annotation -----------------------------------------

  #hard-code column name so it will fit with the scBaptism naming conventions
  column_name = "scb_CIA_labels"


  sce_query <- CIA::CIA_classify(data = sce_query,
                                 signatures_input = markers_list,
                                 similarity_threshold = similarity_threshold,
                                 column_name = column_name,
                                 n_cpus = n_cpus)


  # return sce with new annotation -----------------------------

  #adding annotation column to SCE is done by CIA_classify already


  if (verbose) message("CIA annotation done")
  return(sce_query)

}
