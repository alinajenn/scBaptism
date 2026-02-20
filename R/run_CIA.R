#' run_CIA
#'
#' @param sce_query SCE to be annotated
#' @param markers_list List of marker genes
#' @param similarity_theshold threshold wether the highest score is significantly higher than others
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
#'
#' # load SCE from iuseisee
#'
#' sce_annotated <- readRDS(
#'   file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE")
#' )
#'
#' sce_query <- run_SingleR(sce_query = sce_annotated, markers_list = sce_annotated$labels_main)
#'
#' # plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_CIA_res")
#'
#'
#'@family marker family
run_CIA <- function(sce_query,
                    markers_list, #markers input (list)
                    similarity_threshold = 0, #is highest score significantly higher? If not then cell=unassigned
                    column_name = "scb_CIA_res", #name of metadata column, default = CIA_Prediction
                    n_cpus = 1, # number of cpu cores, default
                    verbose = FALSE)

{

  # checks ----------------------------------------------------



  # transformation --------------------------------------------

  #input as SCE is accepted

  #signatures_input should be "a list where each element is a vector of gene names with names of the list elements representing signature names"




  # running annotation -----------------------------------------


  sce_query <- CIA::CIA_classify(data = sce_query,
                                 signatures_input = markers_list,
                                 similarity_threshold = similarity_threshold,
                                 column_name = column_name,
                                 n_cpus = n_cpus)


  # return sce with new annotation -----------------------------

  # no input needed, adding annotation to SCE is done by CIA_Classify already

  if (verbose) message("CIA annotation done")
  return(sce_query)

}
