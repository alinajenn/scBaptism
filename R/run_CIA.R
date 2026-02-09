#' run_CIA
#'
#' @param sce_query
#' @param markers_list
#'
#'
#' @returns sce_query : a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom CIA_classify CIA
#' @importFrom
#'
#' @examples
#'
#'@family marker family
run_CIA <- function(sce_query,
                    markers_list, #markers input (list)
                    similarity_threshold = 0, #is highest score significantly higher? If not then cell=unassigned
                    column_name = "scb_CIA_res", #name of metadata column, default = CIA_Prediction
                    n_cpus = 1) # number of cpu cores, default

{

  # checks ----------------------------------------------------



  # transformation ---------------------------------------------

  #input as SCE is accepted

  #signatures_input should be "a list where each element is a vector of gene names with names of the list elements representing signature names"




  # running annotation ------------------------------------------


  sce_query <- CIA::CIA_classify(data = sce_query,
                                 signatures_input = markers_list,
                                 similarity_threshold = similarity_threshold,
                                 column_name = column_name,
                                 n_cpus = n_cpus)   #TODO this parameter needs to be the same everywhere (here, up & in call)



  #score_all_signatures: leave out ncups parameter, otherwise does not work (only about paralell computation)

  # return sce with new annotation ------------------------------

  #no input from me needed, adding annotation to SCE is done by CIA_Classify already



}
