#' scBaptism_majority
#'
#' @param input_vector vector to do the majority voting on
#'
#'
#' @returns sce_query: SCE object with an added metadata column for the majority annotation
#'
#' @export
#'
#'
#'
#' @examples
#'
#' #loading example dataset
#' library(iUSEiSEE)
#'
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'

#to-do handle ties
get_winner <- function(input_vector) {
  rank_table <- table(input_vector)
  winner <- max(rank_table)
  winner_name <- names(rank_table[rank_table == winner])

  return(winner_name)
}

