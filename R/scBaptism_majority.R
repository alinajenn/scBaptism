#' get_winner
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
get_winner <- function(input_row, tie_breaker = "first") {

  #rank all values in that row & get highest rank
  rank_table <- table(input_row)
  winner <- max(rank_table)

  #get all strings that received the maximum score
  max_scorer <- names(rank_table[rank_table == winner])

  # resolve ties according to user choice
  tie_breaker <- match.arg(tie_breaker)

  if (length(max_scorer) == 1) {
    return(max_scorer)                     # clear winner
  } else {
    switch(tie_breaker,
           first  = max_scorer[1],          # deterministic – first in alphabetical order
           na     = NA_character_,    # mark as NA if there is a tie
           concat = paste(sort(max_scorer), collapse = ";")  # patch all labels together
    )
  }


  ###TODO are tie-breakers choosen by user or do we hard-code one?
  # should different kind of tie breakers be handeled differently

  return(max_scorer)
}


#tie breaker percentage function here


  return(ties_percentage)
}

majority_vote <- function(sce_query, anno_columns = NULL, tie_breaker) {

  #get the names of all annotation columns, if user did not provide a vector

  if(is.null(anno_columns)) {

  #names of all colData
  col_names <- names(colData(sce_query))

  #get all annotation columns starting with "scb" (scBaptism)
  anno_cols_pattern <- as.character("^scb")

  anno_columns <- grep(anno_cols_pattern, col_names, ignore.case = TRUE, value = TRUE)
  }

  #initialize data.frame with first entry

  #######################columns need to be properly named
  sce_df <- as.data.frame(SummarizedExperiment::colData(sce_query)[[anno_columns[1]]])
  #delete first column from anno_columns, so it does not get added twice
  #sce_df <- sce_df[, -c(1)]
  anno_columns <- anno_columns[-1]

  #get all annotation from colData as a dataframe

  for (item in anno_columns) {
    sce_df <- cbind(sce_df, as.data.frame(colData(sce_query)[[item]]))
  }


  #perform the majority vote
  sce_df <- sce_df %>%
    dplyr::mutate(
      majority_labels = apply(sce_df, 1,  #1 means we apply to rows
                            get_winner, tie_breaker = "first")
    )

                          #all_of() selects the columns with this name (can be a vector of strings)
                          # select(.) means we use the data that is being piped into it

  # write result back into the SCE object

  SummarizedExperiment::colData(sce_query)$scb_majority_labels <- sce_df$majority_labels

  return(sce_query)

}
