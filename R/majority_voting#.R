#' .get_winner
#'
#' @param input_vector vector to do the majority voting on (in our case a row from the annotation data)
#' @param tie_breaker choose what happens with a tie between the annotations, options are "na" and "first", default is "na"
#'
#' @returns max_scorer the name of the most occurring string from the input vector
#'
#'
#' @noRd
#'
#'
.get_winner <- function(input_row, tie_breaker = "na") {

  #rank all values in that row & get highest rank
  rank_table <- table(input_row)
  winner <- max(rank_table)

  #get all strings that received the maximum score
  max_scorer <- names(rank_table[rank_table == winner])

  # resolve ties according to user choice


  if (length(max_scorer) == 1) {
    return(max_scorer)
  } else {
    switch(tie_breaker,
           first  = max_scorer[1],
           na     = NA_character_
           #concat = paste(sort(max_scorer), collapse = ";")  # patch all labels together
    )
  }

  return(max_scorer)
}



#' majority_vote
#'
#' @param sce_query SingleCellExperiment object with multiple annotations
#' @param anno_columns Vector with names of annotation columns, if NULL all columns starting with "scb" will be used
#' @param tie_breaker choose what happens with a tie between the annotations, options are "na" and "first", default is "na"
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr mutate
#'
#' @returns sce_query: SCE object with an added metadata column for the majority annotation
#'
#' @export
#'
#' @examples
#'
#' # loading example dataset, will also serve as a reference
#' library(iUSEiSEE)
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' # use scBaptism to get several single cell annotations
#' sce_annotated <- run_scBaptism(sce_query = sce_annotated, anno_methods = c("clustifyr", "scPred", "scClassify"), reference = sce_annotate, ref_labs = "labels_main")
#'
#' # perform a majority vote on all annotations obtained through scBaptism, with choosing the first element in case of a tie
#' sce_annotated <- majority_vote(sce_query = sce_annotated, anno_columns = NULL, tie_breaker = "first")
#'
#' # plot the majority vote annotation
#'
#'
#'
majority_vote <- function(sce_query, anno_columns = NULL, tie_breaker = "na") {

  #get the names of all annotation columns, if user did not provide a vector

  if(is.null(anno_columns)) {

  #names of all colData
  col_names <- names(SummarizedExperiment::colData(sce_query))

  #get all annotation columns starting with "scb" (scBaptism)
  anno_cols_pattern <- as.character("^scb")

  anno_columns <- grep(anno_cols_pattern, col_names, ignore.case = TRUE, value = TRUE)
  }

  #initialize data.frame with first entry

  sce_df <- as.data.frame(SummarizedExperiment::colData(sce_query)[[anno_columns[1]]])

  #delete first column from anno_columns, so it does not get added twice
  anno_columns_rest <- anno_columns[-1]

  #get all annotation from colData as a dataframe

  for (item in anno_columns_rest) {
    sce_df <- cbind(sce_df, as.data.frame(SummarizedExperiment::colData(sce_query)[[item]]))
  }

  #name the columns
  colnames(sce_df) <- anno_columns


  #perform the majority vote
  sce_df <- sce_df %>%
    dplyr::mutate(
      majority_labels = apply(sce_df, 1,  # 1 means we apply to rows
                            .get_winner, tie_breaker = tie_breaker)
    )

                          #all_of() selects the columns with this name (can be a vector of strings)
                          # select(.) means we use the data that is being piped into it
  #since the resulting column is now a list, convert it into a character vector

  sce_df$majority_labels <- vapply(sce_df$majority_labels, `[`, character(1), 1)

  # write result back into the SCE object

  SummarizedExperiment::colData(sce_query)$scb_majority_labels <- sce_df$majority_labels

  return(sce_query)

}
