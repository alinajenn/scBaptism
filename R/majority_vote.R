#' .get_winner
#'
#' @param input_vector vector to do the majority voting on (in our case a row from the annotation data)
#' @param tie_breaker choose what happens with a tie between the annotations, options are "na", "concat", "LCA" and "first", default is "concat"
#' @param cl_graph graph mapping out cell relationships, needed when user selects "LCA" as tie breaker
#'
#' @importFrom ontoProc getOnto
#' @importFrom igraph make_graph
#'
#' @returns max_scorer the name of the most occurring string from the input vector
#'
#'
#' @noRd
#'
#'
.get_winner <- function(input_row, tie_breaker = "concat", cl_graph = NULL) {

  #rank all values in that row & get highest rank
  rank_table <- table(input_row)
  winner <- max(rank_table)

  #get all strings that received the maximum score
  max_scorer <- names(rank_table[rank_table == winner])

  # resolve ties according to user choice


  if (length(max_scorer) == 1) {
    return(max_scorer)
  } else {
    max_scorer <- switch(EXPR = tie_breaker,
           'first' = max_scorer[1],
           'concat' = paste(sort(max_scorer), collapse = "|"),
           'LCA' = .tiebreaker_LCA(max_scorer, cl_graph),
           'na' = NA_character_

           #sort sorts the vector in a specific order
    )
    return(max_scorer)
  }

}


#' .tiebreaker_LCA
#'
#' @param candidates Input into the
#' @param graph graph mapping out the relationships of the cells
#'
#' @importFrom rols OlsSearch olsSearch
#' @importFrom ontoProc findCommonAncestors
#' @importFrom igraph V
#' @importFrom stats na.omit
#'
#' @returns sce_query: SCE object with an added metadata column for the majority annotation
#'
#' @noRd
#'
#'
#'
.tiebreaker_LCA <- function(candidates, graph) {
  suppressWarnings({
    if (is.null(candidates) || !is.character(candidates) || length(candidates) == 0) {
      return(NULL)
    }

    cand_cl <- c()

    for (cand in candidates) {
      if (is.na(cand) || nchar(trimws(cand)) == 0) next

      tryCatch({
        res  <- OlsSearch(cand, ontology = "CL")
        hits <- olsSearch(res, all = FALSE)
        obo_ids <- hits@response[["obo_id"]]
        obo_ids <- obo_ids[grepl("^CL[:_]", obo_ids)]
        if (length(obo_ids) == 0) next
        cand_cl <- c(cand_cl, obo_ids[1])
      }, error = function(e) NULL)
    }

    cand_cl <- unique(na.omit(cand_cl))
    cand_cl <- cand_cl[cand_cl %in% V(graph)$name]

    if (length(cand_cl) < 2) return(NULL)

    lca <- tryCatch(
      findCommonAncestors(cand_cl, g = graph),
      error = function(e) NULL
    )

    if (is.null(lca) || length(lca@rownames) == 0 || length(lca@rownames[[1]]) == 0) {
      return(NULL)
    }

    lca_id <- lca@rownames[[1]]
    if (is.na(lca_id) || nchar(trimws(lca_id)) == 0) return(NULL)

    lca_label <- tryCatch({
      lca_res <- OlsSearch(lca_id, ontology = "CL")
      lca_hit <- olsSearch(lca_res, all = FALSE)
      lca_hit@response[["label"]][[1]]
    }, error = function(e) NULL)

    return(lca_label)
  })
}




#' majority_vote
#'
#' @param sce_query SingleCellExperiment object with multiple annotations
#' @param anno_columns Vector with names of annotation columns, if NULL all columns starting with "scb" will be used
#' @param tie_breaker choose what happens with a tie between the annotations, options are "concat", "LCA", "na" and "first"
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
#' sce_annotated <- readRDS(file = system.file("datasets",
#'                                             "sce_pbmc3k.RDS",
#'                                             package = "iUSEiSEE"))
#'
#' # use scBaptism to get several single cell annotations
#' sce_annotated <- run_scBaptism(sce_query = sce_annotated,
#'                                anno_methods = c("CelliDref", "SingleR", "scClassify"),
#'                                reference = sce_annotated,
#'                                ref_labs = "labels_main")
#'
#' # perform a majority vote on all annotations obtained through scBaptism
#' # with choosing the first annotation in case of a tie
#' sce_annotated <- majority_vote(sce_query = sce_annotated,
#'                                anno_columns = NULL,
#'                                tie_breaker = "first")
#'
#' # plot the majority vote annotation
#' scater::plotTSNE(sce_annotated, color_by = "scb_majority_labels")
#'
#'
#'
majority_vote <- function(sce_query, anno_columns = NULL, tie_breaker) {

  if(tie_breaker == "LCA") {
    cl_g <- ontoProc::getOnto()

    parents <- cl_g$parents
    self <- rep(names(parents), lengths(parents))

    cl_graph <- igraph::make_graph(rbind(unlist(parents), self))
  }

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
      majority_labels = apply(sce_df, 1, .get_winner, tie_breaker = tie_breaker, cl_graph = cl_graph)
      )


  # transform the majority annotation column from list to character vector
  # this is useful for visualizations later
  sce_df$majority_labels <- vapply(sce_df$majority_labels, `[`, character(1), 1)

  # write result back into the SCE object

  SummarizedExperiment::colData(sce_query)$scb_majority_labels <- sce_df$majority_labels

  return(sce_query)

}
