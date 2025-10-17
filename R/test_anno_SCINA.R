#' Title
#'
#' @param query
#' @param reference
#' @param ...
#'
#' @returns A SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom SCINA SCINA
#' @importFrom SingleCellExperiment colData
#'
#' @examples
#'
#' # TODO
run_SCINA <- function(sce_query, reference, ...) {  #... for the default values

  #"transfrom" (aka select) query from SCE into a normalized expression matrix
  # just call counts? Usually they should be called counts (or normcounts, logcounts), but maybe check for it?
  # normalize counts if not available?

  #transform reference into list containing signature vectors

  transf_query <- sce_query # with some sauce on it

  anno_res <- SCINA::SCINA(transf_query, referlibrence, ...) #this way or load complete package?
  #import should be handled by namespace?


  #add output (vector with celltype predictions) to input SCE
  colData(query)$SCINA_res <- anno_res$cell_labels
  #probabilities are also available as output


  return(sce_out)

}
