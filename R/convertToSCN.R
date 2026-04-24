#' convertToSCN
#'
#' extract sampTab and expDat sce object into regular S3 objects
#' @param sce_object SCE object to be converted to SCN format
#' @param exp_type type of assay data in the SCE (eg "counts", "logcounts", "normcounts" etc)
#'
#' @importFrom SingleCellExperiment counts logcounts normcounts
#' colData
#'
#' @export
#'
#' @return list with sampTab (the metadata) and expDat (the expression matrix),
#' the names of the rows (aka cells) will be stored in the sample_name column
convertToSCN <- function(sce_object, exp_type = "counts"){
  #extract metadata
  sampTab = as.data.frame(colData(sce_object, internal = TRUE))
  sampTab$sample_name = rownames(sampTab)

  #extract expression matrix
  if(exp_type == "counts"){
    expDat = counts(sce_object)
  }

  if(exp_type == "normcounts"){
    expDat = normcounts(sce_object)
  }

  if(exp_type == "logcounts"){
    expDat = logcounts(sce_object)
  }

  return(list(sampTab = sampTab, expDat = expDat))

}
