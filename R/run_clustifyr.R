#' run_clustifyr
#'
#' @param sce_query SCE to be annotated
#' @param reference SCE object that acts as a reference
#' @param ref_labs List of gene labels or column from the references colData
#' @param clusters name of clusters in the query
#' @param query_genes vector of genes of interest to compare. If NULL, then common genes between the expr_mat and ref_mat will be used for comparision.
#' @param n_genes number of genes limit for Seurat variable genes, by default 1000, set to 0 to use all variable genes (generally not recommended)
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param n_perm number of permutations, set to 0 by default
#' @param compute_method method(s) for computing similarity scores
#' @param pseudobulk_method method used for summarizing clusters, options are mean (default), median, truncate (10% truncated mean), or trimean, max, min
#' @param rm0	consider 0 as missing data, recommended for per_cell
#' @param obj_out whether to output object instead of cor matrix
#' @param seurat_out output cor matrix or called seurat object (deprecated, use obj_out instead)
#' @param vec_out only output a result vector in the same order as metadata
#' @param rename_prefix prefix to add to type and r column names
#' @param clustifyr_threshold identity calling minimum correlation score threshold, only used when obj_out = TRUE
#' @param low_threshold_cell option to remove clusters with too few cells
#' @param exclude_genes a vector of gene names to throw out of query
#' @param if_log input data is natural log, averaging will be done on unlogged data
#' @param organism for GO term analysis, organism name: human - 'hsapiens', mouse - 'mmusculus'
#' @param plot_name name for saved pdf, if NULL then no file is written (default)
#' @param rds_name name for saved rds of rank_diff, if NULL then no file is written (default)
#' @param expand_unassigned test all ref clusters for unassigned results
#' @param use_var_genes if providing a seurat object, use the variable genes (stored in seurat_object@var.genes) as the query_genes.
#' @param dr	stored dimension reduction
#' @param verbose display message after annotation is finished
#'
#'
#' @returns sce_query a SingleCellExperiment object, with the extra info on the
#' annotated cells
#'
#' @export
#'
#' @importFrom clustifyr clustify
#' @importFrom clustifyr object_ref
#' @importFrom SummarizedExperiment colData
#'
#'
#' @examples
#'
#' library(iUSEiSEE)
#' library(Seurat)
#' library(dplyr)
#'
#' #load SCE from iUSEiSEE
#'
#' sce_annotated <- readRDS(file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))
#'
#' #run the annotation
#' sce_annotated <- run_clustifyr(sce_query = sce_annotated,
#'                                reference = sce_annotated,
#'                                ref_labs = "labels_main")
#'
#' #plot the existing annotation with scater(t-SNE)
#' scater::plotTSNE(sce_annotated, color_by = "scb_clustifyr_labels")
#'
#'
#'@family reference-based family
#'
run_clustifyr <- function(sce_query,
                          reference,
                          ref_labs,
                          clusters,
                          query_genes = NULL,
                          per_cell = FALSE,
                          n_perm = 0,
                          compute_method = "spearman",
                          pseudobulk_method = "mean",
                          dr = "umap",
                          obj_out = TRUE,
                          seurat_out = obj_out,
                          vec_out = FALSE,
                          clustifyr_threshold = "auto",
                          rm0 = FALSE,
                          rename_prefix = NULL,
                          exclude_genes = c(),
                          metadata = NULL,
                          organism = "hsapiens",
                          plot_name = NULL,
                          rds_name = NULL,
                          expand_unassigned = FALSE,
                          verbose = FALSE
                          )

{

  # checks ----------------------------------------------------------------


  # transformation --------------------------------------------------------

  #input is fine as SCE

  #bulid the reference with object_ref() function from clustifyr

  transf_ref <- clustifyr::object_ref(
    input = reference,
    cluster_col = ref_labs # they call it cluster_col, but ask for cell identities in the tutorial
  )


  # running annotation-----------------------------------------------------

  sce_query <- clustifyr::clustify(
    input = sce_query,
    ref_mat = transf_ref,
    cluster_col = clusters, #actual cluster information
    per_cell = per_cell,
    n_perm = n_perm,
    compute_method = compute_method,
    pseudobulk_method = pseudobulk_method,
    dr = dr,
    obj_out = obj_out,
    seurat_out = seurat_out,
    vec_out = vec_out,
    threshold = clustifyr_threshold,
    rm0 = rm0,
    rename_prefix = rename_prefix,
    exclude_genes = exclude_genes,
    metadata = metadata,
    organism = organism,
    plot_name = plot_name,
    rds_name = rds_name ,
    expand_unassigned = expand_unassigned,
    verbose = verbose
    )


  # return input SCE with new annotation-----------------------------------


  #gets added automatically to the sce_query, but need to rename column

  names(SummarizedExperiment::colData(sce_query))[which(names(SummarizedExperiment::colData(sce_query))=="type")]="scb_clustifyr_labels"

  #rename column for correlation coefficient

  names(SummarizedExperiment::colData(sce_query))[which(names(SummarizedExperiment::colData(sce_query))=="r")]="scb_clustifyr_corr_coef"

  #message
  if(verbose) message("clustifyr annotation done")

  return(sce_query)

}
