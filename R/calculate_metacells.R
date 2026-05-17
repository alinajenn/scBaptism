# metacellify function
# make MCs before running any annotation

calculate_metacells <- function(sce_query){  ##give options for parameter customization (gamma, mode)

  #get the log normalized expression matrix

  query_matrix <- SummarizedExperiment::assay(sce_query, "logcounts")

  #run the metacell calculations, metacell identity for each single cell

  mc_identity <- SuperCell::SCimplify(X = query_matrix)

  #get average gene expression for each metacell

  mc_matrix <- SuperCell::supercell_GE(ge = query_matrix,
                                       groups = mc_identity$membership)

  #create new SCE
  ############could be replaced with supercell_2_sce

  sce_metacells <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = mc_matrix))

  #calculate logcounts, does this make sense or is it doppelt-gemoppelt?
  #Problem: I need counts & logcounts in order to run_scBaptism

  #sce_metacells <- scuttle::logNormCounts(sce_metacells)

  #add cellnames (form scClassify)
  colnames(sce_metacells) <- paste0("cell_", seq_len(ncol(sce_metacells)))

  #clusters (for clustifyr)
  library(scran)

  g <- buildSNNGraph(sce_metacells, k=10, use.dimred = 'PCA')
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(sce_metacells) <- factor(clust)

  names(SummarizedExperiment::colData(sce_metacells))[which(names(SummarizedExperiment::colData(sce_metacells))=="label")]="clusters"

  #run PCA
  sce_metacells <- scater::runPCA(sce_metacells)

  #calcuate tSNE (for plotting later)
  ##################possibly add other dim reductions
  sce_metacells <- scater::runTSNE(sce_metacells)

  return(sce_metacells)

  ####################add in cell names?

}

#Leas example

# sce_counts <- assay(sce, "logcounts", cell.annotation = sce$labels_main)
#
# super <- SCimplify(sce_counts, gamma = 20)
#
# super_counts <- supercell_GE(sce_counts, super$membership, mode="average")
#
# sce_super <- supercell_2_sce(
#   SC.GE = super_counts,
#   SC = super
# )


######################old
#add information back to the query SCE
#SummarizedExperiment::colData(sce_query)$metacells <- MC_identity$membership
#return(sce_query)
