#metacellify funtion

calculate_metacells <- function(sce_query){

  #get the log normalized expression matrix

  query_matrix <- SummarizedExperiment::assay(sce_query, "logcounts")

  #run the metacell calculations, metacell identity for each single cell

  mc_identity <- SuperCell::SCimplify(X = query_matrix)

  ######################old
  #add information back to the query SCE
  #SummarizedExperiment::colData(sce_query)$metacells <- MC_identity$membership
  #return(sce_query)

  #get average expression for each metacell

  mc_matrix <- SuperCell::supercell_GE(ge = query_matrix,
                                       groups = mc_identity$membership)

  #create new SCE

  sce_metacells <- SingleCellExperiment(assays = list(counts = mc_matrix))

  #calculate logcounts

  sce_metacells <- scuttle::logNormCounts(sce_metacells)

  #calcuate tSNE
  ##################possibly add other dim reductions
  sce_metacells <- scater::runTSNE(sce_metacells)

  return(sce_metacells)

  ####################add in cell names?

}
