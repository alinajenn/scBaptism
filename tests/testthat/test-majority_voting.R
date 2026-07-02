test_that("majority_vote function works", {
  # run scBaptism to annotate with SingleR and CIA
  mv_result <- majority_vote(sce_query = sce_annotated,
                              anno_columns = c("labels_main", "labels_fine"),
                              tie_breaker = "first")

  expect_s4_class(mv_result, "SingleCellExperiment")
  expect_true(all(c("scb_majority_labels") %in% names(colData(mv_result))))

})
