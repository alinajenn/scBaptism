test_that("metacell function works", {
  # run scBaptism to annotate with SingleR and CIA
  mc_ten_result <- calculate_metacells(sce_query = sce_annotated,
                               annotation_col = NULL,
                               gamma = 10,
                               n.pc = 10)

  expect_s4_class(mc_ten_result, "SingleCellExperiment")
  expect_true(ncol(mc_ten_result) == 264)

  mc_anno_result <- calculate_metacells(sce_query = sce_annotated,
                               annotation_col = "labels_main",
                               gamma = 30,
                               n.pc = 10)

  expect_true("annotation" %in% names(colData(mc_anno_result)))

})

