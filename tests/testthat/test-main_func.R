test_that("main function works", {
  # run scBaptism to annotate with SingleR and CIA
  anno_result <- run_scBaptism(sce_query = sce_annotated,
                               anno_methods = c("SingleR", "CIA"),
                               markers_list = markers_lists,
                               reference = sce_annotated,
                               ref_labs = "labels_main")

  expect_s4_class(anno_result, "SingleCellExperiment")
  expect_true(all(c("scb_SingleR_labels", "scb_CIA_labels") %in% names(colData(anno_result))))

  expect_message(run_scBaptism(sce_query = sce_annotated,
                               anno_methods = c("CIA"),
                               markers_list = markers_lists,
                               reference = sce_annotated,
                               ref_labs = "labels_main",
                               verbose = TRUE), regexp = "CIA annotation done")

  anno_result_allmethods <- run_scBaptism(sce_query = sce_annotated,
                                          anno_methods = c("SingleR", "scmap",
                                                           "CIA", "scClassify"),
                                          markers_list = markers_lists,
                                          reference = sce_annotated,
                                          ref_labs = "labels_main")
  expect_s4_class(anno_result_allmethods, "SingleCellExperiment")

})



test_that("triggering diverse errors", {
  expect_error(run_scBaptism(sce_query = sce_annotated,
                             anno_methods = c("SingleR"),
                             markers_list = sce_annotated,
                             reference = sce_annotated,
                             ref_labs = "labels_main"),
               regexp = "markers_list")

  expect_error(run_scBaptism(sce_query = sce_annotated,
                             anno_methods = c("CIA"),
                             markers_list = markers_lists,
                             reference = markers_lists,
                             ref_labs = "labels_main"),
               regexp = "SingleCellExperiment")

})



