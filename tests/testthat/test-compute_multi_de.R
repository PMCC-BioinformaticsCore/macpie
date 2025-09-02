test_that("Results are a non-empty list in the slot tools of a Seurat object", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  control_samples <- "DMSO_0"
  treatment_samples <- testdata$combined_id[!grepl("DMSO", testdata$combined_id)][1:10]
  data <- compute_multi_de(testdata, treatment_samples, control_samples, num_cores = 1, method = "edgeR")
  expect_true(length(data) > 0 && ncol(data@tools[["diff_exprs"]][[1]]) == 6)
})
