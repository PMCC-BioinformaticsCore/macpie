test_that("plot_distance returns a pheatmap", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  p <- plot_distance(testdata, group_by = "combined_id", treatment = "DMSO_0")
  expect_s3_class(p, "pheatmap")
  
})
