test_that("plot_qc_metrics returns a ggplot", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  qc_stats <- compute_qc_metrics(testdata, group_by = "combined_id", order_by = "median")
  p <- plot_qc_metrics(qc_stats, group_by = "combined_id", metrics = "sd_value")
  expect_s3_class(p, "gg")
})
