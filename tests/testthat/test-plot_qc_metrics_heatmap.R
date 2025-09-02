test_that("plot_qc_metrics_heatmap returns a ggplot", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  qc_stats <- compute_qc_metrics(testdata, group_by = "combined_id", order_by = "median")
  p <- plot_qc_metrics_heatmap(stats_summary = qc_stats$stats_summary, group_by = "combined_id")
  expect_s3_class(p, "gg")
})
