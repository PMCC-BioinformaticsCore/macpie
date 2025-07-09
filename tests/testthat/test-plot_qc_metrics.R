test_that("plot_qc_metrics returns a ggplot", {
  data("mini_mac")
  qc_stats <- compute_qc_metrics(mini_mac, group_by = "combined_id", order_by = "median")
  p <- plot_qc_metrics(qc_stats, group_by = "combined_id", metrics = "sd_value")
  expect_s3_class(p, "gg")
})
