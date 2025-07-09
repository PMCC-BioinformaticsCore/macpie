test_that("plot_qc_metrics_heatmap returns a ggplot", {
  data("mini_mac")
  qc_stats <- compute_qc_metrics(mini_mac, group_by = "combined_id", order_by = "median")
  p <- plot_qc_metrics_heatmap(stats_summary = qc_stats$stats_summary, group_by = "combined_id")
  expect_s3_class(p, "gg")
})
