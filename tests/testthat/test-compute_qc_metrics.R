test_that("a list of a data frame and a plot is returned", {
  data("mini_mac")
  results <- compute_qc_metrics(
    data = mini_mac, group_by = "combined_id", order_by = "median")
  
  expect_type(results, "list")
  expect_named(results, c("stats_summary", "plot"))
  expect_s3_class(results$plot, "gg")
  
})
