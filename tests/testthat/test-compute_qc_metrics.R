test_that("a list of a data frame and a plot is returned", {
  data("mini_mac")
  results <- compute_qc_metrics(
    data = mini_mac, group_by = "combined_id", order_by = "median")
  
  expect_type(results, "list")
  expect_named(results, c("stats_summary", "plot"))
  expect_s3_class(results$plot, "gg")
  
})

test_that("compute_qc_metrics handles valid order_by column in metadata", {
  data(mini_mac, package = "macpie")
  expect_no_error(
    compute_qc_metrics(mini_mac, group_by = "combined_id", order_by = "Row")
  )
})

test_that("compute_qc_metrics output stats_summary includes expected columns", {
  data(mini_mac, package = "macpie")
  result <- compute_qc_metrics(data = mini_mac, group_by = "combined_id", order_by = "median")
  expect_true(all(c("sd_value", "mad_value", "group_median", "z_score", "IQR") %in% colnames(result$stats_summary)))
})