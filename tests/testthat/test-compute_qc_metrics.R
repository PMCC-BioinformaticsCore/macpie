test_that("a list of a data frame and a plot is returned", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  results <- compute_qc_metrics(
    data = testdata, group_by = "combined_id", order_by = "median")
  
  expect_type(results, "list")
  expect_named(results, c("stats_summary", "plot"))
  expect_s3_class(results$plot, "gg")
  
})

test_that("compute_qc_metrics handles valid order_by column in metadata", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  expect_no_error(
    compute_qc_metrics(testdata, group_by = "combined_id", order_by = "median")
  )
})

test_that("compute_qc_metrics output stats_summary includes expected columns", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  result <- compute_qc_metrics(data = testdata, group_by = "combined_id", order_by = "median")
  expect_true(all(c("sd_value", "mad_value", "group_median", "z_score", "IQR") %in% colnames(result$stats_summary)))
})
