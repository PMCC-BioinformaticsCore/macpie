test_that("plot_multi_de returns a pheatmap object", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  p <- plot_multi_de(testdata,group_by = "combined_id",
               value = "log2fC", p_value_cutoff = 0.01, direction="up",
               n_genes = 10, control = "DMSO_0", by="fc")
  expect_s3_class(p, "pheatmap")

})


test_that("plot_multi_de works with value = 'lcpm'", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  p <- plot_multi_de(testdata, group_by = "combined_id",
                     value = "lcpm", p_value_cutoff = 0.01, direction = "up",
                     n_genes = 5, control = "DMSO_0", by = "fc")
  expect_s3_class(p, "pheatmap")
})

test_that("plot_multi_de works with value = 'metric'", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  p <- plot_multi_de(testdata, group_by = "combined_id",
                     value = "metric", p_value_cutoff = 0.01, direction = "up",
                     n_genes = 5, control = "DMSO_0", by = "fc")
  expect_s3_class(p, "pheatmap")
})

test_that("plot_multi_de throws error if direction is invalid", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  expect_error(
    plot_multi_de(testdata, direction = "up-regulated"),
    "Value of the direction"
  )
})