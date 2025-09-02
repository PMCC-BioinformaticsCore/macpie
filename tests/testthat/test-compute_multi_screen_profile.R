test_that("Results are a non-empty dataframe in seurat slot tools", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  testdata <- compute_multi_screen_profile(testdata, target = "Staurosporine_10",
                                  n_genes_profile = 20, direction = "up", num_cores = 2)
  expect_true(length(testdata) > 0 && nrow(testdata@tools[["screen_profile"]]) > 0)
})
