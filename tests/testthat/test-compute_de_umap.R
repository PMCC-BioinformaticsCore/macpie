test_that("Results are a non-empty dataframe in seurat slot tools", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  testdata <- compute_de_umap(testdata)
  expect_true(length(testdata) > 0 && nrow(testdata@reductions[["umap_de"]]) > 0)
})
