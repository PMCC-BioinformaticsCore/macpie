test_that("Results are a non-empty dataframe in seurat slot tools", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  data("genesets")
  testdata <- compute_multi_enrichr(testdata, genesets = genesets)
  expect_true(length(testdata) > 0 && length(testdata@tools[["pathway_enrichment"]]) > 0)
})
