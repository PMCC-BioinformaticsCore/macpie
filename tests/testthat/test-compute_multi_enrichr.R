test_that("Results are a non-empty dataframe in seurat slot tools", {
  data("mini_mac")
  data("genesets")
  mini_mac <- compute_multi_enrichr(mini_mac, genesets = genesets)
  expect_true(length(mini_mac) > 0 && length(mini_mac@tools[["pathway_enrichment"]]) > 0)
})
