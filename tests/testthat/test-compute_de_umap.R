test_that("Results are a non-empty dataframe in seurat slot tools", {
  data("mini_mac")
  mini_mac <- compute_de_umap(mini_mac)
  expect_true(length(mini_mac) > 0 && nrow(mini_mac@tools[["umap_de"]]) > 0)
})
