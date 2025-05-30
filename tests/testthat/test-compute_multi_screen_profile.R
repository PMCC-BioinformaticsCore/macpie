test_that("Results are a non-empty dataframe in seurat slot tools", {
  data("mini_mac")
  mini_mac <- compute_multi_screen_profile(mini_mac, target = "Staurosporine_10",
                                  n_genes_profile = 20, direction = "up", num_cores = 2)
  expect_true(length(mini_mac) > 0 && nrow(mini_mac@tools[["screen_profile"]]) > 0)
})
