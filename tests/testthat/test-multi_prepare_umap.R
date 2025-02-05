test_that("Results are a non-empty dataframe in seurat slot tools", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  mac <- multi_prepare_umap(mac)
  expect_true(length(mac) > 0 && nrow(mac@tools[["umap_de"]]) > 0)
})
