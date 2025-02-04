test_that("Results are a non-empty dataframe in seurat slot tools", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  mac <- multi_screen_profile(mac, target = "Staurosporine_10",
                                  n_genes_profile = 100, direction = "up", num_cores = 2)
  expect_true(length(mac) > 0 && nrow(mac@tools[["screen_profile"]]) > 0)
})
