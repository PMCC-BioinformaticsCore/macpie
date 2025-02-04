test_that("Results are a non-empty data frame", {
  file_path <- system.file("extdata", "PMMSq033/de_screen.Rds", package = "macpie")
  de_list <- readRDS(file_path)
  fgsea_results <- multi_screen_profile(de_list, target = "Staurosporine_10",
                                  n_genes_profile = 100, direction = "up", num_cores = 2)
  expect_equal(class(fgsea_results),c("data.table","data.frame"))
})
