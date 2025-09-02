test_that("plot_multi_screen_profile returns a ggplot", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  testdata <- compute_multi_screen_profile(testdata, target = "Staurosporine_10", n_genes_profile = 100)
  p <- plot_multi_screen_profile(testdata, color_by = "Sample_type")
  expect_s3_class(p, "gg")
})
