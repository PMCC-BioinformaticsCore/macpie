test_that("test plot_gene_ranks returns a ggplot", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  p <- plot_gene_ranks(testdata,group_by = "combined_id",
  samples = "Staurosporine_10", scale = TRUE)
  expect_s3_class(p, "gg")

})

