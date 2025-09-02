test_that("test plot_counts returns a ggplot ", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  genes <- testdata@tools$diff_exprs$Staurosporine_10$gene[1:6]
  p <- plot_counts(testdata, genes = genes, group_by = "combined_id",
  treatment_samples = "Staurosporine_10",
  control_samples = "DMSO_0",
  normalisation = "clr")
  expect_s3_class(p, "gg")
})
