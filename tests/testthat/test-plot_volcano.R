test_that("ggplot is produced", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"
  top_table <- compute_single_de(testdata, treatment_samples, control_samples, method = "limma_voom")
  p <- plot_volcano(top_table)
  expect_true(inherits(p, "ggplot"))
})
