test_that("Produces non-empty dataframe with a specific top gene.", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"

  #perform differential expression
  top_table<-compute_single_de(testdata, treatment_samples, control_samples,method = "limma_voom")
  expect_s3_class(top_table, "data.frame")
  expect_true(all(c("gene", "log2FC", "metric", "p_value", "p_value_adj") %in% colnames(top_table)))
  expect_true(nrow(top_table) > 0 & top_table[top_table$gene == "AHCY","log2FC"]>0)
})

