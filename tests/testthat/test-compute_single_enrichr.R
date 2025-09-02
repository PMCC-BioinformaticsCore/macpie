test_that("Results are a non-empty data frame", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"
  top_table <- compute_single_de(testdata, treatment_samples, control_samples, method = "limma_voom")
  top_genes <- top_table$gene[top_table$p_value_adj<0.01]
  pe <- compute_single_enrichr(genes = top_genes, db = "MSigDB_Hallmark_2020", species = "human")
  expect_true(length(pe) > 0 && nrow(pe) > 0)
})
