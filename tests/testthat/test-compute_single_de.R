test_that("Produces non-empty dataframe with a specific top gene.", {
  data("mini_mac")
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"

  #perform differential expression
  top_table<-compute_single_de(mini_mac, treatment_samples, control_samples,method = "limma_voom")
  expect_true(nrow(top_table) > 0 & top_table[top_table$gene == "RAB5C","log2FC"]>0)
})
