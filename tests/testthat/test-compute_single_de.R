test_that("Produces non-empty dataframe with a specific top gene.", {

  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"

  #perform differential expression
  top_table<-compute_single_de(mac, treatment_samples, control_samples,method = "limma_voom")
  expect_true(nrow(top_table) > 0 & top_table[top_table$gene == "ANAPC16","log2FC"]>0)
})
