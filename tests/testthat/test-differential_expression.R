test_that("Produces non-empty dataframe with a specific top gene.", {

  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"

  #perform differential expression
  top_table<-differential_expression(mac, treatment_samples, control_samples,method = "limma_voom")
  expect_true(nrow(top_table) > 0 & top_table["ANAPC16","logFC"]>0)
})
