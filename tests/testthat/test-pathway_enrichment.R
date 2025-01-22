test_that("Results are a non-empty data frame", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"
  top_table <- differential_expression(mac, treatment_samples, control_samples, method = "limma_voom")
  top_genes <- top_table$gene[top_table$p_value_adj<0.01]
  pe <- pathway_enrichment(genes = top_genes, db = "MSigDB_Hallmark_2020", species = "human")
  expect_true(length(pe) > 0 && nrow(pe) > 0)
})
