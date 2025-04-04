test_that("Results are a non-empty dataframe in seurat slot tools", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  file_path <- system.file("extdata", "PMMSq033/pathways.Rds", package = "macpie")
  genesets <- readRDS(file_path)
  mac <- compute_multi_enrichr(mac, genesets = genesets)
  expect_true(length(mac) > 0 && length(mac@tools[["pathway_enrichment"]]) > 0)
})
