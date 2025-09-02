test_that("find_clusters_de_umap returns as a seurat object with relevent columns", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  testdata <- compute_de_umap(testdata)
  cluster_testdata <- find_clusters_de_umap(testdata)
  expect_s4_class(cluster_testdata, "Seurat")
  expect_true(all(c("cluster", "UMAPde_1", "UMAPde_2") %in% colnames(cluster_testdata@meta.data)))

})
