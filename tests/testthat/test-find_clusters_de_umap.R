test_that("find_clusters_de_umap returns as a seurat object with relevent columns", {
  data(mini_mac)
  mini_mac <- compute_de_umap(mini_mac)
  cluster_mini_mac <- find_clusters_de_umap(mini_mac)
  expect_s4_class(cluster_mini_mac, "Seurat")
  expect_true(all(c("cluster", "UMAPde_1", "UMAPde_2") %in% colnames(cluster_mini_mac@meta.data)))

})
