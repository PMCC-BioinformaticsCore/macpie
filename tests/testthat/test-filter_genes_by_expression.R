test_that("multiplication works", {
  data(mini_mac)
  filtered_mini_mac <- filter_genes_by_expression(mini_mac,
                                         group_by = "combined_id",
                                         min_counts = 10, min_samples = 2)
  
  expect_true(nrow(filtered_mini_mac) < nrow(mini_mac))
  expect_s4_class(filtered_mini_mac, "Seurat")
})
