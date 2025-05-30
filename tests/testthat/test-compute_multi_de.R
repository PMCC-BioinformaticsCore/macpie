test_that("Results are a non-empty list in the slot tools of a Seurat object", {
  data("mini_mac")
  control_samples <- "DMSO_0"
  treatment_samples <- mini_mac$combined_id[!grepl("DMSO", mini_mac$combined_id)][1:10]
  data <- compute_multi_de(mini_mac, treatment_samples, control_samples, num_cores = 1, method = "edgeR")
  expect_true(length(data) > 0 && ncol(data@tools[["diff_exprs"]][[1]]) == 6)
})
