test_that("Results are a non-empty list in the slot tools of a Seurat object", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  data <- readRDS(file_path)[50:150]
  control_samples <- "DMSO_0"
  treatment_samples <- data$combined_id[!grepl("DMSO", data$combined_id)]
  data <- compute_multi_de(data, treatment_samples, control_samples, num_cores = 1, method = "edgeR")
  expect_true(length(data) > 0 && ncol(data@tools[["diff_exprs"]][[1]]) == 6)
})
