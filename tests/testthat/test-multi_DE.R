test_that("Results are a non-empty list", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)[50:150]
  control_samples<-"DMSO_0"
  treatment_samples <- mac %>%
    select(combined_id) %>%
    filter(!grepl("DMSO", combined_id)) %>%
    pull() %>%
    unique()
  mac <- multi_DE(mac, treatment_samples, control_samples, num_cores = 2, method = "edgeR")
  expect_true(length(mac) > 0 && ncol(mac@tools[["diff_exprs"]][[1]]) == 6)
})
