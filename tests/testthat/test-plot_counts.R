test_that("test plot_counts returns a ggplot ", {
  data(mini_mac)
  genes <- mini_mac@tools$diff_exprs$Staurosporine_10$gene[1:6]
  p <- plot_counts(mini_mac, genes = genes, group_by = "combined_id",
  treatment_samples = "Staurosporine_10",
  control_samples = "DMSO_0",
  normalisation = "clr")
  expect_s3_class(p, "gg")
})
