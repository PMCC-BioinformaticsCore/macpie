test_that("test plot_gene_ranks returns a ggplot", {
  data("mini_mac")
  p <- plot_gene_ranks(mini_mac,group_by = "combined_id",
  samples = "Staurosporine_10", scale_y_log = TRUE)
  expect_s3_class(p, "gg")

})

