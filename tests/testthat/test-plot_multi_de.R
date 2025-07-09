test_that("plot_multi_de returns a pheatmap object", {
  data("mini_mac")
  p <- plot_multi_de(mini_mac,group_by = "combined_id",
               value = "log2fC", p_value_cutoff = 0.01, direction="up",
               n_genes = 10, control = "DMSO_0", by="fc")
  expect_s3_class(p, "pheatmap")

})


test_that("plot_multi_de works with value = 'lcpm'", {
  data("mini_mac")
  p <- plot_multi_de(mini_mac, group_by = "combined_id",
                     value = "lcpm", p_value_cutoff = 0.01, direction = "up",
                     n_genes = 5, control = "DMSO_0", by = "fc")
  expect_s3_class(p, "pheatmap")
})

test_that("plot_multi_de works with value = 'metric'", {
  data("mini_mac")
  p <- plot_multi_de(mini_mac, group_by = "combined_id",
                     value = "metric", p_value_cutoff = 0.01, direction = "up",
                     n_genes = 5, control = "DMSO_0", by = "fc")
  expect_s3_class(p, "pheatmap")
})

test_that("plot_multi_de throws error if direction is invalid", {
  data("mini_mac")
  expect_error(
    plot_multi_de(mini_mac, direction = "up-regulated"),
    "Value of the direction"
  )
})