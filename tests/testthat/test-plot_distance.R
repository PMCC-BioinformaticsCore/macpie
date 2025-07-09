test_that("plot_distance returns a pheatmap", {
  data(mini_mac)
  p <- plot_distance(mini_mac, group_by = "combined_id", treatment = "DMSO_0")
  expect_s3_class(p, "pheatmap")
  
})
