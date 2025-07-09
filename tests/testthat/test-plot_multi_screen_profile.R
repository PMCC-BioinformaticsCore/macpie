test_that("plot_multi_screen_profile returns a ggplot", {
  data(mini_mac)
  mini_mac <- compute_multi_screen_profile(mini_mac, target = "Staurosporine_10")
  p <- plot_multi_screen_profile(mini_mac, color_by = "Sample_type")
  expect_s3_class(p, "gg")
})
