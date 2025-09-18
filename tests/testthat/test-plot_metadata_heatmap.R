test_that("ggplot is produced", {
  metadata_file_path <- system.file("extdata", "PMMSq033_metadata.csv", package = "macpie")
  p <- plot_metadata_heatmap(metadata_file=metadata_file_path)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})
