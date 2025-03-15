test_that("ggplot is produced", {
  metadata_file_path <- system.file("extdata", "PMMSq033/PMMSq033_metadata.csv", package = "macpie")
  p <- plot_metadata_heatmap(metadata_file=metadata_file_path)
  expect_equal(class(p),c("patchwork", "gg", "ggplot"))
})
