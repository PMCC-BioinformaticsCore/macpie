test_that("Produces non-empty dataframe", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  df <- fetch_normalised_counts(mac[,1:10])
  expect_true(nrow(df) > 0)
})
