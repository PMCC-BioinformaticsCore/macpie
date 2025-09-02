test_that("Produces non-empty dataframe", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  df <- compute_normalised_counts(testdata)
  expect_true(nrow(df) > 0)
})
