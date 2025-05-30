test_that("Produces non-empty dataframe", {
  data("mini_mac")
  df <- compute_normalised_counts(mini_mac)
  expect_true(nrow(df) > 0)
})
