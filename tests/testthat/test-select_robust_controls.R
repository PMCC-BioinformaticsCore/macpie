test_that("multiplication works", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  expect_error(select_robust_controls(testdata, samples=NULL, orig_ident="PLATE"))
  res <- select_robust_controls(testdata, samples="DMSO_0", orig_ident="testdata", make_plots = FALSE)
  expect_true(is.list(res))
})
