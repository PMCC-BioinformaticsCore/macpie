test_that("ggplot is produced", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  p<-plot_mds(testdata,"Treatment_1")
  expect_true(inherits(p, "ggplot"))
})
