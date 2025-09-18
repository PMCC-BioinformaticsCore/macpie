test_that("ggplot is produced", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  rows<-testdata[["Row"]] %>%
    pull()
  p<-plot_rle(data = testdata, labels = rows, log = TRUE)
  expect_true(inherits(p, "ggplot"))
})
