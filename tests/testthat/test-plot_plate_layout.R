test_that("ggplot is produced", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  p<-plot_plate_layout(testdata,"nCount_RNA","Treatment_1")
  expect_equal(class(p),c("gg","ggplot"))
})
