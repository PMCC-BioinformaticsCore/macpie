test_that("ggplot is produced", {
  rds_file<-system.file("/extdata/PMMSq033/PMMSq033.rds", package = "macpie")
  expect_true(file.exists(rds_file)) # Check if the file exists
  test_data<-readRDS(rds_file)
  rows<-test_data[["Row"]] %>%
    pull()
  p<-plot_rle(data = test_data, labels = rows, log = TRUE)
  expect_equal(class(p),c("gg","ggplot"))
})
