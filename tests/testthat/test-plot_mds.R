test_that("ggplot is produced", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  p<-plot_mds(mac[,1:50],"Treatment_1")
  expect_equal(class(p),c("gg","ggplot"))
})
