test_that("ggplot is produced", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  p<-plot_plate_layout(mac,"nCount_RNA","Treatment_1")
  expect_equal(class(p),c("gg","ggplot"))
})
