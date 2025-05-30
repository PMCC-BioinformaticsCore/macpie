test_that("ggplot is produced", {
  data("mini_mac")
  p<-plot_plate_layout(mini_mac,"nCount_RNA","Treatment_1")
  expect_equal(class(p),c("gg","ggplot"))
})
