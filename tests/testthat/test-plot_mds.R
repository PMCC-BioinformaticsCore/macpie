test_that("ggplot is produced", {
  data("mini_mac")
  p<-plot_mds(mini_mac,"Treatment_1")
  expect_equal(class(p),c("gg","ggplot"))
})
