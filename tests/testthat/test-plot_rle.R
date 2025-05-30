test_that("ggplot is produced", {
  data("mini_mac")
  rows<-mini_mac[["Row"]] %>%
    pull()
  p<-plot_rle(data = mini_mac, labels = rows, log = TRUE)
  expect_equal(class(p),c("gg","ggplot"))
})
