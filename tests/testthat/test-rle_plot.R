test_that("ggplot is produced", {
  rds_file<-system.file("/extdata/PMMSq033/PMMSq033.rds", package = "macpie")
  mac<-readRDS(rds_file)
  count_matrix<-as.matrix(mac@assays$RNA$counts)
  colnames(count_matrix)<-mac$Well_ID
  p<-rle_plot(count_matrix = count_matrix, id = mac$Well_ID, feature = mac$Row, logged=FALSE)
  expect_equal(class(p),c("gg","ggplot"))
})
