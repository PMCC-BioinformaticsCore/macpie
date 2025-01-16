test_that("ggplot is produced", {
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  mac <- readRDS(file_path)
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"
  top_table <- differential_expression(mac, treatment_samples, control_samples, method = "limma_voom")
  p <- plot_volcano(top_table)
  expect_equal(class(p),c("gg","ggplot"))
})
