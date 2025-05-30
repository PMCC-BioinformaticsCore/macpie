test_that("ggplot is produced", {
  data("mini_mac")
  treatment_samples="Staurosporine_0.1"
  control_samples<-"DMSO_0"
  top_table <- compute_single_de(mini_mac, treatment_samples, control_samples, method = "limma_voom")
  p <- plot_volcano(top_table)
  expect_equal(class(p),c("gg","ggplot"))
})
