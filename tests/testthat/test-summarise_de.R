test_that("summarise_de returns a tibble object", {
  data("mini_mac")
  top_table <- mini_mac@tools$diff_exprs$Staurosporine_10
  result_1 <- summarise_de(top_table, lfc_threshold = 1, padj_threshold = 0.05)
  result_2 <- summarise_de(mini_mac, lfc_threshold = 1,
  padj_threshold = 0.01, multi=TRUE)
  
  expect_s3_class(result_1, "tbl_df")
  expect_s3_class(result_2, "tbl_df")
})
