test_that("summarise_de returns a tibble object", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  top_table <- testdata@tools$diff_exprs$Staurosporine_10
  result_1 <- summarise_de(top_table, lfc_threshold = 1, padj_threshold = 0.05)
  result_2 <- summarise_de(testdata, lfc_threshold = 1,
  padj_threshold = 0.01, multi=TRUE)
  
  expect_s3_class(result_1, "tbl_df")
  expect_s3_class(result_2, "tbl_df")
})
