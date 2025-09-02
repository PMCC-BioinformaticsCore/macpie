test_that("multiplication works", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  filtered_testdata <- filter_genes_by_expression(testdata,
                                         group_by = "combined_id",
                                         min_counts = 10, min_samples = 2)
  
  expect_true(nrow(filtered_testdata) < nrow(testdata))
  expect_s4_class(filtered_testdata, "Seurat")
})
