test_that("Results are a non-empty list", {
  genesets <- download_geneset("human", "MSigDB_Hallmark_2020")
  expect_true(length(genesets) > 0 && inherits(genesets, "list"))
})
