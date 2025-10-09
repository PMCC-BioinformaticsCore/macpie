test_that("subsample_genes works", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  expect_error(subsample_genes(testdata,
                  ngene = 1000,
                  gselect = "random", 
                  seed = 1))
  expect_error(subsample_genes(testdata,
                  ngene = 100,
                  gselect = "min", 
                  seed = 1))
  subsampled_genes <- subsample_genes(testdata,
                                     ngene = 50,
                                     gselect = "random", 
                                     seed = 1)
  expect_s4_class(subsampled_genes,"Seurat")
})
