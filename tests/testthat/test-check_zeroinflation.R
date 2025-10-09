test_that("chech_zeroinflation works", {
  # Load test dataset from tests/testthat/testdata
  testdata <- get_testdata()
  #only match one group should fail the minimum requirement of two groups
  expect_error(check_zeroinflation(data = testdata, group_by = "combined_id",
                                    samples = c("DrugA_10", "DMSO_0"),
                                    cutoffs = c(0.1, 0.2)))
  #cutoffs out of range
  expect_error(check_zeroinflation(data = testdata, group_by = "combined_id",
                                   samples = c("Luminespib_10", "DMSO_0"),
                                    cutoffs = 1.2))
  # expect a list as output
  res <- check_zeroinflation(data = testdata, group_by = "combined_id",
                             samples = c("Luminespib_10", "DMSO_0"),
                             cutoffs = c(0.1, 0.2))
  expect_type(res, "list")
}) 

