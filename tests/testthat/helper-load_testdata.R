# tests/testthat/helper-load_testdata.R
get_testdata <- function(name = c("testdata")) {
  name <- match.arg(name)
  path <- testthat::test_path("testdata", paste0(name, ".rds"))
  readRDS(path)
}
