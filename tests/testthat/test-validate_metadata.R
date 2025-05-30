test_that("Non-empty output", {
  file_path <- system.file("/extdata/PMMSq033_metadata.csv", package = "macpie")
  metadata <- read_metadata(file_path)
  metadt_qc <- validate_metadata(metadata)
  expect_true(nrow(metadt_qc) > 0)
})
