test_that("Checking metadata works", {
  file_path <- system.file("/extdata/PMMSq033_metadata.csv", package = "macpie")
  result <- read_metadata(file_path)
  test_columns <- all(colnames(result) %in% c("Plate_ID", "Well_ID", "Row", "Column", "Species",
                                              "Cell_type", "Model_type", "Time", "Unit",
                                              "Treatment_1", "Concentration_1", "Unit_1",
                                              "Sample_type", "Barcode", "Project"))
  expect_true(test_columns)
})
