test_that("read_metadata reads valid CSV correctly", {
  # Create a temporary CSV file with correct columns
  tmp_csv <- tempfile(fileext = ".csv")
  metadata <- tibble::tibble(
    Plate_ID = "Plate1",
    Well_ID = "A01",
    Row = "A",
    Column = 1,
    Species = "human",
    Treatment_1 = "DMSO",
    Concentration_1 = 1,
    Sample_type = "control",
    Barcode = "ABC123"
  )
  readr::write_csv(metadata, tmp_csv)
  
  result <- read_metadata(tmp_csv)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("Plate_ID", "Barcode") %in% colnames(result)))
})

test_that("read_metadata returns NULL when required columns are missing", {
  tmp_csv <- tempfile(fileext = ".csv")
  incomplete_metadata <- tibble::tibble(
    Plate_ID = "Plate1",
    Well_ID = "A01"
    # Missing many required columns
  )
  readr::write_csv(incomplete_metadata, tmp_csv)
  result <- read_metadata(tmp_csv)
  
  expect_null(result)
})

test_that("read_metadata throws error on nonexistent file", {
  expect_error(read_metadata("nonexistent.csv"), "does not exist")
})

test_that("read_metadata rejects unsupported file formats", {
  tmp_txt <- tempfile(fileext = ".txt")
  writeLines("Hello world", con = tmp_txt)
  expect_error(read_metadata(tmp_txt), "Unsupported file format")
})