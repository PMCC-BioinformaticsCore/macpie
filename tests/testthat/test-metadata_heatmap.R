library(testthat)
library(ggplot2)



test_that("metadata_heatmap handles missing Well_ID column", {
  temp_file <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      Plate_ID = c("P1", "P1"),
      Compound_ID = c("C1", "C2")
    ),
    temp_file, row.names = FALSE
  )

  # Expect an error if Well_ID is missing
  expect_error(metadata_heatmap(temp_file), "The metadata is missing a 'Well_ID' column.")
})

test_that("metadata_heatmap handles invalid Well_ID format", {
  temp_file <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      Well_ID = c("INVALID1", "INVALID2"),
      Plate_ID = c("P1", "P1"),
      Compound_ID = c("C1", "C2")
    ),
    temp_file, row.names = FALSE
  )

  # Expect an error for invalid Well_ID format
  expect_error(metadata_heatmap(temp_file), "The 'Well_ID' column contains invalid entries.")
})

test_that("metadata_heatmap handles missing Plate_ID column", {
  temp_file <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      Well_ID = c("A01", "B01"),
      Compound_ID = c("C1", "C2")
    ),
    temp_file, row.names = FALSE
  )

  # Expect an error if Plate_ID is missing
  expect_error(metadata_heatmap(temp_file), "The metadata is missing a 'Plate_ID' column.")
})

test_that("metadata_heatmap handles missing annotation columns", {
  temp_file <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      Well_ID = c("A01", "B01"),
      Plate_ID = c("P1", "P1")
    ),
    temp_file, row.names = FALSE
  )

  # Expect an error if no annotation columns are present
  expect_error(metadata_heatmap(temp_file), "No annotation columns found to generate heatmaps.")
})

test_that("metadata_heatmap generates and saves a plot", {
  temp_file <- tempfile(fileext = ".csv")
  output_file <- tempfile(fileext = ".jpg")
  write.csv(
    data.frame(
      Well_ID = c("A01", "A02", "B01", "B02"),
      Plate_ID = c("P1", "P1", "P1", "P1"),
      Compound_ID = c("C1", "C2", "C3", "C4")
    ),
    temp_file, row.names = FALSE
  )

  # Run the function with output_file specified
  metadata_heatmap(temp_file, output_file = output_file)

  # Reconstruct the actual file path used by the function
  expected_file <- paste0(
    sub(".(png|jpg|pdf)$", "", output_file),
    "_Plate_P1.",
    tools::file_ext(output_file)
  )

  # Check that the file was created
  expect_true(file.exists(expected_file))

  # Clean up
  unlink(expected_file)
})


test_that("metadata_heatmap produces a ggplot object", {
  # Create a temporary CSV file for testing
  temp_file <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      Well_ID = c("A01", "A02", "B01", "B02"),
      Plate_ID = c("P1", "P1", "P1", "P1"),
      Compound_ID = c("C1", "C2", "C3", "C4")
    ),
    temp_file, row.names = FALSE
  )

  # Capture the output
  output <- metadata_heatmap(temp_file)

  # Test that the output is a ggplot object
  expect_s3_class(output, "gg")
  expect_s3_class(output, "ggplot")
})
