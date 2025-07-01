test_that("Compute_smiles returns correct SMILES and handles edge cases", {
  skip_if_not_installed("webchem")
  skip_if_offline()
  
  # Create mock input tibble (like a tidyseurat object)
  mock_data <- tibble::tibble(
    Treatment_1 = c("Aspirin", "Caffeine", "NonExistentCompound_123")
  )
  
  # Run the function
  result <- compute_smiles(mock_data, compound_column = "Treatment_1")
  
  # Check columns
  expect_true("smiles" %in% colnames(result))
  
  # Check that known SMILES are returned for valid compounds
  expect_true(!is.na(result$smiles[result$Treatment_1 == "Aspirin"]))
  expect_true(!is.na(result$smiles[result$Treatment_1 == "Caffeine"]))
  
  # Check that unknown compound returns NA
  expect_true(is.na(result$smiles[result$Treatment_1 == "NonExistentCompound_123"]))
  
  # Check that number of rows is unchanged
  expect_equal(nrow(result), nrow(mock_data))
})