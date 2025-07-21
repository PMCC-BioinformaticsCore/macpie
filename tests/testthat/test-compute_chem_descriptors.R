test_that("multiplication works", {
  mock_data <- tibble::tibble(
    Treatment = c("Aspirin", "Caffeine", "NonExistentCompound_123")
  )
  result <- compute_smiles(mock_data, compound_column = "Treatment" )
  data <- compute_chem_descriptors(result, 
     compound_column = "Treatment", 
     treatment_ids = mock_data$Treatment, 
     descriptors = "org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor")
  expect_true("Fsp3" %in% colnames(data))
})
