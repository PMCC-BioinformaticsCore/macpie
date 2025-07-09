test_that("multiplication works", {
  mock_data <- tibble::tibble(
    Treatment_1 = c("Aspirin", "Caffeine", "NonExistentCompound_123")
  )
  result <- compute_smiles(mock_data,compound_column = "Treatment_1")
  data <- compute_chem_descriptors(result, descriptors =
  "org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor")
  expect_true("Fsp3" %in% colnames(data))
})
