test_that("aggregate_by_de correctly collapses replicates", {
  # Load test dataset from macpie if available
  file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
  skip_if(file_path == "", "Test file not found in macpie package.")
  
  mac <- readRDS(file_path)
  
  # Check if DE data exists
  expect_true("diff_exprs" %in% names(mac@tools), info = "diff_exprs should exist in tools slot")
  
  # Run function
  collapsed_mac <- aggregate_by_de(mac)
  
  # Should now have 1 observation per combined_id
  n_conditions <- length(unique(mac$combined_id))
  expect_equal(ncol(collapsed_mac), n_conditions-1)
  
  # Assay should be 'DE'
  expect_equal(DefaultAssay(collapsed_mac), "DE")
  
  # Metadata should not contain NA for shared metadata columns
  metadata_cols <- colnames(mac@meta.data)
  shared_cols <- metadata_cols[metadata_cols %in% colnames(collapsed_mac@meta.data)]
  non_na_props <- sapply(shared_cols, function(col) {
    mean(!is.na(collapsed_mac@meta.data[[col]]))
  })
  expect_true(mean(non_na_props) > 0.5, info = "Most metadata columns should be retained after collapsing")
})
