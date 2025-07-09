test_that("convert_human_to_mouse_ids returns a list of mouse gene symbols", {
  name_list <- convert_human_to_mouse(c("BRCA1", "TRAF1", "MYBL1"))
  expect_equal(
    name_list,
    c("Brca1", "Traf1", "Mybl1")
  )
})
