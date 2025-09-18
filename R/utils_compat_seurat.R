# Internal: robust accessor across SeuratObject v4/v5
.get_assay_data <- function(obj, assay = "RNA", what = c("counts", "data", "scale.data")) {
  what <- match.arg(what)
  if (utils::packageVersion("SeuratObject") < "5.0.0") {
    SeuratObject::GetAssayData(obj, assay = assay, slot  = what)
  } else {
    SeuratObject::GetAssayData(obj, assay = assay, layer = what)
  }
}
