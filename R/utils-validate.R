#' Validate macpie input object
#' Accepts Seurat, tidyseurat, or plain data.frame/tibble (easily extendable).
#' @keywords internal
validate_macpie_input <- function(data, required = NULL) {
  ok_class <- inherits(data, "Seurat") ||
    inherits(data, "tidyseurat") ||
    inherits(data, "data.frame") # tibble inherits data.frame
  
  if (!ok_class) {
    stop("Input must be a Seurat, tidyseurat, or data.frame/tibble.", call. = FALSE)
  }

}
