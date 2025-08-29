#' A tiny demo dataset for quick-start examples
#'
#' `demo` is a small Seurat object derived from the package's internal
#' example dataset `mini_mac`, reduced to a small numberof genes
#' for fast examples in the quick start in README.
#'
#' @format A Seurat object with:
#' \describe{
#'   \item{assays$RNA}{counts and data matrices with ~200 genes × ~48 cells}
#'   \item{meta.data}{columns such as \code{Sample_type}, \code{Treatment_1}, \code{Concentration_1}}
#' }
#' @details The object is created in \code{data-raw/02-build-new_demo.R} and
#'   saved with \code{usethis::use_data(new_demo, compress = "xz")}.
#' @seealso \code{\link{mini_mac}}
#' @examples
#' data(demo)
#' demo
"demo"
