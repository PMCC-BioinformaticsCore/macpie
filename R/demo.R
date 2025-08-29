#' A tiny demo dataset for quick-start examples
#'
#' `demo` is a small Seurat object derived from the package's internal
#' example dataset `mini_mac`, reduced to a small number of genes
#' for fast examples in the quick start in README.
#'
#' @format A Seurat object with:
#' \describe{
#'   \item{assays$RNA}{counts and data matrices with ~200 genes}
#'   \item{meta.data}{columns such as \code{Sample_type}, \code{Treatment_1}, \code{Concentration_1}}
#' }
#' @details The object is created in \code{data-raw/demo.R} and
#'   saved with \code{usethis::use_data(demo, compress = "xz")}.
#' @seealso \code{\link{mini_mac}}
#' @usage data(demo)
"demo"
