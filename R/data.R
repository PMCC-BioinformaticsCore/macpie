#' A small example MACseq plate
#'
#' A 384-cell × 300-gene subset of the PMMSq033 dataset, included for
#' runnable examples and fast CRAN checks.
#'
#' @format A \code{\link[Seurat]{Seurat}} object with one assay:
#'   \describe{
#'     \item{RNA}{counts slot only, no normalized data, with metadata}
#'   }
#' @usage data(mini_mac)
#' @source Prepared by subsetting \code{PMMSq033.rds} in \code{data/}
"mini_mac"
