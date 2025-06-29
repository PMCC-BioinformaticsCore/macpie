% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_rle.R
\name{plot_rle}
\alias{plot_rle}
\title{Create an RLE Plot}
\usage{
plot_rle(
  data,
  barcodes = NULL,
  label_column = NULL,
  labels = NULL,
  log = TRUE,
  batch = NULL,
  normalisation = NULL,
  spikes = NULL
)
}
\arguments{
\item{data}{A tidyseurat object merged with metadata. Must contain columns
"Well_ID", "Row", "Column".}

\item{barcodes}{A vector of sample barcodes corresponding to
Cells(seurat_object).}

\item{label_column}{A metadata column name to group the barcodes.}

\item{labels}{A vector of labels of the same length as 'barcodes" to group
the barcodes.}

\item{log}{A logical value indicating whether data should be log-transformed.
Defaults to \code{TRUE}.}

\item{batch}{Either empty, a single value, or a vector corresponding to the
number of samples.}

\item{normalisation}{One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
"edgeR", "RUVg", "RUVs", "RUVr", "limma_voom", "zinb". If empty, defaults to raw.}

\item{spikes}{List of genes to use as spike controls in RUVg}
}
\value{
A ggplot object representing the RLE plot.
}
\description{
This function generates a Relative Log Expression (RLE) plot for visualizing
the distribution of expression data after normalization or log transformation.
}
\details{
The function performs the following steps:
\itemize{
\item Ensures integrity of input data
\item Log-transforms the data
\item Computes the RLE by subtracting the row medians from each value.
\item Creates a boxplot using ggplot2 to visualize the distribution of RLE values.
}
}
\examples{
data("mini_mac")
p <- plot_rle(mini_mac, label_column = "Row")

}
