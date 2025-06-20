% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_mds.R
\name{plot_mds}
\alias{plot_mds}
\title{Plot MDS dimensionality reduction}
\usage{
plot_mds(
  data = NULL,
  group_by = NULL,
  label = NULL,
  max_overlaps = NULL,
  n_labels = 50
)
}
\arguments{
\item{data}{A tidyseurat object merged with metadata. Must contain columns
"Well_ID", "Row", "Column"}

\item{group_by}{A string specifying which column in data will be used to
color the samples.}

\item{label}{A string specifying which column in data will be used to
label a sample.}

\item{max_overlaps}{Maximum number of overlaps for ggrepel}

\item{n_labels}{An integer specifying number of labels to show, based on PC1 and PC2 extremeness. Default set to 50}
}
\value{
ggplot object
}
\description{
This function uses limma's plot MDS to visualise groupoing of the data
points. The function is used with gene.selection = "common" to increase the
speed.
}
\examples{
data("mini_mac")
p <- plot_mds(mini_mac, group_by = "Sample_type", label = "Sample_type", 
n_labels = 10, max_overlaps = 100)
}
