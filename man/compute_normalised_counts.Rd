% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_normalised_counts.R
\name{compute_normalised_counts}
\alias{compute_normalised_counts}
\title{Retrieve normalised counts of MAC-seq data}
\usage{
compute_normalised_counts(
  data = NULL,
  method = NULL,
  batch = NULL,
  k = NULL,
  spikes = NULL,
  max_counts = NULL
)
}
\arguments{
\item{data}{A tidyseurat object merged with metadata. Must contain columns
"Well_ID", "Row", "Column"}

\item{method}{One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
"edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"}

\item{batch}{Either empty, a single value, or a vector corresponding to the
number of samples}

\item{k}{Parameter k for RUVSeq methods, check RUVSeq tutorial}

\item{spikes}{List of genes to use as spike controls}

\item{max_counts}{Maximum count for a gene across all samples}
}
\value{
Data frame of normalised counts
}
\description{
This function retrieves  counts from a number of methods that are available
for normalisation, with the default of limma-voomLmFit.
}
\examples{
data(mini_mac)
compute_normalised_counts(mini_mac)
}
