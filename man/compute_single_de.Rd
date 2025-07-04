% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_single_de.R
\name{compute_single_de}
\alias{compute_single_de}
\title{Retrieve normalised counts of MAC-seq data}
\usage{
compute_single_de(
  data = NULL,
  treatment_samples = NULL,
  control_samples = NULL,
  method = NULL,
  batch = 1,
  k = 2,
  spikes = NULL
)
}
\arguments{
\item{data}{A tidyseurat object merged with metadata. Must contain columns
"Well_ID", "Row", "Column"}

\item{treatment_samples}{Value in the column "combined_id" representing replicates of treatment samples in the data}

\item{control_samples}{Value in the column "combined_id"  representing replicates of control samples in the data}

\item{method}{One of "Seurat_wilcox", "DESeq2", "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"}

\item{batch}{Either empty, a single value, or a vector corresponding to the
number of samples}

\item{k}{Parameter k for RUVSeq methods, check RUVSeq tutorial}

\item{spikes}{List of genes to use as spike controls}
}
\value{
Data frame of DE counts
}
\description{
This function retrieves  counts from a number of methods that are available
for normalisation, with the default of limma-voomLmFit.
}
\examples{
data(mini_mac)
treatment_samples="Staurosporine_0.1"
control_samples<-"DMSO_0"
top_table <- compute_single_de(mini_mac, treatment_samples, control_samples, method = "limma_voom")
}
