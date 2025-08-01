% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_counts.R
\name{plot_counts}
\alias{plot_counts}
\title{Generate a box plot to show gene expression (CPM)}
\usage{
plot_counts(
  data = NULL,
  genes = NULL,
  group_by = NULL,
  treatment_samples = NULL,
  control_samples = NULL,
  color_by = NULL,
  normalisation = NULL,
  batch = 1
)
}
\arguments{
\item{data}{A tidyseurat object merged with metadata. Must contain columns
"Well_ID", "Row", "Column".}

\item{genes}{Genes to be plotted}

\item{group_by}{A column that specifies the treatment group in the input data}

\item{treatment_samples}{Value in the column "combined_id" representing replicates of treatment samples in the data}

\item{control_samples}{Value in the column "combined_id"  representing replicates of control samples in the data}

\item{color_by}{A column that specifies the group coloring}

\item{normalisation}{One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
"edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"}

\item{batch}{To indicate patch factor}
}
\value{
a ggplot2 object
}
\description{
This is the function to generate a box plot to show CPM levels of DE genes
among selected treatment samples and control samples.
}
\examples{
data(mini_mac)
genes <- mini_mac@tools$diff_exprs$Staurosporine_10$gene[1:6]
p <- plot_counts(mini_mac, genes = genes, group_by = "combined_id", 
treatment_samples = "Staurosporine_10", 
control_samples = "DMSO_0",
normalisation = "clr")
}
