% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_multi_de.R
\name{plot_multi_de}
\alias{plot_multi_de}
\title{Generate heatmap of DE genes from multiple treatments
This is the function to generate a heatmap of DE genes from running compute_multi_DE
that shared by more than one treatment group. There are a few options available
to help you to extract shared DE genes.}
\usage{
plot_multi_de(
  data = NULL,
  group_by = NULL,
  value = NULL,
  p_value_cutoff = 0.01,
  direction = "both",
  n_genes = 10,
  control = "DMSO_0",
  by = "fc",
  gene_list = NULL
)
}
\arguments{
\item{data}{A tidyseurat object merged with metadata. Must contain columns
"Well_ID", "Row", "Column".}

\item{group_by}{A column that specifies the treatment group in the input data}

\item{value}{To use logCPM, log2FC or metric in the heatmap}

\item{p_value_cutoff}{Cutoff for adjusted p-value (column p_value_adj), default 0.01}

\item{direction}{Direction to select up or down regulated genes or in both directions}

\item{n_genes}{Top n genes to be extracted from each treatment comparison}

\item{control}{The control group to be included in the final heatmap, usually DMSO_0}

\item{by}{Extract top n genes by either absolute fold change or by adjusted p-value}

\item{gene_list}{External list of genes to plot the heatmap on}
}
\value{
a pheatmap object
}
\description{
Generate heatmap of DE genes from multiple treatments
This is the function to generate a heatmap of DE genes from running compute_multi_DE
that shared by more than one treatment group. There are a few options available
to help you to extract shared DE genes.
}
\examples{
data("mini_mac")
plot_multi_de(mini_mac,group_by = "combined_id", 
             value = "log2fC", p_value_cutoff = 0.01, direction="up", 
             n_genes = 10, control = "DMSO_0", by="fc")
}
